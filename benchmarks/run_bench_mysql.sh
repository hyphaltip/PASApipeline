#!/bin/bash
# run_bench_mysql.sh <run_dir> <cpu>
#
# Same benchmark as run_bench.sh, but against a MySQL/MariaDB backend
# spun up fresh in a Singularity container per run -- mirrors the
# production pattern in Fungi_BFD's nextflow/funannotate.nf
# (FUNANNOTATE_TRAIN / RNASEQ_PREPARE processes): a per-job mariadb.sif
# instance with its own datadir, port, and PASA conf.txt, torn down after.
#
# Unlike production (which rsyncs a pre-populated datadir with an
# already-granted user), this initializes a FRESH, isolated datadir via
# mysql_install_db each run, then grants a dedicated pasa_bench user --
# no dependency on any account-specific state.
#
# Usage:
#   ./run_bench_mysql.sh run_mysql_baseline 4
set -euo pipefail

FUNANNOTATE_ENV="/opt/linux/rocky/8.x/x86_64/pkgs/funannotate/1.9/.pixi/envs/default"
export PATH="$FUNANNOTATE_ENV/bin:$PATH"

MARIADB_SIF="/bigdata/stajichlab/shared/lib/mariadb/mariadb.sif"

RUN_DIR="$1"
CPU="${2:-4}"

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
RUN_PATH="$REPO_ROOT/$RUN_DIR"

rm -rf "$RUN_PATH"
cp -r "$REPO_ROOT/sample_data" "$RUN_PATH"
cd "$RUN_PATH"

module load singularity

INSTANCE_NAME="pasa_bench_mysql_$$"
MYSQL_SCRATCH="$RUN_PATH/mysql_db"
mkdir -p "$MYSQL_SCRATCH/db" "$MYSQL_SCRATCH/conf"

PORT=$(shuf -i 13000-13999 -n1)
MYHOSTNAME=$(hostname -s)
DB_USER="pasa_bench"
DB_PASS="pasa_bench_pw"
DB_NAME="pasa_bench_sample"

echo "=== Initializing fresh MariaDB datadir at $MYSQL_SCRATCH/db ==="
singularity exec -B "$MYSQL_SCRATCH/db:/var/lib/mysql" "$MARIADB_SIF" \
  mysql_install_db --datadir=/var/lib/mysql --auth-root-authentication-method=normal >/dev/null

cp /rhome/jstajich/.pasa/pasa_conf/my.cnf "$MYSQL_SCRATCH/conf/my.cnf"
sed -i "s/^port = .*/port = $PORT/" "$MYSQL_SCRATCH/conf/my.cnf"
sed -i "s/^bind-address.*/bind-address = 0.0.0.0/" "$MYSQL_SCRATCH/conf/my.cnf" || true

stop_mysqldb() { singularity instance stop "$INSTANCE_NAME" 2>/dev/null || true; }
trap "stop_mysqldb; exit 130" SIGHUP SIGINT SIGTERM
trap "stop_mysqldb" EXIT

echo "=== Starting MariaDB instance '$INSTANCE_NAME' on port $PORT ==="
export SINGULARITY_BINDPATH="$MYSQL_SCRATCH/db:/var/lib/mysql,$MYSQL_SCRATCH/conf/my.cnf:/etc/mysql/my.cnf"
singularity instance start --writable-tmpfs "$MARIADB_SIF" "$INSTANCE_NAME" /usr/bin/mysqld_safe

echo "=== Waiting for MariaDB to accept connections ==="
for i in $(seq 1 30); do
  if singularity exec instance://"$INSTANCE_NAME" mysqladmin --socket=/var/run/mysqld/mysqld.sock ping >/dev/null 2>&1; then
    echo "MariaDB is up after ${i}s"
    break
  fi
  sleep 1
  if [ "$i" -eq 30 ]; then
    echo "ERROR: MariaDB did not come up in time" >&2
    exit 1
  fi
done

echo "=== Creating dedicated user + database ==="
singularity exec instance://"$INSTANCE_NAME" mysql --socket=/var/run/mysqld/mysqld.sock -u root <<SQL
CREATE DATABASE IF NOT EXISTS \`$DB_NAME\`;
CREATE USER IF NOT EXISTS '$DB_USER'@'%' IDENTIFIED BY '$DB_PASS';
GRANT ALL PRIVILEGES ON \`$DB_NAME\`.* TO '$DB_USER'@'%';
FLUSH PRIVILEGES;
SQL

# PASA's own DBI connect wants the DATABASE parsed to look like a MySQL
# reference (configure_db_driver: unqualified-path = mysql), and connects
# to MYSQLSERVER=host:port read from conf.txt via --PASACONF.
PASACONF="$MYSQL_SCRATCH/conf/pasa-local.config.txt"
cat > "$PASACONF" <<CONF
MYSQLSERVER=${MYHOSTNAME}:${PORT}
MYSQL_RW_USER=${DB_USER}
MYSQL_RW_PASSWORD=${DB_PASS}
MYSQL_RO_USER=${DB_USER}
MYSQL_RO_PASSWORD=${DB_PASS}
PASA_ADMIN_EMAIL=bench@localhost
PASA_ADMIN_DB=PASA2_admin
USE_PASA_DB_SETUP_HOOK=false
BASE_PASA_URL=http://localhost/cgi-bin/
HOOK_PERL_LIBS=__PASAHOME__/SAMPLE_HOOKS
HOOK_GENE_STRUCTURE_UPDATER=GFF3::GFF3_annot_updater::get_updater_obj
CONF
export PASACONF

# Point the DATABASE value at the mysql db name (unqualified -> mysql driver)
sed -i "s#^DATABASE=.*#DATABASE=$DB_NAME#" mysql.confs/alignAssembly.config mysql.confs/annotCompare.config

echo "=== Running PASA align-assembly benchmark in $RUN_DIR against MySQL (CPU=$CPU) ==="
/usr/bin/time -v ./__run_sample_pipeline.pl \
  --align_assembly_config mysql.confs/alignAssembly.config \
  --annot_compare_config mysql.confs/annotCompare.config \
  --TRANSDECODER --ALIGNERS minimap2 --CPU "$CPU" --just_align_assembly \
  --stringent_alignment_overlap 30 \
  > pasa_run.stdout.log 2> pasa_run.log

echo "=== Done. Parsing per-command timing from pasa_run.log ==="
python3 "$REPO_ROOT/benchmarks/parse_bench_log.py" pasa_run.log
