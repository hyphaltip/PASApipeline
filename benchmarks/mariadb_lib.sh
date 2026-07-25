#!/bin/bash
# mariadb_lib.sh -- shared helpers for bringing a completed PASA run's
# MariaDB datadir back up in a Singularity instance.
#
# Source this; it defines:
#   stage_datadir <training_dir> <workdir>   -> copies datadir/conf/genome
#   mariadb_start <workdir>                  -> sets PORT, DB_USER, DB_PASS
#   mariadb_stop                             -> stops the instance
#
# Reads/sets: MARIADB_SIF, INSTANCE_NAME, PORT, DB_USER, DB_PASS, MYHOSTNAME
#
# Two things make this fiddlier than run_bench_mysql.sh's fresh-datadir case:
#  * the production datadir's root account has a password we do not hold, and
#    its user grant is scoped to whichever node that run landed on
#    ('jstajich'@'r22'), so we must add a '%'-scoped grant via a first pass
#    with skip-grant-tables set in my.cnf ([mysqld] -- passing it as a
#    mysqld_safe argument through `singularity instance start` does NOT reach
#    mysqld);
#  * the credentials themselves come from the run's own preserved
#    pasa_conf file rather than being minted.

MARIADB_SIF="${MARIADB_SIF:-/bigdata/stajichlab/shared/lib/mariadb/mariadb.sif}"

stage_datadir() {
    local training_root="$1" workdir="$2"
    local src="$training_root/training"

    [ -d "$src/mysql_db/db" ] || { echo "ERROR: no datadir at $src/mysql_db/db" >&2; return 1; }
    [ -s "$src/genome.fasta" ] || { echo "ERROR: no genome at $src/genome.fasta" >&2; return 1; }

    echo "=== Staging scratch copy into $workdir ==="
    rm -rf "$workdir"
    mkdir -p "$workdir/mysql_db"
    cp -rL "$src/mysql_db/db"   "$workdir/mysql_db/db"
    cp -rL "$src/mysql_db/conf" "$workdir/mysql_db/conf"
    cp -L  "$src/genome.fasta"  "$workdir/genome.fasta"
    chmod -R u+w "$workdir"
}

_wait_for_mysql() {
    local i
    for i in $(seq 1 60); do
        if singularity exec instance://"$INSTANCE_NAME" \
             mysqladmin --socket=/var/run/mysqld/mysqld.sock ping >/dev/null 2>&1; then
            echo "MariaDB up after ${i}s"; return 0
        fi
        sleep 1
    done
    echo "ERROR: MariaDB did not start" >&2; return 1
}

mariadb_stop() { singularity instance stop "$INSTANCE_NAME" 2>/dev/null || true; }

mariadb_start() {
    local workdir="$1" db_name="$2"

    INSTANCE_NAME="pasa_sweep_$$"
    PORT=$(shuf -i 14000-14999 -n1)
    MYHOSTNAME=$(hostname -s)

    local src_pasaconf
    src_pasaconf=$(ls "$workdir"/mysql_db/conf/pasa-local*.config.txt 2>/dev/null | head -1)
    [ -n "$src_pasaconf" ] || { echo "ERROR: no preserved pasa-local*.config.txt" >&2; return 1; }
    DB_USER=$(sed -n 's/^MYSQL_RW_USER=//p'     "$src_pasaconf" | head -1)
    DB_PASS=$(sed -n 's/^MYSQL_RW_PASSWORD=//p' "$src_pasaconf" | head -1)
    [ -n "$DB_USER" ] && [ -n "$DB_PASS" ] || { echo "ERROR: no credentials in $src_pasaconf" >&2; return 1; }
    echo "Reusing preserved credentials for user '$DB_USER'"

    sed -i "s/^port = .*/port = $PORT/" "$workdir/mysql_db/conf/my.cnf"
    sed -i "s/^user = .*/user = $USER/" "$workdir/mysql_db/conf/my.cnf"
    sed -i "s/^bind-address.*/bind-address = 0.0.0.0/" "$workdir/mysql_db/conf/my.cnf" || true

    export SINGULARITY_BINDPATH="$workdir/mysql_db/db:/var/lib/mysql,$workdir/mysql_db/conf/my.cnf:/etc/mysql/my.cnf"

    echo "=== Pass 1: skip-grant-tables, adding host-independent grant ==="
    sed -i "/^\[mysqld\]/a skip-grant-tables" "$workdir/mysql_db/conf/my.cnf"
    singularity instance start --writable-tmpfs "$MARIADB_SIF" "$INSTANCE_NAME" /usr/bin/mysqld_safe
    _wait_for_mysql || return 1

    singularity exec instance://"$INSTANCE_NAME" \
      mysql --socket=/var/run/mysqld/mysqld.sock -u root <<SQL
FLUSH PRIVILEGES;
CREATE USER IF NOT EXISTS '$DB_USER'@'%' IDENTIFIED BY '$DB_PASS';
ALTER USER '$DB_USER'@'%' IDENTIFIED BY '$DB_PASS';
GRANT ALL PRIVILEGES ON \`$db_name\`.* TO '$DB_USER'@'%';
FLUSH PRIVILEGES;
SQL

    mariadb_stop
    sleep 5

    echo "=== Pass 2: restarting normally on port $PORT ==="
    sed -i "/^skip-grant-tables$/d" "$workdir/mysql_db/conf/my.cnf"
    singularity instance start --writable-tmpfs "$MARIADB_SIF" "$INSTANCE_NAME" /usr/bin/mysqld_safe
    _wait_for_mysql || return 1
}
