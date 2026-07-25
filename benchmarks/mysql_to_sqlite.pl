#!/usr/bin/env perl
#
# mysql_to_sqlite.pl --server host:port --db <mysql_db> --user U --password P \
#                    --sqlite /path/to/out.sqlite [--schema <sqliteschema>]
#
# Copy a populated PASA MySQL/MariaDB database into a fresh SQLite database so
# the same dataset can be benchmarked on both backends. The two schemas in
# schema/ declare an identical set of 27 tables, so this is a straight
# table-by-table row copy -- no name or column mapping needed.
#
# Rows are inserted inside a single transaction per table with a prepared
# statement, which is the difference between minutes and hours for the
# alignment table.
use strict;
use warnings;
use DBI;
use Getopt::Long;
use FindBin;

my ($server, $mysql_db, $user, $password, $sqlite_path, $schema, $help);
my $batch = 10_000;

GetOptions(
    'server=s'   => \$server,
    'db=s'       => \$mysql_db,
    'user=s'     => \$user,
    'password=s' => \$password,
    'sqlite=s'   => \$sqlite_path,
    'schema=s'   => \$schema,
    'batch=i'    => \$batch,
    'help|h'     => \$help,
) or die "bad options\n";

$schema ||= "$FindBin::Bin/../schema/cdna_alignment_sqliteschema";

if ($help || !($server && $mysql_db && defined $user && defined $password && $sqlite_path)) {
    die "usage: $0 --server host:port --db DB --user U --password P --sqlite OUT.sqlite [--schema F]\n";
}
die "ERROR: $sqlite_path already exists\n" if -e $sqlite_path;
die "ERROR: no schema at $schema\n" unless -s $schema;

my ($host, $port) = split(/:/, $server);
$port ||= 3306;

print STDERR "-- connecting to MySQL $mysql_db at $host:$port\n";
my $my_dbh = DBI->connect(
    "dbi:mysql:database=$mysql_db;host=$host;port=$port", $user, $password,
    { RaiseError => 1, AutoCommit => 1, mysql_use_result => 1 },
) or die "cannot connect to mysql: $DBI::errstr";

print STDERR "-- creating SQLite db at $sqlite_path\n";
my $lite_dbh = DBI->connect("dbi:SQLite:dbname=$sqlite_path", "", "",
    { RaiseError => 1, AutoCommit => 1 }) or die "cannot create sqlite: $DBI::errstr";

# Bulk-load pragmas. journal/synchronous are restored to the pipeline's normal
# settings at the end -- these are load-time only, not what gets benchmarked.
$lite_dbh->do("PRAGMA journal_mode = OFF");
$lite_dbh->do("PRAGMA synchronous = OFF");
$lite_dbh->do("PRAGMA cache_size = -200000");

print STDERR "-- applying schema $schema\n";
{
    open(my $fh, '<', $schema) or die "cannot read $schema: $!";
    my $sql = do { local $/; <$fh> };
    close $fh;
    # split on ';' at statement end; the PASA sqlite schema has no embedded ';'
    for my $stmt (split /;\s*\n/, $sql) {
        next unless $stmt =~ /\S/;
        next if $stmt =~ /^\s*--/;
        eval { $lite_dbh->do($stmt) };
        warn "  warn: statement failed: $@\n  ($stmt)\n" if $@;
    }
}

my @tables = do {
    open(my $fh, '<', $schema) or die $!;
    my @t;
    while (<$fh>) { push @t, $1 if /create\s+table\s+(\w+)/i; }
    close $fh;
    @t;
};
printf STDERR "-- copying %d tables\n", scalar @tables;

my $grand = 0;
for my $table (@tables) {
    # column order from the SQLite side, so INSERT matches the target layout
    my $cols_sth = $lite_dbh->prepare("PRAGMA table_info($table)");
    $cols_sth->execute;
    my @cols;
    while (my $r = $cols_sth->fetchrow_hashref) { push @cols, $r->{name}; }
    unless (@cols) { warn "  skip $table: not present in sqlite db\n"; next; }

    my $collist = join(",", map { "`$_`" } @cols);
    my $sel = eval { $my_dbh->prepare("select $collist from `$table`") };
    if ($@) { warn "  skip $table: $@"; next; }
    eval { $sel->execute };
    if ($@) { warn "  skip $table: $@"; next; }

    # The sqlite schema seeds some tables (URL_templates, URL_var_names) with
    # default rows via INSERT statements. Copying MySQL's rows on top of those
    # collides on UNIQUE constraints, so clear first -- the goal is for SQLite
    # to mirror MySQL exactly, not to merge with the schema's defaults.
    $lite_dbh->do("delete from $table");

    my $placeholders = join(",", ("?") x @cols);
    my $ins = $lite_dbh->prepare(
        "insert into $table (" . join(",", @cols) . ") values ($placeholders)");

    my $n = 0;
    $lite_dbh->begin_work;
    while (my $row = $sel->fetchrow_arrayref) {
        $ins->execute(@$row);
        if (++$n % $batch == 0) {
            $lite_dbh->commit;
            $lite_dbh->begin_work;
            print STDERR "\r  $table: $n rows";
        }
    }
    $lite_dbh->commit;
    $grand += $n;
    printf STDERR "\r  %-32s %10d rows\n", $table, $n;
}

# Leave the file in the state the pipeline expects to run against.
$lite_dbh->do("PRAGMA journal_mode = WAL");
$lite_dbh->do("PRAGMA synchronous = NORMAL");
$lite_dbh->do("ANALYZE");

$my_dbh->disconnect;
$lite_dbh->disconnect;
printf STDERR "-- done: %d rows total -> %s\n", $grand, $sqlite_path;
