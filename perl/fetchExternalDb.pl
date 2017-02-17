#!/usr/bin/perl

# fetchExternalDb.pl
# A Perl script/wrapper which fetches the necessary external databases (GO and Proteome
# databases) in order to have a more complete and functional version of the browser. The
# Proteome databases are required for mm9
# Author      : Panagiotis Moulos (moulos@fleming.gr)
# Created     : 12 - 09 - 2012 (dd - mm - yyyy)
# Last Update : 14 - 09 - 2012 (dd - mm - yyyy)
# Version     : 1.0

# INFO: Information about what is going on in this script can be found in:
#       1. /media/HD2/kent/src/hg/makeDb/doc/uniProt/sp120323.txt
#       2. /media/HD2/kent/src/hg/makeDb/doc/proteins/120806.txt
#       3. The Gene Ontology site where we can download raw MySQL dumps
# INFO: The proteome database is too hard to build... We are dumping it directly from
#       UCSC... Try not to do this very often
# TODO: Remove the constant paths from within this file... Make another YAML at some point...
#       I am just not doing it now because I will be so confused...
# TODO: Better organize code... The fetch_ and build_ functions are almost identical

use strict;
use IO::File;
use Getopt::Long;
use File::Temp;
use File::Spec;
use File::Basename;
use File::Path qw(make_path remove_tree);
use Archive::Extract;
use LWP::Simple;
use DBI;

use constant MAXCORES => 6;
use constant PATH_TO_SPTODB => "/media/HD2/kenthome/bin/x86_64/spToDb";
use constant PATH_TO_SPDBSC => "/media/HD2/kent/src/hg/proteins/spToDb/spToDb.sql";
use constant PATH_TO_SPTOVS => "/media/HD2/kenthome/bin/x86_64/spDbAddVarSplice";
use constant NAME_OF_CURPRO => "proteins111004";
use constant NAME_OF_CURUNI => "uniProt";
use constant NAME_OF_CURGOT => "go";

# Make sure output is unbuffered
select(STDOUT);
$|=1;

# This script must be run as root
die "\n\nYou must run this script as root or sudoer!\n\n" if ($> != 0);

# On Ctrl-C, do cleanup
$SIG{INT} = \&catch_cleanup;

# Set defaults
our $scriptname = "fetchExternalDb.pl";
our @dbs;           # Which databases to fetch? Can be go, proteome, uniprot
our $source;        # Local dir or download
our @dbdata;        # Username and password for the DB to avoid hardcoding
our $ncore = 1; # Number of cores for parallel MySQL build
our $log = 0;       # Keep log?
our $silent = 0;    # Display verbose messages
our $help = 0;      # Help?

# Global variables
our $tmpdir = File::Temp->newdir(); # Temp directory for this session
our $logfilehandle = &open_log_file if ($log); # Log file if requested

# Check for the presence of YAML, required!!!
&try_module("Net::FTP");

# Advertise
&advertise;

# Check inputs
&check_inputs;

# Record progress...
my $date = &now;
disp("$date - Started...\n");

# Do the job...
if (&smatch("uniprot",@dbs))
{
    my $tmpuni = &fetch_uniprot;
    &build_uniprot($tmpuni);
    &grant_access("gbuser","uniProt");
}
if (&smatch("proteome",@dbs))
{
    my $tmpome = &fetch_proteome;
    &build_proteome($tmpome);
    &grant_access("gbuser","proteome");
}
if (&smatch("go",@dbs))
{
    my $tmpgo = &fetch_go;
    &build_go($tmpgo);
    &grant_access("gbuser","go");
}

$date = &now;
disp("\n$date - Finished!\n");

# Close the log file if requested
&close_log_file if ($log);

# Process inputs
sub check_inputs
{
    my $stop;
    GetOptions("database|d=s{,}" => \@dbs,
               "source|i=s" => \$source,
               "dbdata|p=s{,}" => \@dbdata,
               "ncore|n=i" => \$ncore,
               "log|l" => \$log,
               "silent|s" => \$silent,
               "help|h" => \$help);
    # Check if the required arguments are set
    if ($help)
    {
        &program_usage;
        exit;
    }
    $stop .= "--- Please provide either a local data source or the word download ---" if (!$source);
    $stop .= "--- Please provide database connection data ---\n" if (!@dbdata);
    $stop .= "--- --dbdata should be consisted of two strings! ---\n"
        if (@dbdata && $#dbdata+1 != 2);
    if (!@dbs)
    {
        disp("Databases to fetch not given... Defaults to all (proteome, uniprot, go)...");
    }
    else
    {
        foreach my $db (@dbs)
        {
            if ($db !~ m/(proteome)|(uniprot)|(go)/i)
            {
                $stop = "--- The databases must consist only of the strings proteome, uniprot, go ---";
                last;
            }
        }
    }
    if ($stop)
    {
            print "\n$stop\n";
            print "Type perl $scriptname --help for help in usage.\n\n";
            exit;
    }
    if ($ncore > 1)
    {
        eval "&tryModule(\"Parallel::Loops\")";
        if ($@)
        {
            disp("Module Parallel::Loops not found, proceeding with one core...");
            $ncore = 1;
        }
        else { use Parallel::Loops; }
        if ($ncore > MAXCORES)
        {
            my $c = MAXCORES;
            disp("The maximum number of cores allowed is $c...");
            $ncore = MAXCORES;
        }
    }
    else { $ncore = 1; }
    disp("Script log will be written in fcdb_YYYYMMDDHHMMSS.log...") if ($log);
}

sub fetch_uniprot
{
    my $tmpsql = File::Spec->catfile($tmpdir,"uniprot.sql");
    if ($source =~ m/download/i)
    {
        my $cmd = "mysqldump --host=genome-mysql.cse.ucsc.edu --user=genome --lock-tables=false --verbose ".NAME_OF_CURUNI." > $tmpsql";
        my $cmdBad = system($cmd);
        die "\nSomething went wrong with mysqldump!!!\n" if ($cmdBad);
    }
    else
    {
        (-d $source) ? ($tmpsql = File::Spec->catfile($source,"uniprot.sql")) :
            (die "\nThe source directory is not existent! Please use \"download as source\"");
    }
    return($tmpsql);
}

sub build_uniprot
{
    disp("Creating the uniprot database... Long procedure...");
    &create_db("uniProt");
    my $cmd = "mysql -u$dbdata[0] -p$dbdata[1] uniProt < $_[0]";
    my $bad = system($cmd);
    die "\nSomething went wrong with uniprot building!!!\n" if ($bad);
}

sub fetch_proteome
{
    my $tmpsql = File::Spec->catfile($tmpdir,"proteome.sql");
    if ($source =~ m/download/i)
    {
        my $cmd = "mysqldump --host=genome-mysql.cse.ucsc.edu --user=genome --lock-tables=false --verbose ".NAME_OF_CURPRO." > $tmpsql";
        my $cmdBad = system($cmd);
        die "\nSomething went wrong with mysqldump!!!\n" if ($cmdBad);
    }
    else
    {
        (-d $source) ? ($tmpsql = File::Spec->catfile($source,"proteome.sql")) :
            (die "\nThe source directory is not existent! Please use \"download as source\"");
    }
    return($tmpsql);
}

sub build_proteome
{
    disp("Creating the proteome database... Long procedure...");
    &create_db("proteome");
    my $cmd = "mysql -u$dbdata[0] -p$dbdata[1] proteome < $_[0]";
    my $bad = system($cmd);
    die "\nSomething went wrong with proteome building!!!\n" if ($bad);
}

sub fetch_go
{
    my $tmpsql = File::Spec->catfile($tmpdir,"go.sql");
    if ($source =~ m/download/i)
    {
        my $cmd = "mysqldump --host=genome-mysql.cse.ucsc.edu --user=genome --lock-tables=false --verbose ".NAME_OF_CURGOT." > $tmpsql";
        my $cmdBad = system($cmd);
        die "\nSomething went wrong with mysqldump!!!\n" if ($cmdBad);
    }
    else
    {
        (-d $source) ? ($tmpsql = File::Spec->catfile($source,"go.sql")) :
            (die "\nThe source directory is not existent! Please use \"download as source\"");
    }
    return($tmpsql);
}

sub build_go
{
    disp("Creating the go database... Long procedure...");
    &create_db("go");
    my $cmd = "mysql -u$dbdata[0] -p$dbdata[1] go < $_[0]";
    my $bad = system($cmd);
    die "\nSomething went wrong with go building!!!\n" if ($bad);
}

sub load_txt_table
{
    my ($file,$table,$database) = @_;
    disp("Loading table $table in database $database from file $file...");
    my $conn = &open_connection($database);
    my $q = "LOAD DATA INFILE '$file' INTO TABLE $table FIELDS TERMINATED BY '\t' LINES TERMINATED BY '\n'";
    $conn->do($q);
    &close_connection($conn);
}

sub load_sql_table
{
    my ($file,$database) = @_;
    disp("Loading table in database $database from file $file...");
    my $cmd = "mysql -u$dbdata[0] -p$dbdata[1] -v $database < $file";
    my $bad = system($cmd);
    die "\nSomething went wrong with loading of file $file in $database!\n" if ($bad);
}

sub create_db
{
    my $db = shift @_;
    my $crq = "CREATE DATABASE `$db` DEFAULT CHARACTER SET utf8 COLLATE utf8_general_ci;";
    my $drq = "DROP DATABASE `$db`";
    my $conn = DBI->connect("dbi:mysql:;host=localhost;port=3306",$dbdata[0],$dbdata[1]);   
    if (&check_existence($db)) # Drop if already exists
    {
        $conn->do($drq);
        $conn->do($crq);
    }
    else { $conn->do($crq); }
    &close_connection($conn);
}

sub check_existence
{
    my $dbcheck = shift @_;
    my $out = 1;
    my $conn = DBI->connect("dbi:mysql:database=information_schema;host=localhost;port=3306",$dbdata[0],$dbdata[1]);
    my $query = "SELECT `SCHEMA_NAME` FROM `SCHEMATA` WHERE `SCHEMA_NAME` = \"$dbcheck\"";
    my $sth = $conn->prepare($query);
    $sth->execute();
    $out = 0 if (!$sth->rows());
    $sth->finish();
    &close_connection($conn);
    return($out);
}

sub grant_access
{
    my ($u,$d) = @_;
    my $conn = &open_connection;
    my $q = "GRANT SELECT, INSERT, UPDATE, DELETE, INDEX, CREATE, DROP, ALTER, ".
             "CREATE TEMPORARY TABLES ON $d.* TO '$u'\@'localhost'; FLUSH PRIVILEGES;";
    $conn->do($q);
    &close_connection($conn);
}

sub open_connection
{   
    my $database = shift @_;
    my $hostname = "localhost";
    my $conn = DBI->connect("dbi:mysql:database=$database;host=$hostname;port=3306",$dbdata[0],$dbdata[1]);
    return $conn;
}

sub close_connection
{ 
    my $conn = shift @_;
    $conn->disconnect();
}

sub open_log_file
{
    my $logname = "fcdb_".&now.".log";
    my $lfh = IO::File->new();
    $lfh->open(">$logname");
    return($lfh);
}

sub close_log_file { $logfilehandle->close; }

sub now
{
    my $format = shift @_;
    $format = "human" if (!$format);
    my ($sec,$min,$hour,$day,$month,$year) = localtime(time);
    $year += 1900;
    $month++;
    $month = "0".$month if (length($month)==1);
    $day = "0".$day if (length($day)==1);
    $hour = "0".$hour if (length($hour)==1);
    $min = "0".$min if (length($min)==1);
    $sec = "0".$sec if (length($sec)==1);
    ($format ne "machine") ? (return($day."/".$month."/".$year." ".$hour.":".$min.":".$sec)) :
    (return($year.$month.$day.$hour.$min.$sec));
}

sub smatch
{
    my ($self,$s,@a) = @_;
    if (grep(/^$s$/,@a))
    {
        return(1);
    }
    return(0);
}

sub catch_cleanup 
{
    print STDERR "\nCatching ctrl-C, cleaning temporary files!\n";
    &cleanup;
    die;
}

sub cleanup 
{
    remove_tree($tmpdir);
}

sub advertise
{
    disp($scriptname);
    disp("Local UCSC Genome Browser installation helper for external databases...\nCopyright: Panagiotis Moulos (moulos\@fleming.gr)\n");
}

sub disp
{
    print "\n@_" if (!$silent);
    print $logfilehandle "\n@_" if (!$silent && $log);
}

sub try_module
{
    my $module = shift @_;
    eval "require $module";
    if ($@)
    {
        my $killer = "Module $module is required to continue with the execution. If you are in\n". 
                     "Windows and you have ActiveState Perl installed, use the Package Manager\n".
                     "to get the module. If you are under Linux, log in as a super user (or use\n".
                     "sudo under Ubuntu) and type \"perl -MCPAN -e shell\" (you will possibly have\n".
                     "to answer some questions). After this type \"install $module\" to install\n".
                     "the module. If you don't know how to install the module, contact your\n".
                     "system administrator.";
        die "\n$killer\n\n";
    }
}

sub init_local_wb
{
    my $z = " "x$_[0];
    print "|$z|";
    print "\b"x($_[0]+1);
}

sub end_local_wb { print "\n"; }

sub program_usage 
{
    # The look sucks here but it is actually good in the command line
    my $usagetext = << "END";
    
$scriptname
Local UCSC Genome Browser installation helper for external databases.

Author : Panagiotis Moulos (moulos\@fleming.gr)

Main usage
$scriptname [OPTIONS]

--- Required ---
  --source|i        Data source for loading the external databases. It
            should be either a directory just with the raw downloaded
            data files, or the word "download".
            It should be provided given the --build switch.
  --dbdata|c        Connection data for the local database. It should be a
            vector of length two containing a username and a password.
--- Optional ---
  --database|d      Which one(s) of the three supported databases to
            download and build. Can be proteome, uniprot, go. For example
            --database proteome go.
  --silent|s        Use this option if you want to turn informative 
            messages off.
  --help|h      Display this help text.
    
This program fetches several external databases for the optimal function
of the UCSC Genome Browser.

END
    print $usagetext;
    exit;
}


# OTHER VERSIONS OF FETCHING UNIPROT... THEY MIGHT BE USEFUL AS CODE CHUNKS

sub fetch_uniprot_2
{
    my ($f,$ftp,$tmpuni,$tmpzip,$fsize,$cmd,$tpath,$ok);
    if ($source =~ m/download/i) # Case where we download from HTTP
    {
        use Net::FTP;

        my @filestofetch = (
            "/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz",
            "/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.dat.gz",
            "/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot_varsplic.fasta.gz"
        );

        $ftp = Net::FTP->new("ftp.expasy.org",Debug => 0) or die "\nCannot connect to FTP: $@","\n";
        $ftp->login("anonymous") or die "\nCannot login ",$ftp->message,"\n";
        $ftp->binary; # Saves a lot of trouble...

        $tmpuni = File::Spec->catdir($tmpdir,"uniprot");
        foreach $f (@filestofetch)
        {
            disp("Getting file $f from expasy server...");
            &init_local_wb(50);
            $tmpzip = File::Spec->catfile($tmpuni,$f);
            $fsize = $ftp->size($f);
            $ftp->hash(\*STDOUT,$fsize/50);
            $ftp->get($f,$tmpzip) or die "\nFTP get failed ",$ftp->message,"\n";
            &end_local_wb;
        }
    }
    else
    {
        (-d $source) ? ($tmpuni = $source) :
            (die "\nThe source directory is not existent! Please use \"download as source\"");
    }
    return($tmpuni);
}

sub build_uniprot_2
{
    # Build the text files
    my $tmpuni = $_[0];
    disp("Creating the uniProt database tables...");
    my $tabpath = make_path(File::Spec->catdir($tmpuni,"tabfiles"));
    my $cmd = "zcat $tmpuni/*.dat.gz | ".PATH_TO_SPTODB." stdin $tabpath";
    my $ok = system($cmd);

    # Create the database tables
    &create_db("uniProt");
    &create_uniprot_tables;
    &load_uniprot_tables($tabpath);
    &patch_uniprot_tables($tmpuni);
}

sub create_uniprot_tables
{
    my $conn = &open_connection("uniProt");
    # Suck in the tables file from kent source
    open(BIG,PATH_TO_SPDBSC);
    my $bigquery = do { local $/;  <BIG> };
    close(BIG);
    $conn->do($bigquery);
    &close_connection($conn);
}

sub load_uniprot_tables
{
    opendir(DIR,$_[0]) or die "\nUnable to open $_[0]: $!\n";
    my @files = grep { !/^\.\.?$/ } readdir(DIR);
    my ($file,$table);

    if ($ncore == 1)
    {
        foreach $file (@files)
        {
            $table = basename($file);
            $table =~ s/\.txt(\.gz)?$//;
            &load_txt_table($file,$table,"uniProt");
        }
    }
    else
    {
        my $pl = Parallel::Loops->new($ncore);
        disp("Parallely loading uniprot tables...");        
        $pl->foreach(\@files, sub {
            my $t = basename($_);
            $t =~ s/\.txt(\.gz)?$//;
            &load_txt_table($_,$t,"uniProt");
        });
    }
}

sub patch_uniprot_tables
{
    my $patch = File::Spec->catfile($_[0],"uniprot_sprot_varsplic.fasta.gz");
    my $cmd = "cat $patch | ".PATH_TO_SPTOVS." uniProt stdin .";
    my $ok = system($cmd);
    disp("Whooops, something went wrong with loading splice variants!") if ($ok);
}

sub load_go_tables
{
    opendir(DIR,$_[0]) or die "\nUnable to open $_[0]: $!\n";
    my @files = grep { !/^\.\.?$/ } readdir(DIR);
    my $file;

    if ($ncore == 1)
    {
        foreach $file (@files)
        {
            &load_sql_table($file,"go");
        }
    }
    else
    {
        my $pl = Parallel::Loops->new($ncore);
        disp("Parallely loading go tables...");     
        $pl->foreach(\@files, sub {
            &load_sql_table($_,"go");
        });
    }
}
