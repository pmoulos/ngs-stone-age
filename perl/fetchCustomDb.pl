#!/usr/bin/perl

# fetchCustomDb.pl
# A Perl script/wrapper which fetches the necessary tables to build a customized database
# for the UCSC Genome Browser, which is not minimal but not full either. The customized
# data are provided in a YAML file which is constructed after the user research with UCSC
# Table Browser to determine which tables are necessary for each feature of the desired
# organisms. This script supposes that a minimal installation of the browser exists and
# does NOT load the tables. Use loadDb.sh from kent scripts for that.
# Author      : Panagiotis Moulos (moulos@fleming.gr)
# Created     : 05 - 09 - 2012 (dd - mm - yyyy)
# Last Update : 07 - 02 - 2013 (dd - mm - yyyy)
# Version     : 1.1

# TODO: Support other ways of downloading, HTTP (LWP::Simple) and FTP.
# TODO: Hard-code the YAML text so a file can be created on the fly

# In version 1.1, having a hash of default tables is deprecated due to difficulties in maintaining
# both a hard-coded version of the default tables and a YAML file. Now, a YAML file is required or
# written directly when not given, used and the deleted. Added an rsync process for genbank data.

use strict;
use IO::File;
use Getopt::Long;
use File::Temp;
use File::Spec;
use File::Basename;
use File::Path qw(make_path remove_tree);
use List::Util qw(min max sum);
use Archive::Extract;
use LWP::Simple;
use DBI;

use constant MAXCORES => 6;
             
# Make sure output is unbuffered
select(STDOUT);
$|=1;

# This script must be run as root
die "\n\nYou must run this script as root or sudoer!\n\n" if ($> != 0);

# On Ctrl-C, do cleanup
$SIG{INT} = \&catch_cleanup;

# Set defaults
our $scriptname = "fetchCustomDb.pl";
our $paramfile; # YAML parameters file with several options
our $dryrun;        # Run rsync with --dry-run?
our $build = 0; # Build the downloaded tables?
our @reqorgs;       # Download some of the organisms in the configuration file?
our @gbdata;        # Download gbdb, database, genbank data?
our @dbdata;        # Username and password for the DB to avoid hardcoding
our $ncore = 1; # Number of cores for parallel MySQL build
our $log = 0;       # Keep log?
our $silent = 0;    # Display verbose messages
our $help = 0;      # Help?

# Global variables
our $phref; # Parameters hash
our $tmpdir = File::Temp->newdir(); # Temp directory for this session
our $logfilehandle = &open_log_file if ($log); # Log file if requested

# Check for the presence of YAML, required!!!
&try_module("YAML");
&try_module("File::Rsync");

# Advertise
&advertise;

# Check inputs
&check_inputs;

# Record progress...
my $date = &now;
disp("$date - Started...\n");

# Read the parameters file or load the defaults if not given
($paramfile) ? ($phref = &read_param_file) : ($phref = &load_default_params);
&fix_tables;

# Create the list of files in the temporary file
my ($goldentablelist,$gbdbtablelist) = &create_table_list;

# Initialize rsync object
my $rsyncGolden = File::Rsync->new({
    "archive" => 1,
    "compress" => 1,
    "partial" => 1,
    "recursive" => 1,
    "progress" => 1,
    "stats" => 1,
    "files-from" => $goldentablelist,
    "verbose" => 1,
    "human-readable" => 1,
    "dry-run" => $dryrun ? 1 : 0,
    "outfun" => \&outnow,
    "errfun" => \&errnow
});

my $rsyncGbdb = File::Rsync->new({
    "archive" => 1,
    "compress" => 1,
    "partial" => 1,
    "recursive" => 1,
    "progress" => 1,
    "stats" => 1,
    "files-from" => $gbdbtablelist,
    "verbose" => 1,
    "human-readable" => 1,
    "dry-run" =>  $dryrun ? 1 : 0,
    "outfun" => \&outnow,
    "errfun" => \&errnow
});

my $rsyncGenbank = File::Rsync->new({
    "archive" => 1,
    "compress" => 1,
    "partial" => 1,
    "recursive" => 1,
    "progress" => 1,
    "stats" => 1,
    "verbose" => 1,
    "human-readable" => 1,
    "dry-run" =>  $dryrun ? 1 : 0,
    "outfun" => \&outnow,
    "errfun" => \&errnow
});

# Do the job!
print "\n";
$rsyncGolden->exec({
        src => $phref->{"download"}.$phref->{"goldenserver"},
        dest => $phref->{"goldenlocal"}
    }) or die "\n".join("",@{$rsyncGolden->err})."\n";
#print "\n".join("",@{$rsyncGolden->out})."\n";
$rsyncGbdb->exec({  
        src => $phref->{"download"}.$phref->{"gbdbserver"},
        dest => $phref->{"gbdblocal"}
    }) or die "\n".join("",@{$rsyncGbdb->err})."\n";
#print "\n".join("",@{$rsyncGbdb->out})."\n";
$rsyncGenbank->exec({   
        src => $phref->{"download"}.$phref->{"genbankserver"},
        dest => $phref->{"genbanklocal"}
    }) or die "\n".join("",@{$rsyncGbdb->err})."\n";
#print "\n".join("",@{$rsyncGenbank->out})."\n";

if ($build)
{
    disp("Building the downloaded databases... It will take some time, please wait...");
    &build_dbs;
}

$date = &now;
disp("\n$date - Finished!\n");

# Close the log file if requested
&close_log_file if ($log);

# Process inputs
sub check_inputs
{
    my $stop;
    GetOptions("param|p=s" => \$paramfile,
               "dry|d" => \$dryrun,
               "build|b" => \$build,
               "someorg|m=s{,}" => \@reqorgs,
               "somedata|e=s{,}" => \@gbdata,
               "dbdata|c=s{,}" => \@dbdata,
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
    $stop .= "--- Please provide database connection data ---\n" if (!@dbdata && $build);
    $stop .= "--- --dbdata should be consisted of two strings! ---\n"
        if (@dbdata && $#dbdata+1 != 2 && $build);
    if ($stop)
    {
            print "\n$stop\n";
            print "Type perl $scriptname --help for help in usage.\n\n";
            exit;
    }
    if ($ncore)
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
    if (@reqorgs)
    {
        my @tmporgs = @reqorgs;
        my @validorgs = ("hg18","hg19","mm9","mm10","dm3");
        @reqorgs=();
        for (my $i=0; $i<@tmporgs; $i++)
        {
            $tmporgs[$i] = lc($tmporgs[$i]);
        }
        foreach my $rorg (@tmporgs)
        {
            if (!&smatch($rorg,@validorgs))
            {
                disp("Skipping unsupported organism: $rorg...");
                next;
            }
            else { push(@reqorgs,$rorg); }
        }
    }
    if (@gbdata)
    {
        my @tmpdata = @gbdata;
        my @validata = ("gbdb","table","genbank");
        @gbdata=();
        for (my $i=0; $i<@tmpdata; $i++)
        {
            $tmpdata[$i] = lc($tmpdata[$i]);
        }
        foreach my $gdat (@tmpdata)
        {
            if (!&smatch($gdat,@validata))
            {
                disp("Skipping unsupported data type: $gdat...");
                next;
            }
            else { push(@gbdata,$gdat); }
        }
    }
    disp("Parameter file not given! Will try to load defaults...") if (!$paramfile);
    disp("Script log will be written in fcdb_YYYYMMDDHHMMSS.log...") if ($log);
}

sub read_param_file
{
    use YAML qw(LoadFile Dump);
    my $pfh;
    eval
    {
        open($pfh,"<",$paramfile);
        $phref = LoadFile($pfh);
        close($pfh);
    };
    if ($@)
    {
        disp("Bad parameter file! Will try to load defaults...");
        $phref = &load_default_params();
    }
    return($phref);
}

sub create_table_list
{
    my ($org,$table);
    my $fgolden = File::Spec->catfile($tmpdir,"goldenbtables.txt");
    my $fgbdb = File::Spec->catfile($tmpdir,"gbdbtables.txt");
    #my $fgolden = "/media/HD2/scripts/goldentables.txt";
    #my $fgbdb = "/media/HD2/scripts/gbdbtables.txt";
    my $fhgolden = IO::File->new();
    my $fhgbdb = IO::File->new();
    $fhgolden->open(">$fgolden");
    $fhgbdb->open(">$fgbdb");
    foreach $org (keys(%{$phref->{"org"}}))
    {
        next if (@reqorgs && (!&smatch($org,@reqorgs)));
        foreach $table (@{$phref->{"org"}->{$org}->{"tables"}})
        {
            print $fhgolden $org.$phref->{"dbprefix"}."/".$table.$phref->{"dbsuffix"}."\n";
            print $fhgolden $org.$phref->{"dbprefix"}."/".$table.$phref->{"txsuffix"}."\n";
            #print $fhgolden $phref->{"goldenserver"}."/".$org.$phref->{"dbprefix"}."/".$table.$phref->{"dbsuffix"}."\n";
            #print $fhgolden $phref->{"goldenserver"}."/".$org.$phref->{"dbprefix"}."/".$table.$phref->{"txsuffix"}."\n";
        }
        foreach $table (@{$phref->{"org"}->{$org}->{"gbdb"}})
        {
            print $fhgbdb $org."/".$table."\n";
            #print $fhgbdb $phref->{"gbdbserver"}."/".$org."/".$table."\n";
        }
    }
    $fhgolden->close;
    $fhgbdb->close;
    return($fgolden,$fgbdb);
}

sub fetch_chrom_info
{
    my ($org,$url,$tmpzip,$ae,$ok,$unzipped,$line);
    my @cols;
    my %chrominfo;
    
    foreach $org (keys(%{$phref->{"splits"}}))
    {
        $url = "http://hgdownload.cse.ucsc.edu/goldenPath/".$org."/database/chromInfo.txt.gz";
        $tmpzip = File::Spec->catfile($tmpdir,$org."chromInfo.txt.gz");
        disp("Downloading chromosome info file for $org...");
        my $tmp = getstore($url,$tmpzip);
        disp("Uncompressing...");
        $ae = Archive::Extract->new(archive => $tmpzip);
        $ok = $ae->extract(to => $tmpdir) or die "\n",$ae->error,"\n";
        $unzipped = File::Spec->catfile($tmpdir,$ae->files->[0]);

        my $fh = IO::File->new();
        $fh->open("< $unzipped");
        while($line = <$fh>)
        {
            $line =~ s/\r|\n$//g;
            @cols = split(/\t/,$line);
            push(@{$chrominfo{$org}},$cols[0]);
        }
        $fh->close;
    }

    return(\%chrominfo);
}

sub fix_tables
{
    my ($org,$expr,$elem,$chr);
    my $chrominfo = &fetch_chrom_info;

    # Fix the per chromosome split ones according to the given problematic tables
    if (defined($phref->{"splits"}))
    {
        foreach $org (keys(%{$phref->{"splits"}}))
        {
            my (@tofix,@toadd);
            $expr = join("|",@{$phref->{"splits"}->{$org}->{"tables"}});
            @tofix = grep { $_ =~ /$expr/ } @{$phref->{"org"}->{$org}->{"tables"}};
            @{$phref->{"org"}->{$org}->{"tables"}} = grep { $_ !~ /$expr/ } @{$phref->{"org"}->{$org}->{"tables"}};

            foreach $elem (@tofix)
            {
                foreach $chr (@{$chrominfo->{$org}})
                {
                    push(@toadd,$chr."_".$elem);
                }
            }
            
            push(@{$phref->{"org"}->{$org}->{"tables"}},@toadd);
        }
    }

    # Remove the possible duplicates in the tables list
    foreach $org (keys(%{$phref->{"org"}}))
    {
        my %uh = &unique(@{$phref->{"org"}->{$org}->{"tables"}});
        if (scalar keys(%uh) < sum(values(%uh)))
        {
            @{$phref->{"org"}->{$org}->{"tables"}} = keys(%uh);
        }
    }
}

sub load_default_params
{
    # Code to write YAML file directly
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

sub unique
{
    my @list = @_;
    my (%seen,$item);
    foreach $item (@list) 
    {
        $seen{$item}++;
    }
    return(%seen);
}

sub catch_cleanup 
{
    print STDERR "\nCatching Ctrl-C, cleaning temporary files!\n";
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
    disp("Local UCSC Genome Browser installation helper... Copyright: Panagiotis Moulos (moulos\@fleming.gr)\n");
}

sub disp
{
    print "\n@_" if (!$silent);
    print $logfilehandle "\n@_" if (!$silent && $log);
}

sub outnow
{
    print $_[0];
    print $logfilehandle $_[0] if ($log);
}

sub errnow
{
    print $_[0];
    print $logfilehandle $_[0] if ($log);
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

sub build_dbs
{
    my ($db,$cmd,$notok);
    &drop_dbs; # Drop them, else loadDb.sh does not run!
    if ($ncore == 1)
    {
        foreach $db (keys(%{$phref->{"org"}}))
        {
            $cmd = "sh ".$phref->{"loadscript"}." ".$phref->{"environment"}." ".$db;
            $notok = system($cmd);
            if ($notok) { disp("Whooops, something went wrong while loading $db..."); }
        }
    }
    else
    {
        my %notoks;
        my @orgs = keys(%{$phref->{"org"}});
        my $pl = Parallel::Loops->new($ncore);
        $pl->share(\%notoks);
        disp("Parallely loading databases...");     
        $pl->foreach(\@orgs, sub {
            my $cmdlocal = "sh ".$phref->{"loadscript"}." ".$phref->{"environment"}." ".$_;
            $notoks{$_} = system($cmdlocal);
        });
    }
}

sub drop_dbs
{
    foreach my $db (keys(%{$phref->{"org"}}))
    {
        &drop_db($db);
    }
}

sub drop_db
{
    my $db = shift @_;
    my $drq = "DROP DATABASE `$db`";
    my $conn = DBI->connect("dbi:mysql:;host=localhost;port=3306",$dbdata[0],$dbdata[1]);
    $conn->do($drq) if (&check_existence($db)); # Drop if already exists
    $conn->disconnect();
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
    $conn->disconect();
    return($out);
}

sub program_usage 
{
    # The look sucks here but it is actually good in the command line
    my $usagetext = << "END";
    
$scriptname
Local UCSC Genome Browser installation helper.

Author : Panagiotis Moulos (moulos\@fleming.gr)

Main usage
$scriptname [OPTIONS]

--- Optional ---
  --param|p     A YAML parameters file containing the locations of the
            remote and local UCSC directories and which tables to fetch. See
            the fetchcustom.yml for a straightforward commented example.
  --dry|d       This switch defines whether to run rsync with --dry-run to
            determine the size of the files to be downloaded.
  --build|b     (Re)Build the databases from the downloaded tables.
  --dbdata|c        Connection data for the local database. It should be a
            vector of length two containing a username and a password.
            It should be provided given the --build switch.
  --silent|s        Use this option if you want to turn informative 
            messages off.
  --help|h      Display this help text.
    
This program fetches custom tables from the UCSC remote server and puts them
in the proper directories in the local installation. It does NOT build the
tables. You have to run the Kent loadDb.sh script for that. See the documentation
in kent source. This helper also supposes the PRE-EXISTING minimal installation
of UCSC Genome Browser.

END
    print $usagetext;
    exit;
}

## Template
#foreach $org (keys(%{$phref->{"org"}}))
#{
    ## Template for golden path downloads
    #$rsync->exec({
        #src => $phref->{"download"}.$phref->{"goldenserver"}."/$org".$phref->{"dbprefix"},
        #dest => $phref->{"goldenlocal"}."/$org".$phref->{"dbprefix"},
    #});
    ## Template for gbdb downloads
    #$rsync->exec({
        #src => $phref->{"download"}.$phref->{"gbdbserver"}."/$org",
        #dest => $phref->{"gbdblocal"}."/$org".$phref->{"dbprefix"},
    #});
#}
