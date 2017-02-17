#!/usr/bin/perl

# Description later... Works also with duplicates
#
# Author      : Panagiotis Moulos (pmoulos@eie.gr)
# Created     : 09 - 02 - 2009 (dd - mm - yyyy)
# Last Update : 16 - 02 - 2009 (dd - mm - yyyy)
# Version     : 1.0

use strict;
use Getopt::Long;
use File::Basename;

# Make sure output is unbuffered
select(STDOUT);
$|=1;

# Set defaults
our $scriptname = "easyintersect.pl";
our $fileA;		   	# Input file A 
our $fileB;			# Input file B
our $ucA = 0;		# Unique ID column for file A
our $ucB = 0;		# Unique ID column for file B
our @out;		 	# Output files
our $header = 0;	# Files contain header?
our $silent = 0;  	# Display verbose messages
our $help = 0;	   	# Help?

# Check inputs
&checkInputs;

my (%allA,%allB,%union,%isect);
my ($e,@curr);

# Suck everything in
disp("Reading $fileA...");
open (INA,$fileA) or die "\nThe file $fileA does not exist!\n";;
my $headerA = <INA> if ($header);
while (my $line = <INA>)
{
	$line =~ s/\r|\n$//g;
	@curr = split(/\t/,$line);
	#$allA{$curr[$ucA]} = join("\t",@curr[0..$uc-1])."\t".join("\t",@curr[$uc+1..$#curr]);
	$allA{$curr[$ucA]} = join("\t",@curr);
}
close(INA);
disp("Reading $fileB...");
open (INB,$fileB) or die "\nThe file $fileB does not exist!\n";;
my $headerB = <INB> if ($header);
while (my $line = <INB>)
{
	$line =~ s/\r|\n$//g;
	@curr = split(/\t/,$line);
	#$allB{$curr[$ucB]} = join("\t",@curr[0..$uc-1])."\t".join("\t",@curr[$uc+1..$#curr]);
	$allB{$curr[$ucB]} = join("\t",@curr);
}
close(INB);

# Union, intersection and only sets
foreach $e (keys(%allA))
{ 
	$union{$e} = $allA{$e}; 
}

foreach $e (keys(%allB)) 
{
    if ($union{$e}) 
    { 
    	$isect{$e} = $allA{$e}."\t".$allB{$e}; # Include all elements
    }
    else
	{
		$union{$e} = $allB{$e}; # Add from B
    }
    if ($isect{$e})
    {
    	delete $allA{$e};
    	delete $allB{$e};
    }
}

# Fix headers for proper output for intersection if they exist
my $fheader;
if ($header)
{
	#my @hB = split(/\t/,$headerB);
	#my $fheader = join("\t",@hB[0..$uc-1])."\t".join("\t",@hB[$uc+1..$#hB]);
	my $th = $headerA;
	$th =~ s/\r|\n$//g;
	$fheader = "$th\t$headerB";
}

# Write output
my ($bA,$bB);
$bA = fileparse($fileA,'\..*?');
$bB = fileparse($fileB,'\..*?');
foreach my $opt (@out)
{
	printHash($fileA,$fileB,$headerA,"UNION",%union) if ($opt eq "union");
	printHash($fileA,$fileB,$fheader,"INTERSECTION",%isect) if ($opt eq "intersection");
	printHash($fileA,$fileB,$headerA,"ONLY_$bA",%allA) if ($opt eq "onlyA");
	printHash($fileA,$fileB,$headerB,"ONLY_$bB",%allB) if ($opt eq "onlyB");
}
disp("Finished!\n");


# Process inputs
sub checkInputs
{
    my $stop;
    GetOptions("inA|a=s" => \$fileA,
    		   "inB|b=s" => \$fileB,
    		   "uniA|u=i" => \$ucA,
    		   "uniB|v=i" => \$ucB,
    		   "output|o=s{,}" => \@out,
    		   "header|d" => \$header,
    		   "silent|s" => \$silent,
    		   "help|h" => \$help);
    # Check if the required arguments are set
    if ($help)
    {
    	&programUsage;
    	exit;
    }
    $stop .= "--- Please specify input files ---\n" if (!($fileA || $fileB));
    if ($stop)
    {
            print "\n$stop\n";
            print "Type perl $scriptname --help for help in usage.\n\n";
            exit;
    }
    if ($ucA && $ucB)
    {
    	$ucA--;	
    	$ucB--;
    }
    if ($ucA && !$ucB)
    {
    	disp("Unique ID column for fileB not given... assuming the same as fileA...");
    }
    elsif (!$ucA && $ucB)
    {
    	disp("Unique ID column for fileA not given... assuming the same as fileB...");
    }
    # Check proper output format
    if (@out)
    {
		foreach my $c (@out)
		{
			if ($c ne "union" && $c ne "intersection" && $c ne "onlyA" && $c ne "onlyB")
			{
				my $msg = "WARNING! --output options should be one or more of \"union\", \"intersection\",".
				          " \"onlyA\", \"onlyB\"\n".
						  "Using default (\"intersection\")...";
				disp($msg);
				@out = ("intersection");
			}
		}
	}
	else
	{
		disp("Output file type not given... Using default (intersection)") ;
		@out = ("intersection");
	}
}

sub createOutputFile
{
	my $inA = shift @_;
	my $inB = shift @_;
	my $type = shift @_;
	my ($baseA,$dirA,$extA) = fileparse($inA,'\..*?');
	my $baseB = fileparse($inB,'\..*?');
	return($dirA.$baseA."_".$baseB."_".$type.$extA);
}

sub printHash
{
	my $infileA = shift @_;
	my $infileB = shift @_;
	my $he = shift @_;
	my $otype = shift @_;
	my %hash = @_;
	my $outfilename = &createOutputFile($infileA,$infileB,$otype);
	disp("Writing output in $outfilename...");
	open(OUTPUT,">$outfilename");
	print OUTPUT "$he" if ($he);
	foreach my $ck (sort keys(%hash))
	{
		print OUTPUT "$hash{$ck}\n";
	}
	close(OUTPUT);
}

sub disp
{
	print "\n@_" if (!$silent);
}

sub programUsage 
{
	# The look sucks here but it is actually good in the command line
	my $usagetext = << "END";
	
$scriptname
A perl program to perform set operations between any kind of delimited files as
long as they contain a column with an ID for each line. The ID does not have to
be unique as long as two lines with the same ID have exactly the same content.
Generally, the files should have the same structure (e.g. both have header lines,
have the same number of columns) and the ID column MUST be the same for both files,
else the behaviour of the script can be unpredictable. ATTENTION! This script does
NOT perform overlap in genomic regions, it operates based only on file element IDs.
Output can consist of up to four different files containing the union, intersection
of the two files as well as elements being only in fileA or fileB. The intersection
file contains the rest of the columns for both files.

Author : Panagiotis Moulos (pmoulos\@eie.gr)

Main usage
$scriptname --inA fileA --inB fileB [OPTIONS]

--- Required ---
  --inA|a  file  First input file
  --inB|b  file  Second input file
--- Optional ---
  --uniA|u		Use this option to provide the script the column number
			in fileA that contains element IDs. Defaults to the first
			column.
  --uniB|v		Use this option to provide the script the column number
			in fileB that contains element IDs. Defaults to the first
			column.
  --output|o		Use this option to determine which intersection output
			filetypes you wish to retrieve.	Possible choices are:
   				"union" for retrieving the union of the two files.
   				"intersection" for retrieving the union of the two
   				files.
				"onlyA" for retrieving elements only in fileA.
				"onlyB" for retrieving elements only in fileB.
  --header|d		Use this option if you have a header line in
			your input files. It will also be written in your output
			files. Should be the same! Defaults to no header.
  --silent|s		Use this option if you want to turn informative 
  			messages off.
  --help|h		Display this help text.

The main output of the program is up to four files in their original format containing
any additional data columns.

END
	print $usagetext;
	exit;
}
