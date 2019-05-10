#!/usr/bin/perl

# Program to shuffle the order of the sequences in a FASTA file.
#
# Author      : Panagiotis Moulos (pmoulos@eie.gr)
# Created     : 22 - 04 - 2009 (dd - mm - yyyy)
# Last Update : XX - XX - XXXX (dd - mm - yyyy)
# Version     : 1.0

use strict;
use Getopt::Long;
use POSIX qw(floor
 			 ceil);

# Make sure output is unbuffered
select(STDOUT);
$|=1;

# Set defaults
our $scriptname = "fasta2tab.pl";
our @input;         # Input file(s)
our @output;        # Output file(s)
our $addheader = 0; # Add a header line to the output?
our $waitbar = 0;   # Waitbar?
our $silent = 0;    # Do not display any messages
our $help = 0;      # Help?

# Check inputs
&checkInputs;

my ($i,$j);
my @nlines;
if ($waitbar)
{
	for ($i=0; $i<@input; $i++)
	{
		$nlines[$i] = &countLines($input[$i]);
	}
}

for ($i=0; $i<@input; $i++)
{
	# Quick and dirty checking about FASTA file
	open(FCHK,$input[$i]) or die "\nThe file $input[$i] does not exist!\n";
	if (<FCHK> !~ /^>/)
	{
		disp("File $input[$i] does not appear to be a FASTA file. Exiting...");
		exit;
	}
	close(FCHK);
	
	my $currid;
	disp("Converting file $input[$i]");
	open(INPUT,$input[$i]) or die "\nThe file $input[$i] does not exist!\n";
	&waitbarInit if ($waitbar);
	$output[$i] = &createOutputFile($input[$i]) if (!$output[$i]);
	open(OUTPUT,">$output[$i]");
	print OUTPUT "SeqID\tSequence\n" if ($addheader);
	while (my $line = <INPUT>)
	{
		&waitbarUpdate($.,$nlines[$i]) if ($waitbar);
		$line =~ s/\r|\n$//g;
		if ($line =~ /^>/)
		{
			print OUTPUT "\n" if ($currid);
			$currid = $line;
			$currid =~ s/^>//g;
			print OUTPUT "$currid\t";
		}
		elsif ($currid)
		{
			print OUTPUT "$line";
		}
	}
	print OUTPUT "\n"; # It can confuse tools like wc if we don't do that...	
	close(INPUT);
	close(OUTPUT);
}
disp("Done!\n");

# Process inputs
sub checkInputs
{
    my $stop;
    GetOptions("input|i=s{,}" => \@input,
    		   "output|o=s{,}" => \@output,
    		   "addheader|d" => \$addheader,
    		   "waitbar|b" => \$waitbar,
    		   "silent|s" => \$silent,
    		   "help|h" => \$help);
    # Check if the required arguments are set
    if ($help)
    {
    	&programUsage;
    	exit;
    }
    $stop .= "--- Please specify input file(s) ---\n" if (!@input);
    if ($stop)
    {
            print "\n$stop\n";
            print "Type perl $scriptname --help for help in usage.\n\n";
            exit;
    }
    if (@output)
    {
    	my $inlen = @input;
    	my $outlen = @output;
    	if ($inlen != $outlen)
    	{
    		disp("Output filenames not as many as input files... Output filenames will be auto-generated");
			@output = ();
		}
	}
}

sub createOutputFile
{
	use File::Basename;
	my $in = shift @_;
	my ($base,$dir) = fileparse($in,'\..*?');
	return($dir.$base."_Tab.txt");
}

sub waitbarInit
{
	my $initlen = $_[0];
	$initlen = 50 if ($initlen eq "");
	my $printlen = ' 'x$initlen;
	print "\nProgress\n";
	print "|$printlen|\n";
	print("|");
}

sub waitbarUpdate
{
	my $curr = $_[0];
	my $tot = $_[1];
	my $waitbarlen = $_[2];
	$waitbarlen = 50 if ($waitbarlen eq "");
	my $step;
	if ($tot > $waitbarlen)
	{
		$step = ceil($tot/$waitbarlen);
		print "#" if ($curr%$step == 0);
	}
	else
	{
		$step = floor($waitbarlen/$tot);
		print "#" x $step;
	}
	if ($curr == $tot)
	{
		my $rem = $waitbarlen - floor($tot/$step);
		($rem != 0) ? (print "#" x $rem."|\n") : print "|\n";
	}
}

sub countLines
{
	open(IN,$_[0]) or die "\nThe file $_[0] does not exist!\n\n";
	my $totlines=0;
	$totlines += tr/\n/\n/ while sysread(IN,$_,2**16);
	close(IN);
	return $totlines;
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
A perl program to convert a FASTA file to a tabular sequence file. The FASTA ID
will be the unique sequence ID column for each sequence

Author : Panagiotis Moulos (pmoulos\@eie.gr)

Main usage
$scriptname --input file(s) [OPTIONS]

--- Required ---
  --input|i  file(s)  Input tab delimited text file(s)
--- Optional ---
  --output|o		Output filename(s). If multiple input is given,
  			names given should be as many as the input files.
  --addheader|d		Use this option if you want a header line with
			column names in the output file. Defaults to no header.
  --waitbar|w		Use this option if you wish to display a simple
			progress bar while selecting and printing the final files.			
  --silent|s		Use this option if you want to turn informative 
  			messages off.
  --help|h		Display this help text.

The main output of the program is up to four files in their original format containing
any additional data columns.

END
	print $usagetext;
	exit;
}
