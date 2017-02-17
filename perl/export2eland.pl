#!/usr/bin/perl

# Program to convert extended eland _export.txt fike to standard eland
#
# Author      : Panagiotis Moulos (pmoulos@eie.gr)
# Created     : 10 - 11 - 2009 (dd - mm - yyyy)
# Last Update : XX - XX - XXXX (dd - mm - yyyy)
# Version     : 1.0

use strict;
use Getopt::Long;
use File::Basename;
use Switch;
use POSIX qw(floor
 			 ceil);

# Make sure output is unbuffered
select(STDOUT);
$|=1;

# Set defaults
our $scriptname = "export2eland.pl";
our @input;        # Input file(s)
our $mismatch = 1; # Allowed mismatches
our $filQC = 0;    # Filter QC
our $filNM = 0;    # Filter NM
our $filO = 0; 	   # Filter other non-matcing codes
our $waitbar = 0;  # Use waitbar
our $silent = 0;   # Do not display any messages
our $help = 0;     # Help?

# Check inputs
&checkInputs;

my ($i,$j);
my @nlines;
if ($waitbar)
{
	for ($i=0; $i<@input; $i++)
	{
		disp("Counting lines for file $input[$i]");
		$nlines[$i] = &countLines($input[$i]);
	}
}

for ($i=0; $i<@input; $i++)
{
	my (@output,@conts,@outline);
	my ($id,$seq,$match,$exmatch,$e1match,$e2match,$chr,$start,$strand,$N,$sub1,$sub2);
	my $miscount;
	disp("Converting file $input[$i]");
	open(INPUT,$input[$i]) or die "\nThe file $input[$i] does not exist!\n";
	&waitbarInit if ($waitbar);
	$output[$i] = &createOutputFile($input[$i]);
	open(OUTPUT,">$output[$i]");
	
	while (my $line = <INPUT>)
	{
		&waitbarUpdate($.,$nlines[$i]) if ($waitbar);
		$line =~ s/\r|\n$//g;
		@conts = split(/\t/,$line);
		next if ($conts[10] eq "QC" && $conts[21] eq "N");
		next if (!$filQC && $conts[10] eq "QC");
		next if (!$filNM && $conts[10] eq "NM");
		next if (!$filO && $conts[10] !~ 'chr[0-9XY]*');
		$miscount = ($conts[14] =~ tr/[ACGTN]//);
		next if ($miscount > $mismatch);
		$match = 0;
		$id = join(":",@conts[0..5]);
		$seq = $conts[8];
		$match = "NM" if ($filO && $conts[10] !~ 'chr[0-9XY]*');
		$match = "QC" if ($filQC && $conts[10] eq "QC");
		$match = "NM" if ($filNM && $conts[10] eq "NM");
		
		if (!$match)
		{
			switch ($miscount)
			{
				case 0 
				{ 
					$match = "U0";
					($exmatch,$e1match,$e2match) = (1,0,0);
				}
				case 1 
				{
					$match = "U1";
					($exmatch,$e1match,$e2match) = (0,1,0);
				}
				case 2
				{
					$match = "U2";
					($exmatch,$e1match,$e2match) = (0,0,1);
				}
				else 
				{
					$match = "RM";
					($exmatch,$e1match,$e2match) = (0,0,0);
				}
			}
		}
		
		$chr = $conts[10];
		$start = $conts[12];
		$strand = $conts[13];
		($N,$sub1,$sub2) = ("NA","NA","NA");
		
		@outline = ($id,$seq,$match,$exmatch,$e1match,$e2match,$chr,$start,$strand,$N,$sub1,$sub2);
		print OUTPUT join("\t",@outline),"\n";
	}
	close(INPUT);
	close(OUTPUT);
}
disp("Done!\n");

# Process inputs
sub checkInputs
{
    my $stop;
    GetOptions("input|i=s{,}" => \@input,
    		   "mismatch|m=i" => \$mismatch,
    		   "nofilterQC|q" => \$filQC,
    		   "nofilterNM|n" => \$filNM,
    		   "nofilterO|o" => \$filO,
    		   "waitbar|w" => \$waitbar,
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
}

sub createOutputFile
{
	use File::Basename;
	my $in = shift @_;
	my ($base,$dir,$ext) = fileparse($in,'\..*?');
	($base =~ /export/) ? ($base =~ s/export/eland_result/) : ($base = $base."_eland_result");
	return($dir.$base.$ext);
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
A perl program to convert a _export.txt GERALD output to standard non-paired
end read ELAND output with some assumptions.

Author : Panagiotis Moulos (pmoulos\@eie.gr)

Main usage
$scriptname --input file(s) [OPTIONS]

--- Required ---
  --input|i  file(s)  Input tab delimited text file(s)
--- Optional ---
  --mismatch|m		The number of allowed mismatches.
  --nofilterQC|q		Include quality filtered excluded reads.
  --nofilterNM|n		Include non-matched excluded reads.
  --nofilterO|o		Include excluded reads for other reasons.
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


#if (!($filQC || $filNM))
#{
	#$match = "NM" if (!$filO && $conts[10] !~ 'chr[0-9XY]*');
	#$match = $conts[10] if ($conts[10] eq "NM" || $conts[10] eq "QC");
#}
