#!/usr/bin/perl

# Program to shuffle the order of the sequences in a FASTA file.
#
# Author      : Panagiotis Moulos (pmoulos@eie.gr)
# Created     : 15 - 04 - 2009 (dd - mm - yyyy)
# Last Update : 22 - 04 - 2009 (dd - mm - yyyy)
# Version     : 1.0

use strict;
use Getopt::Long;
use POSIX qw(floor
 			 ceil);

# Make sure output is unbuffered
select(STDOUT);
$|=1;

# Set defaults
our $scriptname = "tab2fasta.pl";
our @input;       # Input file(s)
our @output;      # Output file(s)
our @idseq;       # Columns that contain ID to become FASTA ID anda sequence
our $width = 50;  # Fasta line width
our $header = 0;  # Tabular file has header?
our $waitbar = 0; # Waitbar?
our $silent = 0;  # Do not display any messages
our $help = 0;    # Help?

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
	my @conts;
	disp("Converting file $input[$i]");
	open(INPUT,$input[$i]) or die "\nThe file $input[$i] does not exist!\n";
	&waitbarInit if ($waitbar);
	$output[$i] = &createOutputFile($input[$i]) if (!$output[$i]);
	open(OUTPUT,">$output[$i]");
	my $line = <INPUT> if ($header);
	while (my $line = <INPUT>)
	{
		&waitbarUpdate($.,$nlines[$i]) if ($waitbar);
		$line =~ s/\r|\n$//g;
		@conts = split(/\t/,$line);
		print OUTPUT ">$conts[$idseq[0]]\n";
		my $seq = $conts[$idseq[1]];
		while ($seq)
		{
			my $olin = substr($seq,0,$width,"");
			print OUTPUT "$olin\n";
		}
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
    		   "output|o=s{,}" => \@output,
    		   "idseq|q=i{,}" => \@idseq,
    		   "width|w=i" => \$width,
    		   "header|d" => \$header,
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
    if (@idseq)
    {
    	if (!$idseq[1])
    	{
    		disp("ID or sequence columns not given... Using defaults (1,2)...");
    		@idseq = (0,1);
		}
		else
		{
    		$idseq[0]--;
    		$idseq[1]--;
		}	
    }
    else
    {
    	disp("ID and sequence columns not given... Using defaults (1,2)...");
		@idseq = (0,1);
	}
}

sub createOutputFile
{
	use File::Basename;
	my $in = shift @_;
	my ($base,$dir) = fileparse($in,'\..*?');
	return($dir.$base."_FASTA.fa");
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
A perl program to convert a tabular file with several columns among which
one is a line unique ID and one has a sequence to a FASTA file. The unique
ID will be the FASTA ID line for each sequence

Author : Panagiotis Moulos (pmoulos\@eie.gr)

Main usage
$scriptname --input file(s) [OPTIONS]

--- Required ---
  --input|i  file(s)  Input tab delimited text file(s)
--- Optional ---
  --output|o		Output filename(s). If multiple input is given, names
  				given should be as many as the input files. 
  --idseq|q		A vector with two values, the first one should be the
  			column in the input file containing unique IDs and the
  			other should be the column containing the sequence. Defaults
  			to (1,2). Example: --idseq 1 4
  --width|w		An integer denoting the desired width of sequence lines
  			in FASTA output. Deaults to 50.
  --header|d		Use this option if you have a header line in
			your input files. It will not be written in your output
			FASTA files. Defaults to no header.
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


#my @seq = split(undef,$conts[$idseq[1]]);
#my $tim = int(@seq/$width);
#for ($j=0; $j<$tim; $j++)
#{
	#my $olin = join("",splice(@seq,0,$width));
	#print OUTPUT "$olin\n";
#}
#if ($seq[0])
#{
	#my $olin = join("",@seq);
	#print OUTPUT "$olin\n";
#}
