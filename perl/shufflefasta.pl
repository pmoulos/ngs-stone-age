#!/usr/bin/perl

# Program to shuffle the order of the sequences in a FASTA file.
#
# Author      : Panagiotis Moulos (pmoulos@eie.gr)
# Created     : 18 - 03 - 2009 (dd - mm - yyyy)
# Last Update : 19 - 03 - 2009 (dd - mm - yyyy)
# Version     : 1.0

use strict;
use File::Basename;

# Make sure output is unbuffered
select(STDOUT);
$|=1;

# Set defaults
my $input = shift @ARGV;
my $output = shift @ARGV;
my $randexist = 0;  # Does Math::Random exist?

die "\nSyntax: (perl) shufflefasta.pl INPUT_FILE OUTPUT_FILE\n\n" if (!$input || !$output);

# Check if the module Math::Random exists, required for random number generation
eval
{
	require Math::Random;
};
if ($@)
{
	my $warn = "Module Math::Random would speed up the procedure. Proceeding with conventional\n".
    		   "perl rand() function (may take up a lot of memory if you have large files...\n";
	print "\n$warn\n";
    $randexist = 0;
}
else
{
	eval "use Math::Random qw(random_permutation random_permuted_index)";
    $randexist = 1;
}

if (!$output)
{
	my ($base,$dir,$ext) = fileparse($input,'\..*?');
	$output = $dir.$base."_SHUFFLED".$ext;
}

print "\nShuffling...";
my ($currid,%file);
open(INPUT,$input)  or die "\nThe file $input does not exist!\n";;
while (my $line = <INPUT>)
{
	$line =~ s/\r|\n$//g;
	if ($line =~ /^>/)
	{
		$currid = $line;
	}
	elsif ($currid)
	{
		push(@{$file{$currid}},$line);
	}
}	

my @shuffkeys = keys(%file);
if ($randexist)
{
	@shuffkeys = random_permutation(@shuffkeys);
	open(OUTPUT,">$output");
	foreach my $ce (@shuffkeys)
	{
		print OUTPUT "$ce\n";
		print OUTPUT join ("\n",@{$file{$ce}}),"\n";
	}
}
else
{
	open(OUTPUT,">$output");
	my $choice = splice(@shuffkeys,rand(@shuffkeys),1);
	print OUTPUT "$choice\n";
	print OUTPUT join ("\n",@{$file{$choice}}),"\n";
}
print "\nDone!...\n";
