#!/usr/bin/perl

# Description later...
#
# Author      : Panagiotis Moulos (pmoulos@eie.gr)
# Created     : 08 - 06 - 2009 (dd - mm - yyyy)
# Last Update : 24 - 06 - 2009 (dd - mm - yyyy)
# Version     : 1.0
#

use strict;
use File::Basename;
use Getopt::Long;
 			 
# Make sure output is unbuffered
select(STDOUT);
$|=1;

# Set defaults
our $scriptname = "peakmotifocc.pl";
our $peaks;		 # The peak-region-whatever file (BED format)
our @scanres;	 # Scan results (GFF or BED) for that specific peak file
our $idcol;		 # peakid column (default 4) 
our $type;	     # GFF or BED scan results, empty for autodetection
our $mode;	     # Search mode: if peak IDs and gff use them, else intersectbed.pl
our $npass = 3;  # Number of passes for dynamic binary search
our $header = 0; # Peak files have headers?
our $silent = 0; # Display verbose messages
our $help = 0;	 # Help?

# Check inputs
&checkInputs;

if (!$type)
{
	disp("Trying to auto-detect scan result type (gff or bed) and adjust mode option...");
	open(DETECT,$scanres[0]) or die "\nThe file $scanres[0] does not exist!\n";
	my @check = split(/\t/,<DETECT>);
    my $len = @check;
	if (($len == 3 || $len == 6) && $check[0] =~ /chr(\d+|X|Y)/ && $check[1] =~ /\d+/ && $check[2] =~ /\d+/)
	{
		disp("Scan result files appear to be BED files. Proceeding...");
		$type = "bed";
		$mode = "search";
	}
	elsif ($len >= 9 && $check[3] =~ /\d+/ && $check[4] =~ /\d+/)
	{
		disp("Scan result files appear to be GFF files. Proceeding...");
		$type = "gff";
		$mode = "keys";
	}
	else { die "\nCannot recognize scan result input files!\n"; }
}

my $windows = 0; # Windows or Linux?
my ($i,$j); # General indices
my $opeaks = $peaks; # Keep original filenames
my @oscanres = @scanres; # Keep original filenames

if ($type eq "bed" && $mode eq "search") # We have to sort the files
{
	# Simulate try and catch
	if ($^O !~ /MSWin/) # Case of linux, easy sorting
	{
		# Sort peak files
		if ($header)
		{
			disp("Sorting file $peaks...");
			my $wid = "$$";
			my $sortcmd = "awk 'NR==1; NR > 1 {print \$0 | \" sort -k1,1 -k2g,2\"}' $peaks > /tmp/tempreg.in$wid";
			`$sortcmd`;
			$peaks = "/tmp/tempreg.in$wid";
		}
		else
		{	
			disp("Sorting file $peaks...");
			`sort -k1,1 -k2g,2 $peaks > /tmp/tempreg.in$$ `;
			$peaks = "/tmp/tempreg.in$$";
		}
		# Sort scan result bed files
		for ($i=0; $i<@scanres;$i++)
		{
			disp("Sorting file $scanres[$i]...");
			`sort -k1,1 -k2g,2 $scanres[$i] > /tmp/tempbed$i.in$$ `;
			$scanres[$i] = "/tmp/tempbed$i.in$$";
		}
	}
	else # We are in Windows... package required
	{
		eval
		{
			require File::Sort;
		};
		if ($@)
		{
			my $killer = "Module File::Sort is required to continue with file sorting. If you\n".
						 "have ActiveState Perl installed, use the Package Manager to get the\n".
						 "module. If you don't know how to install the package, sort the files\n".
						 "using another tool like Excel or contact your system administrator.";
			die "\n$killer\n\n";
		}
		else # Else sort the bed files according to chromosome and start position
		{
			if ($header) # Cannot find a solution for direct sorting in Windows... sorry :-(
			{
				my $dmsg = "Module File::Sort can't sort a file with a header line without possible\n".
						   "messing up data. Please sort files outside $scriptname first (e.g. using\n".
						   "Excel or something similar.";
				die "\n$dmsg\n\n";
			}
			eval "use File::Sort qw(sort_file)"; # Like this or interpreter complains
			# Sort peak/region/whatever file
			disp("Sorting region file $peaks...");
			sort_file(
			{
				I => $peaks,
				o => "tempreg.tmp",
				k => ['1,1','2n,2'],
				t => "\t"
			});
			$peaks = "tempreg.tmp";
			# Sort scan result files
			for ($i=0; $i<@scanres;$i++)
			{
				disp("Sorting region file $scanres[$i]...");
				sort_file(
				{
					I => $scanres[$i],
					o => "tempbed$i.tmp",
					k => ['1,1','2n,2'],
					t => "\t"
				});
				$scanres[$i] = "tempbed$i.tmp";
			}
		}
		$windows = 1;
	}
}

# Start doing job, easy if search mode is keys
if ($mode eq "keys")
{	
	# A help variable
	my @temp;
	foreach my $s (@scanres)
	{
		push(@temp,"-");
	}
	# Suck in peak regions
	my (%peakrecs,%hitrecs,@conts);
	disp("Reading primary region file $opeaks...");
	open(INPEAK,"$peaks");
	my $peakhead = <INPEAK> if ($header);
	$peakhead =~ s/\r|\n$//g;
	while (my $line = <INPEAK>)
	{
		$line =~ s/\r|\n$//g;
		@conts = split(/\t/,$line);
		$peakrecs{$conts[$idcol]} = $line;
		@{$hitrecs{$conts[$idcol]}} = @temp;
	}
	close(INPEAK);
	for ($j=0; $j<@scanres; $j++)
	{
		disp("Processing secondary region file $oscanres[$j]...");
		open(INSCAN,"$scanres[$j]");
		if ($type eq "gff")
		{
			while (my $line = <INSCAN>)
			{
				$line =~ s/\r|\n$//g;
				my @scans = split(/\t/,$line);
				${$hitrecs{$scans[0]}}[$j] = "+" if ($peakrecs{$scans[0]});
			}
		}
		elsif ($type eq "bed")
		{
			# Who knows...
			my @chk = split(/\t/,<INSCAN>);
			die "\nScan result BED file $scanres[$j] not in proper format for \"keys\" search mode"
				if (!$chk[3]);
			sysseek(INSCAN,0,0);
			while(my $line = <INSCAN>)
			{
				$line =~ s/\r|\n$//g;
				my @scans = split(/\t/,$line);
				${$hitrecs{$scans[3]}}[$j] = "+" if ($peakrecs{$scans[3]});
			}
		}
		close(INSCAN);
	}
	my $output = &createOutputFile($opeaks);
	disp("Writing output in $output...");
	open(OUTPUT,">$output");
	# Deal with header... output must have one because of hits in each file
	if ($header)
	{
		print OUTPUT "$peakhead\t",join("\t",@oscanres),"\n";
	}
	else
	{
		my $outhead;
		for (my $k=0; $k<@conts; $k++)
		{
			$outhead .= "\Column$k\t";
		}
		print OUTPUT "$outhead",join("\t",@oscanres),"\n";
	}
	foreach my $ke (keys(%peakrecs))
	{
		print OUTPUT "$peakrecs{$ke}\t",join("\t",@{$hitrecs{$ke}}),"\n";
	}
	close(OUTPUT);
}
elsif ($mode eq "search") # Binary search on sorted files
{	
	# A help variable
	my @temp;
	foreach my $s (@scanres)
	{
		push(@temp,"-");
	}
	# Suck in peak regions
	my (%peakrecs,%hitrecs,@conts);
	disp("Reading primary region file $opeaks...");
	open(INPEAK,"$peaks");
	my $peakhead = <INPEAK> if ($header);
	$peakhead =~ s/\r|\n$//g;
	while(my $line = <INPEAK>)
	{
		next if ($line =~/^chrM/);
		next if ($line =~/rand/);
		$line =~ s/\r|\n$//g;
		my ($chr,@rest) = split(/\t/,$line);
		push(@{$peakrecs{$chr}},join("\t",@rest));
		push(@{$hitrecs{$chr}},join("\t",@temp));
	}
	close(INPEAK);
	for ($j=0; $j<@scanres; $j++)
	{
		disp("Processing secondary region file $oscanres[$j]...");
		my ($line,$cchr,$n,$bsr,$ci,@crest,@currvals);
		open(INSCAN,"$scanres[$j]");
		while(my $line = <INSCAN>)
		{
			next if ($line =~/^chrM/);
			next if ($line =~/rand/);
			$line =~ s/\r|\n$//g;
			($cchr,@crest) = split(/\t/,$line);
			
			@currvals = @{$peakrecs{$cchr}} if ($peakrecs{$cchr});
			$n = 0;
			while ($n < $npass) 
			{
				($ci,$bsr) = &binSearchAny($crest[0],$crest[1],@currvals);
				#print "\n",join("\n",@{$hitrecs{$cchr}}),"\n";
				if ($bsr) # Found in overlap
				{
					my @sf = split(/\t/,${$hitrecs{$cchr}}[$ci]);
					$sf[$j] = "+" ;
					${$hitrecs{$cchr}}[$ci] = join("\t",@sf);
				}
				if ($ci)
				{
					splice(@currvals,$ci,1);
					$n++;
				} else { last; }
			}
		}
		close(INSCAN);
	}
	my $output = &createOutputFile($opeaks);
	disp("Writing output in $output...");
	open(OUTPUT,">$output");
	# Deal with header... output must have one because of hits in each file
	if ($header)
	{
		print OUTPUT "$peakhead\t",join("\t",@oscanres),"\n";
	}
	else
	{
		my $outhead;
		for (my $k=0; $k<@conts; $k++)
		{
			$outhead .= "\Column$k\t";
		}
		print OUTPUT "$outhead",join("\t",@oscanres),"\n";
	}
	foreach my $outchr (sort(keys(%peakrecs)))
	{
		my @a = @{$hitrecs{$outchr}};
		for (my $ind=0; $ind<@a; $ind++)
		{
			print OUTPUT "$outchr\t${$peakrecs{$outchr}}[$ind]\t${$hitrecs{$outchr}}[$ind]\n";
		}
	}
	close(OUTPUT);
}

# Remove garbage
(!$windows) ? (`rm /tmp/temp*.in$$ `) : (`del /f /q temp*.tmp `) if ($type eq "bed" && $mode eq "search");
disp("Finished!\n\n");
	

# Process inputs
sub checkInputs
{
    my $stop;
    GetOptions("input|i=s" => \$peaks,
    		   "motif|r=s{,}" => \@scanres,
    		   "idcol|c=i" => \$idcol,
    		   "type|t=s" => \$type,
    		   "mode|m=s" => \$mode,
    		   "pass|n=i" => \$npass,
    		   "header|d" => \$header,
    		   "silent|s" => \$silent,
    		   "help|h" => \$help);
    # Check if the required arguments are set
    if ($help)
    {
    	&programUsage;
    	exit;
    }
    $stop .= "--- Please specify input region file ---\n" if (!$peaks);
    $stop .= "--- Please specify scan result file(s) ---\n" if (!@scanres);
    $stop .= "--- --mode option should be one of \"keys\" or \"search\" (case sensitive) ---\n" 
    	if ($mode && $mode ne "keys" && $mode ne "search");
    if ($stop)
    {
            print "\n$stop\n";
            print "Type perl $scriptname --help for help in usage.\n\n";
            exit;
    }
    if ($npass < 0)
	{
		disp("--pass option should be >0... Using default (3)");
		$npass = 3;
	}
    if (!$mode)
	{
		disp("Search mode not given. Using binary search mode...");
		$mode = "search";
	}
	# Check range id columns
	if (!$idcol)
	{
		disp("ID column not given. Search mode will be \"search\"... ");
		$idcol = 3;
		$mode = "search";
	} else { $idcol--; }
	disp("Scan result type not given. Will try to autodetect...") if (!$type);
}

sub binSearchAny
{
	my $start = shift @_;
	my $end = shift @_;
	my @areas = @_;
	my ($ind,$currstart,$currend);
	my ($l,$u) = (0,$#areas);
	$u = 1 if ($#areas == 0); # Kavourmadies...
	while ($l <= $u)
	{
		$ind = int(($l + $u)/2);
		($currstart,$currend) = split(/\t/,$areas[$ind]);
		if (($currstart >= $start && $currend <= $end) ||
		    ($currstart <= $start && $currend >= $end) ||
			($currstart < $start && $currend < $end && $start < $currend) ||
			($currstart > $start && $currend > $end && $end > $currstart))
		{
			return($ind,$areas[$ind]);
		}
		else
		{
			$u = $ind - 1 if ($end <= $currstart);
            $l = $ind + 1 if ($start >= $currend);
		}
	}
	return(0);
}

sub createOutputFile
{
	my $in = shift @_;
	my ($base,$dir,$ext) = fileparse($in,'\..*?');
	return($dir.$base."_HITS".$ext);
}

sub disp
{
	print "\n@_" if (!$silent);
}

sub programUsage
{
	my $usagetext = << "END";

$scriptname
A perl program which given a region file in BED format (can also contain unique region
IDs) and one or more BED or GFF files (gff as output from p53scan program) will create
an output file containing the regions (e.g. peaks) of the initial input file and as many
columns appended as the other series of files. Each column is named with the name of the
second series of input files and contains a "+" if the area in the bed/gff file falls
within the region under interrogation or "-" otherwise. For example, if the input file
contains a set of peaks and the user wants to see in which peaks the set of motifs used
for motif scanning is found and where is not found. The main input consists of one 
peak/whatever region BED file and is unique and also of a set of BED/GFF files which are
derived from a motif scan algorithm for several motifs. In the case of GFF the program
has been tested with GFF output from p53scan. Also, in the case of GFF the first column
should contain IDs same as (derived from same dataset, i.e. scanning of the region file
for several motifs) the region IDs in the input region file. Other options include the
type of input motif scan results, the search mode to produce the results etc. Attention
should be paid to the fact that motif scan results should be derived from the same file
(or partial) as the input file, or else no hits might be found!

Author : Panagiotis Moulos (pmoulos\@eie.gr)

Main usage
$scriptname --input file --motif file(s) [OPTIONS]

--- Required ---
  --input|i  file  A BED file with peak/whatever regions. Can contain unique ID.
  --motif|r  file  Motif occurences preferably from scanthemall.pl in bed or gff format
                   (can be from p53scan or pwmscan also)
--- Optional ---
  --idcol|c		Use this option to provide the script the column number
			in the peak/region file that contains unique region ID. It  
			should consist of 1 number (e.g. --idcol 1). If not given,
			non-existence is assumed and the program will run on "search"
			mode.
  --type|t		Use this option to tell the script what kind of motif scan
  			result output file you have. It should be one of "bed" or "gff".
  			If it is a gff from p53scan for example and peak/region IDs have
  			been used to identify sequence hits, the script can run in "keys"
  			mode given the fact that these IDs are the same as in the region
  			file. If type is "bed", --idcol is given and the 4th column of
  			the BED file that contains the scan hits have IDs corresponding
  			to the region file, "keys" mode can be used, else you have to use
  			"search" mode. If type is not given, program will try to detect.
  			If auto-detection fails, program will exit.
  --mode|m		Use this option to provide the script the search mode that
  			will be used for assigning hits to specific regions. It can be
  			"keys" if both files have the same region IDs or "search" if
  			there are no keys, or if the motif scan results are BED files
  			without bedline IDs (4th column). "keys" mode works faster as
  			there is direct mapping. "search" uses a binary search algorithm
  			and the input files are sorted first. Option should be given
  			for example as --mode search. Default is "search".
  --pass|n		Use this option to supply the program the number of times
			that the dynamic binary search algorithm will search regions
			from the motif scan results for each region of the input file.
			One pass returnsat maximum one hit because if the algorithm 
			finds a hit from input, it will stop. The number of passes 
			determines how many times the algorithm will search for 
			overlapping regions. Use a larger number of passes for when
			regions from motif scan results are likely to overlap many 
			regions from input (e.g. when multiple motif occurences are
			likely to exist in a region. It defaults to 3.
  --header|d		Use this option if you have a header line in your
			peak/region input file. It will also be written in your 
			output files. Defaults to no header.
  --silent|s		Use this option if you want to turn informative 
  			messages off.
  --help|h		Display this help text.

The main output of the program is the region file as is but with several columns (one
for each motif) appended to it. Output file name will be the same but with "_HITS" 
appended to it.

END
	print $usagetext;
	exit;
}

sub printH
{
	my %h = @_;
	foreach my $k (keys(%h))
	{
		print "\n$k\t",join("\t",@{$h{$k}});
	}
}
