#!/usr/bin/perl

use strict;
use Getopt::Long;
use File::Basename;
use File::Temp;
use File::Spec;
             
# Make sure output is unbuffered
select(STDOUT);
$|=1;

# Set defaults
our $scriptname = "hyperassignpeaks.pl";
our @input; # The peak file(s) (BED format)
our $region; # Significant regions to be assigned to peaks (BED format)
our $background; # Background regions to be assigned to peaks (BED format)
our @span; # Upstream and downstream  (default +/-100k from TSS)
our $where; # Promoter or coding?
our @sbcols; # Columns containing unique sig/back region ID and strand (default (4,5))
our @pcols; # Columns containing unique peak region ID and mode (default (4,5))
our @expcols; # Columns containing gene expression
our @out; # Output filetype
our $test; # Default hypergeometric test
our $pval; # Hypergeometric p-value cutoff (default 0.05)
our $redun; # Redundancy level of assignment
our $nomode; # If no peak mode given, center or boundary?
our $log; # Keep log?
our $silent; # Display verbose messages
our $help; # Help?

# Check inputs
&check_inputs;

# Bavard
my $date = &now;
&disp("$date - Started...\n");
if ($where eq "promoter")
{
	&disp("Genomic span from TSS : ($span[0],$span[1])");
}
else
{
	&disp("Genomic span from TSS : ($span[0], TES)");
}
&disp("Assignment redundancy level : $redun");
&disp("No mode presence action : $nomode");
&disp("Over-representation test : $test");
&disp("p-value threshold : $pval") if ($test ne "none");
&disp("Chosen output(s) : ",join("\, ",@out),"\n");

# General indices and some info about gff output
my ($i,$j,$k);
my $gffreq = my $pdataout = 0;
foreach my $o (@out)
{
	$gffreq = 1 if ($o =~ /gff/);
	$pdataout = 1 if ($o =~ /peakdata|bed/);
}

# Some intialization
my (%sigID,%sigStart,%sigEnd);
my (%backID,%backStart,%backEnd);
my (@all,$chr,$start,$end,$id,$strand,$expr,$line);
my $lensig = my $lenback = 0;
my @lenpeak;

# Variable for matrix generation
my (%hasPeak,%hasExpression);
tie %hasPeak, "Tie::IxHash::Easy" if (&smatch("matrix-number",@out)
	|| &smatch("matrix-peaks",@out) || &smatch("matrix-presence",@out));

# Suck in significant region file
open (REG,$region) or die "\nThe file $region does not exist!\n";;
&disp("Reading region file $region...");
$line = <REG>;
my $reghead = &decide_header($line);
if (!$reghead)
{
	seek(REG,0,0);
}
else
{
	if (@expcols)
	{
		my @ht = split("\t",$reghead);
		$hasExpression{"header"} = join("\t",@ht[@expcols]);
	}
}
while ($line = <REG>)
{
	next if ($line =~/^chrM/);
	next if ($line =~/rand/);
	$line =~ s/\r|\n$//g;
	@all = split(/\t/,$line);
	$chr = $all[0];
	$start = $all[1];
	$end = $all[2];
	$id = $all[$sbcols[0]];
	$strand = $all[$sbcols[1]];
	$expr = join("\t",@all[@expcols]) if (@expcols);
	#if ($strand == 1 || $strand eq "+" || $strand eq "F")
	if ($strand =~ m/^[+1F]$/)
	{
		push(@{$sigStart{$chr}},$start);
		push(@{$sigEnd{$chr}},$end);
	}
	#elsif ($strand == -1 || $strand eq "-" || $strand eq "R")
	elsif ($strand =~ m/^(-1)|[-R]$/)
	{
		push(@{$sigStart{$chr}},$end);
		push(@{$sigEnd{$chr}},$start);
	}
	else # Some self-defense
	{
		&disp("Improper strand format... Skipping line $. from $region");
		next;
	}
	push(@{$sigID{$chr}},$id);
	$lensig++;

	# Initiate the hash to keep record of peaks for a peak matrix file generation
	if (&smatch("matrix-number",@out) || &smatch("matrix-peaks",@out)
		|| &smatch("matrix-presence",@out))
	{
		for ($i=0; $i<@input; $i++)
		{
			$hasPeak{$id}{basename($input[$i])} = ();
		}
		$hasExpression{"data"}{$id} = $expr;
	}
}
close(REG);

# Suck in background region file if test is to be performed
if ($test ne "none")
{
	open (BACK,$background) or die "\nThe file $background does not exist!\n";
	&disp("Reading background file $background...");
	$line = <BACK>;
	my $backhead = &decide_header($line);
	seek(BACK,0,0) if (!$backhead);
	while ($line = <BACK>)
	{
		next if ($line =~/^chrM/);
		next if ($line =~/rand/);
		$line =~ s/\r|\n$//g;
		@all = split(/\t/,$line);
		$chr = $all[0];
		$start = $all[1];
		$end = $all[2];
		$id = $all[$sbcols[0]];
		$strand = $all[$sbcols[1]];
		if ($strand == 1 || $strand eq "+" || $strand eq "F")
		{
			push(@{$backStart{$chr}},$start);
			push(@{$backEnd{$chr}},$end);
		}
		elsif ($strand == -1 || $strand eq "-" || $strand eq "R")
		{
			push(@{$backStart{$chr}},$end);
			push(@{$backEnd{$chr}},$start);
		}
		else # Some self-defense
		{
			disp("Improper strand format... Skipping line $. from $background");
			next;
		}
		push(@{$backID{$chr}},$id);
		$lenback++;
	}
	close(BACK);
}

# Some self-defense... Check uniqueness of gene IDs...
die "\nThe region IDs in regions file $region are not unique! Exiting...\n" 
	if (!&check_unique(\%sigID));
die "\nThe region IDs in background file $background are not unique! Exiting...\n" 
	if (!&check_unique(\%backID) && $test ne "none");

# Read and process peak files
for ($i=0; $i<@input; $i++)
{
	&disp("\nReading file $input[$i]...");
	
	my (%peakID,%peakMode,%peakScore,%countSig,%countBack,%allPeakData,
		%countGenesSig,%countGenesBack);
	my ($pstart,$pend,%peakStarts,%peakEnds) if ($gffreq);
	my (@pall,$pchr,$pmode,$pid,$phead,$pscore);
	my (@starts,@ends,@modes);
	my ($limit1,$limit2);
	my (@peakchrs,@peakids,@peakscores,@geneids);
	my (@sigallpeakids,@backallpeakids,@sigassgenes,@backassgenes);
	my ($currchr,$currdist,$currpeak,$currgene,$currout,@currgenes);
	my (%genePeaksSig,%peaksGenesSig,%distPeakBased,%distGeneBased,
		%peakIndex,%geneIndex);
	my ($p,$cp,@elems);
	my %finalPeaks;
	
	open(INPUT,$input[$i]) or die "\nThe file $input[$i] does not exist!\n";
	$line = <INPUT>;
	$phead = &decide_header($line);
	seek(INPUT,0,0) if (!$phead);
	while ($line = <INPUT>)
	{
		$line =~ s/\r|\n$//g;
		@pall = split(/\t/,$line);
		$pchr = $pall[0];
		$pstart = $pall[1];
		$pend = $pall[2];
		$pid = $pall[$pcols[0]];
		if ($pcols[1]) {
			$pmode = $pall[$pcols[1]];
		}
		else {
			if ($nomode eq "center") {
				$pmode = $pstart + &round(($pend - $pstart)/2);
			}
			elsif ($nomode eq "boundary" || $nomode eq "closest"
				|| $nomode eq "furthest") {
				$pmode = [$pstart,$pend];
			}
		}
		$pscore = $pall[$pcols[2]] if ($pcols[2]);
		push(@{$peakID{$pchr}},$pid);
		push(@{$peakMode{$pchr}},$pmode);
		push(@{$peakScore{$pchr}},$pscore) if ($pscore);
		if ($gffreq) # Starts and ends required in this case for GFF files
		{
			push(@{$peakStarts{$pchr}},$pstart);
			push(@{$peakEnds{$pchr}},$pend);
		}
		$allPeakData{$pid} = join("\t",@pall) if ($pdataout);
		$lenpeak[$i]++;
	}
	close(INPUT);
	
	&disp("Processing file $input[$i]...");
	die "\nThe region IDs in file $input[$i] are not unique! Exiting...\n" 
		if (!&check_unique(\%peakID));
	
	@peakchrs = keys(%peakID);
	
	# Do stuff with siginificant file
	if ($test ne "none")
	{
		&disp("Associating query regions with subject regions in foreground ".
			"and background region files :");
		&disp("Subject region file : $region");
		&disp("Background region file : $background");
	}
	else
	{
		&disp("Associating query regions with subject regions in ".
			"region file : $region");
	}
	foreach $currchr (@peakchrs)
	{
		&disp("Queries at $currchr...");
		
		@modes = @{$peakMode{$currchr}};
		@peakids = @{$peakID{$currchr}};
		@peakscores = @{$peakScore{$currchr}} if $pcols[2];
		
		if ($sigID{$currchr}) # Could not have subjects at a specific chromosome
		{
			@geneids = @{$sigID{$currchr}};
			@starts = @{$sigStart{$currchr}};
			@ends = @{$sigEnd{$currchr}};
		
			for ($j=0; $j<@starts; $j++)
			{
				for ($k=0; $k<@modes; $k++)
				{
					$currdist = &dist($starts[$j],$ends[$j],$modes[$k]);
					if ($where eq "promoter")
					{
						$limit1 = $span[0];
						$limit2 = $span[1];
					}
					elsif ($where eq "coding")
					{
						$limit1 = $span[0];
						$limit2 = abs($ends[$j]-$starts[$j]+$limit1);
					}
					elsif ($where eq "downtes") {
						$limit1 = $ends[$j];
						$limit2 = $ends[$j] + $span[0];
					}
					if ($currdist > $limit1 && $currdist < $limit2)
					{
						push(@{$genePeaksSig{$currchr}{$geneids[$j]}},
							$peakids[$k]);
						push(@{$peaksGenesSig{$currchr}{$peakids[$k]}},
							$geneids[$j]);
						push(@{$distGeneBased{$currchr}{$geneids[$j]}},
							$currdist);
						push(@{$distPeakBased{$currchr}{$peakids[$k]}},
							$currdist);
						push(@{$peakIndex{$currchr}{$peakids[$k]}},$k);
						push(@{$geneIndex{$currchr}{$geneids[$j]}},$j);
						push(@sigallpeakids,$peakids[$k]);
						push(@sigassgenes,$geneids[$j]);
						(!$pcols[2]) ? 
						(push(@{$hasPeak{$geneids[$j]}{basename($input[$i])}},
							$peakids[$k]."_".$currdist)) :
						(push(@{$hasPeak{$geneids[$j]}{basename($input[$i])}},
							$peakids[$k]."_".$currdist."_".$peakscores[$k]));
					}
				}
			}
			# if non-redundant
			if ($redun eq "genecentric")
			{
				my (%geneDist,%peakInd);
				my (@pGenes,@pDists,@pInds,@gPeaks,@gDists,@gInds,
					@sortedGenesByDist,@remaining,@matrixElems);
				my ($index,$pe,$di,$in,$hs,$hsi,$redPeak,$g);
				my @redPeaks = keys(%{$peaksGenesSig{$currchr}});
				foreach $redPeak (@redPeaks)
				{
					#print "\n----------  $redPeak";
					@pGenes = @{$peaksGenesSig{$currchr}{$redPeak}};
					@pDists = @{$distPeakBased{$currchr}{$redPeak}};
					@pInds = @{$peakIndex{$currchr}{$redPeak}};
					@geneDist{@pGenes} = @pDists;
					@peakInd{@pGenes} = @pInds;
					@sortedGenesByDist = sort { 
						abs($geneDist{$a}) <=> abs($geneDist{$b}) 
					} keys(%geneDist);
					# We need the first and based on this we adjust the rest of 
					# the hashes
					#print "\nF: ",$sortedGenesByDist[0];
					@{$peaksGenesSig{$currchr}{$redPeak}} = 
						($sortedGenesByDist[0]);
					@{$distPeakBased{$currchr}{$redPeak}} = 
						($geneDist{$sortedGenesByDist[0]});
					@{$peakIndex{$currchr}{$redPeak}} = 
						$peakInd{$sortedGenesByDist[0]};
					# Then we need to kick out the rest of the sorted genes from
					# gene based hashes
					@remaining = @sortedGenesByDist[1..$#sortedGenesByDist];
					if (@remaining)
					{
						#print "\nR: ",join(", ",@remaining);
						foreach $g (@remaining)
						{
							#print "\n$g\n";
							@gPeaks = @{$genePeaksSig{$currchr}{$g}};
							#print join(", ",@gPeaks);
							@gDists = @{$distGeneBased{$currchr}{$g}};
							@gInds = @{$geneIndex{$currchr}{$g}};
							$index = 0;
							$index++ until $gPeaks[$index] eq $redPeak;
							$pe = splice(@gPeaks,$index,1);
							$di = splice(@gDists,$index,1);
							$in = splice(@gInds,$index,1);
							@{$genePeaksSig{$currchr}{$g}} = @gPeaks;
							@{$distGeneBased{$currchr}{$g}} = @gDists;
							@{$geneIndex{$currchr}{$g}} = @gInds;                       
							$hs = $pe."_".$di;
							#print "\nP: $hs";
							@matrixElems = @{$hasPeak{$g}{basename($input[$i])}};
							#print "\nD: ",join(", ",@matrixElems);
							if (scalar @matrixElems == 1)
							{
								#print "\nF: ",join(", ",@matrixElems),"\n";
								@{$hasPeak{$g}{basename($input[$i])}} = ();
							}
							else
							{
								$hsi = grep { 
									${$hasPeak{$g}{basename($input[$i])}}[$_] =~ m/$hs/ 
								} 0..$#matrixElems;
								my @z = splice(@matrixElems,$hsi,1);
								#print "\nF: ",join(", ",@z),"\n";
								@{$hasPeak{$g}{basename($input[$i])}} = 
									@matrixElems;
							}
						}
					}
					# Empty the temp hashes
					%geneDist = ();
					%peakInd = ();
				}
			}
			if ($redun eq "peakcentric") { }
		}
		
		# A function must be added to correct the above hashes for gene-centric
		# redundancy
		
		# If no test performed, no need to do anything with background file
		if ($test ne "none") 
		{
			# Could not have genes at chromosome (unlikely for background...)
			if ($backID{$currchr})
			{
				@geneids = @{$backID{$currchr}};
				@starts = @{$backStart{$currchr}};
				@ends = @{$backEnd{$currchr}};
			
				for ($j=0; $j<@starts; $j++)
				{
					for ($k=0; $k<@modes; $k++)
					{
						$currdist = &dist($starts[$j],$ends[$j],$modes[$k]);
						if ($currdist > $span[0] && $currdist < $span[1])
						{
							push(@backallpeakids,$peakids[$k]);
							push(@backassgenes,$geneids[$j]);
						}
					}
				}
			}
		}
	}
	
	# Get peak counts in significant and background file
	%countSig = &unique(@sigallpeakids);
	%countBack = &unique(@backallpeakids) if ($test ne "none");
	%countGenesSig = &unique(@sigassgenes);
	%countGenesBack = &unique(@backassgenes) if ($test ne "none");
	
	# Run hypergeometric test
	if ($test eq "hypgeom")
	{
		&disp("Running hypergeometric test for each assigned query region..."); 
		@elems = keys(%countSig);
		for ($j=0; $j<@elems; $j++)
		{
			$p = abs(1 - hypergeom_cdf($lensig,$lenback - $lensig,
				$countBack{$elems[$j]},$countSig{$elems[$j]}));
			($p*@elems > 1) ? ($cp = 1) : ($cp = $p*@elems); # Bonferroni type correction
			$finalPeaks{$elems[$j]} = 
				"$p\t$cp\t$countSig{$elems[$j]}/$countBack{$elems[$j]}" 
					if ($p < $pval);
		}
	}
	elsif ($test eq "chi2")
	{
		&disp("Running chi-square test for each assigned peak..."); 
		@elems = keys(%countSig);
		for ($j=0; $j<@elems; $j++)
		{
			$p = &chisquarecont($countSig{$elems[$j]},
				$lensig - $countSig{$elems[$j]},$countBack{$elems[$j]},
				$lenback - $countBack{$elems[$j]});
			($p*@elems > 1) ? ($cp = 1) : ($cp = $p*@elems); # Bonferroni type correction
			$finalPeaks{$elems[$j]} = 
				"$p\t$cp\t$countSig{$elems[$j]}/$countBack{$elems[$j]}" 
					if ($p < $pval);
		}
	}
	elsif ($test eq "none")
	{
		&disp("No statistical testing performed...");   
		@elems = keys(%countSig);
		for ($j=0; $j<@elems; $j++)
		{
			$p = "NA";
			$cp = "NA";
			$finalPeaks{$elems[$j]} = "$p\t$cp\t$countSig{$elems[$j]}";
		}
	}
	
	my $sp = keys(%finalPeaks);
	my $ap = @elems;
	my $bp = keys(%countBack) if ($test ne "none");
	my $sag = keys(%countGenesSig);
	my $sbg = keys(%countGenesBack) if ($test ne "none");
	if ($test ne "none")
	{
		&disp("\nAssigned peaks in significant list : $ap out ".
			"of $lenpeak[$i] peaks in $sag out of $lensig genes");
		&disp("Assigned peaks in background : $bp out of $lenpeak[$i] peaks ".
			"in $sbg out of $lenback genes");
		&disp("Over-represented at p-value<$pval : $sp out of $ap\n");
	}
	else
	{
		&disp("\nAssigned peaks in gene list : $ap out ".
			"of $lenpeak[$i] peaks in $sag out of $lensig genes");
	}
	
	# Free some memory...
	($sp,$ap,$bp,$sag,$sbg,@sigassgenes,@backassgenes,
		%countGenesSig,%countGenesBack) = 
	(undef,undef,undef,undef,undef,undef,undef,undef,undef);
	
	# Construct output
	foreach my $opt (@out)
	{
		if ($opt eq "stats")
		{
			my $outfile = &create_output_file($input[$i],$opt);
			my $co;
			&disp("Writing output in $outfile...");
			open(OUTPUT,">$outfile");
			foreach $co (sort(keys(%finalPeaks)))
			{
				print OUTPUT "$co\t$finalPeaks{$co}\n";
			}
		}
		if ($opt =~ /gff-peak/)
		{ 
			my $outfile = &create_output_file($input[$i],$opt);
			&disp("Writing output in $outfile...");
			open(OUTPUT,">$outfile");
			print OUTPUT "Chromosome\tProgram\tFeature\tStart\tEnd\tp-value\t".
				"Strand\tFrame\tGeneID\tPeakID\tCorrected p-value\t".
				"Enrichment\tDistance\n" if ($opt eq "gff-peak-db");
			foreach $currchr (@peakchrs)
			{
				foreach $currpeak (keys(%finalPeaks))
				{   
					if ($peaksGenesSig{$currchr}{$currpeak})
					{
						@currgenes = @{$peaksGenesSig{$currchr}{$currpeak}};
						my ($q,$cq,$r) = split(/\t/,$finalPeaks{$currpeak});
						my $pos = 0;
						foreach $currgene (@currgenes)
						{   
							$currout = "$currchr\tHyperAssignPeaks\tPeak\t".
									   "${$peakStarts{$currchr}}[${$peakIndex{$currchr}{$currpeak}}[$pos]]\t".
									   "${$peakEnds{$currchr}}[${$peakIndex{$currchr}{$currpeak}}[$pos]]\t".
									   "$q\t+\t\.\t$currgene\t$currpeak\t$cq\t$r\t".
									   "${$distPeakBased{$currchr}{$currpeak}}[$pos]\n";
							print OUTPUT "$currout";
							$pos++;
						}
					}
				}
			}
			close(OUTPUT);
		}
		if ($opt =~ /gff-gene/)
		{
			my $outfile = &create_output_file($input[$i],$opt);
			&disp("Writing output in $outfile...");
			open(OUTPUT,">$outfile");
			print OUTPUT "Chromosome\tProgram\tFeature\tStart\tEnd\t".
				"p-value\tStrand\tFrame\tPeakID\tGeneID\tCorrected p-value\t".
				"Enrichment\tDistance\n" if ($opt eq "gff-gene-db");
			foreach $currchr (@peakchrs)
			{
				foreach $currgene (keys(%{$genePeaksSig{$currchr}}))
				{   
					my $pos = 0;
					foreach $currpeak (@{$genePeaksSig{$currchr}{$currgene}})
					{
						if ($finalPeaks{$currpeak})
						{
							my ($q,$cq,$r) = split(/\t/,$finalPeaks{$currpeak});
							my $s = ${$sigStart{$currchr}}[${$geneIndex{$currchr}{$currgene}}[$pos]];
							my $e = ${$sigEnd{$currchr}}[${$geneIndex{$currchr}{$currgene}}[$pos]];
							my $st = "+";
							if ($s > $e)
							{
								($s,$e) = ($e,$s); # Swap them
								$st = "-";
							}
							$currout = "$currchr\tHyperAssignPeaks\tGene\t$s\t$e\t$q\t$st\t\.\t$currpeak\t".
									   "$currgene\t$cq\t$r\t${$distGeneBased{$currchr}{$currgene}}[$pos]\n";
							print OUTPUT "$currout";
						}
						$pos++;
					}
				}
			}
			close(OUTPUT);
		}
		if ($opt eq "peak")
		{
			my $outfile = &create_output_file($input[$i],$opt);
			&disp("Writing output in $outfile...");
			open(OUTPUT,">$outfile");
			foreach $currchr (@peakchrs)
			{
				foreach $currpeak (keys(%finalPeaks))
				{   
					if ($peaksGenesSig{$currchr}{$currpeak})
					{
						@currgenes = @{$peaksGenesSig{$currchr}{$currpeak}};            
						print OUTPUT "$currpeak\t",join("\ ",@currgenes),"\n";
					}
				}
			}
			close(OUTPUT);
		}
		if ($opt eq "gene")
		{
			my $outfile = &create_output_file($input[$i],$opt);
			&disp("Writing output in $outfile...");
			open(OUTPUT,">$outfile");
			foreach $currchr (@peakchrs)
			{
				foreach $currgene (keys(%{$genePeaksSig{$currchr}}))
				{   
					my @outpeaks;
					foreach $currpeak (@{$genePeaksSig{$currchr}{$currgene}})
					{
						push(@outpeaks,$currpeak) if ($finalPeaks{$currpeak});
					}
					print OUTPUT "$currgene\t",join("\ ",@outpeaks),"\n" if (@outpeaks);
				}
			}
			close(OUTPUT);
		}
		if ($opt eq "pretty-peak")
		{
			my $outfile = &create_output_file($input[$i],$opt);
			&disp("Writing output in $outfile...");
			open(OUTPUT,">$outfile");
			print OUTPUT "PeakID/GeneID\tp-value/Distance\tBonferroni p-value/Strand\tEnrichment\n\n";
			foreach $currchr (@peakchrs)
			{
				foreach $currpeak (keys(%finalPeaks))
				{   
					if ($peaksGenesSig{$currchr}{$currpeak})
					{
						print OUTPUT "$currpeak\t$finalPeaks{$currpeak}\n";
						@currgenes = @{$peaksGenesSig{$currchr}{$currpeak}};
						my $pos = 0;
						foreach $currgene (@currgenes)
						{   
							my $ct = "+";
							my $cd = ${$distPeakBased{$currchr}{$currpeak}}[$pos];
							$ct = "-" if ($cd > 0);
							print OUTPUT "$currgene\t$cd\t$ct\t\n";
							$pos++;
						}
						print OUTPUT "\n";
					}
				}
			}
			close(OUTPUT);
		}
		if ($opt eq "pretty-gene")
		{
			my $outfile = &create_output_file($input[$i],$opt);
			&disp("Writing output in $outfile...");
			open(OUTPUT,">$outfile");
			print OUTPUT "GeneID/PeakID\tStrand/p-value\tBonferroni p-value\tEnrichment\tDistance\n\n";
			foreach $currchr (@peakchrs)
			{
				foreach $currgene (keys(%{$genePeaksSig{$currchr}}))
				{   
					my $s = ${$sigStart{$currchr}}[${$geneIndex{$currchr}{$currgene}}[0]];
					my $e = ${$sigEnd{$currchr}}[${$geneIndex{$currchr}{$currgene}}[0]];
					my $st = "+";
					if ($s > $e)
					{
						($s,$e) = ($e,$s); # Swap them
						$st = "-";
					}
					my (@outpeaks,@dis);
					my $pos = 0;
					foreach $currpeak (@{$genePeaksSig{$currchr}{$currgene}})
					{
						if ($finalPeaks{$currpeak})
						{
							push(@outpeaks,$currpeak);
							push(@dis,${$distGeneBased{$currchr}{$currgene}}[$pos]);        
						}
						$pos++;
					}
					if (@outpeaks)
					{
						print OUTPUT "$currgene\t$st\t\t\t\n";
						for (my $cpo=0; $cpo<@outpeaks; $cpo++)
						{
							print OUTPUT "$outpeaks[$cpo]\t$finalPeaks{$outpeaks[$cpo]}\t$dis[$cpo]\n";
						}
						print OUTPUT "\n";
					}
				}
			}
			close(OUTPUT);
		}
		if ($opt eq "peakdata")
		{
			my $outfile = &create_output_file($input[$i],$opt);
			&disp("Writing output in $outfile...");
			open(OUTPUT,">$outfile");
			print OUTPUT "$phead\n" if ($phead);
			my @fkeys = keys(%finalPeaks);
			foreach $currpeak (@fkeys)
			{
				print OUTPUT "$allPeakData{$currpeak}\n";
			}
			close(OUTPUT);
		}
		if ($opt eq "bed")
		{
			my $outfile = &create_output_file($input[$i],$opt);
			&disp("Writing output in $outfile...");
			open(OUTPUT,">$outfile");
			my @fkeys = keys(%finalPeaks);
			foreach $currpeak (@fkeys)
			{
				my @spl = split(/\t/,$allPeakData{$currpeak});
				print OUTPUT $spl[0]."\t".$spl[1]."\t".$spl[2]."\t".$spl[$pcols[0]]."\t".$spl[$pcols[1]]."\t.\n";
			}
			close(OUTPUT);
		}
		print_gene_or_peak($input[$i],"all-peak",%peaksGenesSig) 
			if ($opt eq "all-peak");
		print_gene_or_peak($input[$i],"all-gene",%genePeaksSig) 
			if ($opt eq "all-gene");
	}
}

&print_matrix(\%hasPeak,"matrix-number",\%hasExpression) 
	if (&smatch("matrix-number",@out));
&print_matrix(\%hasPeak,"matrix-presence",\%hasExpression) 
	if (&smatch("matrix-presence",@out));
&print_matrix(\%hasPeak,"matrix-peaks",\%hasExpression) 
	if (&smatch("matrix-peaks",@out));

$date = &now;
&disp("$date - Finished!\n\n");

sub hypergeom_pdf
{   
    my ($n,$m,$N,$i) = @_;
    my $loghyp1 = logfact($m) + logfact($n) + logfact($N) + logfact($m+$n-$N);
    my $loghyp2 = logfact($i) + logfact($n-$i) + logfact($m+$i-$N) + 
		logfact($N-$i) + logfact($m+$n);
    return(exp($loghyp1 - $loghyp2));
}

sub hypergeom_cdf
{  
    my ($n,$m,$N,$i) = @_; 
    my @lessthan = (0..$i);
    my $cum = 0; #;-)
    foreach my $j (@lessthan)
    {
        $cum += hypergeom_pdf($n,$m,$N,$j);
    }
    return($cum);
}

sub logfact 
{
    my $a = $_[0];
    return gammaln($a + 1.0);
}

sub gammaln 
{
    my $xx = $_[0];
    my $eps = 1; 
    my @cof = (76.18009172947146, -86.50532032941677,
               24.01409824083091, -1.231739572450155,
               0.12086509738661e-2, -0.5395239384953e-5);
    my $y = my $x = $xx;
    my $tmp = $x + 5.5;
    $tmp -= ($x + .5) * log($tmp);
    my $ser = 1.000000000190015;
    for my $j (0..5) 
    {
        $ser += $cof[$j]/++$y;
    }
    ($x = $eps) if ($x == 0);
    return(-$tmp + log(2.5066282746310005*$ser/$x));
}

sub dist
{
    my ($s,$e,$m) = @_;
    my $d;
    my ($d1,$d2);
    if (ref $m eq "ARRAY") { # Boundary mode
		my @bs = @{$m};
		if ($s < $e) # + strand
		{
			$d1 = $bs[0] - $s;
			$d2 = $bs[1] - $s;
		}
		else # - strand, reversed while reading file
		{
			$d1 = $s - $bs[0];
			$d2 = $s - $bs[1];
		}
		if ($nomode eq "boundary" || $nomode eq "closest")
		{
			$d = ($d1 < $d2) ? ($d1) : ($d2);
		}
		elsif ($nomode eq "furthest")
		{
			$d = ($d1 < $d2) ? ($d2) : ($d1);
		}
	}
	else { # Center mode
		if ($s < $e) # + strand
		{
			#(($m > $s) && ($m < $e)) ? ($d = 0) : ($d = $m - $s);
			$d = $m - $s;
		}
		else # - strand, reversed while reading file
		{
			#(($m < $s) && ($m > $e)) ? ($d = 0) : ($d = $s - $m);
			$d = $s - $m;
		}
	}
    return($d);
}

sub chisquarecont
{
    my ($a,$b,$c,$d) = @_;
    my $tot = $a + $b + $c + $d;
    my $A = (($a + $b)*($a + $c))/$tot;
    my $B = (($a + $b)*($b + $d))/$tot;
    my $C = (($c + $d)*($a + $c))/$tot;
    my $D = (($c + $d)*($b + $d))/$tot;
    
    # Chi square statistic
    my $x2 = (($a - $A)**2)/$A + (($b - $B)**2)/$B + (($c - $C)**2)/$C +
             (($d - $D)**2)/$D;

    return(1-chdtr(1,$x2));
}

sub print_gene_or_peak
{
    my ($infile,$otype,%inhash) = @_;
    my ($outchr,$ind,$outhash);
    my (@k,@v);
    my $outfilename = &create_output_file($infile,$otype);
    &disp("Writing output in $outfilename...");
    my @chrs = keys(%inhash);
    open(OUTPUT,">$outfilename");
    foreach $outchr (sort(@chrs))
    {
        $outhash = $inhash{$outchr};
        @k = keys(%$outhash);
        @v = values(%$outhash);
        for ($ind=0; $ind<@k; $ind++)
        {
            print OUTPUT "$k[$ind]\t";
            print OUTPUT join("\ ",@{$v[$ind]});
            print OUTPUT "\n";
        }
    }
    close(OUTPUT);
}

sub print_matrix
{
    my ($inhash,$type,$exprhash) = @_;
    my ($row,$column,$colhash);
    my $outfilename = &create_output_file($input[0],$type);
    &disp("Writing output in $outfilename...");
    my @rows = keys(%$inhash);
    open(OUTPUT,">$outfilename");
    my $headhash = $inhash->{$rows[0]};
    my @headers = keys(%$headhash);
    (!defined($exprhash->{"header"})) ? 
		(print OUTPUT "GeneID\t",join("\t",@headers),"\n") :
		(print OUTPUT "GeneID\t",join("\t",@headers),"\t",
			$exprhash->{"header"},"\n");
    if ($type =~ m/number/i)
    {
        foreach $row (@rows)
        {
            $colhash = $inhash->{$row};
            my @tmp = keys(%$colhash);
            my @v;
            foreach $column (@tmp)
            {
                (defined($colhash->{$column}) 
					&& scalar @{$colhash->{$column}} > 0) ? 
					(push(@v,scalar(@{$colhash->{$column}}))) : (push(@v,0));
            }
            (!$exprhash->{"data"}->{$row}) ? 
				(print OUTPUT "$row\t",join("\t",@v),"\n") :
				(print OUTPUT "$row\t",join("\t",@v),"\t",
					$exprhash->{"data"}->{$row},"\n");
        }
    }
    elsif ($type =~ m/presence/i)
    {
        foreach $row (@rows)
        {
            $colhash = $inhash->{$row};
            my @tmp = keys(%$colhash);
            my @v;
            foreach $column (@tmp)
            {
                (defined($colhash->{$column}) 
					&& scalar @{$colhash->{$column}} > 0) ? 
					(push(@v,"+")) : (push(@v,"-"));
            }
            (!$exprhash->{"data"}->{$row}) ? 
				(print OUTPUT "$row\t",join("\t",@v),"\n") :
				(print OUTPUT "$row\t",join("\t",@v),"\t",
					$exprhash->{"data"}->{$row},"\n");
        }
    }
    elsif ($type =~ m/peaks/i)
    {
        foreach $row (@rows)
        {
            $colhash = $inhash->{$row};
            my @tmp = keys(%$colhash);
            my @v;
            foreach $column (@tmp)
            {
                (defined($colhash->{$column}) 
					&& scalar @{$colhash->{$column}} > 0) ? 
					(push(@v,join("; ",@{$colhash->{$column}}))) : 
					(push(@v,"NP"));
            }
            (!$exprhash->{"data"}->{$row}) ? 
				(print OUTPUT "$row\t",join("\t",@v),"\n") :
				(print OUTPUT "$row\t",join("\t",@v),"\t",
					$exprhash->{"data"}->{$row},"\n");
        }
    }
    close(OUTPUT);
}

sub create_output_file
{
    my ($in,$type) = @_;
    my $ext;
    my ($base,$dir) = fileparse($in,'\.[^.]*');
    if ($type =~/gff/)
    {
        ($type =~/db/) ? ($ext = ".txt") : ($ext = ".gff");
    }
    elsif ($type =~/bed/) { $ext = ".bed" }
    else { $ext = ".txt" }
    if ($type =~ /matrix/)
    {
        return($dir."gene-peak-$type.txt");
    }
    if ($type =~ /bed/)
    {
        return($dir.$base."_ASSIGNED".$ext);
    }
    else
    {
        return($dir.$base."_ASSIGNED_".$type.$ext);
    }
}

sub check_unique
{
    my %h = %{$_[0]};
    my @vals = values(%h);
    my %ch = &unique(@vals);
    (scalar @vals == scalar keys(%ch)) ? (return(1)) : (return(0));
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

# Just parse parameters, checking is now performed by the module
sub check_inputs
{
    GetOptions("input|i=s{,}" => \@input,
               "region|r=s" => \$region,
               "background|b=s" => \$background,
               "span|n=i{,}" => \@span,
               "where|w=s" => \$where,
               "idstrand|c=i{,}" => \@sbcols,
               "idmode|m=i{,}" => \@pcols,
               "expression|e=i{,}" => \@expcols,
               "test|t=s" => \$test,
               "pvalue|p=f" => \$pval,
               "redundancy|d=s" => \$redun,
               "nomode|a=s" => \$nomode,
               "outformat|o=s{,}" => \@out,
               "log|l" => \$log,
               "silent|s" => \$silent,
               "help|h" => \$help);
    if ($help)
    {
        &program_usage;
        exit;
    }
    
    # Check fatal
    my $stop;
    $stop .= "--- Please specify input file(s) ---\n" if (!@input);
    $stop .= "--- Please specify region file(s) ---\n" if (!$region);
    
    if ($stop) {
        &disp("$stop\n");
        &disp("Type perl $scriptname --help for help in usage.\n\n");
        exit;
    }
    
    # Check required packages
    my $status;
    if ($test eq "chi2")
    {
        $status = eval { &try_module("Math::Cephes") };
        if ($status)
        {
            &disp("Module Math::Cephes is required to perform the chi-square ",
				"test! Using default (none)...");
            $test = "none";
        }
        else { use Math::Cephes; }
    }
    if (@out && &smatch("matrix",@out))
    {
        $status = eval { &try_module("Tie::IxHash::Easy") };
        if ($status)
        {
            &disp("Module Tie::IxHash::Easy is required for one or more of ",
				"the selected outputs! Using default (gff-peak)...");
            @out = ("gff-peak");
        }
        else { use Tie::IxHash::Easy; }
    }
    
    # Check the rest
    # Check statistical test
    if ($test) {
		if ($test ne "hypgeom" && $test ne "chi2" && $test ne "none" 
			&& $test ne "auto")
		{
			&disp("test parameter should be one of \"hypgeom\", \"chi2\",",
				"\"none\" or \"auto\"! Using default (none)...");
			$test = "none";
		}
	}
    else # A test must be set
    {
        &disp("--test parameter not given! Using default (none)...");
        $test = "none";
    }
    # Check if span given
    if (!@span)
    {
        &disp("Search range from region start points (e.g. TSS) not given! ",
			"Using defaults (-10kbp,10kbp or -10kb)");
        @span = (-10000,10000);
        if ($where eq "coding")
        {
            @span = (-10000);
        }
    }
    # Check if where given
    if ($where ne "promoter" && $where ne "coding")
    {
        &disp("The where parameter must be one of \"promoter\" or \"coding\"! ",
			"Using default (promoter)...");
        $where = "promoter";
        if (@span)
        {           
            $span[1] = 10000;
        }
    }
    else
    {
        if ($where eq "promoter")
        {
            if (!$span[1])
            {
                &disp("The where parameter is \"promoter\" but the second ",
					"\"span\" argument is not given! Using default (10000)...");
                @span = ($span[0]);
            }
        }
        if ($where eq "coding")
        {
            if (@span)
            {
                @span = ($span[0]);
            }
        }
    }
    # Check if id and strand columns given for sig/back files
    if (!@sbcols)
    {
        &disp("Unique ID and strand columns for region and background files ",
			"not given! Using defaults as from BED format (4,6)...");
        @sbcols = (3,5);
    }
    else # Proper perl indexing
    {
        $sbcols[0]--;
        $sbcols[1]--;
    }
    # Check if id and mode columns given for peak files
    if (!@pcols)
    {
        &disp("Unique ID and mode columns for query region files not given!",
			"Using default ID as from BED format (4) and");
        &disp("query modes will be determined by the nomode parameter...");
        @pcols = (3);
    }
    else # Proper perl indexing
    {
        $pcols[0]--;
        if (!$pcols[1])
        {
            &disp("Input query regions modes not given! Action will be",
				"determined by the nomode parameter...");
        }
        else
        {
            $pcols[1]--;
        }
        if (!$pcols[2])
        {
            &disp("Additional query regions scores not given! They will not",
				"be reported...");
        }
        else
        {
            $pcols[2]--;
        }
    }
    # Check if expression column is given in genes file and set proper Perl indexing
    if (@expcols)
    {
        my $l = scalar @expcols;
        for (my $i=0; $i<$l; $i++)
        {
             $expcols[$i]--;
        }
    }
    # Check redundancy
    if ($redun)
    {
        if ($redun ne "all" && $redun ne "genecentric" 
            && $redun ne "peakcentric")
        {
            my $msg = "WARNING! redundancy parameter should be one of".
				"\"all\", \"genecentric\" or \"peakcentric\"!...".
				" Using default (all)...";
            &disp($msg);
            $redun = "all";
        }
    }
    else
    {
        &disp("Redundancy preference not given! Using default (all)");
        $redun = "all";
    }
    # Check what to do in no mode
    if ($nomode)
    {
        if ($nomode ne "center" && $nomode ne "boundary" 
			&& $nomode ne "closest" && $nomode ne "furthest")
        {
            my $msg = "WARNING! nomode parameter should be one of ".
				"\"center\", \"boundary\", \"closest\" or \"furthest\"! ".
				"Using default (center)...";
            &disp($msg);
            $nomode = "center";
        }
    }
    else
    {
        &disp("No mode preference not given! Using default (center)");
        $nomode = "center";
    }
    # Check proper output format
    if (@out)
    {
        foreach my $c (@out)
        {
            if ($c ne "bed" && $c ne "stats" && $c ne "gff-peak" 
				&& $c ne "gff-gene" && $c ne "peak" &&  $c ne "gene" 
                && $c ne "all-peak" && $c ne "all-gene" && $c ne "pretty-peak" 
                && $c ne "pretty-gene" && $c ne "gff-peak-db" 
                && $c ne "gff-gene-db" && $c ne "peakdata" 
                && $c ne "matrix-number" && $c ne "matrix-presence" 
                && $c ne "matrix-peaks")
            {
                my $msg = "WARNING! outformat parameter options should be one ".
					"or more of \"bed\", \"gff-peak\", \"gff-gene\",\n".
                    "\"peak\", \"gene\", \"all-peak\", \"all-gene\", ".
                    "\"pretty-peak\", \"pretty-gene\", \"gff-peak-db\", \n".
					"\"gff-gene-db\", \"peakdata\", \"stats\", ".
					"\"matrix-number\", \"matrix-presence\" or ".
                    "\"matrix-peaks\"! \nUsing default (\"gff-peak\")...";
                &disp($msg);
                @out = ("gff-peak");
            }
        }
    }
    else
    {
        &disp("Output format not given! Using default (gff-peak)");
        @out = ("gff-peak");
    }
}

sub try_module {
    my ($module,@fun) = @_;
    eval "require $module";
    if ($@) {
        my $killer = "Module $module is required to continue with the\n". 
            "execution. If you are in Windows and you have ActiveState Perl\n".
            "installed, use the Package Manager to get the module. If you \n".
            "are under Linux, log in as a super user (or use sudo under\n".
            "Ubuntu) and type \"perl -MCPAN -e shell\" (you will possibly\n".
            "have to answer some questions). After this type \"install \n".
            "$module\" to install the module. If you don't know how to\n".
            "install the module, contact your system administrator.\n";
        die "\n$killer\n\n";
    }
    else {
        if (@fun) {
            my $funs = join(" ",@fun);
            eval "use $module qw($funs)";
        }
        else { eval "use $module"; }
    }
}

sub now
{
    my ($self,$format) = @_;
    $format = "human" if (!$format);
    my ($sec,$min,$hour,$day,$month,$year) = localtime(time);
    $year += 1900;
    $month++;
    $month = "0".$month if (length($month)==1);
    $day = "0".$day if (length($day)==1);
    $hour = "0".$hour if (length($hour)==1);
    $min = "0".$min if (length($min)==1);
    $sec = "0".$sec if (length($sec)==1);
    ($format ne "machine") ? 
    (return($day."/".$month."/".$year." ".$hour.":".$min.":".$sec)) :
    (return($year.$month.$day.$hour.$min.$sec));
}

sub decide_header
{
    my $line= $_[0];
    $line =~ s/\r|\n$//g;
    my @cols = split(/\t/,$line);
    if ($cols[0] =~ m/^chr/ && $cols[1] =~ m/\d+/ && $cols[2] =~ m/\d+/)
    {
        return(0); # Does not contain a header, is proper bed line
    }
    else
    {
        return($line);
    }
}

sub smatch {
    my ($s,@a) = @_;
    if (grep(/^$s$/,@a)) {
        return(1);
    }
    return(0);
}

sub round {
    my ($number,$digits) = @_;
    if (!$digits) {
        return int($number + .5*($number <=> 0));
    }
    else {
        return(sprintf("%.".$digits."f",$number));
    }
}

sub disp {
    my @msg = @_;
    print STDERR "\n@msg" if (!$silent);
}

sub program_usage 
{
    # The look sucks here but it is actually good in the command line
    my $usagetext = << "END";
    
$scriptname
A perl program to assign ChIP-Seq peaks to a set of regulated genes based on
the gene-peak distances between the genes of the regulated set and the gene-
peak distances between the genes of a background set (e.g. the whole genome).
In order to find peaks associatevely enriched in the set of the regulated
genes as compared to the background set, the program uses the hypergeometric
test on the population of peaks at a p-value threshold. There can be multiple
output format, all of them containing the gene-peak associations in different
formats. The hypergeometric method is NOT verified and there are other methods
out there that may perform better. Use at your own risk. The tools works very
nicely to calculate peak-gene distances with a very nice output.

Author : Panagiotis Moulos (moulos\@fleming.gr)

Main usage
$scriptname --input input(s) --region regfile --background background [OPTIONS]

--- Required ---
  --input|i  file(s)    Peak BED file(s) containing a column with a
            UNIQUE peak ID and a column with the peak mode (the
            point with the highest tag pile-up) or a location that
            the user thinks as the point of the peak from which 
            the distance to the genes will be calculated.
  --region|r  file  A BED file with the set of regulated genes,
            containing a column with a UNIQUE gene (or whatever
            region) ID and a column with the gene strand.
--- Optional ---
  --background|b        A BED file with the set of background 
            genes, containing a column with a UNIQUE gene (or 
            whatever region) ID and a column with the gene strand.
            Required when running a statistical test.
  --span|d      Use this option to set the genomic span (distance
            upstream and downstream from TSS) into which the program
            will look for any peaks. It should be two values (e.g.
            --span -50000 50000) and defaults to (-100000,100000).
  --where|w         Use this parameter to tell Assign.pm whether to 
            check for query regions. When  "promoter", it will check 
            for queries upstream and downstream of subjects according 
            to --span. If "coding", the second --span argument is 
            ignored and it automatically becomes 
            subject_end - subject_start + span_1 to check for
            queries inside the subject regions. If "downtes", the 
            second --span argument is ignored and it automatically 
            becomes subject_end + span_1 to check for presence
            e.g. downstream of transcriptional end sites. The option 
            names "promoter", "coding" and "downtes" are indicative. 
            Queries and subjects can be any genomic regions of interest.
  --idstrand|t      The columns in BOTH the gene files where their
            unique IDs and strands are. You should provide two values
            (e.g. --idstrand 4 5) where the first denotes the unique
            ID column and the second the strand column. It defaults
            to (4,5).
  --idmode|m        The columns in the peak files where their unique
            IDs and modes are. You should provide two values (e.g
            --idmode 4 5) where the first denotes the unique ID
            column and the second the mode column. It defaults to (4,5).
            Optionally, you can provide three values, where the 3rd 
            represents a peak score if available. This will be reported
            when using the "matrix-peak" output. The values must be 
            provided strictly with the following order: id column, 
            mode column, score column.
  --test|t      What over-representation statistical test to perform.
            Can be one of hypgeom for hypergeometric test, chi2 for
            chi-square test, auto for automatic selection and none for
            no testing. Defaults to hypgeom.
  --pvalue|p        The hypergeometric test p-value threshold. It
            defaults to 0.05.
  --redundancy|d        The reundancy level when assigning peaks to genes.
            It can be "genecentric" for assigning multiple peaks to one
            gene, "peakcentric" to allow one peak to many genes or "all"
            to allow a multi-to-multi assignment (default). The option
            "peakcentric" is not yet implemented.
  --nomode|a        If query region modes are not provided, then the
            anchor points for determining distance from region areas are
            calculated from the center of the query regions ("center",
            default), from the closest boundary ("closest"), from the
            furthest boundary ("furthest") or  automatically 
            ("boundary") where the boundary closest to the start of the 
            region are will be used.
  --outformat|o     Use this option to determine which output format
            filetype(s) you wish to retrieve.   Possible choices are:
                "stats" for retrieving the significantly associated 
                peaks with their p-values, Bonferroni corrected p-values
                and enrichment ratios.
                "gff-peak" for retrieving a peak-based gff file which
                contains additional columns with peak ids, distances
                and enrichment ratios. The score column is the p-value.
                "gff-gene" for similar to "gff-peak" but gene-based.
                "gff-peak-db" for same as "gff-peak" but with a header,
                suitable for incorporating to a database.
                "gff-gene-db" for same as "gff-gene" but with a header,
                suitable for incorporating to a database.
                "peak" for a simple file which contains the significantly
                associated peak IDs in the first column and a list of 
                associated genes in the second column.
                "gene" for similar to "peak" but gene based.
                "all-peak" for a simple file which contains ALL (based on
                distance) associated peak IDs in the first column and a 
                list of associated genes in the second column.
                "all-gene" for similar to "all-peak" but gene based.
                "pretty-peak" for retrieving a more human-readable format
                "bed" for retrieving a 6-column BED file suitable for a 
                genome browser without additional data.
                quite self-explicating (please see output).
                "pretty-gene" similar to "pretty-peak" but gene-based
                (please see output).
                "peakdata" for retrieving only the assigned peaks from
                the original peak file.
                "matrix-number" to retrieve a spreadsheet-like file where 
                rows correspond to genes (or the --region file) and columns
                correspond to peak files. The cell (i,j) contains the number
                of peaks in peak file j assigned to gene i.
                "matrix-presence" to retrieve a spreadsheet-like file where 
                rows correspond to genes (or the --region file) and columns
                correspond to peak files. The cell (i,j) contains "+" if peak
                in peak file j assigned to gene i, "-" otherwise.
                "matrix-peaks" to retrieve a spreadsheet-like file where 
                rows correspond to genes (or the --region file) and columns
                correspond to peak files. The cell (i,j) contains the peaks
                in peak file j assigned to gene i, "NP" otherwise.
            Example: --outformat stats gff-peak pretty-gene matrix
  --source|u        Use this option to set the online data source in
            the case of selecting one of the prefefined region templates
            with --region. Can be one of "ucsc", "refseq" or "ensembl".
            Defaults to "ensembl".
  --gversion|g      Use this option to set the version of the genome
            for data to be downloaded. It can be "hg19", "hg18" for human
          "mm10", "mm9" for mouse, "rn5" for rat, "dm3" for fruitfly and
          "danrer7" for zebrafish.
  --expression|e      An array of column numbers which may contain
            expression (or other custom) values in the region file.
  --log|l       Output a log file. It can be a file name or empty for
                auto-generation.
  --silent|s        Use this option if you want to turn informative 
            messages off.
  --help|h      Display this help text.
    
The main output of the program is up to nine files with information on gene-peak
association.

END
    print $usagetext;
    exit;
}
