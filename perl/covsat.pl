#!/usr/bin/perl -w

# covsat.pl
# A Perl script to create coverage saturation plots for deep sequencing experiments using
# BEDTools and R
#
# Author      : Panagiotis Moulos (moulos@fleming.gr)
# Created     : 09 - 07 - 2012 (dd - mm - yyyy)
# Last Update : 27 - 10 - 2012 (dd - mm - yyyy)
# Version     : 1.0

# TODO: Fix the streaming exons

use strict;
use IO::File;
use Getopt::Long;
use File::Temp;
use File::Spec;
use File::Basename;
use Switch;
use DBI;
use POSIX qw(floor ceil);

use constant GPATH => "/opt/NGSTools/BEDTools/genomes";
use constant MAXCORES => 6;
use constant REMOTE_HOST => "genome-mysql.cse.ucsc.edu";
use constant REMOTE_USER => "genome";
use constant DISTANT => 50000;

# Make sure output is unbuffered
select(STDOUT);
$|=1;

# On Ctrl-C, do cleanup
$SIG{INT} = \&catch_cleanup;

# Set defaults
our $scriptname = "covsat.pl";
our @input; 			# Input BED files (BED supported right now)
our $type;				# Genome, features or custom coverage
our $genome;			# The genome to be used
our $mappable;		# The mappable genome fraction
our $resolution;		# The resolution (density) of the plot
our $step;				# The step that the data are splt to (e.g. 500kbp)
our $opath;			# Path to generate output plots
our $ncore;			# How many cores to use for calculation
our $bam;				# BAM files as input?
our $figtype = "png";	# Output figure type
our @dbdata;			# Database connection credentials (typically, gbuser)
our $silent = 0;		# Display verbose messages
our $help = 0;			# Help?

# Check inputs
&check_inputs;

# Define the genome file
my $gfile;
switch ($genome)
{
	case "hg18" { $gfile = File::Spec->catfile(GPATH,"human.hg18.genome") }
	case "hg19" { $gfile = File::Spec->catfile(GPATH,"human.hg19.genome") }	
	case "mm8" { $gfile = File::Spec->catfile(GPATH,"mouse.mm8.genome") }
	case "mm9" { $gfile = File::Spec->catfile(GPATH,"mouse.mm9.genome") }
	else { $gfile = $genome; }
}

# Create and link a temporary directory for file operations
our $tmpdir = File::Temp->newdir();

# Record progress
my $date = &now;
disp("$date - Started...");

# If bam files are provided, they must be converted to bed first in tmpdir (ouch!!!)
@input = &bam_to_bed(@input) if ($bam);

# Process bed file(s)
my %covhash;
our %tcov;
foreach my $f (@input) 
{ 
	my ($nreads,$victim,$cur,$x,$y,$nam,$Rscript,$RSfile,$ROfile,$Tfile,$fh,$ff);
	my (@xx,@yy);
	%covhash = ();
	%tcov = ();
	
	disp("Preprocessing file $f...");
	disp("  Counting number of reads for file $f...");
	$nreads = &count_lines($f);
	disp("  Shuffling file $f...");
	$ff = basename($f);
	$victim = File::Spec->catfile($tmpdir,"$ff.shf");
	`shuf $f > $victim`;
	$cur = 0;
	
	if ($ncore == 1)
	{		
		if ($resolution) 
		{
			$step = ceil($nreads/$resolution);
			disp("  Calculating coverage saturation using a resolution of $resolution points in the coverage plot...");
		}
		else { disp("  Calculating coverage saturation using a step of $step for $nreads reads..."); }
		while ($cur < $nreads)
		{
			$cur += $step;
			disp("    Current depth: $cur");
			$tcov{$cur} = File::Spec->catfile($victim.".cov.".$cur);
			switch ($type)
			{
				case /genome/
				{
					`head -$cur $victim | sort -k1,1 -k2g,2 -u | genomeCoverageBed -i stdin -g $gfile -max 50 | grep genome > $tcov{$cur}`;
					$covhash{$cur} = &calc_genome_coverage($tcov{$cur},$mappable);
				}
				case /features/
				{
					# Probably some function call here
				}
				case /custom/
				{
					`head -$cur $victim | sort -k1,1 -k2g,2 -u | coverageBed -a stdin -b $genome -hist | grep all > $tcov{$cur}`;
					$covhash{$cur} = &calc_custom_coverage($tcov{$cur});
				}
			}
		}
	}
	else
	{
		my @depth;
		my $pl = Parallel::Loops->new($ncore);
		$pl->share(\%covhash);
		
		if ($resolution)
		{
			$step = ceil($nreads/$resolution);
			disp("  Parallely calculating coverage saturation using a resolution of $resolution points in the coverage plot...");
		}
		else { disp("  Parallely calculating coverage saturation using a step of $step for $nreads reads..."); }
		while ($cur < $nreads)
		{
			$cur += $step;
			push(@depth,$cur);
			$tcov{$cur} = File::Spec->catfile($victim.".cov.".$cur);
		}	
		$pl->share(\%tcov);		
		
		$pl->foreach(\@depth, sub {
			disp("    Current depth: $_");
			if ($type =~ m/genome/)
			{
					`head -$_ $victim | sort -k1,1 -k2g,2 -u | genomeCoverageBed -i stdin -g $gfile -max 50 | grep genome > $tcov{$_}`;
					$covhash{$_} = &calc_genome_coverage($tcov{$cur},$mappable);
			}
			elsif ($type =~ m/features/) { }
			elsif ($type =~ m/custom/)
			{
				`head -$_ $victim | sort -k1,1 -k2g,2 -u | coverageBed -a stdin -b $genome -hist | grep all > $tcov{$cur}`;
				$covhash{$_} = &calc_custom_coverage($tcov{$cur});
			}
		});	
	}
	
	# Figures and texts
	foreach (sort { $a <=> $b } keys(%covhash))
	{
		push(@xx,$_);
		push(@yy,$covhash{$_});
	}
	$x = "c(".join(",",@xx).")";
	$y = "c(".join(",",@yy).")";
	$nam = File::Spec->catfile($opath,"$ff.cov.$figtype");
	$Rscript = &write_R($x,$y,$figtype,$nam);
	
	$RSfile = File::Spec->catfile($tmpdir,"$ff.cov.R");
	$ROfile = File::Spec->catfile($tmpdir,"$ff.cov.Rout");
	$Tfile = File::Spec->catfile($opath,"$ff.cov.txt");
	
	# Make the figure
	$fh = IO::File->new();
	$fh->open(">$RSfile");
	print $fh $Rscript;
	$fh->close;
	
	#`R CMD BATCH --vanilla $RSfile $ROfile`;
	`Rscript --vanilla $RSfile`;
	
	# Write the file
	$fh = IO::File->new();
	$fh->open(">$Tfile");
	print $fh "Depth\tCoverage\n";
	foreach (sort { $a <=> $b } keys(%covhash))
	{
		print $fh "$_\t$covhash{$_}\n";
	}
	$fh->close;
}

$date = &now;
disp("$date - Finished!\n");

#&printHash(\%covhash);


# DB AND MAIN ARITHMETIC FUNCTIONS

sub get_tssh
{
	my $source = lc($_[0]);
	my $tssh; # Hash of hashes: {chr}->{strand}->{gene}->{start} most efficient memory usage
	switch($source)
	{
		case /ensembl/
		{
			$tssh = &get_ensembl_tssh;
		}
		case /refseq/
		{
			$tssh = &get_refseq_tssh;
		}
		case /ucsc/
		{
			$tssh = &get_ucsc_tssh;
		}
	}
	return($tssh);
}

sub get_ensembl_tssh
{
	my %etssh;
	my ($data,$oldStart,$newStart,$newEnd);
	my $conn = &open_connection($genome);
	my $q = "SELECT `chrom`, `strand`, `txStart`, `txEnd`, `name2` FROM `ensGene`";
	my $sth = $conn->prepare($q);
	$sth->execute();
	#while ($data = $conn->selectrow_hashref($q)) # There is something causing infinite loop
	while ($data = $sth->fetchrow_hashref())
	{
		if (!$etssh{$data->{"chrom"}}{$data->{"strand"}}{$data->{"name2"}})
		{
			($data->{"strand"} == "+") ?
			$etssh{$data->{"chrom"}}{$data->{"strand"}}{$data->{"name2"}} = {$data->{"txStart"}} :
			$etssh{$data->{"chrom"}}{$data->{"strand"}}{$data->{"name2"}} = {$data->{"txEnd"}};
		}
		else
		{
			$oldStart = $etssh{$data->{"chrom"}}{$data->{"strand"}}{$data->{"name2"}};
			$newStart = $data->{"txStart"};
			$newEnd = $data->{"txEnd"};
			if ($data->{"strand"} == "+")
			{
				$etssh{$data->{"chrom"}}{$data->{"strand"}}{$data->{"name2"}} = $newStart if ($newStart < $oldStart);
			}
			elsif ($data->{"strand"} == "-")
			{
				$etssh{$data->{"chrom"}}{$data->{"strand"}}{$data->{"name2"}} = $newEnd if ($newEnd > $oldStart);
			}
		}
	}
	&close_connection($conn);
	return (\%etssh);
}

sub get_refseq_tssh
{
	my %rtssh;
	my ($data,$oldStart,$newStart,$newEnd);
	my $conn = &open_connection($genome);
	my $q = "SELECT `geneName`, `chrom`, `strand`, `txStart`, `txEnd` FROM `refFlat`";
	my $sth = $conn->prepare($q);
	$sth->execute();
	#while ($data = $conn->selectrow_hashref($q))
	while ($data = $sth->fetchrow_hashref())
	{
		if (!$rtssh{$data->{"chrom"}}{$data->{"strand"}}{$data->{"geneName"}})
		{
			($data->{"strand"} == "+") ?
			$rtssh{$data->{"chrom"}}{$data->{"strand"}}{$data->{"geneName"}} = {$data->{"txStart"}} :
			$rtssh{$data->{"chrom"}}{$data->{"strand"}}{$data->{"geneName"}} = {$data->{"txEnd"}};
		}
		else
		{
			$oldStart = $rtssh{$data->{"chrom"}}{$data->{"strand"}}{$data->{"geneName"}};
			$newStart = $data->{"txStart"};
			$newEnd = $data->{"txEnd"};
			if ($data->{"strand"} == "+")
			{
				$rtssh{$data->{"chrom"}}{$data->{"strand"}}{$data->{"geneName"}} = $newStart if ($newStart < $oldStart);
			}
			elsif ($data->{"strand"} == "-")
			{
				$rtssh{$data->{"chrom"}}{$data->{"strand"}}{$data->{"geneName"}} = $newEnd if ($newEnd > $oldStart);
			}
		}
	}
	&close_connection($conn);
	return (\%rtssh);
}

sub get_ucsc_tssh
{
	my %utssh;
	my ($data,$oldStart,$newStart,$newEnd);
	my $conn = &open_connection($genome);
	my $q = "SELECT knownCanonical.chrom AS chrom, `chromStart`, `chromEnd`, `strand`, `transcript` ".
			"FROM `knownCanonical` INNER JOIN `knownGene` ON knownCanonical.transcript=knownGene.name";
	my $sth = $conn->prepare($q);
	$sth->execute();
	#while ($data = $conn->selectrow_hashref($q))
	while ($data = $sth->fetchrow_hashref())
	{
		if (!$utssh{$data->{"chrom"}}{$data->{"strand"}}{$data->{"transcript"}})
		{
			($data->{"strand"} == "+") ?
			$utssh{$data->{"chrom"}}{$data->{"strand"}}{$data->{"transcript"}} = {$data->{"chromStart"}} :
			$utssh{$data->{"chrom"}}{$data->{"strand"}}{$data->{"transcript"}} = {$data->{"chromEnd"}};
		}
		else
		{
			$oldStart = $utssh{$data->{"chrom"}}{$data->{"strand"}}{$data->{"transcript"}};
			$newStart = $data->{"chromStart"};
			$newEnd = $data->{"chromEnd"};
			if ($data->{"strand"} == "+")
			{
				$utssh{$data->{"chrom"}}{$data->{"strand"}}{$data->{"transcript"}} = $newStart if ($newStart < $oldStart);
			}
			elsif ($data->{"strand"} == "-")
			{
				$utssh{$data->{"chrom"}}{$data->{"strand"}}{$data->{"transcript"}} = $newEnd if ($newEnd > $oldStart);
			}
		}
	}
	&close_connection($conn);
	return (\%utssh);
}

sub get_exonh
{
	my $source = lc($_[0]);
	my $exonh; # Hash of hashes: {chr}->{strand}->{exon}->{start\tend} most efficient memory usage
	switch($source)
	{
		case /ensembl/
		{
			$exonh = &get_ensembl_exonh;
		}
		case /refseq/
		{
			$exonh = &get_refseq_exonh;
		}
		case /ucsc/
		{
			$exonh = &get_ucsc_exonh;
		}
	}
	return($exonh);
}

sub get_ensembl_exonh
{
	my %eexonh;
	my ($data,$chr,$name,$strand,$i,$k);
	my (@starts,@ends,@coords,@ses,@ustarts,@uends);
	my $conn = &open_connection($genome);
	my $q = "SELECT `chrom`, `strand`, `exonStarts`, `exonEnds`, `name2`  FROM `ensGene`";
	my $sth = $conn->prepare($q);
	$sth->execute();
	while ($data = $sth->fetchrow_hashref())
	{
		push(@{$eexonh{$data->{"chrom"}}{$data->{"strand"}}{$data->{"name2"}}{"start"}},split(",",$data->{"exonStarts"}));
		push(@{$eexonh{$data->{"chrom"}}{$data->{"strand"}}{$data->{"name2"}}{"end"}},split(",",$data->{"exonEnds"}));	
	}
	foreach $chr (keys(%eexonh))
	{
		foreach $strand (keys(%{$eexonh{$chr}}))
		{ 
			foreach $name (keys(%{$eexonh{$chr}{$strand}}))
			{
				@starts = @{$eexonh{$chr}{$strand}{$name}{"start"}};
				@ends = @{$eexonh{$chr}{$strand}{$name}{"end"}};
				for ($i=0; $i<@starts; $i++)
				{
					push(@coords,$starts[$i]."-".$ends[$i]);
				}
				my %u = &unique(@coords);
				foreach my $k (keys(%u))
				{
					@ses = split("-",$k);
					push(@ustarts,$ses[0]);
					push(@uends,$ses[1]);
				}
				@{$eexonh{$chr}{$strand}{$name}{"start"}} = @ustarts;
				@{$eexonh{$chr}{$strand}{$name}{"end"}} = @uends;
			}
		}
	}
	&close_connection($conn);
	return (\%eexonh);
}

sub get_refseq_exonh
{
	my %rexonh;
	my ($data,$chr,$name,$strand,$i,$k);
	my (@starts,@ends,@coords,@ses,@ustarts,@uends);
	my $conn = &open_connection($genome);
	my $q = "SELECT `geneName`, `chrom`, `strand`, `exonStarts`, `exonEnds` FROM `refFlat`";
	my $sth = $conn->prepare($q);
	$sth->execute();
	while ($data = $sth->fetchrow_hashref())
	{
		push(@{$rexonh{$data->{"chrom"}}{$data->{"strand"}}{$data->{"geneName"}}{"start"}},split(",",$data->{"exonStarts"}));
		push(@{$rexonh{$data->{"chrom"}}{$data->{"strand"}}{$data->{"geneName"}}{"end"}},split(",",$data->{"exonEnds"}));	
	}
	foreach $chr (keys(%rexonh))
	{
		foreach $strand (keys(%{$rexonh{$chr}}))
		{ 
			foreach $name (keys(%{$rexonh{$chr}{$strand}}))
			{
				@starts = @{$rexonh{$chr}{$strand}{$name}{"start"}};
				@ends = @{$rexonh{$chr}{$strand}{$name}{"end"}};
				for ($i=0; $i<@starts; $i++)
				{
					push(@coords,$starts[$i]."-".$ends[$i]);
				}
				my %u = &unique(@coords);
				foreach my $k (keys(%u))
				{
					@ses = split("-",$k);
					push(@ustarts,$ses[0]);
					push(@uends,$ses[1]);
				}
				@{$rexonh{$chr}{$strand}{$name}{"start"}} = @ustarts;
				@{$rexonh{$chr}{$strand}{$name}{"end"}} = @uends;
			}
		}
	}
	&close_connection($conn);
	return (\%rexonh);
}

sub get_ucsc_exonh
{
	my %uexonh;
	my ($data,$chr,$name,$strand,$i,$k);
	my (@starts,@ends,@coords,@ses,@ustarts,@uends);
	my $conn = &open_connection($genome);
	my $q = "SELECT `name`, knownGene.chrom, knownGene.strand, knownGene.exonStarts, knownGene.exonEnds ".
			"FROM `knownGene` INNER JOIN `knownCanonical` ON knownGene.name=knownCanonical.transcript";
	my $sth = $conn->prepare($q);
	$sth->execute();
	while ($data = $sth->fetchrow_hashref())
	{
		push(@{$uexonh{$data->{"chrom"}}{$data->{"strand"}}{$data->{"name"}}{"start"}},split(",",$data->{"exonStarts"}));
		push(@{$uexonh{$data->{"chrom"}}{$data->{"strand"}}{$data->{"name"}}{"end"}},split(",",$data->{"exonEnds"}));	
	}
	foreach $chr (keys(%uexonh))
	{
		foreach $strand (keys(%{$uexonh{$chr}}))
		{ 
			foreach $name (keys(%{$uexonh{$chr}{$strand}}))
			{
				@starts = @{$uexonh{$chr}{$strand}{$name}{"start"}};
				@ends = @{$uexonh{$chr}{$strand}{$name}{"end"}};
				for ($i=0; $i<@starts; $i++)
				{
					push(@coords,$starts[$i]."-".$ends[$i]);
				}
				my %u = &unique(@coords);
				foreach my $k (keys(%u))
				{
					@ses = split("-",$k);
					push(@ustarts,$ses[0]);
					push(@uends,$ses[1]);
				}
				@{$uexonh{$chr}{$strand}{$name}{"start"}} = @ustarts;
				@{$uexonh{$chr}{$strand}{$name}{"end"}} = @uends;
			}
		}
	}
	&close_connection($conn);
	return (\%uexonh);
}

sub stream_ensembl_exonh
{
	my %sexonh;
	my ($data,$chr,$name,$strand,$i,$k);
	my (@starts,@ends,@coords,@ses,@ustarts,@uends);
	my $conn = &open_connection($genome);
	my $q = "SELECT `chrom`, `strand`, `exonStarts`, `exonEnds`, `name2`  FROM `ensGene` ORDER BY  `name2`,`chrom`";
	my $sth = $conn->prepare($q);
	$sth->execute();
	while ($data = $sth->fetchrow_hashref())
	{
		push(@{$sexonh{$data->{"chrom"}}{$data->{"strand"}}{$data->{"name2"}}{"start"}},split(",",$data->{"exonStarts"}));
		push(@{$sexonh{$data->{"chrom"}}{$data->{"strand"}}{$data->{"name2"}}{"end"}},split(",",$data->{"exonEnds"}));
	}
	foreach $chr (keys(%sexonh))
	{
		foreach $strand (keys(%{$sexonh{$chr}}))
		{ 
			foreach $name (keys(%{$sexonh{$chr}{$strand}}))
			{
				@starts = @{$sexonh{$chr}{$strand}{$name}{"start"}};
				@ends = @{$sexonh{$chr}{$strand}{$name}{"end"}};
				for ($i=0; $i<@starts; $i++)
				{
					push(@coords,$starts[$i]."-".$ends[$i]);
				}
				my %u = &unique(@coords);
				foreach my $k (keys(%u))
				{
					@ses = split("-",$k);
					push(@ustarts,$ses[0]);
					push(@uends,$ses[1]);
				}
				@{$sexonh{$chr}{$strand}{$name}{"start"}} = @ustarts;
				@{$sexonh{$chr}{$strand}{$name}{"end"}} = @uends;
			}
		}
	}
	&close_connection($conn);
	return (\%sexonh);
}

#sub get_refseq_exonh
#{
	#my %rexonh;
	#my ($data,$chr,$name,$strand,$i,$k);
	#my (@starts,@ends,@coords,@ses,@ustarts,@uends);
	#my $conn = &open_connection($genome);
	#my $q = "SELECT `geneName`, `chrom`, `strand`, `exonStarts`, `exonEnds` FROM `refFlat`";
	#my $sth = $conn->prepare($q);
	#$sth->execute();
	#while ($data = $sth->fetchrow_hashref())
	#{
		#push(@{$rexonh{$data->{"chrom"}}{$data->{"strand"}}{$data->{"geneName"}}{"start"}},split(",",$data->{"exonStarts"}));
		#push(@{$rexonh{$data->{"chrom"}}{$data->{"strand"}}{$data->{"geneName"}}{"end"}},split(",",$data->{"exonEnds"}));	
	#}
	#foreach $chr (keys(%rexonh))
	#{
		#foreach $strand (keys(%{$rexonh{$chr}}))
		#{ 
			#foreach $name (keys(%{$rexonh{$chr}{$strand}}))
			#{
				#@starts = @{$rexonh{$chr}{$strand}{$name}{"start"}};
				#@ends = @{$rexonh{$chr}{$strand}{$name}{"end"}};
				#for ($i=0; $i<@starts; $i++)
				#{
					#push(@coords,$starts[$i]."-".$ends[$i]);
				#}
				#my %u = &unique(@coords);
				#foreach my $k (keys(%u))
				#{
					#@ses = split("-",$k);
					#push(@ustarts,$ses[0]);
					#push(@uends,$ses[1]);
				#}
				#@{$rexonh{$chr}{$strand}{$name}{"start"}} = @ustarts;
				#@{$rexonh{$chr}{$strand}{$name}{"end"}} = @uends;
			#}
		#}
	#}
	#&close_connection($conn);
	#return (\%rexonh);
#}

#sub get_ucsc_exonh
#{
	#my %uexonh;
	#my ($data,$chr,$name,$strand,$i,$k);
	#my (@starts,@ends,@coords,@ses,@ustarts,@uends);
	#my $conn = &open_connection($genome);
	#my $q = "SELECT `name`, knownGene.chrom, knownGene.strand, knownGene.exonStarts, knownGene.exonEnds ".
			#"FROM `knownGene` INNER JOIN `knownCanonical` ON knownGene.name=knownCanonical.transcript";
	#my $sth = $conn->prepare($q);
	#$sth->execute();
	#while ($data = $sth->fetchrow_hashref())
	#{
		#push(@{$uexonh{$data->{"chrom"}}{$data->{"strand"}}{$data->{"name"}}{"start"}},split(",",$data->{"exonStarts"}));
		#push(@{$uexonh{$data->{"chrom"}}{$data->{"strand"}}{$data->{"name"}}{"end"}},split(",",$data->{"exonEnds"}));	
	#}
	#foreach $chr (keys(%uexonh))
	#{
		#foreach $strand (keys(%{$uexonh{$chr}}))
		#{ 
			#foreach $name (keys(%{$uexonh{$chr}{$strand}}))
			#{
				#@starts = @{$uexonh{$chr}{$strand}{$name}{"start"}};
				#@ends = @{$uexonh{$chr}{$strand}{$name}{"end"}};
				#for ($i=0; $i<@starts; $i++)
				#{
					#push(@coords,$starts[$i]."-".$ends[$i]);
				#}
				#my %u = &unique(@coords);
				#foreach my $k (keys(%u))
				#{
					#@ses = split("-",$k);
					#push(@ustarts,$ses[0]);
					#push(@uends,$ses[1]);
				#}
				#@{$uexonh{$chr}{$strand}{$name}{"start"}} = @ustarts;
				#@{$uexonh{$chr}{$strand}{$name}{"end"}} = @uends;
			#}
		#}
	#}
	#&close_connection($conn);
	#return (\%uexonh);
#}

sub tssh_to_bed
{
	# This function takes a tssh HoH and converts it to a temporary bed file to be used with
	# the coverage functions of the bed tools
	# left and right MUST be in absolute numbers, e.g. for 1kb upstream, left must be -1000
	# while for 5kb downstream, right must be +5000
	# so for 5-10kb upstream, (left,right)=(-10000,-5000)
	my ($left,$right,$tssh) = @_;
	my ($chr,$strand,$id,$start,$end);
	my $rand = int(rand(1000000)) + 1000000; # Safety switch for very quick MySQL queries
	my $tmpbed = File::Spec->catfile($tmpdir,&now("machine")."_$rand.bed");
	my $fh = IO::File->new();
	$fh->open(">$tmpbed");
	foreach $chr (sort keys(%{$tssh}))
	{
		next if ($chr =~ m/chrM|rand/);
		foreach $strand (keys(%{$tssh->{$chr}}))
		{
			foreach $id (sort { $tssh->{$chr}->{$strand}->{$a} <=> $tssh->{$chr}->{$strand}->{$b} } keys(%{$tssh->{$chr}->{$strand}}))
			{
				if ($strand eq "+")
				{
					$start = $tssh->{$chr}->{$strand}->{$id} + $left;
					$end = $tssh->{$chr}->{$strand}->{$id} + $right - 1;
				}
				elsif ($strand eq "-")
				{
					$start = $tssh->{$chr}->{$strand}->{$id} - $right + 1;
					$end = $tssh->{$chr}->{$strand}->{$id} - $left;
				}
				print $fh "$chr\t$start\t$end\t$id\t0\t$strand";
			}
		}
	}
	$fh->close;
	return($tmpbed);
}

sub exonh_to_bed
{
	my $exonh = $_[0];
	my (@starts,@ends,@coords,@ses,@ustarts,@uends);
	my ($chr,$id,$strand,$i,$k);
	my $tmpbed = File::Spec->catfile($tmpdir,&now("machine").".bed");
	my $fh = IO::File->new();
	$fh->open(">$tmpbed");
	foreach $chr (sort keys(%{$exonh}))
	{
		next if ($chr =~ m/chrM|rand/);
		foreach $strand (keys(%{$exonh->{$chr}}))
		{
			foreach $id (keys(%{$exonh->{$chr}->{$strand}}))
			{
				@starts = @{$exonh->{$chr}->{$strand}->{$id}->{"start"}};
				@ends = @{$exonh->{$chr}->{$strand}->{$id}->{"end"}};
				for ($i=0; $i<@starts; $i++)
				{
					print $fh "$chr\t$starts[$i]\t$ends[$i]\t$id"."_"."$i\t0\t$strand\n";
				}
			}
		}
	}
	$fh->close;
	return($tmpbed);
}

sub bam_to_bed
{
	my @in = @_;
	my @out;
	foreach my $f (@in)
	{
		my $base = fileparse($f,'\.[^.]*');
		my $bed = File::Spec->catfile($tmpdir,$base.".bed");
		disp("Converting BAM file $f to BED file $bed...");
		`samtools view -b $f | bamToBed -i stdin > $bed`;
		push(@out,$bed);
	}
	return(@out);
}

sub write_bed_line
{
	my $fh = shift @_;
	my @line = @_;
	print $fh join("\t",@line),"\n";
}

sub open_connection
{   
	my $database = shift @_;
	my ($hostname,$conn);
	if (&check_existence($database))
	{
		$hostname = "localhost";
		$conn = DBI->connect("dbi:mysql:database=$database;host=$hostname;port=3306",$dbdata[0],$dbdata[1]);
	}
	else # Connect to the public MySQL host at UCSC
	{
		$hostname = REMOTE_HOST;
		$conn = DBI->connect("dbi:mysql:database=$database;host=$hostname;port=3306",REMOTE_USER);
	}
    return $conn;
}

sub close_connection
{ 
    my $conn = shift @_;
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
	&close_connection($conn);
	return($out);
}

sub calc_genome_coverage 
{ 
	my ($in,$mappable) = @_;
	my $mapped = 0;
	my $effective;
	
	open(GENCOV,$in);
	my $line = <GENCOV>; # Skip first line, 0 coverage
	my @cols = split(/\t/,$line);
	my $gsize = $cols[3]; # But grab the total genome size
	while ($line = <GENCOV>)
	{
		@cols = split(/\t/,$line);
		$mapped += $cols[2];
	}
	close(GENCOV);
	
	($mappable <= 1) ? ($effective = $mappable*$gsize) : ($effective = $mappable);
	return(sprintf("%.5f", $mapped/$effective));
}

sub calc_custom_coverage 
{ 
	my $in = shift @_;
	my $mapped = 0;
	
	open(GENCOV,$in);
	my $line = <GENCOV>; # Skip first line, 0 coverage
	my @cols = split(/\t/,$line);
	my $size = $cols[3]; # But grab the total custom space size
	while ($line = <GENCOV>)
	{
		@cols = split(/\t/,$line);
		$mapped += $cols[2];
	}
	close(GENCOV);
	
	return(sprintf("%.5f", $mapped/$size));
}

sub write_R
{
	my ($x,$y,$t,$n) = @_;
	
	my $script = 
		"$t(\"$n\")\n".	
		"par(mai=c(0.8,0.8,0.5,0.5))\n".
		"plot($x,$y,col=\"skyblue3\",xaxt=\"n\",yaxt=\"n\",type=\"l\",lwd=2,xlab=\"\",ylab=\"\")\n".
		"points($x,$y,col=\"red3\",pch=20,cex=1.3)\n".
		"axis(1,at=$x,labels=$x/1e+6,cex.axis=0.8,font=2,padj=-1,tcl=-0.3)\n".
		"axis(2,cex.axis=0.8,padj=0.5)\n".
		"title(main=\"Coverage plot\",cex.main=1)\n".
		"mtext(side=1,text=expression(bold(paste(\"Number of reads (x 10\"^\"6\",\")\",sep=\"\"))),line=2,cex=0.9,font=2)\n".
		"mtext(side=2,text=\"Coverage saturation (%)\",line=2,cex=0.9,font=2)\n".
		"grid()\n".
		"dev.off()\n";
	return($script); 
}

# OTHER AUX FUNCTIONS

# Process inputs
sub check_inputs
{
    my $stop;
    GetOptions("input|i=s{,}" => \@input,
    		   "genome|g=s" => \$genome,
    		   "resolution|r=i" => \$resolution,
    		   "step|t=i" => \$step,
    		   "outpath|o=s" => \$opath,
    		   "figtype|f=s" => \$figtype,
    		   "ncore|n=i" => \$ncore,
    		   "type|y=s" => \$type,
    		   "bam|b" => \$bam,
    		   "silent|s" => \$silent,
    		   "help|h" => \$help);
    # Check if the required arguments are set
    if ($help)
    {
    	&program_usage;
    	exit;
    }
    $stop .= "--- Please specify input BED files ---\n" if (!@input);
    $stop .= "--- Please specify the coverage type to be calculated ---\n" if (!$type);
    #$stop .= "--- Please specify the genome under investigation ---\n" if (!$genome);
    if ($type)
    {
		$type = lc($type);
		switch ($type)
		{
			case /genome/
			{
				if (!$mappable && $genome)
				{
					switch ($genome)
					{
						case /hg18|hg19/ { $mappable = 2700000000 }
						case /mm8|mm9/ { $mappable = 1865500000 }
						case /ce/ { $mappable = 90000000 }
						case /dm3/ { $mappable = 120000000 }
						else 
						{ 
							$stop .= "--- Selected genome not supported, please provide a fraction of the mappable genome ---\n"; 
						}
					}
				}
				else
				{
					$stop .= "--- Mappable genome fraction must be a number between 0 and 1 (not 0) ---\n"
						if ($mappable <=0 || $mappable > 1);
				}
			}
			case /features/
			{
				$stop .= "--- Database connection credentials should be a vector of 2 elements ---\n"
					if (!@dbdata || scalar(@dbdata) != 2);
				$stop .= "--- Selected genome $genome not supported! ---\n"
					if (!$genome || ($genome ne "hg18" && $genome ne "hg19" && $genome ne "mm9" && $genome ne "dm3"));
			}
			case /custom/
			{
				$stop .= "--- Please provide a features BED file to be used for coverage calculation ---\n"
					if (!$genome);
			}
		}
	}
	if ($stop)
    {
            print "\n$stop\n";
            print "Type perl $scriptname --help for help in usage.\n\n";
            exit;
    }
    if ($ncore)
    {
		eval "&try_module(\"Parallel::Loops\")";
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
    if ($resolution && $resolution <= 1)
    {
		disp("The resolution must be a positive integer greater than 1. Using default (10)...");
		$resolution = 10;
	}
	if ($step && $step <= 1)
    {
		disp("The step must be a positive integer greater than 1. Using default (500000)...");
		$step = 5e+5;
	}
	if ($resolution && $step)
	{
		disp("Only one of --step or --resolution can be specified! Using step (500000)...");
		$resolution = 0;
		$step = 5e+5;
	}
	if (!$resolution && !$step)
	{
		disp("No --step or --resolution specified. Using default step (500000)...");
		$resolution = 0;
		$step = 5e+5;
	}
	if ($figtype ne "png" && $figtype ne "pdf")
	{
		disp("--type must be one of png or pdf! Using png...");
		$figtype = "png";
	}
	if (!$opath)
	{
		use Cwd;
		disp("No output path specified, the plot(s) will be generated in the same folder as the input file(s)...");
		$opath = getcwd;
	}
}

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

sub count_lines
{
	open(IN,$_[0]) or die "\nThe file $_[0] does not exist!\n\n";
	my $totlines=0;
	$totlines += tr/\n/\n/ while sysread(IN,$_,2**16);
	close(IN);
	return $totlines;
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

sub cleanup 
{
	for my $temp (values(%tcov)) 
	{
		unlink($temp) if (-e $temp);
	}
}

sub catch_cleanup 
{
	print STDERR "Catching ctrl-C, cleaning temporary files!";
	cleanup;
	die;
}

sub disp
{
	print "\n@_" if (!$silent);
}

sub program_usage 
{
	# The look sucks here but it is actually good in the command line
	my $usagetext = << "END";
	
$scriptname
A Perl script to create coverage saturation plots for deep sequencing experiments

Author: Panagiotis Moulos (moulos\@fleming.gr)

Main usage
$scriptname --input file(s) --genome genome [OPTIONS]

--- Required ---
  --input|i	Input BED file (BED supported so far).
  --genome|g	The genome under investigation. Can be one of hg18, hg19, mm8, mm9 OR
			a genome file (tab delimited chromosomal lengths) for another prefered
			organism.
  --mappable|m	The mappable genome fraction. It must be provided if --genome is not human
			(hg18, hg19), mouse (mm8, mm9), c.elegans (ce) or fruitfly (dm). It must be a 
			number between 0 and 1 which denotes fraction. The effective genome size for 
			human and mouse is derived from the MACS software manual.
  --type|y		The coverage type to be calculated. Can be one of "genome" for
					total genome coverage, "features" for coverage to be calculated
					over a set of features (to add more...) or "custom" where genome
					is a custom bed file with custom features (e.g. lincRNA locations).
--- Optional ---
  --resolution|r		The coverage graph resolution (how many points should be plotted),
			Default: 10. Note that only one of -r or -s is acceptable.
  --step|t		The coverage generation step, default: 5000000.
  --output|o		Path to output.
  --ncore|n		The number of cores to use on a multiprocessor machine (default:3, max:6).
  --silent|s		Use this option if you want to turn informative 
			messages off.
  --help|h		Display this help text.
	
This program builds the species and genes tables in the KUPKB_Vis database.

END
	print $usagetext;
	exit;
}

# DEBUG
sub printHash
{
	my %h = %{$_[0]};
	print "\n-------------------- Begin hash contents --------------------\n";
	foreach my $k (keys(%h))
	{
		print "$k\t$h{$k}\n";
	}
	print "-------------------- End hash contents --------------------\n";
}
