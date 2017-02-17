#!/usr/bin/perl

# chipseqPQC.pl
# A Perl script/wrapper that uses several tools to perform post-alignment quality control
# for ChIP-Seq experiments given a mapping file between replicates and class names, and a
# YAML file which provides the parameters for every step in a systematic way.
#
# Author      : Panagiotis Moulos (moulos@fleming.gr)
# Created     : 27 - 07 - 2012 (dd - mm - yyyy)
# Last Update : 08 - 08 - 2012 (dd - mm - yyyy)
# Version     : 1.0

# FIXME: Move covsat.pl to /opt and remove perl from the command in iterative_coverage
# FIXME: Change TRUE and FALSE to "yes" and "no" in the parameters file because it is seen
#		 as text and not recognized as boolean, or write a function to do this
# Added local filemap hash variables in subroutines. A little more memory costly but allows
# for later input file subsetting by subsetting the filemap.

use strict;
use IO::File;
use Getopt::Long;
use File::Temp;
use File::Copy;
use File::Spec;
use File::Basename;
use File::Path qw(make_path remove_tree);
use List::Util qw(first max maxstr min minstr shuffle sum);
use Switch;

use Data::Dumper;

use constant TRUE => 1;
use constant FALSE => 0;
			 
# Make sure output is unbuffered
select(STDOUT);
$|=1;

# On Ctrl-C or die, do cleanup
$SIG{INT} = \&catch_cleanup;
#$SIG{__DIE__} = 'catch_die'; # Does not work as expected...

# Set defaults
our $scriptname = "chipseqPQC.pl";
our $mapfile;
our $opath;		 # An output to write all the results inside
our $qcname;	 # A name for the quality control project
our $paramfile;  # YAML parameters file with several options
our $log = 0;		# Keep log?
our $silent = 0; # Display verbose messages
our $help = 0;   # Help?

# Advertise
&advertise;

# Check inputs
&check_inputs;

# Global variables
our $phref;
our $tmpdir = File::Temp->newdir();
our $logfilehandle = &open_log_file if ($log); # Log file if requested

# Check for the presence of YAML, required!!!
&try_module("YAML");
&try_module("Tie::IxHash::Easy");

# Are we in Windows? For future use...
my $winOS = FALSE;
$winOS = TRUE if ($^O =~ /MSWin/);

# Record progress...
my $date = &now;
disp("$date - Started...\n");

# Read the parameters file or load the defaults if not given
($paramfile) ? ($phref = &read_param_file) : ($phref = &load_default_params);

# General variables
my ($RSfile,$ROfile);

# Define the genome file
our $genome = &get_genome;

# Read map file
our ($filemap,$rephash,$ctrlhash) = &read_map_file();

# Determine if we have combined samples in the given files
my $hasCombined = FALSE;
my $hasReplicates = FALSE;
$hasCombined = TRUE if (scalar keys(%$filemap) > 1);
$hasReplicates = TRUE if (scalar keys(%$rephash) != sum(values(%$rephash)));

# If the replicates are not combined, we have to combine them, temporarily
my ($joinfiles,$c,$ss,$cc,$ts,$tc,$dirsam,$dircon);
my (@osamples,@ocontrols);
if (!$hasCombined && $hasReplicates)
{
	$joinfiles = &combine_replicates;
	foreach $c (keys(%$joinfiles))
	{
		#$filemap->{"combined"}->{$c} = $joinfiles->{$c};
		push(@{$filemap->{"combined"}->{$c}},$joinfiles->{$c});
	}
	# Take care also of the control hash...
	foreach $c (keys(%$ctrlhash))
	{
		my (@samples,@controls);
		@osamples = keys(%{$ctrlhash->{$c}});
		@ocontrols = values(%{$ctrlhash->{$c}});

		$dirsam = dirname($osamples[0]);
		$dircon = dirname($ocontrols[0]);

		foreach $ss (@osamples) { push(@samples,basename($ss)); }
		foreach $cc (@ocontrols) { push(@controls,basename($cc)); }
		
		my %u = &unique(@ocontrols);
		if (sum(values(%u)) == scalar(keys(%u)) && scalar(keys(%u)) > 1)
		{
			#$ctrlhash->{$c}->{join(".",keys(%{$ctrlhash->{$c}}))."combined"} = join(".",values(%{$ctrlhash->{$c}}))."combined";
			$ts = File::Spec->catfile($dirsam,join(".",@samples).".combined");
			$tc = File::Spec->catfile($dircon,join(".",@controls).".combined");
			$ctrlhash->{$c}->{$ts} = $tc;
		}
		elsif (sum(values(%u)) != scalar(keys(%u)) && scalar(keys(%u)) == 1)
		{
			$ts = File::Spec->catfile($dirsam,join(".",@samples).".combined");
			$ctrlhash->{$c}->{$ts} = join(".",keys(%u));
		}
		elsif (sum(values(%u)) == scalar(keys(%u)) && scalar(keys(%u)) == 1)
		{
			# Nothing, already there
		}
		elsif (sum(values(%u)) != scalar(keys(%u)) && scalar(keys(%u)) > 1)
		{
			#$ctrlhash->{$c}->{join(".",keys(%{$ctrlhash->{$c}}))} = join(".",keys(%u));
			$ctrlhash->{$c}->{$ts} = join(".",keys(%u));
		}
	}
}

# Step 1: number of unique reads and image...
my $urfile = &count_reads;
($RSfile,$ROfile) = &writeR_color_unreads($urfile);
#`R CMD BATCH --vanilla $RSfile $ROfile`;
`Rscript --vanilla $RSfile`;

# Step 2: several coverages...
my @covfiles = &simple_coverage;
($RSfile,$ROfile) = &writeR_cov_bars(@covfiles);
#`R CMD BATCH --vanilla $RSfile $ROfile`;
`Rscript --vanilla $RSfile`;

# Step 3: iterative coverages
&iterative_coverage;

# Step 4: MACS saturation plots
my @satfiles = &macs_saturation;
($RSfile,$ROfile) = &writeR_macs_sat(@satfiles);
#`R CMD BATCH --vanilla $RSfile $ROfile`;
`Rscript --vanilla $RSfile`;

## Step 5: PCA and others from htSeqTools
disp("\n--- RUNNING QUALITY CONTROL USING htSeqTools... THIS IS A LONG PROCEDURE... ---\n");
my ($classes,$replicates) = &get_classes_and_replicates("both");
my $classRvec = make_R_vector("character",@$classes);
my $repRvec = make_R_vector("character",@$replicates);
($RSfile,$ROfile) = &writeR_htseqtools_qc($classRvec,$repRvec);
#`R CMD BATCH --vanilla $RSfile $ROfile`;
`Rscript --vanilla $RSfile`;

# Step 6: Count correlations, heatmap, hilbert and genomic region wiggle coverage
disp("\n--- RUNNING QUALITY CONTROL USING COUNT AND COVERAGE CORRELATIONS AND HILBERT CURVES... THIS IS A LONG PROCEDURE... ---\n");
($RSfile,$ROfile) = &writeR_cor_qc($classRvec,$repRvec);
#`R CMD BATCH --vanilla $RSfile $ROfile`;
`Rscript --vanilla $RSfile`;

# Step 7: Create peak correlations
my $justfiles = &just_macs;
($RSfile,$ROfile) = &writeR_cor_peaks($justfiles);
#`R CMD BATCH --vanilla $RSfile $ROfile`;
`Rscript --vanilla $RSfile`;

## Step 8: IDR analysis


$date = &now;
disp("\n$date - Finished!\n");

# Close the log file if requested
&close_log_file if ($log);


### MAIN SUBROUTINES ###

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

sub read_map_file
{
	use Tie::IxHash::Easy;
	my $line;
	my @cols;
	my (%filemap,%rephash,%ctrlhash);
	tie %ctrlhash, "Tie::IxHash::Easy";
	my $mfh = IO::File->new();
	$mfh->open("< $mapfile") or die "\nCannot open file $mapfile\n\n";
	$line = <$mfh>; # Skip header
	while ($line = <$mfh>)
	{
		$line =~ s/\r|\n$//g;
		@cols = split(/\t/,$line);
		#$filemap{"type"}{"class"}="filename"; # Structure
		push(@{$filemap{$cols[2]}{$cols[1]}},File::Spec->rel2abs($cols[0]));
		$rephash{$cols[1]}++;
		$ctrlhash{$cols[1]}{File::Spec->rel2abs($cols[0])} = File::Spec->rel2abs($cols[3])
			if ($cols[3] !~ m/(\s*\-\s*)|(NA)|(\t)/);
	}
	$mfh->close;
	return(\%filemap,\%rephash,\%ctrlhash);
}

sub subset_map { }

sub get_classes_and_replicates
{
	my $what = $_[0];
	my %tclmap;
	($_[1]) ? (%tclmap = %{$_[1]}) : (%tclmap = %{$filemap->{"single"}});
	my ($class,$replicate);
	my (@classes,@replicates);
	foreach $class (keys(%tclmap))
	{
		foreach $replicate (@{$tclmap{$class}})
		{
			push(@classes,$class);
			push(@replicates,$replicate);
		}
	}
	if ($what eq "both")
	{
		return(\@classes,\@replicates);
	}
	elsif ($what eq "classes")
	{
		return(@classes);
	}
	elsif ($what eq "replicates")
	{
		return(@replicates);
	}
}

sub get_genome
{
	switch ($phref->{"coverage"}->{"genome"})
	{
		case "hg18" { $genome = File::Spec->catfile($phref->{"general"}->{"genomes"},"human.hg18.genome") }
		case "hg19" { $genome = File::Spec->catfile($phref->{"general"}->{"genomes"},"human.hg19.genome") }	
		case "mm8" { $genome = File::Spec->catfile($phref->{"general"}->{"genomes"},"mouse.mm8.genome") }
		case "mm9" { $genome = File::Spec->catfile($phref->{"general"}->{"genomes"},"mouse.mm9.genome") }
		else { $genome = $phref->{"coverage"}->{"genome"} }
	}
	return($genome);
}

sub combine_replicates
{
	my %tclmap;
	($_[0]) ? (%tclmap = %{$_[0]}) : (%tclmap = %{$filemap->{"single"}});
	my ($dir,$r,$class,$catcom,$joinOK);
	my (@oreplicates,@replicates);
	my %joinfiles;
	foreach $class (keys(%tclmap))
	{
		@oreplicates = @{$tclmap{$class}};
		if (scalar(@oreplicates) > 1)
		{
			# Strip the directory names
			@replicates = ();
			foreach $r (@oreplicates)
			{
				push(@replicates,basename($r));
			}
			disp("Combining all replicates for class $class...");
			$dir = dirname($oreplicates[0]);
			$joinfiles{$class} = File::Spec->catfile($dir,join(".",@replicates).".combined");
			$catcom = "cat ".join(" ",@oreplicates)." | sort -k1,1 -k2g,2 -u > ".$joinfiles{$class};
			$joinOK = system($catcom);
			die "\nWhooops, something went wrong with replicate combining for ".join(" and ",@oreplicates)."!\n\n"
				if ($joinOK);
		}
	}
	return(\%joinfiles);
}

sub count_reads
{
	my %tclmap;
	($_[0]) ? (%tclmap = %{$_[0]}) : (%tclmap = %{$filemap->{"single"}});
	my ($class,$replicate,$n);
	my @replicates;
	my $op = File::Spec->catdir($opath,"count_reads");
	make_path($op);
	my $urfile = File::Spec->catfile($op,"unique_reads.txt");
	my $fh = IO::File->new();
	$fh->open(">$urfile");
	foreach $class (keys(%tclmap))
	{
		foreach $replicate (@{$filemap->{"single"}->{$class}})
		{
			disp("Counting the number of uniquely mapped reads for input file ".basename($replicate)."...");
			$n = &count_lines($replicate);
			print $fh "$replicate\t$n\n";
		}
	}
	$fh->close;
	return($urfile);
}

sub simple_coverage
{
	my %localmap;
	($_[0]) ? (%localmap = %{$_[0]}) : (%localmap = %$filemap);
	my ($type,$class,$file,$covfile);
	my @covfiles;
	foreach $type (keys(%localmap))
	{
		foreach $class (keys(%{$localmap{$type}}))
		{
			foreach $file (@{$localmap{$type}{$class}}) 
			{
				disp("Calculating genomic coverage for input file ".basename($file)."...");
				$covfile = File::Spec->catfile($tmpdir,basename($file).".scov");
				push(@covfiles,$covfile);
				`genomeCoverageBed -i $file -g $genome -max 50 > $covfile`;
			}
		}
	}
	return(@covfiles);
}

sub iterative_coverage
{
	my %localmap;
	($_[0]) ? (%localmap = %{$_[0]}) : (%localmap = %$filemap);
	my ($type,$class,$file,$curimgfile,$curtxtfile,$newimgfile,$newtxtfile,$itOK);
	my $op = File::Spec->catdir($opath,"iterative_coverage");
	make_path($op);
	foreach $type (keys(%localmap))
	{
		foreach $class (keys(%{$localmap{$type}}))
		{
			foreach $file (@{$localmap{$type}{$class}})
			{
				disp("\n--- ITERATIVE COVERAGE PLOT GENERATION FOR FILE ".basename($file)."... ---\n");
				if ($phref->{"coverage"}->{"step"} && !$phref->{"coverage"}->{"resolution"})
				{
					#`perl covsat.pl --input $file --genome $phref->{"coverage"}->{"genome"} --step $phref->{"coverage"}->{"step"} --type $phref->{"general"}->{"imageformat"} --ncore $phref->{"coverage"}->{"ncore"} --opath $op`;
					$itOK = system("perl $phref->{\"coverage\"}->{\"path\"} --input $file --genome $phref->{\"coverage\"}->{\"genome\"} --step $phref->{\"coverage\"}->{\"step\"} --type $phref->{\"general\"}->{\"imageformat\"} --ncore $phref->{\"coverage\"}->{\"ncore\"} --outpath $op");
				}
				elsif (!$phref->{"coverage"}->{"step"} && $phref->{"coverage"}->{"resolution"})
				{
					#`perl covsat.pl --input $file --genome $phref->{"coverage"}->{"genome"} --resolution $phref->{"coverage"}->{"resolution"} --type $phref->{"general"}->{"imageformat"} --ncore $phref->{"coverage"}->{"ncore"} --opath $op`;
					$itOK = system("perl $phref->{\"coverage\"}->{\"path\"} --input $file --genome $phref->{\"coverage\"}->{\"genome\"} --resolution $phref->{\"coverage\"}->{\"resolution\"} --type $phref->{\"general\"}->{\"imageformat\"} --ncore $phref->{\"coverage\"}->{\"ncore\"} --outpath $op");
				}
				else
				{
					#`perl covsat.pl --input $file --genome $phref->{"coverage"}->{"genome"} --step $phref->{"coverage"}->{"step"} --resolution $phref->{"coverage"}->{"resolution"} --type $phref->{"general"}->{"imageformat"} --ncore $phref->{"coverage"}->{"ncore"} --opath $op`;
					$itOK = system("perl $phref->{\"coverage\"}->{\"path\"} --input $file --genome $phref->{\"coverage\"}->{\"genome\"} --step $phref->{\"coverage\"}->{\"step\"} --resolution $phref->{\"coverage\"}->{\"resolution\"} --type $phref->{\"general\"}->{\"imageformat\"} --ncore $phref->{\"coverage\"}->{\"ncore\"} --outpath $op");
				}
				$curimgfile = File::Spec->catfile($opath,$file.".cov.",$phref->{"general"}->{"imageformat"});
				$curtxtfile = File::Spec->catfile($opath,$file.".cov.txt");
				$newimgfile = File::Spec->catfile($op,$file.".cov.",$phref->{"general"}->{"imageformat"});
				$newtxtfile = File::Spec->catfile($op,$file.".cov.txt");
				move($curimgfile,$newimgfile);
				move($curtxtfile,$newtxtfile);
			}
		}
	}
}

sub macs_saturation
{
	my %localctrlhash;
	($_[0]) ? (%localctrlhash = %{$_[0]}) : (%localctrlhash = %$ctrlhash);
	#my %ctrlhash = %{$_[0]};
	my ($class,$sample,$base,$newdiagfile,$newpeakfile,$cmd,$cmdOK);
	my @satfiles;
	my $op = File::Spec->catdir($opath,"macs_saturation");
	make_path($op);
	foreach $class (keys(%localctrlhash))
	{
		foreach $sample (keys(%{$localctrlhash{$class}}))
		{
			disp("\n--- MACS PEAKS SATURATION PLOT FOR FILE ".basename($sample)."... ---\n");
			$base = fileparse($sample,'\.[^.]*');
			$cmd = "$phref->{\"peak\"}->{\"macs\"}->{\"path\"} --treatment=$sample ".
				"--control=$localctrlhash{$class}{$sample} --name=$base ".
				"--mfold=$phref->{\"peak\"}->{\"macs\"}->{\"mfold\"} ".
				"--bw=$phref->{\"peak\"}->{\"macs\"}->{\"bw\"} ".
				"--gsize=$phref->{\"peak\"}->{\"macs\"}->{\"gsize\"} ".
				"--pvalue=$phref->{\"peak\"}->{\"macs\"}->{\"pvalue\"} ".
				"--fe-min=$phref->{\"peak\"}->{\"macs\"}->{\"fe-min\"} ".
				"--fe-max=$phref->{\"peak\"}->{\"macs\"}->{\"fe-max\"} ".
				"--fe-step=$phref->{\"peak\"}->{\"macs\"}->{\"fe-step\"} --diag";
			#$cmdOK = system($cmd);
			$cmdOK = `$cmd`;
			print $logfilehandle "\n$cmdOK" if ($log);
			die "$!" if $?;
			#if (!$cmdOK)
			#{
				$newdiagfile = File::Spec->catfile($op,$base."_diag.xls");
				$newpeakfile = File::Spec->catfile($op,$base."_peaks.xls");
				move($base."_diag.xls",$newdiagfile);
				move($base."_peaks.xls",$newpeakfile);
				unlink($base."_peaks.bed");
				unlink($base."_negative_peaks.xls");
				unlink($base."_summits.bed");
				unlink($base."_model.r");
				push(@satfiles,$newdiagfile);
			#}
			#else
			#{
			#	disp ("Whooops! Something went wrong with MACS peak saturation generation...");
			#}
		}
	}
	return(@satfiles);
}

sub just_macs
{
	my %localctrlhash;
	($_[0]) ? (%localctrlhash = %{$_[0]}) : (%localctrlhash = %$ctrlhash);
	#my %ctrlhash = %{$_[0]};
	my ($class,$sample,$base,$newbedfile,$newpeakfile,$cmd,$cmdOK);
	my %macsfiles;
	my $op = File::Spec->catdir($opath,"peak_correlation");
	make_path($op);
	foreach $class (keys(%localctrlhash))
	{
		foreach $sample (keys(%{$localctrlhash{$class}}))
		{
			if ($sample !~ m/combined$/)
			{
				disp("\n--- MACS DEFAULT PEAK CALLING FOR FILE ".basename($sample)."... ---\n");
				$base = fileparse($sample,'\.[^.]*?');
				$cmd = "$phref->{\"peak\"}->{\"macs\"}->{\"path\"} --treatment=$sample ".
					"--control=$localctrlhash{$class}{$sample} --name=$base ".
					"--mfold=$phref->{\"peak\"}->{\"macs\"}->{\"mfold\"} ".
					"--bw=$phref->{\"peak\"}->{\"macs\"}->{\"bw\"} ".
					"--gsize=$phref->{\"peak\"}->{\"macs\"}->{\"gsize\"} ".
					"--pvalue=$phref->{\"peak\"}->{\"macs\"}->{\"pvalue\"}";
				#$cmdOK = system($cmd);
				$cmdOK = `$cmd`;
				print $logfilehandle "\n$cmdOK" if ($log);
				die "$!" if $?;
				#if (!$cmdOK)
				#{
					$newbedfile = File::Spec->catfile($op,$base."_peaks.bed");
					$newpeakfile = File::Spec->catfile($op,$base."_peaks.xls");
					move($base."_peaks.bed",$newbedfile);
					move($base."_peaks.xls",$newpeakfile);
					unlink($base."_negative_peaks.xls");
					unlink($base."_summits.bed");
					unlink($base."_model.r");
					$macsfiles{$newpeakfile}{"control"} = $localctrlhash{$class}{$sample};
					$macsfiles{$newpeakfile}{"treatment"} = $sample;
				#}
				#else
				#{
				#	disp ("Whooops! Something went wrong with MACS...");
				#}
			}
		}
	}
	return(\%macsfiles);
}

sub writeR_color_unreads
{
	my $f = $_[0];
	my $op = File::Spec->catdir($opath,"count_reads");
	my $o = File::Spec->catfile($op,$qcname.".unique_reads.".$phref->{"general"}->{"imageformat"});
	my $Rscript = 
		"source(file.path(\"$phref->{\"general\"}->{\"rscript\"}\",\"severalQC.R\"))\n".	
		"color.unreads(\"$f\",output=\"$phref->{\"general\"}->{\"imageformat\"}\",fil=\"$o\")\n";
	my $RSfile = File::Spec->catfile($tmpdir,"$qcname.reads.R");
	my $ROfile = File::Spec->catfile($tmpdir,"$qcname.reads.Rout");
	my $fh = IO::File->new();
	$fh->open(">$RSfile");
	print $fh $Rscript;
	$fh->close;
	return($RSfile,$ROfile); 
}

sub writeR_cov_bars
{
	my @f = @_;
	my ($cf,$bf,$imgen,$imchr);
	my $op = File::Spec->catdir($opath,"simple_coverage");
	make_path($op);
	my $Rscript = 
		"source(file.path(\"$phref->{\"general\"}{\"rscript\"}\",\"severalQC.R\"))\n";
	foreach $cf (@f)
	{
		$bf = basename($cf);
		$imgen = File::Spec->catfile($op,$bf.".gen.".$phref->{"general"}->{"imageformat"});
		$imchr = File::Spec->catfile($op,$bf.".chr.".$phref->{"general"}->{"imageformat"});
		#print "\n\n$imgen\n$imchr\n\n";
		$Rscript .= 
			"cov.bar(\"$cf\",what=\"genome\",output=\"$phref->{\"general\"}->{\"imageformat\"}\",fil=\"$imgen\")\n".
			"cov.bar(\"$cf\",what=\"chromosome\",output=\"$phref->{\"general\"}->{\"imageformat\"}\",fil=\"$imchr\")\n";
	}
	my $RSfile = File::Spec->catfile($tmpdir,"$qcname.scov.R");
	my $ROfile = File::Spec->catfile($tmpdir,"$qcname.scov.Rout");
	my $fh = IO::File->new();
	$fh->open(">$RSfile");
	print $fh $Rscript;
	$fh->close;
	return($RSfile,$ROfile); 
}

sub writeR_macs_sat
{
	my @f = @_;
	my ($cf,$img,$base,$dir);
	my $Rscript = 
		"source(file.path(\"$phref->{\"general\"}{\"rscript\"}\",\"severalQC.R\"))\n";
	foreach $cf (@f)
	{
		($base,$dir) = fileparse($cf,'\.[^.]*');
		$base =~ s/_diag//g;
		$img = File::Spec->catfile($dir,$base.".peaksat.".$phref->{"general"}{"imageformat"});
		$Rscript .= 
			"plot.sat.macs(\"$cf\",output=\"$phref->{\"general\"}{\"imageformat\"}\",fil=\"$img\")\n";
	}
	my $RSfile = File::Spec->catfile($tmpdir,"$qcname.peaksat.R");
	my $ROfile = File::Spec->catfile($tmpdir,"$qcname.peaksat.Rout");
	my $fh = IO::File->new();
	$fh->open(">$RSfile");
	print $fh $Rscript;
	$fh->close;
	return($RSfile,$ROfile); 
}

sub writeR_htseqtools_qc
{
	my ($cls,$reps) = @_;
	my @qctypes;
	push(@qctypes,"mds") if $phref->{"htseqtools"}->{"mds"};
	push(@qctypes,"ssd") if $phref->{"htseqtools"}->{"ssd"};
	push(@qctypes,"gini") if $phref->{"htseqtools"}->{"gini"};
	my $qctypesvec = &make_R_vector("character",@qctypes);
	my $op = File::Spec->catdir($opath,"htseqtools_qc");
	make_path($op);
	my $genimg = File::Spec->catfile($op,$qcname.".htseqtools.".$phref->{"general"}{"imageformat"});
	my $Rscript = 
		"source(file.path(\"$phref->{\"general\"}{\"rscript\"}\",\"severalQC.R\"))\n".
		"htseqtools.qc($reps,$cls,qctype=$qctypesvec,n.cores=$phref->{\"htseqtools\"}->{\"ncore\"},".
		"output=\"$phref->{\"general\"}{\"imageformat\"}\",fil=\"$genimg\")\n";
	my $RSfile = File::Spec->catfile($tmpdir,"$qcname.htseqtools.R");
	my $ROfile = File::Spec->catfile($tmpdir,"$qcname.htseqtools.Rout");
	my $fh = IO::File->new();
	$fh->open(">$RSfile");
	print $fh $Rscript;
	$fh->close;
	return($RSfile,$ROfile);
}

sub writeR_cor_qc
{
	my ($cls,$reps) = @_;
	my @qctypes;
	push(@qctypes,"counts") if $phref->{"corqc"}->{"counts"};
	push(@qctypes,"pca") if $phref->{"corqc"}->{"pca"};
	push(@qctypes,"mds") if $phref->{"corqc"}->{"mds"};
	push(@qctypes,"hilbert") if $phref->{"corqc"}->{"hilbert"};
	my $qctypesvec = &make_R_vector("character",@qctypes);
	my $op = File::Spec->catdir($opath,"cor_qc");
	make_path($op);
	my $genimg = File::Spec->catfile($op,$qcname.".corqc.".$phref->{"general"}{"imageformat"});
	my $Rscript = 
		"source(file.path(\"$phref->{\"general\"}{\"rscript\"}\",\"severalQC.R\"))\n".
		"cor.qc($reps,$cls,cor.type=$qctypesvec,org=\"$phref->{\"coverage\"}->{\"genome\"}\",".
		"output=\"$phref->{\"general\"}{\"imageformat\"}\",fil=\"$genimg\")\n";
	my $RSfile = File::Spec->catfile($tmpdir,"$qcname.corqc.R");
	my $ROfile = File::Spec->catfile($tmpdir,"$qcname.corqc.Rout");
	my $fh = IO::File->new();
	$fh->open(">$RSfile");
	print $fh $Rscript;
	$fh->close;
	return($RSfile,$ROfile);
}

sub writeR_cor_peaks
{
	my $mapfile = &create_cor_peaks_map_file($_[0]);
	my @plottypes;
	push(@plottypes,"circle") if $phref->{"corpeaks"}->{"plot"}->{"circle"};
	push(@plottypes,"ellipse") if $phref->{"corpeaks"}->{"plot"}->{"ellipse"};
	push(@plottypes,"heatmap") if $phref->{"corpeaks"}->{"plot"}->{"heatmap"};
	push(@plottypes,"cgram.shade") if $phref->{"corpeaks"}->{"plot"}->{"cgramshade"};
	push(@plottypes,"cgram.pts") if $phref->{"corpeaks"}->{"plot"}->{"cgrampts"};
	push(@plottypes,"simple") if $phref->{"corpeaks"}->{"plot"}->{"simple"};
	push(@plottypes,"image") if $phref->{"corpeaks"}->{"plot"}->{"image"};
	my $plottypesvec = &make_R_vector("character",@plottypes);
	my $op = File::Spec->catdir($opath,"peak_correlation");
	#make_path($op); # Exists
	my $genimg = File::Spec->catfile($op,$qcname.".corpeaks.".$phref->{"general"}{"imageformat"});
	my $Rscript = 
		"source(file.path(\"$phref->{\"general\"}{\"rscript\"}\",\"correlatePeaks.R\"))\n".
		"correlatePeaks(\"$mapfile\",avg.win=$phref->{\"corpeaks\"}{\"avgwin\"},".
		"fdr.cut=$phref->{\"corpeaks\"}{\"fdrcut\"},type=$plottypesvec,".
		"output=\"$phref->{\"general\"}{\"imageformat\"}\",fil=\"$genimg\")\n";
	my $RSfile = File::Spec->catfile($tmpdir,"$qcname.corpeaks.R");
	my $ROfile = File::Spec->catfile($tmpdir,"$qcname.corpeaks.Rout");
	my $fh = IO::File->new();
	$fh->open(">$RSfile");
	print $fh $Rscript;
	$fh->close;
	return($RSfile,$ROfile);
}

sub create_cor_peaks_map_file
{
	my %filemap = %{$_[0]};
	my @files = keys(%filemap);
	my $cp = dirname($files[0]);
	my $mapfile = File::Spec->catfile($cp,"mapfile.txt");
	my $fh = IO::File->new();
	$fh->open(">$mapfile");
	print $fh "filename\ttreatment\tcontrol\n";
	foreach my $f (@files)
	{
		print $fh "$f\t$filemap{$f}{\"treatment\"}\t$filemap{$f}{\"control\"}\n";
	}
	$fh->close;
	return($mapfile);
}

sub make_R_vector
{
	my $type = shift @_;
	my @input = @_;
	my $vec;
	if ($type eq "character") 
	{	
		$vec = "c(\"";
		$vec .= join("\",\"",@input);
		$vec .= "\")";
	}
	elsif ($type eq "numeric") 
	{ 
		$vec = "c(";
		$vec .= join(",",@input);
		$vec .= ")";
	}
	return($vec);
}

### AUXILIARY SUBROUTINES ###

# Process inputs
sub check_inputs
{
    my $stop;
    GetOptions("mapfile|m=s" => \$mapfile,
			   "name|n=s" => \$qcname,
			   "param|p=s" => \$paramfile,
			   "log|l" => \$log,
			   "opath|o=s" => \$opath,
    		   "silent|s" => \$silent,
    		   "help|h" => \$help);
    # Check if the required arguments are set
    if ($help)
    {
    	&program_usage;
    	exit;
    }
    $stop .= "--- Please specify a sample map file ---\n" if (!$mapfile);
    if ($stop)
    {
            print "\n$stop\n";
            print "Type perl $scriptname --help for help in usage.\n\n";
            exit;
    }
    if (!$qcname)
    {
		print "\nQuality control project name not given! Default timestamp (YYYYMMDDHHMMSS) will be used...";
		$qcname = &now("machine");
	}
	if (!$opath)
	{
		use Cwd;
		$opath = File::Spec->catdir(getcwd,"$qcname");
		make_path($opath);
		print "\nNo output path specified, the results will be placed in a folder named like the project in the current path...";
	} else { make_path($opath); }
    disp("Parameter file not given! Will try to load defaults...") if (!$paramfile);
    
}

sub load_default_params
{
	my %h = (
				"general" => {					
					 "rscript" => "/opt/NGSTools/Rlocal",
					 "genomes" => "/opt/NGSTools/BEDTools/genomes",
					 "imageformat" => "png"
				},
				"coverage" => {
					"path" => "/media/HD4/Fleming/dev/covsat.pl",
					"genome" => "hg18",
					"mappable" => "",
					"resolution" => 0,
					"step" => 500000,
					"output" => "",
					"ncore" => 3
				},
				"htseqtools" => {
					"mds" => TRUE,
					"ssd" => TRUE,
					"gini" => TRUE
				},
				"corpeaks" => {
					"plot" => {
						"circle" => TRUE,
						"ellipse" => TRUE,
						"heatmap" => TRUE,
						"cgramshade" => TRUE,
						"cgrampts" => TRUE,
						"simple" => TRUE,
						"image" => TRUE
					},
					"avgwin" => 100,
					"fdrcut" => 1
				},
				"idr" => {
					"caller" => {
						"macs" => {
							"path" => "/usr/bin/macs2",
							"pvalue" => 1e-6
						}
					}
				},
				"peak" => {
					"macs" => {
						"path" => "/usr/bin/macs14",
						"format" => "",
						"petdist" => "",
						"gsize" => "",
						"tsize" => "",
						"bw" => "",
						"pvalue" => "10,30",
						"pvalue" => 1e-6,
						"nolamda" => "",
						"slocal" => "",
						"llocal" => "",
						"on-auto" => "",
						"nomodel" => "",
						"shiftsize" => "",
						"keep-dup" => "",
						"to-large" => "",
						"wig" => "",
						"bdg" => "",
						"single-profile" => "", 
						"space" => "",
						"call-subpeaks" => "",
						"verbose" => "",
						"diag" => "",
						"fe-min" => "",
						"fe-max" => "",
						"fe-step" => ""
					}
				}
		);
	return(\%h);
}

sub init_result_map
{
	my %h = (
		"count_reads" => {
			"file" => "",
			"image" => ""
			},
		"simple_coverage" => {
			"filename" => {
				"genome" => "",
				"chromosome" => ""
				}
			},
		"iterative_coverage" => {
			"filename" => {
				"file" => "",
				"image" => ""
				}
			},
		"macs_saturation" => {
			"filename" => {
				"file" => "",
				"image" => ""
				}
			},
		"htseqtools_qc" => {
			"mds" => "",
			"ssd" => "",
			"gini" => ""
			},
		"cor_qc" => {
			"counts" => "",
			"filename" => {
				"coverage" => {
					"stub_coverage_where" => ""
					},
				"hilbert" => {
					"stub_chromosome" => ""
					}
				}
			},
		#"diffbind" => {
		#	"stub" => ""
		#	},
		"idr" => {
			"stub" => ""
			}
	);
	return(\%h);
}

sub count_lines
{
	open(IN,$_[0]) or die "\nThe file $_[0] does not exist!\n\n";
	my $totlines = 0;
	$totlines += tr/\n/\n/ while sysread(IN,$_,2**16);
	close(IN);
	return $totlines;
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
	print STDERR "\nCatching ctrl-C, cleaning temporary files!";
	&cleanup;
	&close_log_file if ($log);
	die;
}

sub catch_die
{
	&cleanup;
	&close_log_file if ($log);
}

sub cleanup 
{
	remove_tree($tmpdir);
}

sub open_log_file
{
	my $logname = File::Spec->catfile($opath,"chipseqQC_".&now("machine").".log");
	my $cur = &now;
	my $lfh = IO::File->new();
	$lfh->open(">$logname");
	print "\n\nLog file $logname opened in $cur!\n";
	return($lfh);
}

sub close_log_file
{
	my $cur = &now;
	print "\nLog file closed in $cur!\n\n";
	$logfilehandle->close;
}

sub advertise
{
	use Term::ANSIColor;
	print color 'bold yellow on_blue';
	disp($scriptname);
	print color 'bold green';
	disp("Post alignment Quality Control for ChIP-Seq data... Copyright: Panagiotis Moulos (moulos\@fleming.gr)\n");
	print color 'reset';
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

sub program_usage 
{
	# The look sucks here but it is actually good in the command line
	my $usagetext = << "END";

Main usage
$scriptname --param parameter_file.yml [OPTIONS]

--- Required ---
  --mapfile|m		A file with mappings between samples, controls and
			classes. See example.
--- Optional ---
  --param|p		A YAML parameters file containing the locations of the
			sample map file (simple tab delimited, see exmaple), and MACS
			run parameters.
  --name|n		A name for the QC project. If not given, the default
			timestamp will be used.
  --opath|o		Output path for the QC results. If not given, a directory
			with the name --name will be used.
  --silent|s		Use this option if you want to turn informative 
  			messages off.
  --help|h		Display this help text.
	
This program builds the species and genes tables in the KUPKB_Vis database.

END
	print $usagetext;
	exit;
}

# Old implementation of line counting, dependent on OS
#$wccom = "wc -l ".join(" ",@replicates)."| awk 'BEGIN { FS = \" \" } ; { print $2\"\t\"$1 }' > ".$urfile;
#$wcOK = system($wccom);
#die "\nWhooops, something went wrong with line counting for files ".join(", ",@replicates)."!\n\n" if (!$wcOK);
# Old MACS command
#`macs14 --treatment=$sample --control=$localctrlhash{$sample} --name=$base \
#--mfold=$phref->{"peak"}->{"macs"}->{"mfold"} \
#--bw=$phref->{"peak"}->{"macs"}->{"bw"} \
#--gsize=$phref->{"peak"}->{"macs"}->{"gsize"} \
#--pvalue=$phref->{"peak"}->{"macs"}->{"pvalue"} \
#--fe-min=$phref->{"peak"}->{"macs"}->{"fe-min"} \
#--fe-max=$phref->{"peak"}->{"macs"}->{"fe-max"} \
#--fe-step=$phref->{"peak"}->{"macs"}->{"fe-step"} --diag`;
#$scriptname
#Post alignment Quality Control for ChIP-Seq data.
#Author : Panagiotis Moulos (moulos\@fleming.gr)
