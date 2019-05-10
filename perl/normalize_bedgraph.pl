#!/usr/bin/perl

use strict;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use File::Spec;
             
# Make sure output is unbuffered
select(STDOUT);
$|=1;

# Set defaults
our $scriptname = "normalize_bedgraph.pl";
our @bglist;
our @extnorm;
our $sumto;
our $exportfacs;
our $perlonly;
our @output;
our $prerun;
our $prerunlog;
our $ncores;
our $man;
our $log;
our $silent;
our $help;

&check_inputs;

if (!@output) {
	foreach my $f (@bglist) {
		push(@output,&create_output_file($f));
	}
}
elsif (&smatch("stdout",@output)) {
	@output = ();
}

my ($line,$out,$chr,$start,$end,$signal);
my (%wigsum,%normfactor);

if ($perlonly) {
	# Calculate normalization factors or use external (e.g. from edgeR)
	if (@extnorm) {
		for (my $i=0; $i<@bglist; $i++) {
			$normfactor{$bglist[$i]} = $extnorm[$i];
		}
	}
	else {
		if ($ncores==1) {
			for (my $i=0; $i<@bglist; $i++) {
				&disp("Reading ".basename($bglist[$i])."...");
				$wigsum{$bglist[$i]} = 0;
				open(BGIN,$bglist[$i]);
				while ($line = <BGIN>) {
					&disp("  Read $. signals from ".basename($bglist[$i])."...") 
						if ($.%1000000 == 0);
					$line =~ s/\r|\n$//g;
					($chr,$start,$end,$signal) = split(/\t/,$line);
					$wigsum{$bglist[$i]} += ($end - $start)*$signal;
				}
				close(BGIN);
				$normfactor{$bglist[$i]} = $sumto/$wigsum{$bglist[$i]};
			}
		}
		else {
			my $pl1 = Parallel::Loops->new($ncores);
			$pl1->share(\%wigsum,\%normfactor);
			$pl1->foreach(\@bglist, sub {
				&disp("Reading ".basename($_)."...");
				$wigsum{$_} = 0;
				open(BGIN,$_);
				while ($line = <BGIN>)
				{
					$line =~ s/\r|\n$//g;
					($chr,$start,$end,$signal) = split(/\t/,$line);
					$wigsum{$_} += ($end - $start)*$signal;
				}
				close(BGIN);
				$normfactor{$_} = $sumto/$wigsum{$_};
			});
		}
	}

	if (!$prerun && !$prerunlog) {
		# Reparse the bedgraph, rescale and write
		if ($ncores==1) {
			for (my $i=0; $i<@bglist; $i++) {
				if ($output[$i]) {
					&disp("Writing ".basename($output[$i])."...");
					open(OUTPUT,">$output[$i]");
					$out = *OUTPUT;
				}
				else {
					&disp("Writing to stdout...");
					$out = *STDOUT;
				}
				open(BGIN,$bglist[$i]);
				while ($line = <BGIN>) {
					&disp("  Wrote $. normalized signals for "
						.basename($bglist[$i])."...") if ($.%1000000 == 0);
					$line =~ s/\r|\n$//g;
					($chr,$start,$end,$signal) = split(/\t/,$line);
					$signal = &round($signal*$normfactor{$bglist[$i]},2);
					print $out "$chr\t$start\t$end\t$signal\n";
				}
				close(BGIN);
			}
		}
		else {
			my $pl2 = Parallel::Loops->new($ncores);
			$pl2->share(\%wigsum,\%normfactor);
			$pl2->foreach(\@bglist, sub {
				my $outp = &create_output_file($_);
				&disp("Writing ".basename($outp)."...");
				open(OUTPUT,">$out");
				$out = *OUTPUT;
				
				open(BGIN,$_);
				while ($line = <BGIN>)
				{
					$line =~ s/\r|\n$//g;
					($chr,$start,$end,$signal) = split(/\t/,$line);
					$signal = &round($signal*$normfactor{$_},2);
					print $out "$chr\t$start\t$end\t$signal\n";
				}
				close(BGIN);
			});
		}
	}
	else {
		if ($prerunlog) {
			open(LOG,">$prerunlog");
			print LOG "filename\ttotal signal\n";
			foreach my $f (keys(%wigsum)) {
				print LOG basename($f)."\t$wigsum{$f}\n";
			}
			close(LOG);
		}
		else {
			&disp("\nfilename\ttotal signal");
			foreach my $f (keys(%wigsum)) {
				&disp(basename($f)."\t$wigsum{$f}");
			}
			&disp("")
		}
	}
}
else { # Use awk
	# Calculate normalization factors or use external (e.g. from edgeR)
	if (@extnorm) {
		for (my $i=0; $i<@bglist; $i++) {
			$normfactor{$bglist[$i]} = $extnorm[$i];
		}
	}
	else {
		if ($ncores==1) {
			for (my $i=0; $i<@bglist; $i++) {
				&disp("Reading ".basename($bglist[$i])."...");
				$wigsum{$bglist[$i]} = `awk '{sum += (\$3-\$2)*\$4; if (FNR \% 1000000==0) { printf("\\n  Read %d signals from %s",FNR,FILENAME) | "cat 1>&2"}} END {print sum}' $bglist[$i]`;
				$normfactor{$bglist[$i]} = $sumto/$wigsum{$bglist[$i]};
			}
		}
		else {
			my $pl3 = Parallel::Loops->new($ncores);
			$pl3->share(\%wigsum,\%normfactor);
			$pl3->foreach(\@bglist, sub {
				&disp("Reading ".basename($_)."...");
				$wigsum{$_} = `awk '{sum += (\$3-\$2)*\$4; } END {print sum}' $_`;
				$normfactor{$_} = $sumto/$wigsum{$_};
			});
		}
	}

	if (!$prerun && !$prerunlog) {
		# Reparse the bedgraph, rescale and write
		if ($ncores==1) {
			for (my $i=0; $i<@bglist; $i++) {
				if ($output[$i]) {
					&disp("Writing ".basename($output[$i])."...");
					`awk -v var="$output[$i]" '{printf \"%s\\t%s\\t%s\\t%.2f\\n\", \$1,\$2,\$3,\$4*$normfactor{$bglist[$i]}; if (FNR \% 1000000==0) { printf("\\n  Wrote %d signals to %s",FNR,var) | "cat 1>&2"}}' $bglist[$i] > $output[$i]`;
					
				}
				else {
					&disp("Writing to stdout...");
					`awk '{printf \"%s\\t%s\\t%s\\t%.2f\\n\", \$1,\$2,\$3,\$4*$normfactor{$bglist[$i]}}' $bglist[$i]`;
				}
			}
		}
		else {
			my $pl4 = Parallel::Loops->new($ncores);
			$pl4->share(\%wigsum,\%normfactor);
			$pl4->foreach(\@bglist, sub {
				my $outp = &create_output_file($_);
				&disp("Writing ".basename($outp)."...");
				`awk -v var="$outp" '{printf \"%s\\t%s\\t%s\\t%.2f\\n\", \$1,\$2,\$3,\$4*$normfactor{$_};}' $_ > $outp`;
			});
		}
	}
	else {
		if ($prerunlog) {
			open(LOG,">$prerunlog");
			print LOG "filename\ttotal signal\n";
			foreach my $f (keys(%wigsum)) {
				print LOG basename($f)."\t$wigsum{$f}";
			}
			close(LOG);
		}
		else {
			print STDERR "\nfilename\ttotal signal\n";
			foreach my $f (keys(%wigsum)) {
				print STDERR basename($f)."\t$wigsum{$f}";
			}
		}
	}
}

if ($exportfacs) {
	&disp("Exporting normalization factors...");
	open(FACS,">$exportfacs");
	print FACS "file\tnormalization factor\ttotal signal\n";
	foreach my $f (keys(%normfactor)) {
		print FACS basename($f)."\t".&round($normfactor{$f},6)."\t$wigsum{$f}";
		print FACS "\n" if ($perlonly); # For some reason, awk maintains the EOL
	}
	close(FACS);
}

&disp("Finished!\n\n");

# Process inputs
sub check_inputs
{
    GetOptions(
        "input|i=s{,}" => \@bglist,
        "output|o=s{,}" => \@output,
        "extnorm|r=f{,}" => \@extnorm,
        "sumto|s=i" => \$sumto,
        "exportfactors|f=s" => \$exportfacs,
        "perlonly|p" => \$perlonly,
        "prerun|r" => \$prerun,
        "prerunlog|l=s" => \$prerunlog,
        "ncores|c=i" => \$ncores,
        "man|m" => \$man,
        "log|g" => \$log,
        "silent|s" => \$silent,
        "help|h" => \$help
    );
    
    Pod::Usage::pod2usage( -verbose => 1, -exitstatus => 0 ) if ($help);
    Pod::Usage::pod2usage( -exitstatus => 0, -verbose => 2 ) if ($man);
    
    # Check fatal
    my $stop;
    $stop .= "--- Please specify input file(s) ---\n" if (!@bglist);
    
    if ($stop) {
        &disp("$stop\n");
        &disp("Type perl $scriptname --help for help in usage.\n\n");
        exit;
    }

    # Check number of cores and parallel mode
    if ($ncores) {
		if ($ncores>1) {
            my $status = eval { try_module("Parallel::Loops") };
            if ($status) {
                &disp("Module Parallel::Loops is required to use multiple ",
					"cores! Using 1 core...");
                $ncores = 1;
            }
            else { use Parallel::Loops; }
        }
    }
    else {
        $ncores = 1;
    }
    
    # Check signal summarization
    if (!$sumto) {
        &disp("The normalization total signal not given! ",
			"Using default(1000000000)...");
        $sumto = 1000000000;
    }
    if ($sumto !~ m/^-?[1-9]\d*$/ || $sumto == 0) {
        &disp("The signal sum parameter must be a positive or negative ",
			"integer! Using default (1000000000)...");
        $sumto = 1000000000;
    }
    
    # Check presence of external normalization factors
    if (@extnorm) {
        my $e = @extnorm;
        my $b = @bglist;
        if ($e != $b) {
            &disp("The number of external normalization factors given must ",
				"be equal to the number of input files! Ignoring...");
            @extnorm = ();
        }
        if (@extnorm && $sumto) {
            &disp("Normalizing signal sum (sumto) and external normalization ",
				"factors (extnorm) are mutually exclusive! Ignoring sumto...");
            $sumto = 0;
        }
    }
    
    # Check dry run and output files
    if (!$prerun && !$prerunlog) {
        if (@output) {
            if (@output && $output[0] ne "stdout") {
                my $o = @output;
                my $b = @bglist;
                if ($o != $b) {
                    &disp("The number of output files must be equal to the ",
						"number of input files! Autogenerating...");
                    @output = ();
                }
            }
            elsif (!@output) {
                &disp("Output filenames will be autogenerated...");
                @output = ();
            }
        }
        else {
            &disp("Output filenames will be autogenerated...");
            @output = ();
        }
    }

    # Check if we are on Linux for the usage of awk
    if (!$perlonly) {
        if ($^O =~ /MSWin/) { # Windows... bad news...
            &disp("Windows OS detected! Switching to pure ",
				"Perl for file streaming...");
            $perlonly = 1;
        }
    }
}

sub create_output_file {
    my $in = $_[0];
    my ($base,$dir,$ext) = fileparse($in,'\.[^.]*');
    return(File::Spec->catfile($dir,$base."_norm".$ext));
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


__END__

=pod

=head1 NAME

normalize_bedgraph.pl - Normalize bedgraph files to a constant total signal

=head1 SYNOPSIS

normalize_bedgraph.pl --input file1 [file2, file3, ..., filen] [OPTIONS]

Examples:

=over 4

=item perl normalize_bedgraph.pl --input repli_1.bg repli_2.bg --output 
repli_1_norm.bg repli_2_norm.bg --sumto 5000000000 --exportfactors

=item perl normalize_bedgraph.pl --input repli_1.bg repli_2.bg --perlonly 
--extnorm 0.92 1.28

=back

=head1 DESCRIPTION

Normalize bedgraph signal using to a total signal value (similar to RSeqC
normalize_bigwig.py but faster!) or a set of external normalization factors
(e.g. calculated from DESeq or edgeR).

=head1 ARGUMENTS

=over 4

=item input B<(required)>

--input or -i

Input bedgraph file(s). Please be careful as there is checking whether the input
file(s) are indeed bedgraph files. It's ok if they contain more than 4 columns
but the first four must be bedgraph (chromosome, start, end, signal separated by
tabs). Input files need not to be sorted.

=item output B<(optional)>

--output or -o

Output file names. It can be "stdout" for exporting to STDOUT, a set of file
names equal to the number of input files or nothing for the files to be
auto-generated.

=item sumto B<(optional)>

--sumto or -s

Normalize to --sumto total wig signal. Defaults to 1000000000. It is mutually
exclusive with --extnorm with --extnorm in precedence.

=item extnorm B<(optional)>

--extnorm or -e

A set of external normalization factors (e.g. calculated from DESeq or edgeR).
It is mutually exclusive with --sumto with --extnorm in precedence.

=item exportfactors B<(optional)>

--exportfactors or -f

Export normalization factors and signal sums to a file specified by
--exportfactors.

=item perlonly B<(optional)>

--perlonly or -p

Use pure Perl to run the script, otherwise, uses Linux awk. Useful for e.g.
Windows systems but slower.

=item prerun B<(optional)>

--prerun or -r

If this switch is turned on, the script just counts the total wiggle signal in
the input files and prints it on the screen. This is useful in order for example
to determine the total normalization signal (--sumto).

=item prerunlog B<(optional)>

--prerunlog or -l

Writes the output of --prerun to a file specified by --prerunlog. If only the
--prerunlog is specified, --prerun is executed automatically.

=item ncores B<(optional)>

--ncores or -c

Runs the script in parallel mode where one core processes one bedgraph file. It
requires the module Parallel::Loops to be installed.

=item silent B<optional>

--silent or -s

Do not display verbose messages.

=item help B<(optional)>

--help or -h

Display this help text.

=item man B<(optional)>

--man or -m

Display the full manual of the script.

=back

=head1 AUTHOR

Panagiotis Moulos (L<moulos@fleming.gr)>

=cut
