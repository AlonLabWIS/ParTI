#!/usr/bin/perl
use warnings;
use strict;
use Pod::Usage;
use Getopt::Long;

my $nArchetypesMin = 3; # How many archetypes to look for
my $nArchetypesMax = 5; # How many archetypes to look for
my $fracStep = 0.1; # Step size in the fraction of top expressed genes
		  # we'll consider
my $queue = "alon"; # new-all.q or alon
my $tag = "";
my $doSub = 0;
my $verbose = 0;
my $help = 0;
my $man = 0;

GetOptions("verbose!" => \$verbose,
	   "nArchetypesMin|n=i" => \$nArchetypesMin,
	   "nArchetypesMax|x=i"  => \$nArchetypesMax,
	   "fracStep|f=f"  => \$fracStep,
	   "queue|q=s" => \$queue,
	   "tag|t=s" => \$tag,
	   "sub!" => \$doSub,
	   "help|?"  => \$help,
	   man => \$man) or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;

for ( my $nArch = $nArchetypesMin; 
      $nArch <= $nArchetypesMax; 
      $nArch++ ) {
    for ( my $frac = 0; $frac < 1.0; $frac += $fracStep ) {
	my $jobName = sprintf("scanArchFrac_a%d_f%.2f", $nArch, $frac);
	my $mFile = sprintf("%s.m", $jobName);
	my $bFile = sprintf("%s.sh", $jobName);
	# my $outFile = sprintf("%s.out", $jobName);
	my $bOutFile = sprintf("%s.bOut", $jobName);
	my $bErrFile = sprintf("%s.bErr", $jobName);

	# print STDERR "Testing $bOutFile\n";
	if ( -e $bOutFile ) {
		print STDERR "Skipping $jobName because $bOutFile already exists.\n";
		next;
	}
	
	open F, "<../ParTI/TCGAscripts/TCGAscanGeneFrac.m";
	open G, ">".$mFile;
	while (<F>) {
	    $_ =~ s/nArchetypes/$nArch/g;
	    $_ =~ s/myQuantile/$frac/g;
	    print G
	}
	close G;
	close F;
	
	open G, ">".$bFile;
	print G "#BSUB -q $queue"."\n";
	print G "#BSUB -J $jobName"."\n";
	print G '#BSUB -R "select[mem>7000] rusage[mem=7000]"'."\n";
	print G "#BSUB -o $bOutFile"."\n";
	print G "#BSUB -e $bErrFile"."\n";
	# print G "module load matlab/R2012a"."\n";
	print G "hostname"."\n";
	print G "date"."\n";
	# print G sprintf("/opt/MATLAB/R2012a/bin/matlab -nojvm -nodisplay < %s", $mFile)."\n";
	print G sprintf("matlab -nodesktop -nodisplay < %s", $mFile)."\n";
	print G "date"."\n";	
	close G;

	my $cmd = 'bsub < '.$bFile;	
	if ( $doSub ) {
	    `$cmd`;
	} else {
	    print $cmd."\n";
	}

    }
}

__END__

=head1 NAME

template - What it does

=head1 SYNOPSIS

template [options]

Options:
    -help            brief help message
    -man             full documentation

=head1 OPTIONS

=over 8

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

B<This program> will read the given input file(s) and do something
useful with the contents thereof.

=cut
