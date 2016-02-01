#!/usr/bin/perl
use warnings;
use strict;
use Pod::Usage;
use Getopt::Long;

my $nArchFracFile = "TCGA_frac_nArchs.tab"; # file with project code,
				# frac genes and  archetypes to look for
my $doSub = 0;
my $verbose = 0;
my $help = 0;
my $man = 0;

GetOptions("verbose!" => \$verbose,
	   "nArchFracFile|f=s" => \$nArchFracFile,
	   "sub!" => \$doSub,
	   "help|?"  => \$help,
	   man => \$man) or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;

my %projectFrac;
my %projectNArch;
open F, "<".$nArchFracFile;
while (<F>) {
    chomp;
    my @l = split /\t/;
    $projectFrac{$l[0]} = $l[1];
    $projectNArch{$l[0]} = $l[2];
}
close F;

# Now iterate through projects, rewrite fraction and nArch in .m
# script and (submit)
foreach my $project (sort keys %projectFrac) {
    my $jobName = sprintf("ParTI_%s", $project);
    my $mFile = sprintf("%s.m", $jobName);
    
    open F, "<ParTI/TCGAscripts/clusterParTI.m" 
	or die("Could not open cluster matlab ParTI script.\n");
    open G, ">".$project.'_UCSC/'.$mFile;
    while (<F>) {
	$_ =~ s/nArchetypes/$projectNArch{$project}/g;
	$_ =~ s/myQuantile/$projectFrac{$project}/g;
	print G
    }
    close G;
    close F;
    
    my $cmd = 'cd '.$project.'_UCSC && '.
	'/opt/MATLAB/R2012a/bin/matlab -nodesktop -nodisplay < '.
	$mFile;
    if ( $doSub ) {
	print STDERR "===== Now processing $project =====\n";
	`$cmd`;
    } else {
	print $cmd."\n";
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
