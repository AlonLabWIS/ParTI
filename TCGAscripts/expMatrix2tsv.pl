#!/usr/bin/perl
use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;

my $expFile = "";
my $geneListFile = "geneList.list";
my $tsvFile = "expMatrix.tsv";
my $verbose = 0;
my $help = 0;
my $man = 0;

GetOptions("verbose!" => \$verbose,
           "expFile|e=s" => \$expFile,
           "geneList|g=s" => \$geneListFile,
	   "tsvMatrix|c=s" => \$tsvFile,
           "help|?"  => \$help,
           man => \$man) or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;

my @patientIDs = ();
my @genes = ();
my $expression = [];
my $n = 0;

open E, "<".$expFile or die("Could not read from $expFile.\n");
print STDERR "Reading in $expFile...\n" if $verbose;
while (<E>) {
    $_ =~ s/\n$//g;
    my @l = split /\t/;
    if ( @patientIDs == 0 ) {
	shift @l;
	@patientIDs = @l;	
    } else {
	my $thisGene = shift @l;
	push @genes, $thisGene;
	$expression->[@{ $expression }] = \@l;
	# print STDERR "Stored $l[0], $l[1], $l[2], ... in index ".
	#     ( @{ $expression } - 1)." of expression\n" if $verbose;
	# print STDERR "Check: ".
	#     $expression->[@{ $expression } - 1]->[0].", ".
	#     $expression->[@{ $expression } - 1]->[1].", ".
	#     $expression->[@{ $expression } - 1]->[2]."\n";
    }
    $n++;
    # last if $verbose && $n == 3;
}
close E;

print STDERR "Now writing out list of genes to $geneListFile.\n" if $verbose;
open G, ">".$geneListFile or die("Could not write to $geneListFile.\n");
print G join("\n", @genes);
close G;

print STDERR "Now writing out expression data to $tsvFile.\n" if $verbose;
open T, ">".$tsvFile or die("Could not write to $tsvFile.\n");
foreach my $jPatient (1..@patientIDs) {
    print T $patientIDs[$jPatient-1];
    foreach my $iGene (1..@{ $expression }) {
	print T "\t".$expression->[$iGene-1]->[$jPatient-1];
    }
    print T "\n";
}
close T;

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
