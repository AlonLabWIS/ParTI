#!/usr/bin/perl
use strict;
use warnings;
use POSIX;
use Pod::Usage;
use Getopt::Long;

my $inFile = "inMatrix.tsv";
my $outFile = "outMatrix.tsv";
my $value = 1;
my $valFreq = -1;
my $entropyTop = -1;
my $ignoreVal = "";
my $round = 0;
my $verbose = 0;
my $help = 0;
my $man = 0;

GetOptions("verbose!" => \$verbose,
           "in|i=s" => \$inFile,
           "out|o=s" => \$outFile,
	   "value|v=s" => \$value,
	   "valFreq|f=f" => \$valFreq, #in %
	   "entropyTop|e=f" => \$entropyTop, #in %
	   "ignoreVal|g=s" => \$ignoreVal,
	   "round|r!" => \$round,
           "help|?"  => \$help,
           man => \$man) or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;

my @patientIDs = ();
my @genes = ();
my $matrix = [];
my $n = 0;

open E, "<".$inFile or die("Could not read from $inFile.\n");
print STDERR "Reading from $inFile...\n" if $verbose;
while (<E>) {
    chomp;
    my @l = split /\t/;
    if ( @genes == 0 ) {
	shift @l;
	@genes = @l;	
    } else {
	my $thisPatient = shift @l;
	# print STDERR "Got patient $thisPatient\n";
	push @patientIDs, $thisPatient;
	if ( $round ) {
	    foreach my $i (0..$#l) {
		$l[$i] = ( $l[$i] eq "NaN" ? "NaN" : floor($l[$i]) );
	    }
	}
	$matrix->[@{ $matrix }] = \@l;
	# print STDERR "Stored $l[0], $l[1], $l[2], ... in index ".
	#     ( @{ $matrix } - 1)." of matrix\n" if $verbose;
	# print STDERR "Check: ".
	#     $matrix->[@{ $matrix } - 1]->[0].", ".
	#     $matrix->[@{ $matrix } - 1]->[1].", ".
	#     $matrix->[@{ $matrix } - 1]->[2]."\n";
    }
    $n++;
    # last if $verbose && $n == 3;
}
close E;

# print STDERR "Patients: ".join(", ", @patientIDs)."\n" if $verbose;

my @toKeep = 0..$#genes;

if ( $valFreq > 0 ) {
    print STDERR "Keeping only columns for which the value '".$value.
	"' occurs in at least ".$valFreq."% of the records.\n";
    @toKeep = ();
    foreach my $i (0..$#genes) {
	my $nMatches = 0;
	foreach my $j (0..$#patientIDs) {
	    $nMatches++ if $matrix->[$j]->[$i] eq $value;
	}
	push @toKeep, $i if ( $nMatches / @patientIDs ) > ( $valFreq / 100 );
    }
} elsif ( $entropyTop > 0 ) {
    print STDERR "Keeping only the ".$entropyTop."% columns with ".
	"highest entropy\n";
    @toKeep = ();
    my @entropies = ();
    foreach my $i (0..$#genes) {
	$entropies[$i] = 0;

	my %valHash = ();
	foreach my $j (0..$#patientIDs) {
	    my $thisVal = $matrix->[$j]->[$i];
	    if ( ! exists($valHash{$thisVal}) ) {
		$valHash{$thisVal} = 0;
	    }
	    $valHash{$thisVal}++;
	}
	my @counts = values %valHash;
	foreach my $j (0..$#counts) {
	    $counts[$j] = $counts[$j] / @patientIDs;
	    # Counts now contains probabilities.
	    # We now compute base 2 entropy.
	    $entropies[$i] -= ( $counts[$j] == 0 ? 0 : $counts[$j] *
				log($counts[$j])/log(2) );
	}
	# print STDERR "c: ".join(", ", @counts)." -> ".$entropies[$i]."\n";
    }
    # print STDERR "Entropies: ".join(", ", @entropies[0..5]).", ...\n";
    @toKeep = sort { $entropies[$b] <=> $entropies[$a] } 0..$#entropies;
    my $nKept = ceil(@toKeep * $entropyTop / 100);
    @toKeep = @toKeep[0..($nKept-1)];
    print STDERR "Highest entropies: ".
	join(", ", @entropies[@toKeep[0..5]]).", ...\n" if $verbose;
}

print STDERR "Keept ".@toKeep." / ".@genes." features\n".
    "Now booleanizing feature matrix.\n" if $verbose;

my @boolGenes = ();
my @boolMatrix = ();
foreach my $idx (@toKeep) {
    # last if $verbose && $idx > 3;
    my @thisFeature = ();
    my %valHash = ();
    foreach my $j (0..$#patientIDs) {
	# print STDERR "Matrix $j, $idx: ".$matrix->[$j]->[$idx]."\n" if $verbose;
	push @thisFeature, $matrix->[$j]->[$idx];
	$valHash{$matrix->[$j]->[$idx]} = 1;
    }
    my @uniqVals = keys %valHash;

    foreach my $uniqVal (@uniqVals) {
	next if $uniqVal eq $ignoreVal; #Ignore certain values
	push @boolGenes, $genes[$idx]."=".$uniqVal;
	my @thisBoolFeature = ();
	foreach my $x (@thisFeature) {
	    push @thisBoolFeature, ($x eq $uniqVal ? 1 : 0);
	}
	push @boolMatrix, \@thisBoolFeature
    }
}

print STDERR "Now writing out boolean matrix to $outFile.\n" if $verbose;
open G, ">".$outFile or die("Could not write to $outFile.\n");
print G "sampleID\t".join("\t", @boolGenes)."\n";

foreach my $jPatient (0..$#patientIDs) {
    print G $patientIDs[$jPatient];
    foreach my $iGene (0..$#boolGenes) {    
	print G "\t".$boolMatrix[$iGene]->[$jPatient];
    }
    print G "\n";
}
close G;

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
