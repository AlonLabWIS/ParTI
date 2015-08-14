#!/usr/bin/perl
use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;

my $clinicalFile = "clinicalData_reOrdered.tsv";
my $featuresFile = "continuousFeatures.list";
my $tsvFile = "continuousClinicalData_reOrdered.tsv";
my $verbose = 0;
my $help = 0;
my $man = 0;

GetOptions("verbose!" => \$verbose,
	   "clinicalMatrix|c=s" => \$clinicalFile,
           "discreteFeatures|f=s" => \$featuresFile,
	   "tsvMatrix|c=s" => \$tsvFile,
           "help|?"  => \$help,
           man => \$man) or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;

my @continuousFeatures = ();
open F, "<".$featuresFile or die("Could not read from $featuresFile.\n");
while (<F>) {
    chomp;
    push @continuousFeatures, $_;
}
close F;

my $isFirstLine = 1;
my %featToIndex = ();
my @colsToKeep = ();
open F, "<".$clinicalFile or die("Could not read from $clinicalFile.\n");
open G, ">".$tsvFile or die("Could not write to $clinicalFile.\n");
while (<F>) {
    chomp;
    my @l = split /\t/;
    if ( $isFirstLine ) {
	foreach my $i (0..$#l) {
	    $featToIndex{$l[$i]} = $i;
	}
	foreach my $feat (@continuousFeatures) {
	    if ( exists($featToIndex{$feat}) ) {
		push @colsToKeep, $featToIndex{$feat};
	    } else {
		print STDERR "Could not find $feat in $clinicalFile.\n";
	    }
	}
    }
    foreach my $i (0..$#colsToKeep) {
	if ( ! $isFirstLine && ! is_numeric($l[$colsToKeep[$i]]) ) {
	# if ( $l[$colsToKeep[$i]] eq "" || $l[$colsToKeep[$i]] eq "NA" ) {
		$l[$colsToKeep[$i]] = "NaN";
	}
	print G ($i == 0 ? "" : "\t").$l[$colsToKeep[$i]];
    }
    print G "\n";
    $isFirstLine = 0;
}
close G;
close F;

sub getnum {
    use POSIX qw(strtod);
    my $str = shift;
    $str =~ s/^\s+//;
    $str =~ s/\s+$//;
    $! = 0;
    my($num, $unparsed) = strtod($str);
    if (($str eq '') || ($unparsed != 0) || $!) {
	return undef;
    }
    else {
	return $num;
    }
}
sub is_numeric { defined getnum($_[0]) }

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
