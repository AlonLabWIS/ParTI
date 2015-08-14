#!/usr/bin/perl
use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;

## Finds the number of fields in the first row.
## If the number of fields in the second row is one less than in the
## first, add enough empty fields at the end of the line.

my $n = 0;
my $nFields = -1;
my $added = 0;
my @fields = ();

while (<>) {
    $_ =~ s/\n$//g;
    my @l = split /\t/;
    $n++;

    if ( $nFields == -1 ) { 
	$nFields = @l; 
	@fields = @l;
    }

    if ( @l < $nFields ) {
	my $missing = $nFields - @l;
	foreach my $i (1..$missing) {
	    push @l, "NA";
	}
	$added++;
	print STDERR "Added ".$missing.
	    " fields at the end of line $n.\n";
    }
    
    # if ( @l == $nFields - 1 ) {
    # 	push @l, "";
    # 	$added++;
    # }
    if ( @l != $nFields ) {
	print STDERR "\n----- Error -----\n";
	foreach my $i (1..(@l > $nFields ? $nFields : @l)) {
	    print STDERR "Field $i: ".$fields[$i-1]."\n".
		$l[$i-1]."\n";
	}
	die("On line $n, we found ".@l." fields while we ".
	    "expected ".$nFields.". Aborting.\n");
    }
    print join("\t", @l)."\n";
}

if ( $added > 0 ) {
    print STDERR "Padded empty fields at the end of $added lines. ". 
	"Please check the output file visually to verify that ".
	"contents match headers.\n";
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
