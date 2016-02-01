#!/usr/bin/perl
use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;

## Examines every field one by one. Split by '|', remove heading and
## tailing white spaces. Attempts to 'unique' the fields. If more than
## one field remains, replace with a predefined string (Multiple).
my $multCst = '(Multiple)';

sub collapseRecord {
    my $rec = shift(@_);

    my @splitRec = split /\|/, $rec;
    my %uniqRec = ();
    for (my $i=0; $i<@splitRec; $i++) {
	$splitRec[$i] =~ s/^[ ]*//g;
	$splitRec[$i] =~ s/[ ]*$//g;
	$uniqRec{$splitRec[$i]} = 1;
    }
    
    my @uniqRecArray = keys %uniqRec;
    if ( @uniqRecArray == 0 ) {
	return("");
    } elsif ( @uniqRecArray == 1 ) {
	return($uniqRecArray[0]);
    } else {
	return($multCst);
    }
}

if ( collapseRecord("OK") ne "OK" ) { die("test1"); }
if ( collapseRecord("OK | OK") ne "OK" ) { die("test2"); }
if ( collapseRecord("YES | NO") ne $multCst ) { die("test3"); }
if ( collapseRecord("") ne "" ) { die("test4"); }

my $line = 0;
while (<>) {
    $line++;
    my $lastEmpty = 0;

    s/\n$//;
    $lastEmpty = 1 if  /\t$/;
    my @record = split /\t/;
    # print STDERR 'line '.$line.': '.@record."\n";

    foreach (my $i=0; $i<@record; $i++) {
	$record[$i] = collapseRecord($record[$i]);
    }
    # print STDERR 'line '.$line.': '.@record."\n";

    print join("\t", @record).($lastEmpty ? "\t" : "" )."\n";
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
