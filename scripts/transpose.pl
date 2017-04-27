#!/usr/bin/perl
use warnings;
use strict;

#script to transpose files. takes a list of files as an input, transposes each one (doesn't replace it)

my $listfile = shift @ARGV;
chomp $listfile;
open (LIST, "<$listfile") || die "Cannot open $listfile\n";
while (my $file = <LIST>) { #open file list
    chomp $file;
    my @rows;
    my @transposed;
    open (FILE,"$file") || die "cannot open $file\n";
    open (OUTFILE, ">${file}_transposed") or die ("cannot open OUTFILE");
    while(<FILE>) {
        push @rows, [ split /\t/ ];
    }
    foreach my $row (@rows) {
        foreach my $column (0 .. $#{$row}) {
            push(@{$transposed[$column]}, $row->[$column]);
        }
    }

    foreach my $new_row (@transposed) {
        my @currentrow;
        foreach my $new_col (@{$new_row}) {
            my $col = trim($new_col); #get rid of the newlines on the last column
            push(@currentrow, $col);
        }
	my $toprint = join "\t", @currentrow;
        print OUTFILE "$toprint\n";
    }
    close FILE;
    close OUTFILE;
}

#trim leading and trailing whitespace
sub trim {
	my $string = shift;
        $string =~ s/^\s+//;
        $string =~ s/\s+$//;
        return $string;
}
