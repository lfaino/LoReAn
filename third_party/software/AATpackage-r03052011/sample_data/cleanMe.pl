#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;

# ensure we're in the same directory as the script:
chdir "$FindBin::Bin" or die "Error, cannot cd to $FindBin::Bin";


my @keep = qw ( cleanMe.pl 
                runMe.sh 
arab.pep
arab.genomicSeq
arab.cdna
);

my %KEEP = map { + $_ => 1 } @keep;

foreach my $file (<*>) {
	unless ($KEEP{$file}) {
		print STDERR "-removing file $file\n";
		unlink($file);
	}
}

`rm -rf *fasta.chrysalis`;

exit(0);


	
