#!/usr/bin/env perl

use strict;

our $DEBUG = 0;

my $usage = "usage: $0 extFile\n\n";

my $extFile = $ARGV[0] or die $usage;

open (EXT, $extFile) or die "Cannot open $extFile\n";
my $header = <EXT>;

my @chains;
while (<EXT>) {
    unless (/\w/) { next;}
    chomp;
    $_ =~ s/^\s+//; #rm leading whitespace
    ## using var names as in ext.c
    my ($dstart, $dend, $score, $astart, $aend, $orient, $zero1, $zero2, $acc) = split (/\s+/);
    
    push (@chains, 
	  { dstart=> $dstart,
	    dend => $dend,
	    score => $score,
	    astart => $astart,
	    aend => $aend,
	    orient => $orient,
	    zero1 => $zero1,
	    zero2 => $zero2,
	    acc => $acc } );
    
}

## Sort the chains by acc, orient and start

@chains = sort { 
    $a->{acc} cmp $b->{acc}
    ||
    $a->{orient} <=> $b->{orient} 
    ||
	$a->{dstart} <=> $b->{dstart}
} @chains;


my @collapsedChains = shift @chains;

while (@chains) {
    my $prevChain = $collapsedChains[$#collapsedChains];
    
    my ($prev_dstart, 
	$prev_dend, 
	$prev_score, 
	$prev_astart, 
	$prev_aend, 
	$prev_orient, 
	$prev_acc) = ($prevChain->{dstart},
		      $prevChain->{dend},
		      $prevChain->{score},
		      $prevChain->{astart},
		      $prevChain->{aend},
		      $prevChain->{orient},
		      $prevChain->{acc});
    
    my $nextChain = shift @chains;
    my ($next_dstart, 
	$next_dend,
	$next_score,
	$next_astart,
	$next_aend,
	$next_orient,
	$next_acc) = ($nextChain->{dstart},
		      $nextChain->{dend},
		      $nextChain->{score},
		      $nextChain->{astart},
		      $nextChain->{aend},
		      $nextChain->{orient},
		      $nextChain->{acc});
    
    
    print "$next_acc vs $prev_acc, $next_orient vs. $prev_orient, $next_dstart, $prev_dend\t" if $DEBUG;

    if ( $next_acc eq $prev_acc 
	&& 
	 $next_orient == $prev_orient
	 && 
	 $next_dstart <= $prev_dend) {
	
	## merge overlapping entry:
	my @dcoords = sort {$a<=>$b} ($prev_dstart, $prev_dend, $next_dstart, $next_dend);
	my $dstart = shift @dcoords;
	my $dend = pop @dcoords;

	my @acoords = sort {$a<=>$b} ($prev_astart, $prev_aend, $next_astart, $next_aend);
	my $astart = shift @acoords;
	my $aend = pop @acoords;
	
	my @scores = sort {$a<=>$b} ($prev_score, $next_score);
	my $score = pop @scores;
	
	## adjust prevChain contents
	$prevChain->{dstart} = $dstart;
	$prevChain->{dend} = $dend;
	$prevChain->{astart} = $astart;
	$prevChain->{aend} = $aend;
	$prevChain->{score} = $score;

	print "expanding current chain.\n" if $DEBUG;

    } 

    else {

	## terminate previous chain
	push (@collapsedChains, $nextChain);
	print "next chain.\n" if $DEBUG;
    }
}


## Sort collapsed chains by dstart and score:

@collapsedChains = sort { $a->{dstart} <=> $b->{dstart}
			  || 
			      $b->{score} <=> $a->{score} } @collapsedChains;



## output collapsed chains:
print $header;

foreach my $chain (@collapsedChains) {
    my ($dstart, 
	$dend,
	$score,
	$astart,
	$aend,
	$orient,
	$acc, 
	$zero1,
	$zero2 ) = ($chain->{dstart},
		    $chain->{dend},
		    $chain->{score},
		    $chain->{astart},
		    $chain->{aend},
		    $chain->{orient},
		    $chain->{acc},
		    $chain->{zero1},
		    $chain->{zero2}
		    );
    
    printf("%8d %8d %6d %7d %5d %1d %5d %5d %s\n",
	   $dstart, $dend, $score, $astart, $aend, $orient, $zero1, $zero2, $acc);
    
}

exit(0);


