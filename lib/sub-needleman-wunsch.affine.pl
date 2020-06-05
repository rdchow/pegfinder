#!/usr/bin/perl
##############################################################################
# Needleman-Wunsch  Algorithm 

use strict;
use warnings;
use Algorithm::NeedlemanWunsch;

sub needleman_wunsch {
    my ($seq1, $seq2) = @_;
    my $arr1 = [split ("",$seq1)];
    my $arr2 = [split ("",$seq2)];
    my (@a_align, @b_align);

    my $matcher = Algorithm::NeedlemanWunsch->new(\&score_sub);
    $matcher->gap_open_penalty(-2);
    $matcher->gap_extend_penalty(-0.1);
    
    my $score = $matcher->align($arr1,$arr2,
        {   align   => sub {unshift @a_align, $arr1->[shift]; unshift @b_align, $arr2->[shift]},
            shift_a => sub {unshift @a_align, $arr1->[shift]; unshift @b_align,'-'},
            shift_b => sub {unshift @a_align,'-'; unshift @b_align, $arr2->[shift]},
        });

    my $align1 = "@a_align";
    my $align2 = "@b_align";
    $align1 =~ s/(.)\s/$1/seg;
    $align2 =~ s/(.)\s/$1/seg;

    print "$align1\n$align2\n";
   

    #Detect if the 5' or 3' ends of the two sequences are not aligned exactly
    #first look at 5' ends
    my $trimmingStatus5p = 0; 
    # ^ variable that tracks 5' trimming needs. 0: no trimming, 1:trimmed wildtype seq, 2: trimmed edited seq
    my @chars1 = split("",$align1);
    my @chars2 = split("",$align2);

    if ($chars1[0] eq "-"){
        $trimmingStatus5p = 1;
    }
    elsif ($chars2[0] eq "-"){
        $trimmingStatus5p = 2;
    }
    
    #now look at 3' ends
    my $trimmingStatus3p = 0; 
    # ^ variable that tracks 3' trimming needs. 0: no trimming, 1:trimmed wildtype seq, 2: trimmed edited seq
    if ($chars1[-1] eq "-"){
        $trimmingStatus3p = 1;
    }
    elsif ($chars2[-1] eq "-"){
        $trimmingStatus3p = 2;
    }
    my $wtdelcounter = 0; # counts how many gaps are in the wildtype alignment; need for calculating the edit distance
    my $mutdelcounter = 0;
    my %alterhash; #Stores the alterations to generate sequence 2; 1-based coordinates
    for (my $i = 0; $i < length $align1; $i++){     
        if ($chars1[$i] ne $chars2[$i]){
            print $i+1,"\t",$chars1[$i],"\t",$chars2[$i],"\n";
            $alterhash{$i+1} = $chars1[$i].'>'.$chars2[$i];
        }
        if ($chars1[$i] eq "-"){
            $wtdelcounter++;
        }
        if ($chars2[$i] eq "-"){
            $mutdelcounter++;
        }
    }

    print "wildtype deletion: $wtdelcounter\n";
    #determine the most upstream and most downstream positions with alterations
    my $minEditPos;
    my $maxEditPos;
    foreach my $k (sort keys %alterhash){
        if (!defined $minEditPos){
            $minEditPos = $k;
        }
        else {
            if ($k < $minEditPos){
                    $minEditPos = $k;
            }
        }
        if (!defined $maxEditPos){
            $maxEditPos = $k;
        }
        else {
            if ($k > $maxEditPos){
                    $maxEditPos = $k;
            }
        }
    }
    return($align1,$align2,$minEditPos,$maxEditPos,$trimmingStatus5p,$trimmingStatus3p,$wtdelcounter,\%alterhash,$mutdelcounter);


    #create a score matrix
    sub score_sub {
        if (!@_) {
            return -2; # gap penalty
        } 
        return ($_[0] eq $_[1]) ? 1 : -1;
    }
}
1;