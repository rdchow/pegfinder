#!/usr/bin/perl
##############################################################################
# Needleman-Wunsch  Algorithm 

use strict;
use warnings;
sub needleman_wunsch {
    my ($seq1, $seq2) = @_;

    # scoring scheme
    my $MATCH    =  1; # +1 for letters that match
    my $MISMATCH = -1; # -1 for letters that mismatch
    my $GAP      = -1; # -1 for any gap

    # initialization
    my @matrix;
    $matrix[0][0]{score}   = 0;
    $matrix[0][0]{pointer} = "none";
    for(my $j = 1; $j <= length($seq1); $j++) {
        $matrix[0][$j]{score}   = $GAP * $j;
        $matrix[0][$j]{pointer} = "left";
    }
    for (my $i = 1; $i <= length($seq2); $i++) {
        $matrix[$i][0]{score}   = $GAP * $i;
        $matrix[$i][0]{pointer} = "up";
    }

    # fill
    for(my $i = 1; $i <= length($seq2); $i++) {
        for(my $j = 1; $j <= length($seq1); $j++) {
            my ($diagonal_score, $left_score, $up_score);

            # calculate match score
            my $letter1 = substr($seq1, $j-1, 1);
            my $letter2 = substr($seq2, $i-1, 1);                            
            if ($letter1 eq $letter2) {
                $diagonal_score = $matrix[$i-1][$j-1]{score} + $MATCH;
            }
            else {
                $diagonal_score = $matrix[$i-1][$j-1]{score} + $MISMATCH;
            }

            # calculate gap scores
            $up_score   = $matrix[$i-1][$j]{score} + $GAP;
            $left_score = $matrix[$i][$j-1]{score} + $GAP;

            # choose best score
            if ($diagonal_score >= $up_score) {
                if ($diagonal_score >= $left_score) {
                    $matrix[$i][$j]{score}   = $diagonal_score;
                    $matrix[$i][$j]{pointer} = "diagonal";              
                }
                else {
                    $matrix[$i][$j]{score}   = $left_score;
                    $matrix[$i][$j]{pointer} = "left";
                }
            } 
            else {
                if ($up_score >= $left_score) {
                    $matrix[$i][$j]{score}   = $up_score;
                    $matrix[$i][$j]{pointer} = "up";
                }
                else {
                    $matrix[$i][$j]{score}   = $left_score;
                    $matrix[$i][$j]{pointer} = "left";
                }
            }
        }
    }

    # trace-back
    my $align1 = "";
    my $align2 = "";

    # start at last cell of matrix
    my $j = length($seq1);
    my $i = length($seq2);

    while (1) {
        last if $matrix[$i][$j]{pointer} eq "none"; # ends at first cell of matrix

        if ($matrix[$i][$j]{pointer} eq "diagonal") {
            $align1 .= substr($seq1, $j-1, 1);
            $align2 .= substr($seq2, $i-1, 1);
            $i--;
            $j--;
        }
        elsif ($matrix[$i][$j]{pointer} eq "left") {
            $align1 .= substr($seq1, $j-1, 1);
            $align2 .= "-";
            $j--;
        }
        elsif ($matrix[$i][$j]{pointer} eq "up") {
            $align1 .= "-";
            $align2 .= substr($seq2, $i-1, 1);
            $i--;
        }    
    }

    $align1 = reverse $align1; #seq1 alignment
    $align2 = reverse $align2; #seq2 alignment
   
    #Detect if the 5' or 3' ends of the two sequences are not aligned exactly
    #first look at 5' ends
    my $trimmingStatus5p = 0; 
    # ^ variable that tracks 5' trimming needs. 0: no trimming, 1:trimming wildtype seq, 2: trimming edited seq
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
    # ^ variable that tracks 5' trimming needs. 0: no trimming, 1:trimming wildtype seq, 2: trimming edited seq
    if ($chars1[-1] eq "-"){
        $trimmingStatus3p = 1;
    }
    elsif ($chars2[-1] eq "-"){
        $trimmingStatus3p = 2;
    }
    my $wtdelcounter = 0; # counts how many gaps are in the wildtype alignment; need for calculating the edit distance
    #look at the alterations
    my %alterhash; #Stores the alterations to generate sequence 2; 1-based coordinates
    for (my $i = 0; $i < length $align1; $i++){     
        if ($chars1[$i] ne $chars2[$i]){
            print $i+1,"\t",$chars1[$i],"\t",$chars2[$i],"\n";
            $alterhash{$i+1} = $chars1[$i].'>'.$chars2[$i];
        }
        if ($chars1[$i] eq "-"){
            $wtdelcounter++;
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
    return($align1,$align2,$minEditPos,$maxEditPos,$trimmingStatus5p,$trimmingStatus3p,$wtdelcounter);
}
1;