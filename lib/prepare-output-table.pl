#!/usr/bin/perl
use strict;
use warnings;

#$|=2;

sub prep_table_multi { #prepare output table with alternative 3' extensions + oligos
    my ($ranksghashpt,$align2,$minEditPos,$maxEditPos,$wtdelcounter,$maxSgCt,$rgn) = @_;
    my %ranksghash = %$ranksghashpt;
    # sgRNA seq, NA or OT score, cut position, orientation, distance, disrupt, GC % 
    
    #output text
    my $fileout = "RT_len\tRT_seq\tRT_picked\tPBS_len\tPBS_seq\tPBS_picked\t3'_extension_seq\t3'_extension_picked\textensF_oligo\textensR_oligo\tsgRNA_seq\tsgRNA_rank\tsgF_oligo\tsgR_oligo\tsg_Orientation\tsg_Seed/PAM_disrupt\tsg_GC\%\tsg_OnTargetScore\tEnzyme\n";

    # for the top sgRNAs in $ranksghash, identify the candidate RT templates and PBS sequences
    for (my $i = 0; $i < $maxSgCt; $i++){
        my $chosenSG = $ranksghash{$i}[0];
        my $chosenOrientation = $ranksghash{$i}[3];
        my $chosenCutPos = $ranksghash{$i}[2];
        my $gcPctg = $ranksghash{$i}[6];

        #Design the RT templates
        my @rtdata = find_RT($align2,$chosenOrientation,$minEditPos,$maxEditPos,$chosenCutPos,$wtdelcounter);
        my ($rthashpt,$chosenRT,$chosenRTlen,$rttable) = @rtdata;
        my %rthash = %$rthashpt;

        #Design the PBS templates
        my @pbsdata = find_pbs($chosenSG,$gcPctg);
        my ($pbshashpt,$chosenPBS,$chosenPBSlen,$pbstable) = @pbsdata;
        my %pbshash = %$pbshashpt;


        foreach my $k (sort {$a<=>$b} keys %rthash){
            foreach my $j (sort {$a<=>$b} keys %pbshash){
                $fileout .= "$k\t$rthash{$k}\t";

                #RT
                if ($k eq $chosenRTlen){
                    $fileout .= "1\t";
                }
                else {
                    $fileout .= "0\t";
                }
                $fileout .= "$j\t$pbshash{$j}\t";

                #PBS
                if ($j eq $chosenPBSlen){
                    $fileout .= "1\t";
                }
                else{
                    $fileout .= "0\t";
                }
                
                my $ext = $rthash{$k}.$pbshash{$j};
                my $rc_ext = reverse_complement($ext);
                my $rank = $i+1;
                $fileout .= $ext."\t";

                if ($k eq $chosenRTlen && $j eq $chosenPBSlen){
                    $fileout .= "1";
                }
                else {
                    $fileout .= "0";
                }
                $fileout .= "\tgtgc".$ext."\taaaa".$rc_ext."\t".$chosenSG."\t".$rank."\t";
                
                if (substr($chosenSG,0,1) eq "G"){
                    $fileout .= "cacc".$chosenSG."gttttaga\t";
                    $fileout .= "tagctctaaaac".reverse_complement($chosenSG)."\t";
                }
                else {
                    $fileout .= "caccg".$chosenSG."gttttaga\t";
                    $fileout .= "tagctctaaaac".reverse_complement($chosenSG)."\t";      
                }
                $fileout .= $chosenOrientation."\t".$ranksghash{$i}[5]."\t".$ranksghash{$i}[6]."\t".$ranksghash{$i}[1]."\t";
                $fileout .= $rgn."\n";
            }
        }
    }
    return($fileout,1);
}

sub prep_table_chosen{
    my ($rthashpt,$pbshashpt,$chosenRTlen,$chosenPBSlen,$chosenSG,$maxSgCt,$rgn,$chosenOrientation,$gcPctg,$chosenDisrupt) = @_;
    my %rthash = %$rthashpt;
    my %pbshash = %$pbshashpt;
    my $fileout = "RT_len\tRT_seq\tRT_picked\tPBS_len\tPBS_seq\tPBS_picked\t3'_extension_seq\t3'_extension_picked\textensF_oligo\textensR_oligo\tsgRNA_seq\tsgRNA_rank\tsgF_oligo\tsgR_oligo\tsg_Orientation\tsg_Seed/PAM_disrupt\tsg_GC\%\tEnzyme\n";
    
    foreach my $k (sort {$a<=>$b} keys %rthash){
        foreach my $j (sort {$a<=>$b} keys %pbshash){
            $fileout .= "$k\t$rthash{$k}\t";

            #RT
            if ($k eq $chosenRTlen){
                $fileout .= "1\t";
            }
            else {
                $fileout .= "0\t";
            }
            $fileout .= "$j\t$pbshash{$j}\t";

            #PBS
            if ($j eq $chosenPBSlen){
                $fileout .= "1\t";
            }
            else{
                $fileout .= "0\t";
            }
                
            my $ext = $rthash{$k}.$pbshash{$j};
            my $rc_ext = reverse_complement($ext);
            my $rank = "preselected";
            $fileout .= $ext."\t";

            if ($k eq $chosenRTlen && $j eq $chosenPBSlen){
                $fileout .= "1";
            }
            else {
                $fileout .= "0";
            }
            $fileout .= "\tgtgc".$ext."\taaaa".$rc_ext."\t".$chosenSG."\t".$rank."\t";
                
            if (substr($chosenSG,0,1) eq "G"){
                $fileout .= "cacc".$chosenSG."gttttaga\t";
                $fileout .= "tagctctaaaac".reverse_complement($chosenSG)."\t";
            }
            else {
                $fileout .= "caccg".$chosenSG."gttttaga\t";
                $fileout .= "tagctctaaaac".reverse_complement($chosenSG)."\t";      
            }
            $fileout .= $chosenOrientation."\t".$chosenDisrupt."\t".$gcPctg."\t";
            $fileout .= $rgn."\n";
        }
    }
    return($fileout,1);
}
1;