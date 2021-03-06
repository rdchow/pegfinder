#!/usr/bin/perl
use strict;
use warnings;

#$|=2;
sub process_broad_sgRNA { #process the Broad sgRNA file
    my ($data_text,$minEditPos,$maxEditPos,$maxEditDistance,$wtdelcounter) = @_;
    $data_text =~ s/"//g;
    my @rows = split (/[\r\n]+/,$data_text);
    my $head = $rows[0];
    my @headers = split ("\t",$head);

    #save the column index of relevant columns
    my $cutPosCol;
    my $sgRNACol;
    my $orientationCol;
    my $onTargetScoreCol;
    my $offTargetT1Col;
    my $offTargetT2Col;
    my $offTargetT3Col;
    for (my $i = 0; $i < scalar @headers; $i++){
        if ($headers[$i] eq "sgRNA Cut Position (1-based)"){
            $cutPosCol = $i;
        }
        if ($headers[$i] eq "sgRNA Sequence"){
            $sgRNACol = $i;
        }
        if ($headers[$i] eq "Orientation"){
            $orientationCol = $i;
        }
        if ($headers[$i] eq "On-Target Efficacy Score"){
            $onTargetScoreCol = $i;
        }
        if ($headers[$i] eq "# Off-Target Tier I Match Bin I Matches"){
            $offTargetT1Col = $i;
        }
        if ($headers[$i] eq "# Off-Target Tier II Match Bin I Matches"){
            $offTargetT2Col = $i;
        }
        if ($headers[$i] eq "# Off-Target Tier III Match Bin I Matches"){
            $offTargetT3Col = $i;
        }
    }
    #Save the sgRNAs that have the potential to create the desired edit
    #Considering that modifications can range from +0 to +$maxEditDistance of the sgRNA cut site
    my %sghash; # stores the candidate sgRNAs for the editing sgRNA

    for (my $i = 1; $i < scalar @rows; $i++){
       # print $orientationCol,"\n\n";
        my $line = $rows[$i];
        my @lines = split ("\t",$line);
        my $cutPos = $lines[$cutPosCol];
        my $sgRNA = $lines[$sgRNACol];
        my $orientation = $lines[$orientationCol];
        my $onTargetScore = $lines[$onTargetScoreCol];
        my $offTargetT1 = $lines[$offTargetT1Col];
        my $offTargetT2 = $lines[$offTargetT2Col];
        my $offTargetT3 = $lines[$offTargetT3Col];

        #Sense strand sgRNAs:
        if ($orientation eq "sense"){
            print "$cutPos\t$minEditPos\t$maxEditPos\n";
            if ($minEditPos >= $cutPos && $maxEditPos <= ($cutPos+$maxEditDistance) && (($offTargetT1 + $offTargetT2 + $offTargetT3) <= 1)){
                my $distance = $minEditPos-$cutPos;
                $sghash{$sgRNA} = [$onTargetScore,$cutPos,$orientation,$distance];
            }
        }
        #Antisense strand sgRNAs:
        elsif ($orientation eq "antisense"){
            $cutPos--;
            print "$cutPos\t$minEditPos\t$maxEditPos\n";
            if ($minEditPos >= ($cutPos-$maxEditDistance) && $maxEditPos <= ($cutPos+$wtdelcounter) && (($offTargetT1 + $offTargetT2 + $offTargetT3) <= 1)){
                my $distance = ($cutPos-$maxEditPos)+$wtdelcounter; # gaps in the wildtype alignment will push the max edit position closer to the cut site; need to add back in those bases 
                $sghash{$sgRNA} = [$onTargetScore,$cutPos,$orientation,$distance];
            }
        }
    }
    if (scalar keys %sghash > 0){
        print "Found ".scalar(keys %sghash)." candidate sgRNAs that are within $maxEditDistance nt of the edit start.\n";
        # Rank/filter the sgRNAs
        # Strategy of picking:
        # 1) Among sgRNAs with OnTarget >= 0.5, choose sgRNA with shortest distance to edit start. If tied for shortest distance, pick sgRNA with highest OT score
        # 2) If no sgRNAs with OnTarget >= 0.5, choose sgRNA with highest OnTarget score regardless of distance
        my $chosenSG;

        my $OTcheck = 0;
        foreach my $k (keys %sghash){
            if ($sghash{$k}[0] >= 0.5){
                $OTcheck = 1;
            }
        }
        # If there is an sgRNA with OnTarget >= 0.5
        if ($OTcheck == 1){
            my $currentMaxOT = 0;
            my $currentMinDist = 999999999999999;

            foreach my $k (keys %sghash){
                if ($sghash{$k}[0] >= 0.5){ # considering the sgRNAs with OnTarget >= 0.5
                    if ($sghash{$k}[3] < $currentMinDist){ # if shorter distance than previous best sgRNA
                        $chosenSG = $k;
                        $currentMinDist = $sghash{$k}[3];
                        $currentMaxOT = $sghash{$k}[0];
                    }
                    elsif ($sghash{$k}[3] == $currentMinDist){ # if same distance as previous best sgRNA
                        if ($sghash{$k}[0] > $currentMaxOT){ # if current sgRNA has higher OT score
                            $chosenSG = $k;
                            $currentMinDist = $sghash{$k}[3];
                            $currentMaxOT = $sghash{$k}[0];
                        }
                    }
                }
            }
        }

        # If no sgRNA with OnTarget >= 0.5
        elsif ($OTcheck == 0){
            my $currentMaxOT = 0;

            foreach my $k (keys %sghash){
                if ($sghash{$k}[0] > $currentMaxOT){ # if current sgRNA has higher OT score
                    $chosenSG = $k;
                    $currentMaxOT = $sghash{$k}[0];
                }
            }
        }

        #calculate GC content of the chosen sgRNA
        my @sgbases = split ("",$chosenSG);
        my $gcCounter = 0;
        for (my $i = 0; $i < scalar @sgbases; $i++){
            if ($sgbases[$i] eq "G" || $sgbases[$i] eq "C"){
                $gcCounter++;
            }
        }
        my $gcPctg = $gcCounter / scalar (@sgbases) * 100;

        print "Chose $chosenSG as the best sgRNA.\n";
        print "On Target Score: $sghash{$chosenSG}[0]\n";
        print "Cut Position: $sghash{$chosenSG}[1]\n";
        print "Orientation: $sghash{$chosenSG}[2]\n";
        print "Distance To Edit Start: $sghash{$chosenSG}[3]\n";
        print "G/C content: $gcPctg%\n\n";

        my $chosenCutPos = $sghash{$chosenSG}[1];
        my $chosenOrientation = $sghash{$chosenSG}[2];
        my $chosenDistance = $sghash{$chosenSG}[3];   
        return (\%sghash,$chosenSG,$chosenCutPos,$chosenOrientation,$chosenDistance,$gcPctg);
    }
    
    # if no candidate sgRNA found, return value stating none
    else {
        my $val = "none";
        return($val,$val);
    }
}

sub process_crisprscan_sgRNA {
    my ($data_text,$minEditPos,$maxEditPos,$maxEditDistance,$wtdelcounter,$seq1) = @_;
    $data_text =~ s/"//g;
    my @rows = split (/[\r\n]+/,$data_text);
    my $head = $rows[0];
    my @headers = split ("\t",$head);

    #save the column index of relevant columns
    my $sgRNACol;
    my $orientationCol;
    my $onTargetScoreCol;
    my $offTargetCol;
    for (my $i = 0; $i < scalar @headers; $i++){
        if ($headers[$i] eq "seq"){
            $sgRNACol = $i;
        }
        if ($headers[$i] eq "strand"){
            $orientationCol = $i;
        }
        if ($headers[$i] eq "score_doench"){
            $onTargetScoreCol = $i;
        }
        if ($headers[$i] eq "offtarget_number_seed"){
            $offTargetCol = $i;
        }
    }
    #Save the sgRNAs that have the potential to create the desired edit
    #Considering that modifications can range from 0 to +$maxEditDistance of the sgRNA cut site
    my %sghash; # stores the candidate sgRNAs for the editing sgRNA

    for (my $i = 1; $i < scalar @rows; $i++){
       # print $orientationCol,"\n\n";
        my $line = $rows[$i];
        my @lines = split ("\t",$line);
        my $cutPos;
        my $sgRNA = $lines[$sgRNACol];
        $sgRNA = substr($sgRNA,0,20);
        my $orientation = $lines[$orientationCol];
        if ($orientation eq "+"){$orientation = "sense"};
        if ($orientation eq "-"){$orientation = "antisense"};
        my $onTargetScore = $lines[$onTargetScoreCol];
        my $offTarget = $lines[$offTargetCol];
        my $distance;

        #determine the cut position relative to the wildtype sequence (convert to 1-based)
        if ($orientation eq "sense"){
            $seq1 =~ /$sgRNA.GG/;
            $cutPos = $-[0]+18;
            print "$sgRNA \t$cutPos\n";
            $distance = $minEditPos-$cutPos;
            if ($minEditPos >= $cutPos && $maxEditPos <= ($cutPos+$maxEditDistance) && $offTarget == 0){
                $sghash{$sgRNA} = [$onTargetScore,$cutPos,$orientation,$distance];
            }
        }
        elsif ($orientation eq "antisense"){
            my $rc_c_sg = reverse_complement($sgRNA);
            $seq1 =~ /CC.$rc_c_sg/;
            $cutPos = $-[0]+6;
            $distance = ($cutPos-$maxEditPos)+$wtdelcounter; # gaps in the wildtype alignment will push the max edit position closer to the cut site; need to add back in those bases
            if ($minEditPos >= ($cutPos-$maxEditDistance) && $maxEditPos <= ($cutPos+$wtdelcounter) && $offTarget == 0){
                $sghash{$sgRNA} = [$onTargetScore,$cutPos,$orientation,$distance];
            }
        }
    }
    if (scalar keys %sghash > 0){
        print "Found ".scalar(keys %sghash)." candidate sgRNAs that are within $maxEditDistance nt of the edit start.\n";
        # Rank/filter the sgRNAs
        # Strategy of picking:
        # 1) Among sgRNAs with OnTarget >= 25, choose sgRNA with shortest distance to edit start. If tied for shortest distance, pick sgRNA with highest OT score
        # 2) If no sgRNAs with OnTarget >= 25, choose sgRNA with highest OnTarget score regardless of distance
        my $chosenSG;

        my $OTcheck = 0;
        foreach my $k (keys %sghash){
            if ($sghash{$k}[0] >= 25){
                $OTcheck = 1;
            }
        }
        # If there is an sgRNA with OnTarget >= 25
        if ($OTcheck == 1){
            my $currentMaxOT = 0;
            my $currentMinDist = 999999999999999;

            foreach my $k (keys %sghash){
                if ($sghash{$k}[0] >= 25){ # considering the sgRNAs with OnTarget >= 25
                    if ($sghash{$k}[3] < $currentMinDist){ # if shorter distance than previous best sgRNA
                        $chosenSG = $k;
                        $currentMinDist = $sghash{$k}[3];
                        $currentMaxOT = $sghash{$k}[0];
                    }
                    elsif ($sghash{$k}[3] == $currentMinDist){ # if same distance as previous best sgRNA
                        if ($sghash{$k}[0] > $currentMaxOT){ # if current sgRNA has higher OT score
                            $chosenSG = $k;
                            $currentMinDist = $sghash{$k}[3];
                            $currentMaxOT = $sghash{$k}[0];
                        }
                    }
                }
            }
        }

        # If no sgRNA with OnTarget >= 0.5
        elsif ($OTcheck == 0){
            my $currentMaxOT = 0;

            foreach my $k (keys %sghash){
                if ($sghash{$k}[0] > $currentMaxOT){ # if current sgRNA has higher OT score
                    $chosenSG = $k;
                    $currentMaxOT = $sghash{$k}[0];
                }
            }
        }

        #calculate GC content of the chosen sgRNA
        my @sgbases = split ("",$chosenSG);
        my $gcCounter = 0;
        for (my $i = 0; $i < scalar @sgbases; $i++){
            if ($sgbases[$i] eq "G" || $sgbases[$i] eq "C"){
                $gcCounter++;
            }
        }
        my $gcPctg = $gcCounter / scalar (@sgbases) * 100;

        print "Chose $chosenSG as the best sgRNA.\n";
        print "On Target Score: $sghash{$chosenSG}[0]\n";
        print "Cut Position: $sghash{$chosenSG}[1]\n";
        print "Orientation: $sghash{$chosenSG}[2]\n";
        print "Distance To Edit Start: $sghash{$chosenSG}[3]\n";
        print "G/C content: $gcPctg%\n\n";

        my $chosenCutPos = $sghash{$chosenSG}[1];
        my $chosenOrientation = $sghash{$chosenSG}[2];
        my $chosenDistance = $sghash{$chosenSG}[3];   
        return (\%sghash,$chosenSG,$chosenCutPos,$chosenOrientation,$chosenDistance,$gcPctg);
    }
    
    # if no candidate sgRNA found, return value stating none
    else {
        my $val = "none";
        return($val,$val);
    }
}

sub find_broad_nicksgRNA {
    my ($data_text,$minEditPos,$maxEditPos,$maxEditDistance,$chosenCutPos,$chosenOrientation) = @_;
    $data_text =~ s/"//g;
    my @rows = split (/[\r\n]+/,$data_text);

    my $head = $rows[0];
    my @headers = split ("\t",$head);

    #save the column index of relevant columns
    my $cutPosCol;
    my $sgRNACol;
    my $orientationCol;
    my $onTargetScoreCol;
    my $offTargetT1Col;
    my $offTargetT2Col;
    my $offTargetT3Col;
    for (my $i = 0; $i < scalar @headers; $i++){
        if ($headers[$i] eq "sgRNA Cut Position (1-based)"){
            $cutPosCol = $i;
        }
        if ($headers[$i] eq "sgRNA Sequence"){
            $sgRNACol = $i;
        }
        if ($headers[$i] eq "Orientation"){
            $orientationCol = $i;
        }
        if ($headers[$i] eq "On-Target Efficacy Score"){
            $onTargetScoreCol = $i;
        }
        if ($headers[$i] eq "# Off-Target Tier I Match Bin I Matches"){
            $offTargetT1Col = $i;
        }
        if ($headers[$i] eq "# Off-Target Tier II Match Bin I Matches"){
            $offTargetT2Col = $i;
        }
        if ($headers[$i] eq "# Off-Target Tier III Match Bin I Matches"){
            $offTargetT3Col = $i;
        }
    }
    #Save all sgRNAs
    my %sghash; # stores the candidate sgRNAs for the editing sgRNA
    my %allsghash; # stores all identified sgRNAs, used for finding PE3 nick sgRNA
    for (my $i = 1; $i < scalar @rows; $i++){
       # print $orientationCol,"\n\n";
        my $line = $rows[$i];
        my @lines = split ("\t",$line);
        my $cutPos = $lines[$cutPosCol];
        my $sgRNA = $lines[$sgRNACol];
        my $orientation = $lines[$orientationCol];
        my $onTargetScore = $lines[$onTargetScoreCol];
        my $offTargetT1 = $lines[$offTargetT1Col];
        my $offTargetT2 = $lines[$offTargetT2Col];
        my $offTargetT3 = $lines[$offTargetT3Col];
        $allsghash{$sgRNA} = [$onTargetScore,$cutPos,$orientation,$offTargetT1,$offTargetT2,$offTargetT3];
    }

    #Now design the 2nd nicking sgRNA for PE3
    if (scalar keys %allsghash > 0){
        my $chosenNickSG = ""; #final selection for nicking sgRNA for PE3
        my $chosenNickSGPos = ""; #position of the chosen nicking sgRNA for PE
        my $chosenNickSGDist = ""; #distance between PE3 nick and the editing sgRNA cut
        my $chosenNickOrientation = ""; #the orientation of the PE3 secondary nicking sgRNA
        my %nicksghash; #hash for the potential candidates for nicking sgRNA
        foreach my $k (keys %allsghash){ #[$onTargetScore,$cutPos,$orientation,$offTargetT1];
            my $onTargetScore = $allsghash{$k}[0];
            my $cutPos = $allsghash{$k}[1];
            my $orientation = $allsghash{$k}[2];
            my $offTargetT1 = $allsghash{$k}[3];
            my $offTargetT2 = $allsghash{$k}[4];
            my $offTargetT3 = $allsghash{$k}[5];

            #only consider nicking sgRNA if on the opposite strand and 20-100 nt away from the chosen cut site
            if ($chosenOrientation eq "sense"){
                $cutPos--;
                my $dist = $cutPos-$chosenCutPos+1;
                if ($orientation eq "antisense" && (($offTargetT1 + $offTargetT2 + $offTargetT3) <= 1) && (abs($dist)>=20) && (abs($dist) <= 100) ){
                    $nicksghash{$k} = [$onTargetScore,$cutPos,$orientation,$dist];
                }
            }
            elsif ($chosenOrientation eq "antisense"){
                my $dist = $cutPos-$chosenCutPos-1;
                if ($orientation eq "sense" && (($offTargetT1 + $offTargetT2 + $offTargetT3) <= 1) && (abs($dist)>=20) && (abs($dist) <= 100) ){
                    $nicksghash{$k} = [$onTargetScore,$cutPos,$orientation,$dist];
                }
            }
        }
        if (scalar keys %nicksghash > 0){
        #choose the highest On target scoring nicking sgRNA
            my $ctt = 0;
            foreach my $k (sort {$nicksghash{$b}[0] <=> $nicksghash{$a}[0]} keys %nicksghash){
                if ($ctt == 0){
                    $chosenNickSG = $k;
                    $chosenNickSGPos = $nicksghash{$k}[1];
                    $chosenNickOrientation = $nicksghash{$k}[2];
                    $chosenNickSGDist = $nicksghash{$k}[3];
                }
                $ctt++;
            }
                
            return ($chosenNickSG,$chosenNickSGPos,$chosenNickOrientation,$chosenNickSGDist,\%nicksghash);
        }
        
        # if no candidate nicking sgRNA found, return value stating none
        else {
            return("none found","none");
        }
    }
}

sub find_crisprscan_nicksgRNA {
    my ($data_text,$minEditPos,$maxEditPos,$maxEditDistance,$chosenCutPos,$chosenOrientation,$seq1) = @_;
    $data_text =~ s/"//g;
    my @rows = split (/[\r\n]+/,$data_text);

    my $head = $rows[0];
    my @headers = split ("\t",$head);

    #save the column index of relevant columns
    my $sgRNACol;
    my $orientationCol;
    my $onTargetScoreCol;
    my $offTargetCol;
    for (my $i = 0; $i < scalar @headers; $i++){
        if ($headers[$i] eq "seq"){
            $sgRNACol = $i;
        }
        if ($headers[$i] eq "strand"){
            $orientationCol = $i;
        }
        if ($headers[$i] eq "score_doench"){
            $onTargetScoreCol = $i;
        }
        if ($headers[$i] eq "offtarget_number_seed"){
            $offTargetCol = $i;
        }
    }
    #Save all sgRNAs
    my %allsghash; # stores all identified sgRNAs, used for finding PE3 nick sgRNA
    for (my $i = 1; $i < scalar @rows; $i++){
       # print $orientationCol,"\n\n";
        my $line = $rows[$i];
        my @lines = split ("\t",$line);
        my $cutPos;
        my $sgRNA = $lines[$sgRNACol];
        $sgRNA = substr($sgRNA,0,20);
        my $orientation = $lines[$orientationCol];
        if ($orientation eq "+"){$orientation = "sense"};
        if ($orientation eq "-"){$orientation = "antisense"};
        my $onTargetScore = $lines[$onTargetScoreCol];
        my $offTarget = $lines[$offTargetCol];

        #determine the cut position relative to the wildtype sequence (convert to 1-based)
        if ($orientation eq "sense"){
            $seq1 =~ /$sgRNA.GG/;
            $cutPos = $-[0]+18;
        }
        elsif ($orientation eq "antisense"){
            my $rc_c_sg = reverse_complement($sgRNA);
            $seq1 =~ /CC.$rc_c_sg/;
            $cutPos = $-[0]+6;
        }
        $allsghash{$sgRNA} = [$onTargetScore,$cutPos,$orientation,$offTarget];
    }

    #Now design the 2nd nicking sgRNA for PE3
    if (scalar keys %allsghash > 0){
        my $chosenNickSG = ""; #final selection for nicking sgRNA for PE3
        my $chosenNickSGPos = ""; #position of the chosen nicking sgRNA for PE
        my $chosenNickSGDist = ""; #distance between PE3 nick and the editing sgRNA cut
        my $chosenNickOrientation = ""; #the orientation of the PE3 secondary nicking sgRNA
        my %nicksghash; #hash for the potential candidates for nicking sgRNA
        foreach my $k (keys %allsghash){ #[$onTargetScore,$cutPos,$orientation,$offTargetT1];
            my $onTargetScore = $allsghash{$k}[0];
            my $cutPos = $allsghash{$k}[1];
            my $orientation = $allsghash{$k}[2];
            my $offTarget = $allsghash{$k}[3];

            #only consider nicking sgRNA if on the opposite strand and 20-100 nt away from the chosen cut site
            if ($chosenOrientation eq "sense"){
                $cutPos--;
                my $dist = $cutPos-$chosenCutPos+1;
                if ($orientation eq "antisense" && $offTarget == 0 && (abs($dist)>=20) && (abs($dist) <= 100) ){
                    $nicksghash{$k} = [$onTargetScore,$cutPos,$orientation,$dist];
                }
            }
            elsif ($chosenOrientation eq "antisense"){
                my $dist = $cutPos-$chosenCutPos-1;
                if ($orientation eq "sense" && $offTarget == 0 && (abs($dist)>=20) && (abs($dist) <= 100) ){
                    $nicksghash{$k} = [$onTargetScore,$cutPos,$orientation,$dist];
                }
            }
        }
        if (scalar keys %nicksghash > 0){
        #choose the highest On target scoring nicking sgRNA
            my $ctt = 0;
            foreach my $k (sort {$nicksghash{$b}[0] <=> $nicksghash{$a}[0]} keys %nicksghash){
                if ($ctt == 0){
                    $chosenNickSG = $k;
                    $chosenNickSGPos = $nicksghash{$k}[1];
                    $chosenNickOrientation = $nicksghash{$k}[2];
                    $chosenNickSGDist = $nicksghash{$k}[3];
                }
                $ctt++;
            }
                
            return ($chosenNickSG,$chosenNickSGPos,$chosenNickOrientation,$chosenNickSGDist,\%nicksghash);
        }
        
        # if no candidate nicking sgRNA found, return value stating none
        else {
            return("none found","none");
        }
    }
}

#If an sgRNA was preselected, calculate the needed metrics (Cut position, orientation, distance to edit, GC %)
#take the chosen sgRNA sequence, the wildtype DNA sequence, minEditPos, maxEditPos, and maxEditDistance 
sub process_chosen_sgRNA {
    my ($c_sg,$seq1,$minEditPos,$maxEditPos,$maxEditDistance,$wtdelcounter) = @_;
    my $chosenOrientation;
    my $chosenSG = $c_sg;
    my $chosenCutPos;
    my $chosenDistance;
    my $gcPctg;

    #Check if repetitive sequence or palindrome targeted by chosen sgRNA
    my $counter = 0;
    $counter += () = $seq1 =~ /$c_sg.GG/;
    $counter += () = (reverse_complement($seq1) =~ /$c_sg.GG/);
    #print $counter,"\n";
    if ($counter > 1){
        return ("non-unique","non-unique");
    }
    
    else {
        #determine orientation
        if ($seq1 =~ /$c_sg.GG/){
            $chosenOrientation = "sense";
        }
        elsif (reverse_complement($seq1) =~ /$c_sg.GG/){
            $chosenOrientation = "antisense";
        }

        #calculate GC content of the chosen sgRNA
        my @sgbases = split ("",$chosenSG);
        my $gcCounter = 0;
        for (my $i = 0; $i < scalar @sgbases; $i++){
            if ($sgbases[$i] eq "G" || $sgbases[$i] eq "C"){
                $gcCounter++;
            }
        }
        $gcPctg = $gcCounter / scalar (@sgbases) * 100;

        #determine the cut position relative to the wildtype sequence (convert to 1-based)
        if ($chosenOrientation eq "sense"){
            $seq1 =~ /$c_sg.GG/;
            $chosenCutPos = $-[0]+18;
            $chosenDistance = $minEditPos-$chosenCutPos;
        }
        elsif ($chosenOrientation eq "antisense"){
            my $rc_c_sg = reverse_complement($c_sg);
            $seq1 =~ /CC.$rc_c_sg/;
            $chosenCutPos = $-[0]+6;
            $chosenDistance = ($chosenCutPos-$maxEditPos)+$wtdelcounter; # gaps in the wildtype alignment will push the max edit position closer to the cut site; need to add back in those bases
        }

        print "Using $chosenSG as the preselected sgRNA.\n";
        print "Cut Position: $chosenCutPos\n";
        print "Orientation: $chosenOrientation\n";
        print "Distance To Edit Start: $chosenDistance\n";
        print "G/C content: $gcPctg%\n\n";
        my $sgfoundStatus;
        if ($chosenDistance >=0 && $chosenDistance <= $maxEditDistance){
            $sgfoundStatus = 2;
        }
        elsif ($chosenDistance >=0 && $chosenDistance > $maxEditDistance){
            $sgfoundStatus = 1;
        }
        elsif ($chosenDistance <0){
            $sgfoundStatus = 0;
        }
        return ($chosenCutPos,$chosenOrientation,$chosenDistance,$gcPctg,$sgfoundStatus);
    }
}

sub find_choose_sgRNA {
    my ($seq1,$minEditPos,$maxEditPos,$maxEditDistance, $wtdelcounter) = @_;
    my $chosenOrientation;
    my $chosenSG;
    my $chosenCutPos;
    my $chosenDistance;
    my $gcPctg;
    my %sghash;

    #find all candidate spacers in the wildtype sequence
    my @chars = split ("",$seq1);

    #first check sense orientation:
    for (my $i = 21; $i < scalar @chars; $i++){
        my $tempsg;
        if ($chars[$i] eq "G" && $chars[$i+1] eq "G"){
            for (my $k = $i-21; $k <= $i-2; $k++){
                $tempsg .= $chars[$k];
            }
            
            #calculate GC content of the sgRNA
            my @sgbases = split ("",$tempsg);
            my $gcCounter = 0;
            for (my $j = 0; $j < scalar @sgbases; $j++){
                if ($sgbases[$j] eq "G" || $sgbases[$j] eq "C"){
                    $gcCounter++;
                }
            }
            my $tempgcPctg = $gcCounter / scalar (@sgbases) * 100;

            #determine cut position and distance to edit
            $seq1 =~ /$tempsg.GG/;
            my $tempCutPos = $-[0]+18;
            my $tempDistance = $minEditPos-$tempCutPos;

            if ($minEditPos >= $tempCutPos && $maxEditPos <= ($tempCutPos+$maxEditDistance)){
                $sghash{$tempsg} = [$tempsg, "sense", $tempgcPctg,$tempCutPos,$tempDistance];
            }
        }
    }

    #now check antisense orientation:
    for (my $i = 0; $i < scalar @chars-22; $i++){
        my $tempsg;
        if ($chars[$i] eq "C" && $chars[$i+1] eq "C"){
            for (my $k = $i+3; $k <= $i+22; $k++){
                $tempsg .= $chars[$k];
            }
            $tempsg = reverse_complement($tempsg);
            #calculate GC content of the sgRNA
            my @sgbases = split ("",$tempsg);
            my $gcCounter = 0;
            for (my $j = 0; $j < scalar @sgbases; $j++){
                if ($sgbases[$j] eq "G" || $sgbases[$j] eq "C"){
                    $gcCounter++;
                }
            }
            my $tempgcPctg = $gcCounter / scalar (@sgbases) * 100;

            #determine cut position and distance to edit
            my $rc_tempsg = reverse_complement($tempsg);
            $seq1 =~ /CC.$rc_tempsg/;
            my $tempCutPos = $-[0]+6;
            my $tempDistance = ($tempCutPos-$maxEditPos)+$wtdelcounter; # gaps in the wildtype alignment will push the max edit position closer to the cut site; need to add back in those bases
            if ($minEditPos >= ($tempCutPos-$maxEditDistance) && $maxEditPos <= ($tempCutPos+$wtdelcounter)){
                $sghash{$tempsg} = [$tempsg, "antisense", $tempgcPctg,$tempCutPos,$tempDistance];
            }
        }
    }

    # %sghash holds all identified sgRNAs that are candidates for the desired edit.
    # choose a single sgRNA that is closest to the edit site
    my $sgfoundStatus = 0;
    if (scalar keys %sghash > 0){
        foreach my $k (sort { abs($sghash{$a}[4]) <=> abs($sghash{$b}[4]) } keys %sghash){
            if ($sgfoundStatus == 0){
                $chosenSG = $k;
                $chosenCutPos = $sghash{$k}[3];
                $chosenDistance = $sghash{$k}[4];
                $chosenOrientation = $sghash{$k}[1];
                $gcPctg = $sghash{$k}[2];
                $sgfoundStatus++;
            }
        }
        return (\%sghash,$chosenSG,$chosenCutPos,$chosenOrientation,$chosenDistance,$gcPctg);
    }
    # if no candidate sgRNA found, return value stating none
    else {
        my $val = "none";
        return($val,$val);
    }
}


# if no sgRNA designer file uploaded, design nicking sgRNA
sub find_choose_nick_sgRNA{
    my ($seq1,$minEditPos,$maxEditPos,$maxEditDistance,$chosenCutPos,$chosenOrientation) = @_;
    my $chosenNickOrientation;
    my $chosenNickSG;
    my $chosenNickSGPos;
    my $chosenNickSGDist;
    my %nicksghash;

    #find all candidate spacers in the wildtype sequence
    my @chars = split ("",$seq1);

    #if chosen sgRNA is antisense, look for sense sgRNAs for PE3:
    if ($chosenOrientation eq "antisense"){
        for (my $i = 21; $i < scalar @chars; $i++){
            my $tempsg;
            if ($chars[$i] eq "G" && $chars[$i+1] eq "G"){
                for (my $k = $i-21; $k <= $i-2; $k++){
                    $tempsg .= $chars[$k];
                }
                #determine cut position and distance to primary sgRNA cut
                $seq1 =~ /$tempsg.GG/;
                my $tempCutPos = $-[0]+18;
                my $tempDistance = $chosenCutPos-$tempCutPos-1;

                if ((abs($tempDistance)>=20) && (abs($tempDistance) <= 100)){
                    $nicksghash{$tempsg} = [$tempsg, $tempCutPos,"sense", $tempDistance];
                }
            }
        }
    }
    #else, if chosen sgRNA is sense, look for antisense sgRNAs for PE3:
    if ($chosenOrientation eq "sense"){
    #now check antisense orientation:
    for (my $i = 0; $i < scalar @chars-22; $i++){
        my $tempsg;
        if ($chars[$i] eq "C" && $chars[$i+1] eq "C"){
            for (my $k = $i+3; $k <= $i+22; $k++){
                $tempsg .= $chars[$k];
            }
            $tempsg = reverse_complement($tempsg);
        
            #determine cut position and distance to primary sgRNA cut
            my $rc_tempsg = reverse_complement($tempsg);
            $seq1 =~ /CC.$rc_tempsg/;
            my $tempCutPos = $-[0]+6;
            ####TEST LINE###
            my $tempDistance = $tempCutPos-$chosenCutPos+1;
            if ((abs($tempDistance)>=20) && (abs($tempDistance) <= 100) ){
                $nicksghash{$tempsg} = [$tempsg,$tempCutPos,"antisense",$tempDistance];
            }
        }
    }

    if (scalar keys %nicksghash > 0){
        #choose the nicking sgRNA that is closest to 50bp from the primary nick
            my $ctt = 0;
            foreach my $k (sort {abs($nicksghash{$a}[3]-50) <=> abs($nicksghash{$b}[3]-50)} keys %nicksghash){
                if ($ctt == 0){
                    $chosenNickSG = $k;
                    $chosenNickSGPos = $nicksghash{$k}[1];
                    $chosenNickOrientation = $nicksghash{$k}[2];
                    $chosenNickSGDist = $nicksghash{$k}[3];
                }
                $ctt++;
            }
                
            return ($chosenNickSG,$chosenNickSGPos,$chosenNickOrientation,$chosenNickSGDist,\%nicksghash);
        }
        
        # if no candidate nicking sgRNA found, return value stating none
        else {
            return("none found","none");
        }
    }
}


sub reverse_complement {
        my $dna = shift;

        # reverse the DNA sequence
        my $revcomp = reverse($dna);

        # complement the reversed DNA sequence
        $revcomp =~ tr/TCGA/AGCT/;
        return $revcomp;
}

1;
