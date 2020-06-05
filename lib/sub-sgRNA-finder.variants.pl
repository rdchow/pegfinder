#!/usr/bin/perl
use strict;
use warnings;

#$|=2;

#############################
#NG PAM
sub find_NG_sgRNA {
    my ($seq1,$minEditPos,$maxEditPos,$maxEditDistance, $wtdelcounter,$seq2) = @_;
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
        if ($chars[$i] eq "G"){
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
            $seq1 =~ /$tempsg.G/;
            my $tempCutPos = $-[0]+18;
            my $tempDistance = $minEditPos-$tempCutPos;

            if ($minEditPos >= $tempCutPos && $maxEditPos <= ($tempCutPos+$maxEditDistance) && $tempsg !~ /TTTTT/){
                $sghash{$tempsg} = [$tempsg, "sense", $tempgcPctg,$tempCutPos,$tempDistance];
            }
        }
    }

    #now check antisense orientation:
    for (my $i = 0; $i < scalar @chars-21; $i++){
        my $tempsg;
        if ($chars[$i] eq "C"){
            for (my $k = $i+2; $k <= $i+21; $k++){
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
            $seq1 =~ /C.$rc_tempsg/;
            my $tempCutPos = $-[0]+5;
            my $tempDistance = ($tempCutPos-$maxEditPos)+$wtdelcounter; # gaps in the wildtype alignment will push the max edit position closer to the cut site; need to add back in those bases
            if ($minEditPos >= ($tempCutPos-$maxEditDistance) && $maxEditPos <= ($tempCutPos+$wtdelcounter) && $tempsg !~ /TTTTT/){
                $sghash{$tempsg} = [$tempsg, "antisense", $tempgcPctg,$tempCutPos,$tempDistance];
            }
        }
    }

    # %sghash holds all identified sgRNAs that are candidates for the desired edit.
    # choose a single sgRNA that disrupts the sgRNA, or if none, the sgRNA closest to the edit site
    if (scalar keys %sghash > 0){
        foreach my $k (keys %sghash){

            #check if the sgRNA is disrupted in edited sequence
            my $orientation = $sghash{$k}[1];
            my $sgRNA = $k;
            if ($orientation eq "sense"){
                if ($seq2 !~ /$sgRNA.G/){ #if the sgRNA is disrupted in the edited sequence:
                    $sghash{$k}[5] = 1;
                }
                else {
                    $sghash{$k}[5] = 0;
                }
            }
            elsif ($orientation eq "antisense"){
                my $rc_c_sg = reverse_complement($sgRNA);
                if ($seq2 !~ /C.$rc_c_sg/){ #if the sgRNA is disrupted in the edited sequence:
                    $sghash{$k}[5] = 1;
                }
                else {
                    $sghash{$k}[5] = 0;
                }
            }
        }
        
        #generate a top ranked sgRNA hash
        my $counter = 0;
        my %ranksghash;
        
        #sort by disruptive binary, then by distance to edit
        foreach my $k (sort {$sghash{$b}[5] <=> $sghash{$a}[5] || $sghash{$a}[4] <=> $sghash{$b}[4] } keys %sghash){
            if ($counter == 0){ # top rank sgRNA
                $chosenSG = $k;
                $chosenCutPos = $sghash{$k}[3];
                $chosenDistance = $sghash{$k}[4];
                $chosenOrientation = $sghash{$k}[1];
                $gcPctg = $sghash{$k}[2];
            }
            # sgRNA seq,  NA,cut position, orientation, distance, disrupt, GC % 
            $ranksghash{$counter} = [$k,"NA",$sghash{$k}[3],$sghash{$k}[1],$sghash{$k}[4],$sghash{$k}[5],$sghash{$k}[2]];
            $counter++;
            
        }
        return (\%sghash,$chosenSG,$chosenCutPos,$chosenOrientation,$chosenDistance,$gcPctg,\%ranksghash);
    }
    # if no candidate sgRNA found, return value stating none
    else {
        my $val = "none";
        return($val,$val);
    }
}

# if no sgRNA designer file uploaded, design nicking sgRNA; NG PAM
sub find_NG_nick_sgRNA{
    my ($seq1,$minEditPos,$maxEditPos,$maxEditDistance,$chosenCutPos,$chosenOrientation,$minNickDist,$maxNickDist) = @_;
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
            if ($chars[$i] eq "G"){
                for (my $k = $i-21; $k <= $i-2; $k++){
                    $tempsg .= $chars[$k];
                }
                #determine cut position and distance to primary sgRNA cut
                $seq1 =~ /$tempsg.G/;
                my $tempCutPos = $-[0]+18;
                my $tempDistance = $chosenCutPos-$tempCutPos-1;

                if ((abs($tempDistance)>=$minNickDist) && (abs($tempDistance) <= $maxNickDist) && $tempsg !~ /TTTTT/){
                    $nicksghash{$tempsg} = [$tempsg, $tempCutPos,"sense", $tempDistance];
                }
            }
        }
    }

    #else, if chosen sgRNA is sense, look for antisense sgRNAs for PE3:
    if ($chosenOrientation eq "sense"){
    #now check antisense orientation:
        for (my $i = 0; $i < scalar @chars-21; $i++){
            my $tempsg;
            if ($chars[$i] eq "C"){
                for (my $k = $i+2; $k <= $i+21; $k++){
                    $tempsg .= $chars[$k];
                }
                $tempsg = reverse_complement($tempsg);
            
                #determine cut position and distance to primary sgRNA cut
                my $rc_tempsg = reverse_complement($tempsg);
                $seq1 =~ /C.$rc_tempsg/;
                my $tempCutPos = $-[0]+5;

                my $tempDistance = $tempCutPos-$chosenCutPos+1;
                if ((abs($tempDistance)>=$minNickDist) && (abs($tempDistance) <= $maxNickDist) && $tempsg !~ /TTTTT/ ){
                    $nicksghash{$tempsg} = [$tempsg,$tempCutPos,"antisense",$tempDistance];
                }
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

#Find PE3b nicking sgRNA; NG PAM
sub find_NG_pe3b_sgRNA{
    my ($seq1,$minEditPos,$maxEditPos,$maxEditDistance,$chosenCutPos,$chosenOrientation,$seq2,$wtdelcounter,$mutdelcounter) = @_;
    my $chosenNickOrientation;
    my $chosenNickSG;
    my $chosenNickSGPos;
    my $chosenNickSGDist;
    my %nicksghash;

    #Find spacers that are present in seq2 but not in seq1 -- only accessible after the edit
    #Only care about seed region changes (5nt) so just can't match 10N-NGG in seq1

    #find all candidate spacers in the edited sequence
    my @chars = split ("",$seq2);
    my %sghash;

    #if chosen sgRNA is antisense, look for sense sgRNAs for PE3:
    if ($chosenOrientation eq "antisense"){
        for (my $i = 21; $i < scalar @chars; $i++){
            my $tempsg;
            my $seedsg;
            if ($chars[$i] eq "G"){
                for (my $k = $i-21; $k <= $i-2; $k++){
                    $tempsg .= $chars[$k];
                }
                for (my $k = $i-11; $k <= $i-2; $k++){
                    $seedsg .= $chars[$k];
                }
            
                #check if seed sgRNA is absent in seq1
                if ($seq1 !~ /$seedsg.G/ && $tempsg !~ /TTTTT/) {
                    #determine cut position and distance to primary sgRNA cut
                    $seq2 =~ /$tempsg.G/;
                    my $tempCutPos = $-[0]+18;
                    my $tempDistance = $chosenCutPos-$tempCutPos-$mutdelcounter+1;

                    $nicksghash{$tempsg} = [$tempsg, $tempCutPos,"sense", $tempDistance];
                }
            }
        }
    }

    #else, if chosen sgRNA is sense, look for antisense sgRNAs for PE3:
    if ($chosenOrientation eq "sense"){
    #now check antisense orientation:
        for (my $i = 0; $i < scalar @chars-21; $i++){
            my $tempsg;
            my $seedsg;
            if ($chars[$i] eq "C"){
                for (my $k = $i+2; $k <= $i+21; $k++){
                    $tempsg .= $chars[$k];
                }
                for (my $k = $i+2; $k <= $i+11; $k++){
                    $seedsg .= $chars[$k];
                }
                $tempsg = reverse_complement($tempsg);

                #check if seed sgRNA is absent in seq1

                if ($seq1 !~ /C.$seedsg/ && $tempsg !~ /TTTTT/){
                    
                    #determine cut position and distance to primary sgRNA cut
                    my $rc_tempsg = reverse_complement($tempsg);
                    $seq2 =~ /C.$rc_tempsg/;
                    my $tempCutPos = $-[0]+5;

                    my $tempDistance = $tempCutPos-$chosenCutPos+1-$wtdelcounter;
                    $nicksghash{$tempsg} = [$tempsg,$tempCutPos,"antisense",$tempDistance];
                }
            }
        }
    }

    if (scalar keys %nicksghash > 0){
        #choose the nicking sgRNA that is closest to the primary nick
        my $ctt = 0;
        foreach my $k (sort {abs($nicksghash{$a}[3]) <=> abs($nicksghash{$b}[3])} keys %nicksghash){
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

#If an sgRNA was preselected, calculate the needed metrics (Cut position, orientation, distance to edit, GC %)
#take the chosen sgRNA sequence, the wildtype DNA sequence, minEditPos, maxEditPos, and maxEditDistance 
sub process_chosen_NG_sgRNA {
    my ($c_sg,$seq1,$minEditPos,$maxEditPos,$maxEditDistance,$wtdelcounter,$seq2) = @_;
    my $chosenOrientation;
    my $chosenSG = $c_sg;
    my $chosenCutPos;
    my $chosenDistance;
    my $gcPctg;

    $maxEditDistance = 150; #since user prespecified sgRNA, relax distance filter

    #Check if repetitive sequence or palindrome targeted by chosen sgRNA
    my $counter = 0;
    $counter += () = $seq1 =~ /$c_sg.G/;
    $counter += () = (reverse_complement($seq1) =~ /$c_sg.G/);
    #print $counter,"\n";
    if ($counter > 1){
        return ("non-unique","non-unique");
    }
    
    else {
        #determine orientation
        if ($seq1 =~ /$c_sg.G/){
            $chosenOrientation = "sense";
        }
        elsif (reverse_complement($seq1) =~ /$c_sg.G/){
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
            $seq1 =~ /$c_sg.G/;
            $chosenCutPos = $-[0]+18;
            $chosenDistance = $minEditPos-$chosenCutPos;
        }
        elsif ($chosenOrientation eq "antisense"){
            my $rc_c_sg = reverse_complement($c_sg);
            $seq1 =~ /C.$rc_c_sg/;
            $chosenCutPos = $-[0]+5;
            $chosenDistance = ($chosenCutPos-$maxEditPos)+$wtdelcounter; # gaps in the wildtype alignment will push the max edit position closer to the cut site; need to add back in those bases
        }

        #test if chosen sgRNA will disrupt the PAM/seed
        my $chosenDisrupt;
        if ($chosenOrientation eq "sense"){
            if ($seq2 !~ /$c_sg.G/){ #if the sgRNA is disrupted in the edited sequence:
                $chosenDisrupt = 1;
            }
            else {
                $chosenDisrupt = 0;
            }
        }
        elsif ($chosenOrientation eq "antisense"){
            my $rc_c_sg = reverse_complement($c_sg);
            if ($seq2 !~ /C.$rc_c_sg/){ #if the sgRNA is disrupted in the edited sequence:
                $chosenDisrupt = 1;
            }
            else {
                $chosenDisrupt = 0;
            }
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
        return ($chosenCutPos,$chosenOrientation,$chosenDistance,$gcPctg,$sgfoundStatus,$chosenDisrupt);
    }
}

###################
#NNN PAM
sub find_SpRY_sgRNA {
    my ($seq1,$minEditPos,$maxEditPos,$maxEditDistance, $wtdelcounter,$seq2) = @_;
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
        $seq1 =~ /$tempsg/;
        my $tempCutPos = $-[0]+18;
        my $tempDistance = $minEditPos-$tempCutPos;

        if ($minEditPos >= $tempCutPos && $maxEditPos <= ($tempCutPos+$maxEditDistance) && $tempsg !~ /TTTTT/){
            $sghash{$tempsg} = [$tempsg, "sense", $tempgcPctg,$tempCutPos,$tempDistance];
        }
    }

    #now check antisense orientation:
    for (my $i = 0; $i < scalar @chars-21; $i++){
        my $tempsg;
        for (my $k = $i+2; $k <= $i+21; $k++){
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
        $seq1 =~ /$rc_tempsg/;
        my $tempCutPos = $-[0]+3;
        my $tempDistance = ($tempCutPos-$maxEditPos)+$wtdelcounter; # gaps in the wildtype alignment will push the max edit position closer to the cut site; need to add back in those bases
        if ($minEditPos >= ($tempCutPos-$maxEditDistance) && $maxEditPos <= ($tempCutPos+$wtdelcounter) && $tempsg !~ /TTTTT/){
            $sghash{$tempsg} = [$tempsg, "antisense", $tempgcPctg,$tempCutPos,$tempDistance];
        }
    }

    # %sghash holds all identified sgRNAs that are candidates for the desired edit.
    # choose a single sgRNA that disrupts the sgRNA, or if none, the sgRNA closest to the edit site
    if (scalar keys %sghash > 0){
        foreach my $k (keys %sghash){

            #check if the sgRNA is disrupted in edited sequence
            my $orientation = $sghash{$k}[1];
            my $sgRNA = $k;
            if ($orientation eq "sense"){
                if ($seq2 !~ /$sgRNA/){ #if the sgRNA is disrupted in the edited sequence:
                    $sghash{$k}[5] = 1;
                }
                else {
                    $sghash{$k}[5] = 0;
                }
            }
            elsif ($orientation eq "antisense"){
                my $rc_c_sg = reverse_complement($sgRNA);
                if ($seq2 !~ /$rc_c_sg/){ #if the sgRNA is disrupted in the edited sequence:
                    $sghash{$k}[5] = 1;
                }
                else {
                    $sghash{$k}[5] = 0;
                }
            }
        }
        
        #generate a top ranked sgRNA hash
        my $counter = 0;
        my %ranksghash;
        
        #sort by disruptive binary, then by distance to edit
        foreach my $k (sort {$sghash{$b}[5] <=> $sghash{$a}[5] || $sghash{$a}[4] <=> $sghash{$b}[4] } keys %sghash){
            if ($counter == 0){ # top rank sgRNA
                $chosenSG = $k;
                $chosenCutPos = $sghash{$k}[3];
                $chosenDistance = $sghash{$k}[4];
                $chosenOrientation = $sghash{$k}[1];
                $gcPctg = $sghash{$k}[2];
            }
            # sgRNA seq,  NA,cut position, orientation, distance, disrupt, GC % 
            $ranksghash{$counter} = [$k,"NA",$sghash{$k}[3],$sghash{$k}[1],$sghash{$k}[4],$sghash{$k}[5],$sghash{$k}[2]];
            $counter++;
        }
        return (\%sghash,$chosenSG,$chosenCutPos,$chosenOrientation,$chosenDistance,$gcPctg,\%ranksghash);
    }
    # if no candidate sgRNA found, return value stating none
    else {
        my $val = "none";
        return($val,$val);
    }
}

# if no sgRNA designer file uploaded, design nicking sgRNA; NG PAM
sub find_SpRY_nick_sgRNA{
    my ($seq1,$minEditPos,$maxEditPos,$maxEditDistance,$chosenCutPos,$chosenOrientation,$minNickDist,$maxNickDist) = @_;
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
            for (my $k = $i-21; $k <= $i-2; $k++){
                $tempsg .= $chars[$k];
            }
            #determine cut position and distance to primary sgRNA cut
            $seq1 =~ /$tempsg/;
            my $tempCutPos = $-[0]+18;
            my $tempDistance = $chosenCutPos-$tempCutPos-1;

            if ((abs($tempDistance)>=$minNickDist) && (abs($tempDistance) <= $maxNickDist) && $tempsg !~ /TTTTT/){
                $nicksghash{$tempsg} = [$tempsg, $tempCutPos,"sense", $tempDistance];
            }
        }
    }

    #else, if chosen sgRNA is sense, look for antisense sgRNAs for PE3:
    if ($chosenOrientation eq "sense"){
    #now check antisense orientation:
        for (my $i = 0; $i < scalar @chars-21; $i++){
            my $tempsg;
            for (my $k = $i+2; $k <= $i+21; $k++){
                $tempsg .= $chars[$k];
            }
            $tempsg = reverse_complement($tempsg);
            
            #determine cut position and distance to primary sgRNA cut
            my $rc_tempsg = reverse_complement($tempsg);
            $seq1 =~ /$rc_tempsg/;
            my $tempCutPos = $-[0]+3;

            my $tempDistance = $tempCutPos-$chosenCutPos+1;
            if ((abs($tempDistance)>=$minNickDist) && (abs($tempDistance) <= $maxNickDist) && $tempsg !~ /TTTTT/){
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

#Find PE3b nicking sgRNA; NG PAM
sub find_SpRY_pe3b_sgRNA{
    my ($seq1,$minEditPos,$maxEditPos,$maxEditDistance,$chosenCutPos,$chosenOrientation,$seq2,$wtdelcounter,$mutdelcounter) = @_;
    my $chosenNickOrientation;
    my $chosenNickSG;
    my $chosenNickSGPos;
    my $chosenNickSGDist;
    my %nicksghash;

    #Find spacers that are present in seq2 but not in seq1 -- only accessible after the edit
    #Only care about seed region changes (5nt) so just can't match 10N-NGG in seq1

    #find all candidate spacers in the edited sequence
    my @chars = split ("",$seq2);
    my %sghash;

    #if chosen sgRNA is antisense, look for sense sgRNAs for PE3:
    if ($chosenOrientation eq "antisense"){
        for (my $i = 21; $i < scalar @chars; $i++){
            my $tempsg;
            my $seedsg;
            for (my $k = $i-21; $k <= $i-2; $k++){
                $tempsg .= $chars[$k];
            }
            for (my $k = $i-11; $k <= $i-2; $k++){
                $seedsg .= $chars[$k];
            }
            
            #check if seed sgRNA is absent in seq1
            if ($seq1 !~ /$seedsg/ && $tempsg !~ /TTTTT/){
                #determine cut position and distance to primary sgRNA cut
                $seq2 =~ /$tempsg/;
                my $tempCutPos = $-[0]+18;
                my $tempDistance = $chosenCutPos-$tempCutPos-$mutdelcounter+1;

                $nicksghash{$tempsg} = [$tempsg, $tempCutPos,"sense", $tempDistance];
            }
        }
    }

    #else, if chosen sgRNA is sense, look for antisense sgRNAs for PE3:
    if ($chosenOrientation eq "sense"){
    #now check antisense orientation:
        for (my $i = 0; $i < scalar @chars-21; $i++){
            my $tempsg;
            my $seedsg;
            for (my $k = $i+2; $k <= $i+21; $k++){
                $tempsg .= $chars[$k];
            }
            for (my $k = $i+2; $k <= $i+11; $k++){
                $seedsg .= $chars[$k];
            }
            $tempsg = reverse_complement($tempsg);

            #check if seed sgRNA is absent in seq1
            if ($seq1 !~ /$seedsg/ && $tempsg !~ /TTTTT/){
            
                #determine cut position and distance to primary sgRNA cut
                my $rc_tempsg = reverse_complement($tempsg);
                $seq2 =~ /$rc_tempsg/;
                my $tempCutPos = $-[0]+3;

                my $tempDistance = $tempCutPos-$chosenCutPos+1-$wtdelcounter;
                $nicksghash{$tempsg} = [$tempsg,$tempCutPos,"antisense",$tempDistance];
            }
        }
    }

    if (scalar keys %nicksghash > 0){
        #choose the nicking sgRNA that is closest to the primary nick
        my $ctt = 0;
        foreach my $k (sort {abs($nicksghash{$a}[3]) <=> abs($nicksghash{$b}[3])} keys %nicksghash){
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

#If an sgRNA was preselected, calculate the needed metrics (Cut position, orientation, distance to edit, GC %)
#take the chosen sgRNA sequence, the wildtype DNA sequence, minEditPos, maxEditPos, and maxEditDistance 
sub process_chosen_SpRY_sgRNA {
    my ($c_sg,$seq1,$minEditPos,$maxEditPos,$maxEditDistance,$wtdelcounter,$seq2) = @_;
    my $chosenOrientation;
    my $chosenSG = $c_sg;
    my $chosenCutPos;
    my $chosenDistance;
    my $gcPctg;

    $maxEditDistance = 150; #since user prespecified sgRNA, will not restrain distance
    #Check if repetitive sequence or palindrome targeted by chosen sgRNA
    my $counter = 0;
    $counter += () = $seq1 =~ /$c_sg/;
    $counter += () = (reverse_complement($seq1) =~ /$c_sg/);
    #print $counter,"\n";
    if ($counter > 1){
        return ("non-unique","non-unique");
    }
    
    else {
        #determine orientation
        if ($seq1 =~ /$c_sg/){
            $chosenOrientation = "sense";
        }
        elsif (reverse_complement($seq1) =~ /$c_sg/){
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
            $seq1 =~ /$c_sg/;
            $chosenCutPos = $-[0]+18;
            $chosenDistance = $minEditPos-$chosenCutPos;
        }
        elsif ($chosenOrientation eq "antisense"){
            my $rc_c_sg = reverse_complement($c_sg);
            $seq1 =~ /$rc_c_sg/;
            $chosenCutPos = $-[0]+3;
            $chosenDistance = ($chosenCutPos-$maxEditPos)+$wtdelcounter; # gaps in the wildtype alignment will push the max edit position closer to the cut site; need to add back in those bases
        }

        #test if chosen sgRNA will disrupt the PAM/seed
        my $chosenDisrupt;
        if ($chosenOrientation eq "sense"){
            if ($seq2 !~ /$c_sg/){ #if the sgRNA is disrupted in the edited sequence:
                $chosenDisrupt = 1;
            }
            else {
                $chosenDisrupt = 0;
            }
        }
        elsif ($chosenOrientation eq "antisense"){
            my $rc_c_sg = reverse_complement($c_sg);
            if ($seq2 !~ /$rc_c_sg/){ #if the sgRNA is disrupted in the edited sequence:
                $chosenDisrupt = 1;
            }
            else {
                $chosenDisrupt = 0;
            }
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
        return ($chosenCutPos,$chosenOrientation,$chosenDistance,$gcPctg,$sgfoundStatus,$chosenDisrupt);
    }
}

1;
