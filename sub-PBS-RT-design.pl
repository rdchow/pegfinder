#!/usr/bin/perl


##############################################################################
### Using the chosen sgRNA, find an appropriate primer binding site (PBS)

# PBS recommendation: 13 nt, with 40-60% GC content
# Will report all PBS from 8-17 nt length for downstream optimization
# PBS is complementary to the DNA upstream of the cut site
# 5' end of PBS is close to sgRNA cut site, while 3' end of PBS is farther from cut site

# extract the PBS sequences of 8-17 nt, save into %pbshash
sub find_pbs {
    my ($chosenSG, $gcPctg) = @_;
    my %pbshash;
    my $pbstable ='<table style ="width:30%; float = left">';
    $pbstable .= "<tr><th>Length</th><th>PBS_Seq</th></tr>";

    print "Length\tPBS sequence\n";
    for (my $i = 8; $i <= 17; $i++){
        my $pbsSeq = substr($chosenSG,17-$i,$i);
        $pbsSeq = reverse_complement($pbsSeq);
        $pbshash{$i} = $pbsSeq;
        print "$i\t$pbsSeq\n";
        $pbstable .= "<tr><td>$i</td><td>$pbsSeq</td></tr>";
    }

    #Pick a PBS based on GC percentage of the sgRNA, for the sake of having a default option
    #if GC 30-70%, use PBS of 13nt
    my $chosenPBSlen;
    if ($gcPctg >= 30 && $gcPctg <= 70){
        $chosenPBSlen = 13;
    }
    #if GC <30%, use PBS of 14nt
    elsif ($gcPctg < 30){
        $chosenPBSlen = 14;
    }
    #if GC >70%, use PBS of 12nt
    elsif ($gcPctg > 70){
        $chosenPBSlen = 12;
    }
    my $chosenPBS = $pbshash{$chosenPBSlen};
    $pbstable .= "</table>";
    print "\nChose $chosenPBS of length $chosenPBSlen as the PBS.\n\n";
    return(\%pbshash, $chosenPBS, $chosenPBSlen,$pbstable);
}


##############################################################################
# Finally, we design the RT template
sub find_RT {
    my ($align2,$chosenOrientation,$minEditPos,$maxEditPos,$chosenCutPos,$wtdelcounter) = @_;
    my $minimalEditLen; #holds the distance from sgRNA cut to farthest edit, using the alignment coordinates
        
    #pull out the DNA sequence in between the +1 cut site position and $maxEditPos from Sequence 2 alignment
    my @seq2chars = split ("",$align2);
    if ($chosenOrientation eq "sense"){
        $minimalEditLen = $maxEditPos-$chosenCutPos+1; #distance from sgRNA cut to farthest edited base
    }
    elsif ($chosenOrientation eq "antisense"){
        $minimalEditLen = $chosenCutPos-$minEditPos+1+$wtdelcounter;
        #$minimalEditLen = $chosenCutPos-$minEditPos+1-$wtdelcounter; #distance from sgRNA cut to farthest edited base. 
    }
    print "minimal edit length: $minimalEditLen\n";

    #Get the number of deleted bases in $align2, since we need to remove those for calculating the RT template length ($minimalRTLen)
    my $delCounter = 0;
    for (my $i = 0; $i < scalar @seq2chars; $i++){
        if ($seq2chars[$i] eq "-"){
            $delCounter++;
        }
    }
    my $minimalRTLen = $minimalEditLen - $delCounter; #holds the minimal sequence length needed on the RT template to generate the edits
    #print "$minimalRTLen\n";

    #If the chosen sgRNA is on the sense strand:
    my %rthash; # hash that stores the candidate RT templates
    my %allrthash; #hash that stores all candidate RT templates, irrespective of 'C' ending
    my @rtlengths; # array that stores the lengths of the candidate RT templates
    my $rttable = '<table style ="width:50%; float = left">';
    $rttable .= "<tr><th>Length</th><th>TemplateSequence</th><th>Warnings</th>";
    my $rtcounter = 0; #counts how many templates have been found that do not have 'C' ending.

    if ($chosenOrientation eq "sense"){
        #If the minimal RT template is <= 10nt, print all possible 10-16nt RT templates
        if ($minimalEditLen < 10){
            print "The edit distance is < 10 nt ($minimalEditLen nt), with $delCounter deletion base(s). Extracting all 10-16nt RT templates:\n";
            print "Length\tSenseSequence\tTemplate_Seq\tWarnings\n";
            for (my $i = 10; $i <= 16; $i++){
                my $align2Copy = $align2;
                $align2Copy =~ s/-//g;
                my $templateSeq = substr($align2Copy,$chosenCutPos-1,$i); # holds the sequence of the cDNA that will synthesized by the RT
                my $rcTemplate = reverse_complement($templateSeq); # holds the sequence of the RT template, used for cDNA synthesis. Use this for the pegRNA sequences
                print "$i\t$templateSeq\t$rcTemplate\t";
                $rttable .= "<tr><td>$i</td><td>$rcTemplate</td>";
                $allrthash{$i} = $rcTemplate;
                if ($templateSeq =~ /G$/){
                    print "First base of RT template is C; expected to have lower efficiency\n";
                    $rttable .= "<td>First base is 'C'</td></tr>";
                }
                else {
                    $rtcounter++;
                    $rthash{$i} = $rcTemplate; # we will only save RT templates that do not start with C
                    push(@rtlengths,$i);
                    print "\n";
                    $rttable .= "<td></td></tr>";
                }

                #if no RT templates with non-C ending have been found, take all the templates
                if ($rtcounter == 0){
                    foreach my $k (keys %allrthash){
                        $rthash{$k} = $allrthash{$k};
                        push (@rtlengths,$i);
                    }
                }



            }  
        }
        else {
            print "The edit distance is >= 10 nt ($minimalEditLen nt), with $delCounter deletion base(s). Returning all ", ($minimalRTLen+1) , "-", ($minimalRTLen+7), "nt RT templates:\n";
            print "Length\tSenseSequence\tTemplateSequence\tWarnings\n";
            for (my $i = $minimalRTLen+1; $i < $minimalRTLen+7; $i++){
                my $align2Copy = $align2;
                $align2Copy =~ s/-//g;
                my $templateSeq = substr($align2Copy,$chosenCutPos-1,$i); # holds the sequence of the cDNA that will synthesized by the RT
                my $rcTemplate = reverse_complement($templateSeq); # holds the sequence of the RT template, used for cDNA synthesis. Use this for the pegRNA sequences
                print "$i\t$templateSeq\t$rcTemplate\t";
                $rttable .= "<tr><td>$i</td><td>$rcTemplate</td>";
                $allrthash{$i} = $rcTemplate;
                if ($templateSeq =~ /G$/){
                    print "First base of RT template is C; expected to have lower efficiency\n";
                    $rttable .= "<td>First base is 'C'</td></tr>";
                }
                else {
                    $rtcounter++;
                    $rthash{$i} = $rcTemplate; # we will only save RT templates that do not start with C
                    push(@rtlengths,$i);
                    print "\n";
                    $rttable .= "<td></td></tr>";
                }
                #if no RT templates with non-C ending have been found, take all the templates
                if ($rtcounter == 0){
                    foreach my $k (keys %allrthash){
                        $rthash{$k} = $allrthash{$k};
                        push (@rtlengths,$i);
                    }
                }
            } 
        }
        
    }

    #If the chosen sgRNA is on the antisense strand:
    elsif ($chosenOrientation eq "antisense"){
        #If the minimal RT template is <= 10nt, print all possible 10-16nt RT templates
        if ($minimalEditLen < 10){
            print "The edit distance is < 10 nt ($minimalEditLen nt), with $delCounter deletion base(s).\nExtracting all 10-16nt RT templates:\n";
            print "Length\tSenseSequence\tTemplateSequence\tWarnings\n";
            for (my $i = 10; $i <= 16; $i++){
                my $align2Copy = $align2;
                $align2Copy =~ s/-//g;
                my $templateSeq = substr($align2Copy,($chosenCutPos-$delCounter+$wtdelcounter)-$i,$i); # holds the sequence of the cDNA that will synthesized by the RT
                my $rcTemplate = $templateSeq; # holds the sequence of the RT template, used for cDNA synthesis. Use this for the pegRNA sequences
                print "$i\t$templateSeq\t$rcTemplate\t";
                $rttable .= "<tr><td>$i</td><td>$rcTemplate</td>";
                $allrthash{$i} = $rcTemplate;
                if ($templateSeq =~ /^C/){
                    print "First base of RT template is C; expected to have lower efficiency\n";
                    $rttable .= "<td>First base is 'C'</td></tr>";
                }
                else {
                    $rtcounter++;
                    $rthash{$i} = $rcTemplate; # we will only save RT templates that do not start with C
                    push(@rtlengths,$i);
                    print "\n";
                    $rttable .= "<td></td></tr>";
                }
                #if no RT templates with non-C ending have been found, take all the templates
                if ($rtcounter == 0){
                    foreach my $k (keys %allrthash){
                        $rthash{$k} = $allrthash{$k};
                        push (@rtlengths,$i);
                    }
                }
            }  
        }
        else {
            print "The edit distance is >= 10 nt ($minimalEditLen nt), with $delCounter deletion base(s).\nExtracting all ", ($minimalRTLen+1) , "-", ($minimalRTLen+7), "nt RT templates:\n";
            print "Length\tSenseSequence\tTemplateSequence\tWarnings\n";
            for (my $i = $minimalRTLen+1; $i < $minimalRTLen+7; $i++){
                my $align2Copy = $align2;
                $align2Copy =~ s/-//g;
                my $templateSeq = substr($align2Copy,($chosenCutPos-$delCounter+$wtdelcounter)-$i,$i); # holds the sequence of the cDNA that will synthesized by the RT
                my $rcTemplate = $templateSeq; # holds the sequence of the RT template, used for cDNA synthesis. Use this for the pegRNA sequences
                print "$i\t$templateSeq\t$rcTemplate\t";
                $rttable .= "<tr><td>$i</td><td>$rcTemplate</td>";
                $allrthash{$i} = $rcTemplate;
                if ($templateSeq =~ /^C/){
                    print "First base of RT template is C; expected to have lower efficiency\n";
                    $rttable .= "<td>First base is 'C'</td></tr>";
                }
                else {
                    $rtcounter++;
                    $rthash{$i} = $rcTemplate; # we will only save RT templates that do not start with C
                    push(@rtlengths,$i);
                    print "\n";
                    $rttable .= "<td></td></tr>";
                }
                #if no RT templates with non-C ending have been found, take all the templates
                if ($rtcounter == 0){
                    foreach my $k (keys %allrthash){
                        $rthash{$k} = $allrthash{$k};
                        push (@rtlengths,$i);
                    }
                }
            } 
        }
    }
    $rttable .= "</table>";
    #For sake of having a default, use some criteria to pick one RT template
    print "\nChoosing one RT template sequence...\n";
    my $chosenRT;

    #Choose the median size of all candidate RT templates in %rthash
    my $midRTindex;
    if (scalar @rtlengths %2 == 0){ #if even number of candidates, choose the shorter of the two
        $midRTindex = $rtlengths[(scalar(@rtlengths)/2)-1];
    }
    else { #if odd number of candidates, choose the median size
        $midRTindex = $rtlengths[(scalar(@rtlengths)/2)-0.5];
    }
    $chosenRT = $rthash{$midRTindex};
    print "Of the candidate RT templates, chose $chosenRT of length $midRTindex (median size of candidates).\n";
    return(\%rthash,$chosenRT,$midRTindex,$rttable);
}
1;