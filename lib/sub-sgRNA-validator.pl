#!/usr/bin/perl
use strict;
use warnings;
#validate the sgRNA file {
sub validate_sgRNA {
    my ($data_text,$test) = @_;
    $data_text =~ s/"//g;
    my @rows = split (/[\r\n]+/,$data_text);
    my $head = $rows[0];
    my @headers = split ("\t",$head);
    my $cutPosCol;
    my $sgRNACol;
    my $orientationCol;
    my $onTargetScoreCol;
    my $offTargetT1Col;
    my $counter = 0;
    my $boolfile = 0;
    my $crisprscancounter = 0;
    for (my $i = 0; $i < scalar @headers; $i++){
        if ($headers[$i] eq "sgRNA Cut Position (1-based)"){
            $cutPosCol = $i;
            
            $counter++;
        }
        if ($headers[$i] eq "sgRNA Sequence"){
            $sgRNACol = $i;
            $counter++;
        }
        if ($headers[$i] eq "Orientation"){
            $orientationCol = $i;
            $counter++;
        }
        if ($headers[$i] eq "On-Target Efficacy Score"){
            $onTargetScoreCol = $i;
            $counter++;
        }
        if ($headers[$i] eq "# Off-Target Tier I Match Bin I Matches"){
            $offTargetT1Col = $i;
            $counter++;
        }
        if ($headers[$i] eq "score_crisprscan"){
            $crisprscancounter++;
        }
    }
    if ($counter == 5){
        $boolfile = 1; # 1 = Broad results
    }
    if ($crisprscancounter != 0){
        $boolfile = 2; # 2 = CRISPRscan results
    }

    return ($boolfile);
}
1;