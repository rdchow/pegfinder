#!/usr/bin/env perl
use Mojolicious::Lite;
use File::Temp qw/ tempfile tempdir /;
#plugin AutoReload => {};
plugin 'RenderFile';
require "./lib/sub-needleman-wunsch.affine.pl";
require "./lib/sub-sgRNA-finder.general.pl";
require "./lib/sub-sgRNA-validator.pl";
require "./lib/sub-PBS-RT-design.pl";
require "./lib/prepare-output-table.pl";

###kill hypnotoads with --quit,
plugin Config => {
  default => {
    hypnotoad => {
      listen => ["http://*:$ENV{PORT}"]
    }
  }
};

# Upload form in DATA section
get '/' => 'form';

# Multipart upload handler

post '/upload' => sub {
  my $c = shift;
  
  # Process form input # 
  #get the wildtype sequence, make sure it is DNA sequence
  return $c->redirect_to('form') unless my $wt = $c->param('wildtype');
  $wt = uc($wt);
  $wt =~ s/[\r\n]+//g;

  return $c->render(text => 'Wildtype sequence is invalid. Please enter only DNA sequence in plaintext format.', status => 200)
    unless $wt =~ /\A[ACGT]+\z/i;

  #get the edited sequence, make sure it is DNA sequence
  return $c->redirect_to('form') unless my $edit = $c->param('edited');
  $edit = uc($edit);
  $edit =~ s/[\r\n]+//g;
  return $c->render(text => 'Edited/desired sequence is invalid. Please enter only DNA sequence in plaintext format.', status => 200)
    unless $edit =~ /\A[ACGT]+\z/i;

  #check if wt and edit are the same
  if ($wt eq $edit){ return $c->render(text=>'Wildtype and edited sequence are identical. Please reenter the form.', status=>200);}

  #check length of inputs to limit computational costs
  if (length($wt) > 500 || length($edit) > 500){ return $c->render(text=>'Input wildtype or edited sequence is >500 bp. Please reduce the length of the flanking arms around the edit site (recommend 100 bp).', status=>200);}
  
  #Needleman-Wunsch alignment of the two sequences
  my $seq1 = uc($wt);
  my $seq2 = uc($edit);
  my @nw_out = needleman_wunsch($seq1,$seq2);
  my ($align1,$align2,$minEditPos,$maxEditPos,$trimmingStatus5p,$trimmingStatus3p,$wtdelcounter,$alterhashpt,$mutdelcounter) = @nw_out;
  my %alterhash = %$alterhashpt;
  my $align1f = $align1; #formatted $align1 to have 70 chars per line
  my $align2f = $align2; #formmated $align2 to have 70 chars per line

  $align1f =~ s/(.{70})\K(?=.)/<br>/g;
  $align2f =~ s/(.{70})\K(?=.)/<br>/g;
  
  #if alignments are gapped at the 5' or 3' ends, send error message
  if ($trimmingStatus5p != 0 || $trimmingStatus3p != 0){
     return $c->render(text => 'Wildtype and desired sequence are not aligned on the 5\' or 3\' ends.<br>Please trim off the hanging DNA bases and rerun pegFinder (see alignment below). <br><br>If you used sgRNA finder results (Broad or CRISPRScan), please also rerun the sgRNA finders using the revised wildtype sequence (if modified). <br><hr><br>Wildtype:<br>'.$align1f.'<br><br>Edited:<br>'.$align2f.'<br> ', status => 200);
  }

  #if edit is too close to either end, send error message
  if ($minEditPos <= 20 || (length($align1)-$maxEditPos) <= 20){
     return $c->render(text => 'Desired edits are too close to 5\' or 3\' ends (see alignment below).<br>Please include longer arms flanking the desired edits (recommend >100 bp flanks).<br><br>If you used sgRNA finder results (Broad or CRISPRScan), please also rerun the sgRNA finders using the revised wildtype sequence (if modified). <br><hr><br>Wildtype:<br>'.$align1f.'<br><br>Edited:<br>'.$align2f.'<br> ', status => 200);
  }

  #Process the sgRNA input
  #Can either take a file upload, a chosen sgRNA, or nothing
  my $sgfile = $c->param('sgRNA'); # sgRNA file upload
  my $c_sg = $c->param('c_sgRNA'); # preselected sgRNA
  my $pe3Bool; # whether to search for PE3 sgRNAs
  my $rgn = $c->param('enzyme'); # selected enzyme
  my $minNickDist = $c->param('minNickDist'); # minimum nick distance
  my $maxNickDist = $c->param('maxNickDist'); # minimum nick distance
  my $maxSgCt = $c->param('maxSgCt'); # maximum # of sgRNAs to include in design table

  return $c->render(text => 'Minimum nick distance is invalid, please enter a number in plaintext format.', status => 200)
    unless $minNickDist =~ /\A[1234567890]+\z/i;

  return $c->render(text => 'Maximum nick distance is invalid, please enter a number in plaintext format.', status => 200)
    unless $maxNickDist =~ /\A[1234567890]+\z/i;

  return $c->render(text => 'Maximum # of sgRNAs to include in design table is invalid, please enter a number in plaintext format.', status => 200)
    unless $maxSgCt =~ /\A[1234567890]+\z/i;

  if (defined $c->param('PE3cb')){
    $pe3Bool = $c->param('PE3cb');
  }
  else {
    $pe3Bool = 1;
  }

  my ($chosenSG,$chosenCutPos,$chosenOrientation,$chosenDistance,$gcPctg, $chosenDisrupt,%sghash, $sghashpt, $ranksghashpt, %ranksghash,$data_text,$sgfoundStatus, $maxEditDistance, $chosenNickSGPos,$chosenNickSG,$chosenNickSGDist,$nicksghashpt, %nicksghash,$chosenNickOrientation); # declare all the variables
  my ($chosenNickSGPos3b,$chosenNickSG3b,$chosenNickSGDist3b,$nicksghashpt3b, %nicksghash3b,$chosenNickOrientation3b); # declare all the variables

  $maxEditDistance = 150; #default
  
  if ($rgn eq "Cas9-NGG"){
    $maxEditDistance = 150; # only consider sgRNAs if they cut within 150 nt of the desired edit
  }
  elsif ($rgn eq "Cas9-NG"){
    $maxEditDistance = 50; # with greater targeting range, limit max distance
  }
  elsif ($rgn eq "Cas9-SpRY"){
    $maxEditDistance = 20; # with greater targeting range, limit max distance
  }

  my $sgtable = '<table style ="width:60%">';
  my $nicksgtable = '<table style ="width:75%">';

  return $c->render(text => 'Input sequence is too big (>1 kb). Please reduce sequence length.', status => 200)
    if $c->req->is_limit_exceeded; # Check input size

  ###### Case 1: User does not submit sgRNA file, nor a chosen sgRNA ##########
  if ( $c_sg eq "" && $sgfile->asset->slurp eq ""){
    my @sgdata;
    @sgdata = find_choose_sgRNA_general($seq1,$minEditPos,$maxEditPos,$maxEditDistance,$wtdelcounter,$seq2,$rgn);

    #if no sgRNAs were found, report back
    return $c->render(text => 'No candidate sgRNAs found.', status => 200)
      unless ($sgdata[0] ne "none");
      
    #Otherwise, if an sgRNA was found:
    ($sghashpt,$chosenSG,$chosenCutPos,$chosenOrientation,$chosenDistance,$gcPctg,$ranksghashpt) = @sgdata;
    %sghash = %$sghashpt;
    %ranksghash= %$ranksghashpt;
    $sgfoundStatus = 2;
    $sgtable .= "<tr><th>sgRNA_Seq</th><th>CutPosition</th><th>Orientation</th><th>DistanceToEditStart</th><th>Seed/PAM_Disrupt</th><th>Rank</th><th>Chosen</th></tr>";
     # sgRNA seq,  NA,cut position, orientation, distance, disrupt, GC % 
    foreach my $k (sort {$a<=>$b} keys %ranksghash){
      my $numrank = $k+1;
      $sgtable .= "<tr><td>$ranksghash{$k}[0]</td><td>$ranksghash{$k}[2]</td><td>$ranksghash{$k}[3]</td><td>$ranksghash{$k}[4]</td><td>$ranksghash{$k}[5]</td><td>$numrank</td>";
      if ($ranksghash{$k}[0] eq $chosenSG){
        $sgtable .= "<td>X</td></tr>";
      }
      else { 
        $sgtable .= "<td></td></tr>";
      }
    }
    $sgtable .= "</table>";
  
    if ($pe3Bool > 0){ # check mark for finding PE3 secondary guides
      my @nickdata;
      my @pe3bnickdata;
      @nickdata = find_choose_nick_sgRNA_general($seq1, $minEditPos,$maxEditPos,$maxEditDistance,$chosenCutPos,$chosenOrientation,$minNickDist,$maxNickDist,$rgn);
      ($chosenNickSG,$chosenNickSGPos,$chosenNickOrientation,$chosenNickSGDist,$nicksghashpt) = @nickdata;
      @pe3bnickdata = find_pe3b_sgRNA_general($seq1,$minEditPos,$maxEditPos,$maxEditDistance,$chosenCutPos,$chosenOrientation,$seq2,$wtdelcounter,$mutdelcounter,$rgn);
      ($chosenNickSG3b,$chosenNickSGPos3b,$chosenNickOrientation3b,$chosenNickSGDist3b,$nicksghashpt3b) = @pe3bnickdata;
     
      if ($chosenNickSG3b ne "none found"){ # if a PE3b sgRNA was found
        %nicksghash3b = %$nicksghashpt3b;

        $nicksgtable .= "<tr><th>Nicking sgRNA_Seq</th><th>CutPosition</th><th>Orientation</th><th>DistanceTo_pegRNA_nick</th><th>Type</th><th>Chosen</th></tr>";

        foreach my $k (sort {abs($nicksghash3b{$a}[3]) <=> abs($nicksghash3b{$b}[3])} keys %nicksghash3b){
          $nicksgtable .= "<tr><td>$k</tdh><td>$nicksghash3b{$k}[1]</td><td>$nicksghash3b{$k}[2]</td><td>$nicksghash3b{$k}[3]</td><td>PE3b</td>";
          if ($k eq $chosenNickSG3b){
            $nicksgtable .= "<td>X (PE3b)</td></tr>";
          }
          else { 
              $nicksgtable .= "<td></td></tr>";
          }
        }
      }

      if ($chosenNickSG ne "none found"){ # if a PE3 sgRNA was found
        %nicksghash = %$nicksghashpt;
        
        if ($chosenNickSG3b eq "none found"){
          $nicksgtable .= "<tr><th>Nicking sgRNA_Seq</th><th>CutPosition</th><th>Orientation</th><th>DistanceTo_pegRNA_nick</th><th>Type</th><th>Chosen</th></tr>";
        }
        foreach my $k (sort {abs($nicksghash{$a}[3]-50) <=> abs($nicksghash{$b}[3]-50)} keys %nicksghash){
          $nicksgtable .= "<tr><td>$k</tdh><td>$nicksghash{$k}[1]</td><td>$nicksghash{$k}[2]</td><td>$nicksghash{$k}[3]</td><td>PE3</td>";
          if ($k eq $chosenNickSG){
            $nicksgtable .= "<td>X (PE3)</td></tr>";
          }
          else { 
              $nicksgtable .= "<td></td></tr>";
          }
        }
      }
      

      elsif ($chosenNickSG eq "none found" && $chosenNickSG3b eq "none found") {
        $nicksgtable .= '<tr><br>No suitable secondary nicking sgRNA found</tr>';
      }
      $nicksgtable .= "</table>";
    }

  }

  ###### Case 2: sgRNA was pre-chosen: #######
  elsif (defined $c_sg && $c_sg ne ""){
    #if an sgRNA file was uploaded:
    my $boolfile = 0; #tells us if the uploaded sgRNA file is the expected output. 0: invalid/absent, 1: Broad, 2: CRISPRscan
    $data_text = $sgfile->asset->slurp;
    $c_sg =~ s/[\r\n]//g;
    $c_sg = uc ($c_sg);
    $chosenSG = $c_sg;
    $maxSgCt = 1;
    #if sgRNA is not 20nt, or contains non DNA characters
    
    if (length($c_sg) != 20 || $c_sg !~ /\A[ACGT]+\z/i){ 
      return $c->render(text => 'Preselected sgRNA is invalid, please only enter the 20nt spacer sequence.', status => 200);
    }

    #if sgRNA is not in the wildtype sequence:
    elsif ($rgn eq "Cas9-NGG"){
      if ($seq1 !~ /$c_sg.GG/ && reverse_complement($seq1) !~ /$c_sg.GG/ ){
        return $c->render(text => 'Preselected sgRNA is not present in wildtype sequence (using NGG PAM). Please do one of the following: 1) Choose a different sgRNA, 2) Rerun pegFinder without specifying an sgRNA, or 3) or use an sgRNA finder (Broad, CRISPRscan) and upload the results.', status => 200);
      }
    }
    elsif ($rgn eq "Cas9-NG"){
      if ($seq1 !~ /$c_sg.G/ && reverse_complement($seq1) !~ /$c_sg.G/ ){
        return $c->render(text => 'Preselected sgRNA is not present in wildtype sequence (using NG PAM). Please do one of the following: 1) Choose a different sgRNA, or 2) Rerun pegFinder without specifying an sgRNA.', status => 200);
      }
    }
    elsif ($rgn eq "Cas9-SpRY"){
      if ($seq1 !~ /$c_sg/ && reverse_complement($seq1) !~ /$c_sg/ ){
        return $c->render(text => 'Preselected sgRNA is not present in wildtype sequence (using NNN PAM). Please do one of the following: 1) Choose a different sgRNA, or 2) Rerun pegFinder without specifying an sgRNA.', status => 200);
      }
    }

    #if sgRNA is present in wildtype sequence:
    my @sgdata;
    @sgdata = process_chosen_sgRNA_general($c_sg,$seq1,$minEditPos,$maxEditPos,$maxEditDistance,$wtdelcounter,$seq2,$rgn);
       
    #if the chosen sgRNA matches to multiple positions, report back
    return $c->render(text => 'Preselected sgRNA matches to multiple positions in wildtype sequence. Please do one of the following: 1) Choose a different sgRNA, 2) Rerun pegFinder without specifying an sgRNA, or 3) or use an sgRNA finder (Broad, CRISPRscan) and upload the results.', status => 200)
      unless ($sgdata[0] ne "non-unique");

    
    ($chosenCutPos,$chosenOrientation,$chosenDistance,$gcPctg,$sgfoundStatus,$chosenDisrupt) = @sgdata;
    if ($sgfoundStatus == 0){ # if the preselected sgRNA is incompatible
      return $c->render(text=>'Preselected sgRNA is incompatible with desired edit: it cuts 3\' (downstream) of the alterations. Please do one of the following: 1) Choose a different sgRNA, 2) Rerun pegFinder without specifying an sgRNA, or 3) or use an sgRNA finder (Broad, CRISPRscan) and upload the results.');
    }
    elsif ($sgfoundStatus == 1 && $rgn eq "Cas9-NGG"){
        return $c->render(text=>'Preselected sgRNA cuts >'.$maxEditDistance.'nt from the desired edits, and is predicted to be lower efficiency. <br>Please do one of the following: 1) Choose a different sgRNA, 2) Rerun pegFinder without specifying an sgRNA, or 3) or use an sgRNA finder (Broad, CRISPRscan) and upload the results.');
    }
    $sgtable .= "<tr><th>sgRNA_Seq</th><th>CutPosition</th><th>Orientation</th><th>DistanceToEditStart</th><th>Seed/PAM_Disrupt</th><th>Chosen</th></tr>";
    $sgtable .= "<tr><td>$c_sg</td><td>$chosenCutPos</td><td>$chosenOrientation</td><td>$chosenDistance</td><td>$chosenDisrupt</td><td>X</td></tr>";
    $sgtable .= "</table>";

    #evaluate whether we need to find PE3 secondary sgRNAs
    if ($pe3Bool > 0){ # check mark for finding PE3 secondary guides
      if (defined $sgfile && $data_text ne "" && $rgn eq "Cas9-NGG"){ # can only use sgRNA file if NGG PAM
        $boolfile = validate_sgRNA($data_text,"100");#Validate the sgRNA file
      }
      my @nickdata;

      if ($boolfile == 0 || $data_text eq "" || $rgn ne "Cas9-NGG"){ #if the file is invalid or absent, or non-Cas9-NGG was chosen, use the standalone subroutine
        @nickdata = find_choose_nick_sgRNA_general($seq1, $minEditPos,$maxEditPos,$maxEditDistance,$chosenCutPos,$chosenOrientation,$minNickDist,$maxNickDist,$rgn);
      }
      elsif ($boolfile == 1 && $rgn eq "Cas9-NGG") {#if the sgRNA file is a valid Broad file
        @nickdata = find_broad_nicksgRNA($data_text, $minEditPos,$maxEditPos,$maxEditDistance,$chosenCutPos,$chosenOrientation,$minNickDist,$maxNickDist);
      }
      elsif ($boolfile == 2 && $rgn eq "Cas9-NGG") { # if the sgRNA file is a valid CRISPRscan input
        @nickdata = find_crisprscan_nicksgRNA($data_text, $minEditPos,$maxEditPos,$maxEditDistance,$chosenCutPos,$chosenOrientation,$seq1,$minNickDist,$maxNickDist);
      }
      ($chosenNickSG,$chosenNickSGPos,$chosenNickOrientation,$chosenNickSGDist,$nicksghashpt) = @nickdata;
        

      my @pe3bnickdata = find_pe3b_sgRNA_general($seq1,$minEditPos,$maxEditPos,$maxEditDistance,$chosenCutPos,$chosenOrientation,$seq2,$wtdelcounter,$mutdelcounter,$rgn);
      ($chosenNickSG3b,$chosenNickSGPos3b,$chosenNickOrientation3b,$chosenNickSGDist3b,$nicksghashpt3b) = @pe3bnickdata;

      # if a valid sgRNa file was uploaded:
      if ($boolfile == 1 || $boolfile == 2){
        if ($chosenNickSG3b ne "none found"){ # if a PE3b sgRNA was found
          %nicksghash3b = %$nicksghashpt3b;
          $nicksgtable .= "<tr><th>Nicking sgRNA_Seq</th><th>OnTargetScore</th><th>CutPosition</th><th>Orientation</th><th>DistanceTo_pegRNA_nick</th><th>Type</th><th>Chosen</th></tr>";

          foreach my $k (sort {$nicksghash3b{$a}[3] <=> $nicksghash3b{$b}[3] } keys %nicksghash3b){
            $nicksgtable .= "<tr><td>$k</td><td>NA</td><td>$nicksghash3b{$k}[1]</td><td>$nicksghash3b{$k}[2]</td><td>$nicksghash3b{$k}[3]</td><td>PE3b</td>";
            if ($k eq $chosenNickSG3b){
              $nicksgtable .= "<td>X (PE3b)</td></tr>";
            }
            else { 
              $nicksgtable .= "<td></td></tr>";
            }
          }
        }
                  
        if ($chosenNickSG ne "none found"){ # if a PE3 sgRNA was found
          %nicksghash = %$nicksghashpt;
                    
          if ($chosenNickSG3b eq "none found"){
            $nicksgtable .= "<tr><th>Nicking sgRNA_Seq</th><th>OnTargetScore</th><th>CutPosition</th><th>Orientation</th><th>DistanceTo_pegRNA_nick</th><th>Type</th><th>Chosen</th></tr>";
          }
          foreach my $k (sort {$nicksghash{$b}[0] <=> $nicksghash{$a}[0]} keys %nicksghash){
            $nicksgtable .= "<tr><td>$k</tdh><td>$nicksghash{$k}[0]</td><td>$nicksghash{$k}[1]</td><td>$nicksghash{$k}[2]</td><td>$nicksghash{$k}[3]</td><td>PE3</td>";
            if ($k eq $chosenNickSG){
              $nicksgtable .= "<td>X (PE3)</td></tr>";
            }
            else { 
              $nicksgtable .= "<td></td></tr>";
            }
          }
        }

        #if no PE3/PE3b sgRNAs found
        elsif ($chosenNickSG eq "none found" && $chosenNickSG3b eq "none found") {
          $nicksgtable .= '<tr><br>No suitable secondary nicking sgRNA found</tr>';
        }
        $nicksgtable .= "</table>";
      }

      ##### if no sgRNA file was provided:
      else {
        if ($chosenNickSG3b ne "none found"){ # if a PE3b sgRNA was found
          %nicksghash3b = %$nicksghashpt3b;

          $nicksgtable .= "<tr><th>Nicking sgRNA_Seq</th><th>CutPosition</th><th>Orientation</th><th>DistanceTo_pegRNA_nick</th><th>Type</th><th>Chosen</th></tr>";

          foreach my $k (sort {$nicksghash3b{$a}[3] <=> $nicksghash3b{$b}[3]} keys %nicksghash3b){
            $nicksgtable .= "<tr><td>$k</tdh><td>$nicksghash3b{$k}[1]</td><td>$nicksghash3b{$k}[2]</td><td>$nicksghash3b{$k}[3]</td><td>PE3b</td>";
            if ($k eq $chosenNickSG3b){
              $nicksgtable .= "<td>X (PE3b)</td></tr>";
            }
            else { 
              $nicksgtable .= "<td></td></tr>";
            }
          }
        }

        if ($chosenNickSG ne "none found"){ # if a PE3 sgRNA was found
          %nicksghash = %$nicksghashpt;
                
          if ($chosenNickSG3b eq "none found"){
            $nicksgtable .= "<tr><th>Nicking sgRNA_Seq</th><th>CutPosition</th><th>Orientation</th><th>DistanceTo_pegRNA_nick</th><th>Type</th><th>Chosen</th></tr>";
          }
          foreach my $k (sort {abs($nicksghash{$a}[3]-50) <=> abs($nicksghash{$b}[3]-50)} keys %nicksghash){
            $nicksgtable .= "<tr><td>$k</tdh><td>$nicksghash{$k}[1]</td><td>$nicksghash{$k}[2]</td><td>$nicksghash{$k}[3]</td><td>PE3</td>";
            if ($k eq $chosenNickSG){
              $nicksgtable .= "<td>X (PE3)</td></tr>";
            }
            else { 
              $nicksgtable .= "<td></td></tr>";
            }
          }
        }

        # no PE3 or PE3b sgRNAs found
        elsif ($chosenNickSG eq "none found" && $chosenNickSG3b eq "none found") {
          $nicksgtable .= '<tr><br>No suitable secondary nicking sgRNA found</tr>';
        }
        $nicksgtable .= "</table>";
      }
    }
  }

  ###### Case 3: sgRNA was not pre-chosen, and a sgRNA designer file was uploaded (only valid for Cas9-NGG): #######
  elsif ($sgfile ne "" && defined $sgfile){
    my $boolfile = 0; #tells us if the uploaded sgRNA file is valid
    $data_text = $sgfile->asset->slurp;
    my @sgdata;

    #if an sgRNA file was uploaded:
    if (defined $sgfile && $data_text ne ""){
      #Validate the sgRNA file
      $boolfile = validate_sgRNA($data_text,"100");
    }

    if ($boolfile eq "0"){
        $c->render(text => 'sgRNA file is invalid. Please do one of the following: 1) Upload the results file from the Broad sgRNA designer or CRISPRscan on the wildtype sequence, 2) Enter a preselected sgRNA, or 3) Rerun pegFinder and do not upload a results file.', status => 200);
    }
    elsif ($rgn eq "Cas9-NGG"){
      if ($boolfile == 1){ #if valid sgRNA file: Broad 
        @sgdata = process_broad_sgRNA($data_text, $minEditPos,$maxEditPos,$maxEditDistance,$wtdelcounter,$seq1,$seq2);
      }
      elsif ($boolfile == 2){ #if valid sgRNA file: CRISPRscan
        @sgdata = process_crisprscan_sgRNA($data_text, $minEditPos,$maxEditPos,$maxEditDistance,$wtdelcounter,$seq1,$seq2);
      }

      #also find sgRNAs independently to confirm the sgRNA file is derived from the exact wildtype sequence
      my @confirmsgdata = find_choose_sgRNA_general($seq1,$minEditPos,$maxEditPos,$maxEditDistance,$wtdelcounter,$seq2,$rgn);

      #if no sgRNAs were found, report back
      if ($sgdata[0] eq "none"){
        return $c->render(text => 'No candidate sgRNAs from the supplied file could be matched to the input sequence. Please do one of the following: 1) Rerun the Broad sgRNA designer or CRISPRscan on the wildtype sequence, 2) Enter a preselected sgRNA, or 3) Rerun pegFinder and do not upload a results file.', status => 200)
      }
      else { 
        #Otherwise, if an sgRNA was found:
        ($sghashpt,$chosenSG,$chosenCutPos,$chosenOrientation,$chosenDistance,$gcPctg,$ranksghashpt) = @sgdata;
        %sghash = %$sghashpt;
        %ranksghash = %$ranksghashpt;

        #check if the sgRNA file is aligned with the wildtype sequence
        my ($sghashptC,$chosenSGC,$chosenCutPosC,$chosenOrientationC,$chosenDistanceC,$gcPctgC,$ranksghashptC) = @confirmsgdata;
        my %sghashC = %$sghashptC;
        my $counter = 0;
        foreach my $k (keys %sghashC){ #[$tempsg, "sense", $tempgcPctg,$tempCutPos,$tempDistance]; #[$onTargetScore,$cutPos,$orientation,$distance];
          if (defined $sghash{$k}){
            #$counter++;
          
            if (abs($sghash{$k}[1] - $sghashC{$k}[3]) >0 || $sghash{$k}[2] ne $sghashC{$k}[1]){
              $counter++;
            }
          }
        }
        if ($counter != 0){
          return $c->render(text => 'Uploaded sgRNA file is not aligned with the input wildtype sequence. Please do one of the following: 1) Rerun the Broad sgRNA designer or CRISPRscan on the wildtype sequence, 2) Enter a preselected sgRNA, or 3) Rerun pegFinder and do not upload a results file.', status => 200)
        }


        $sgfoundStatus = 2;
        $sgtable .= "<tr><th>sgRNA_Seq</th><th>OnTargetScore</th><th>CutPosition</th><th>Orientation</th><th>DistanceToEditStart</th><th>Seed/PAM_Disrupt</th><th>Rank</th><th>Chosen</th></tr>";
        # sgRNA seq, OT score, cut position, orientation, distance, disrupt, GC % 
        foreach my $k (sort {$a<=>$b} keys %ranksghash){
          my $numrank = $k+1;
          $sgtable .= "<tr><td>$ranksghash{$k}[0]</td><td>$ranksghash{$k}[1]</td><td>$ranksghash{$k}[2]</td><td>$ranksghash{$k}[3]</td><td>$ranksghash{$k}[4]</td><td>$ranksghash{$k}[5]</td><td>$numrank</td>";
          if ($ranksghash{$k}[0] eq $chosenSG){
            $sgtable .= "<td>X</td></tr>";
          }
          else { 
            $sgtable .= "<td></td></tr>";
          }
        }
        $sgtable .= "</table>";

        #find a nicking sgRNA if the box was ticked on
        if ($pe3Bool > 0){
          my @nickdata;

          if ($boolfile == 1) {#if the sgRNA file is a valid Broad file
            @nickdata = find_broad_nicksgRNA($data_text, $minEditPos,$maxEditPos,$maxEditDistance,$chosenCutPos,$chosenOrientation,$minNickDist,$maxNickDist);
          }
          elsif ($boolfile == 2) { # if the sgRNA file is a valid CRISPRscan input
            @nickdata = find_crisprscan_nicksgRNA($data_text, $minEditPos,$maxEditPos,$maxEditDistance,$chosenCutPos,$chosenOrientation,$seq1,$minNickDist,$maxNickDist);
          }
          ($chosenNickSG,$chosenNickSGPos,$chosenNickOrientation,$chosenNickSGDist,$nicksghashpt) = @nickdata;
              
          my @pe3bnickdata = find_pe3b_sgRNA_general($seq1,$minEditPos,$maxEditPos,$maxEditDistance,$chosenCutPos,$chosenOrientation,$seq2,$wtdelcounter,$mutdelcounter,$rgn);
          ($chosenNickSG3b,$chosenNickSGPos3b,$chosenNickOrientation3b,$chosenNickSGDist3b,$nicksghashpt3b) = @pe3bnickdata;

          if ($chosenNickSG3b ne "none found"){ # if a PE3b sgRNA was found
            %nicksghash3b = %$nicksghashpt3b;
            $nicksgtable .= "<tr><th>Nicking sgRNA_Seq</th><th>OnTargetScore</th><th>CutPosition</th><th>Orientation</th><th>DistanceTo_pegRNA_nick</th><th>Type</th><th>Chosen</th></tr>";

            foreach my $k (sort {$nicksghash3b{$a}[3] <=> $nicksghash3b{$b}[3] } keys %nicksghash3b){
              $nicksgtable .= "<tr><td>$k</td><td>NA</td><td>$nicksghash3b{$k}[1]</td><td>$nicksghash3b{$k}[2]</td><td>$nicksghash3b{$k}[3]</td><td>PE3b</td>";
              if ($k eq $chosenNickSG3b){
                $nicksgtable .= "<td>X (PE3b)</td></tr>";
              }
              else { 
                $nicksgtable .= "<td></td></tr>";
              }
            }
          }
                
          if ($chosenNickSG ne "none found"){ # if a PE3 sgRNA was found
            %nicksghash = %$nicksghashpt;
                  
            if ($chosenNickSG3b eq "none found"){
              $nicksgtable .= "<tr><th>Nicking sgRNA_Seq</th><th>OnTargetScore</th><th>CutPosition</th><th>Orientation</th><th>DistanceTo_pegRNA_nick</th><th>Type</th><th>Chosen</th></tr>";
            }
            foreach my $k (sort {$nicksghash{$b}[0] <=> $nicksghash{$a}[0]} keys %nicksghash){
              $nicksgtable .= "<tr><td>$k</tdh><td>$nicksghash{$k}[0]</td><td>$nicksghash{$k}[1]</td><td>$nicksghash{$k}[2]</td><td>$nicksghash{$k}[3]</td><td>PE3</td>";
              if ($k eq $chosenNickSG){
                $nicksgtable .= "<td>X (PE3)</td></tr>";
              }
              else { 
                $nicksgtable .= "<td></td></tr>";
              }
            }
          }
          elsif ($chosenNickSG eq "none found" && $chosenNickSG3b eq "none found") {
            $nicksgtable .= '<tr><br>No suitable secondary nicking sgRNA found</tr>';
          }
          $nicksgtable .= "</table>";
        }
      }
    }
    elsif ($rgn ne "Cas9-NGG"){
      $c->render(text => 'sgRNA result file is only valid with Cas9-NGG. Please do one of the following: 1) Rerun pegFinder with Cas9-NGG as the CRISPR enzyme, 2) Enter a preselected sgRNA, or 3) Rerun pegFinder and do not upload a results file.', status => 200);
    }
  }


##############
  #now that we have a chosen sgRNA, let's proceed
  if ($chosenSG ne "" && $sgfoundStatus == 2){
    #Design the RT templates
    my @rtdata = find_RT($align2,$chosenOrientation,$minEditPos,$maxEditPos,$chosenCutPos,$wtdelcounter);
    my ($rthashpt,$chosenRT,$chosenRTlen,$rttable) = @rtdata;
    my %rthash = %$rthashpt;

    #Design the PBS templates
    my @pbsdata = find_pbs($chosenSG,$gcPctg);
    my ($pbshashpt,$chosenPBS,$chosenPBSlen,$pbstable) = @pbsdata;
    my %pbshash = %$pbshashpt;

    my $extension = $chosenRT.$chosenPBS;
    my $scaffold = "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC";
    
    #########################
    #Compile alternative 3' extensions into temp file for top sgRNAs
    my @tabledata;
    my $fileout;
    if (defined $c_sg && $c_sg ne ""){ # if preselected sgRNA was chosen
      @tabledata = prep_table_chosen($rthashpt,$pbshashpt,$chosenRTlen,$chosenPBSlen,$chosenSG,$maxSgCt,$rgn,$chosenOrientation,$gcPctg,$chosenDisrupt);
      $fileout = $tabledata[0];
    }
    else { 
      @tabledata = prep_table_multi($ranksghashpt,$align2,$minEditPos,$maxEditPos,$wtdelcounter,$maxSgCt,$rgn);
      $fileout = $tabledata[0];
    }
    my $dir = "./downloadFile/";
    my $template = "pegRNAdesignsXXXXXXXX"; # trailing Xs are changed
    my ($fh, $filename) = tempfile( $template, DIR => $dir,SUFFIX => ".txt");
    print $fh $fileout;

    
    #compile the pegRNA and get the oligos
    my $pegRNA = $chosenSG.$scaffold.$chosenRT.$chosenPBS;
    my $oligotable ='<table style ="width:80%; float = left">';
    $oligotable .= "<tr><th>OligoName</th><th>Sequence</th><th>Description</th></tr>";

    if (substr($chosenSG,0,1) eq "G"){
      $oligotable .= '<tr><td>sgF</td><td>cacc'.$chosenSG.'gttttaga</td><td>sgRNA, forward</tr>';
      $oligotable .= '<tr><td>sgR</td><td>tagctctaaaac'.reverse_complement($chosenSG).'</td><td>sgRNA, reverse</td></tr>';
    }
    else {
      $oligotable .= '<tr><td>sgF</td><td>caccg'.$chosenSG.'gttttaga</td><td>sgRNA, forward</tr>';
      $oligotable .= '<tr><td>sgR</td><td>tagctctaaaac'.reverse_complement($chosenSG).'c</td><td>sgRNA, reverse</td></tr>';      
    }
    
    $oligotable .= '<tr><td>scaffF</td><td>GCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCG</td><td>Scaffold, forward (invariant)</td></tr>';
    $oligotable .= '<tr><td>scaffR</td><td>GCACCGACTCGGTGCCACTTTTTCAAGTTGATAACGGACTAGCCTTATTTTAACTTGCTATTTC</td><td>Scaffold, reverse (invariant)</td></tr>';
    $oligotable .= '<tr><td>extensF</td><td>gtgc'.$extension.'</td><td>3\' extension, forward</td></tr>';
    $oligotable .= '<tr><td>extensR</td><td>aaaa'.reverse_complement($extension).'</td><td>3\' extension, reverse</td></tr>';

    #if a PE3 nicking sgRNA was chosen:
    if ($pe3Bool > 0 && (defined $chosenNickSG || defined $chosenNickSG3b)){
      if ($chosenNickSG ne "none found"){
        if (substr($chosenNickSG,0,1) eq "G"){
          $oligotable .= '<tr><td>PE3_sgF</td><td>cacc'.$chosenNickSG.'</td><td>PE3 nick sgRNA, forward</td></tr>';
          $oligotable .= '<tr><td>PE3_sgR</td><td>aaac'.reverse_complement($chosenNickSG).'</td><td>PE3 nick sgRNA, reverse</td></tr>';
        }
        else {
          $oligotable .= '<tr><td>PE3_sgF</td><td>caccg'.$chosenNickSG.'</td><td>PE3 nick sgRNA, forward</td></tr>';
          $oligotable .= '<tr><td>PE3_sgR</td><td>aaac'.reverse_complement($chosenNickSG).'c</td><td>PE3 nick sgRNA, reverse</td></tr>';
        }        
      }
      if ($chosenNickSG3b ne "none found"){
        if (substr($chosenNickSG3b,0,1) eq "G"){
          $oligotable .= '<tr><td>PE3b_sgF</td><td>cacc'.$chosenNickSG3b.'</td><td>PE3b nick sgRNA, forward</td></tr>';
          $oligotable .= '<tr><td>PE3b_sgR</td><td>aaac'.reverse_complement($chosenNickSG3b).'</td><td>PE3b nick sgRNA, reverse</td></tr>';
        }
        else {
          $oligotable .= '<tr><td>PE3b_sgF</td><td>caccg'.$chosenNickSG3b.'</td><td>PE3b nick sgRNA, forward</td></tr>';
          $oligotable .= '<tr><td>PE3b_sgR</td><td>aaac'.reverse_complement($chosenNickSG3b).'c</td><td>PE3b nick sgRNA, reverse</td></tr>';
        }        
      }
    }

    $oligotable .= '</table>';

    #####################
    #Print the output:
    if ($sgfoundStatus == 2){
      if ($pe3Bool == 0){ # did not request PE3 sgRNAs
        $c->render(text=>"
          <head>
          <link rel=\"icon\" href=\"/images/favicon.png\">
          <title>pegFinder: pegRNA designer for Prime Editing</title></head>
          <img src=\"/images/logo4.png\" height=\"100\" width=\"379\" />
          <br><br><hr><br>
          <style>table, th, td {border: 1px solid black;}
          table {border-collapse:collapse;}tr:hover {background-color:#f5f5f5;}
          body {padding-left: 10px; font-family: Arial !important}
          th {background-color: #18a2bf;color: white;}
          th,td{padding-left:10px; padding-right:10px; padding-top: 3px; padding-bottom:3px; text-align: left;}
          </style>

          <b>$rgn</b> chosen as the CRISPR enzyme.<br><br>
          <b><a href=$filename download=\"pegRNAdesigns.txt\" target=\"_blank\">Download full table of candidate 3' extensions & oligo sequences for top sgRNAs.</a></b><br><br>

          <b><u>Recommended selections for pegRNA design</u></b><br>
          
          sgRNA: $chosenSG<br>
          RT template ($chosenRTlen nt): $chosenRT<br>
          PBS ($chosenPBSlen nt): $chosenPBS<br>
          Sense 3' extension: $extension<br>
          Full-length pegRNA: $pegRNA<br><br>

          <p><u><b>Oligos to ligation clone the pegRNA into Addgene #132777 </b></u><br>
          (Note: May need to customize the overhangs for your plasmids.)</p>
          $oligotable <br>

          <hr>
          <p><b><u>Candidate primary sgRNAs</u></b></p>
          $sgtable <br>

          <p><b><u>Candidate RT templates</u></b></p>
          $rttable <br>

          <p><b><u>Candidate PBS sequences</u></b></p>
          $pbstable<br>

          <br><hr>
          <p><b><u>Sequence alignment:</b></u></p>
          Wildtype:<br><tt>$align1f</tt><br><br>Edited:<br><tt>$align2f</tt><br>  

        ");
      }
      elsif (defined $chosenNickSG || defined $chosenNickSG3b){
        if ($chosenNickSG ne "" && $chosenNickSG3b ne ""){ # found both PE3 and PE3b sgRNAs
          $c->render(text=>"
          <head>
          <link rel=\"icon\" href=\"/images/favicon.png\">
          <title>pegFinder: pegRNA designer for Prime Editing</title></head>
          <img src=\"/images/logo4.png\" height=\"100\" width=\"379\" />
          <br><br><hr><br>
          <style>table, th, td {border: 1px solid black;}
          table {border-collapse:collapse;}tr:hover {background-color:#f5f5f5;}
          body {padding-left: 10px; font-family: Arial !important}
          th {background-color: #18a2bf;color: white;}
          th,td{padding-left:10px; padding-right:10px; padding-top: 3px; padding-bottom:3px; text-align: left;}
          </style>

          <b>$rgn</b> chosen as the CRISPR enzyme.<br><br>
          <b><a href=$filename download=\"pegRNAdesigns.txt\" target=\"_blank\">Download full table of candidate 3' extensions & oligo sequences for top sgRNAs.</a></b><br><br>

          <b><u>Recommended selections for pegRNA design</u></b><br>
          sgRNA: $chosenSG<br>
          RT template ($chosenRTlen nt): $chosenRT<br>
          PBS ($chosenPBSlen nt): $chosenPBS<br>
          Sense 3' extension: $extension<br>
          Full-length pegRNA: $pegRNA<br><br>

          PE3 nicking sgRNA: $chosenNickSG<br>
          PE3b nicking sgRNA: $chosenNickSG3b<br><br>

          <p><u><b>Oligos to ligation clone the pegRNA into Addgene #132777 </b></u><br>
          (Note: May need to customize the overhangs for your plasmids.)<br>PE3 nicking sgRNA oligos are designed for standard hU6-sgRNA vectors. <br></p>
          $oligotable <br> <hr>
        
          <p><b><u>Candidate primary editing sgRNAs</u></b></p>
          $sgtable <br>

          <p><b><u>Candidate RT templates</u></b></p>
          $rttable <br>

          <p><b><u>Candidate PBS sequences</b></u></p>
          $pbstable<br>

          <p><b><u>Candidate secondary nicking sgRNAs (for PE3)</b></u></p>
          $nicksgtable<br>

          <br><hr>
          <p><b><u>Sequence alignment:</b></u></p>
          Wildtype:<br><tt>$align1f</tt><br><br>Edited:<br><tt>$align2f</tt><br> 

          ");
        }

        elsif ($chosenNickSG ne "" && $chosenNickSG3b eq ""){ #found PE3 but not PE3b sgRNA
          $c->render(text=>"
          <head>
          <link rel=\"icon\" href=\"/images/favicon.png\">
          <title>pegFinder: pegRNA designer for Prime Editing</title></head>
          <img src=\"/images/logo4.png\" height=\"100\" width=\"379\" />
          <br><br><hr><br>
          <style>table, th, td {border: 1px solid black;}
          table {border-collapse:collapse;}tr:hover {background-color:#f5f5f5;}
          body {padding-left: 10px; font-family: Arial !important}
          th {background-color: #18a2bf;color: white;}
          th,td{padding-left:10px; padding-right:10px; padding-top: 3px; padding-bottom:3px; text-align: left;}
          </style>

          <b>$rgn</b> chosen as the CRISPR enzyme.<br><br>
          <b><a href=$filename download=\"pegRNAdesigns.txt\" target=\"_blank\">Download full table of candidate 3' extensions & oligo sequences for top sgRNAs.</a></b><br><br>

          <b><u>Recommended selections for pegRNA design</u></b><br>
          sgRNA: $chosenSG<br>
          RT template ($chosenRTlen nt): $chosenRT<br>
          PBS ($chosenPBSlen nt): $chosenPBS<br>
          Sense 3' extension: $extension<br>
          Full-length pegRNA: $pegRNA<br><br>

          PE3 nicking sgRNA: $chosenNickSG<br>
          PE3b nicking sgRNA: none found<br><br>

          <p><u><b>Oligos to ligation clone the pegRNA into Addgene #132777 </b></u><br>
          (Note: May need to customize the overhangs for your plasmids.)<br>PE3 nicking sgRNA oligos are designed for standard hU6-sgRNA vectors. <br></p>
          $oligotable <br>

          <hr>

          <p><b><u>Candidate primary editing sgRNAs</u></b></p>
          $sgtable <br>

          <p><b><u>Candidate RT templates</u></b></p>
          $rttable <br>

          <p><b><u>Candidate PBS sequences</b></u></p>
          $pbstable<br>

          <p><b><u>Candidate secondary nicking sgRNAs (for PE3)</b></u></p>
          $nicksgtable<br>

          <br><hr>
          <p><b><u>Sequence alignment:</b></u></p>
          Wildtype:<br><tt>$align1f</tt><br><br>Edited:<br><tt>$align2f</tt><br> 

          ");
        }
        elsif ($chosenNickSG eq "" && $chosenNickSG3b ne ""){#found PE3b but not PE3 sgRNA
          $c->render(text=>"
          <head>
          <link rel=\"icon\" href=\"/images/favicon.png\">
          <title>pegFinder: pegRNA designer for Prime Editing</title></head>
          <img src=\"/images/logo4.png\" height=\"100\" width=\"379\" />
          <br><br><hr><br>
          <style>table, th, td {border: 1px solid black;}
          table {border-collapse:collapse;}tr:hover {background-color:#f5f5f5;}
          body {padding-left: 10px; font-family: Arial !important}
          th {background-color: #18a2bf;color: white;}
          th,td{padding-left:10px; padding-right:10px; padding-top: 3px; padding-bottom:3px; text-align: left;}
          </style>

          <b>$rgn</b> chosen as the CRISPR enzyme.<br><br>
          <b><a href=$filename download=\"pegRNAdesigns.txt\" target=\"_blank\">Download full table of candidate 3' extensions & oligo sequences for top sgRNAs.</a></b><br><br>

          <b><u>Recommended selections for pegRNA design</u></b><br>
          sgRNA: $chosenSG<br>
          RT template ($chosenRTlen nt): $chosenRT<br>
          PBS ($chosenPBSlen nt): $chosenPBS<br>
          Sense 3' extension: $extension<br>
          Full-length pegRNA: $pegRNA<br><br>

          PE3 nicking sgRNA: none found<br>
          PE3b nicking sgRNA: $chosenNickSG3b<br><br>
          
          <p><u><b>Oligos to ligation clone the pegRNA into Addgene #132777 </b></u><br>
          (Note: May need to customize the overhangs for your plasmids.)<br>PE3 nicking sgRNA oligos are designed for standard hU6-sgRNA vectors. <br></p>
          $oligotable <br>

          <hr>

          <p><b><u>Candidate primary editing sgRNAs</u></b></p>
          $sgtable <br>

          <p><b><u>Candidate RT templates</u></b></p>
          $rttable <br>

          <p><b><u>Candidate PBS sequences</b></u></p>
          $pbstable<br>

          <p><b><u>Candidate secondary nicking sgRNAs (for PE3)</b></u></p>
          $nicksgtable<br>

          <br><hr>
          <p><b><u>Sequence alignment:</b></u></p>
          Wildtype:<br><tt>$align1f</tt><br><br>Edited:<br><tt>$align2f</tt><br> 

          ");
        }
      }
      
    }
  }
};

# handling get requests for download file 
get '/downloadFile/:name' => sub{
  my $c = shift;
  my $filename = $c->param('name');
  $c->render_file(filepath => "./downloadFile/$filename.txt",'filename' => 'pegRNAdesigns.txt','cleanup'=> '0');
};

#app->config(hypnotoad => {proxy => 1}};
app->start;
__DATA__

@@ form.html.ep
<!DOCTYPE html>
<html>
  <head>
    <link rel="icon" href="/images/favicon.png">
    <title>pegFinder: pegRNA designer for Prime Editing</title>
    <script src="https://ajax.googleapis.com/ajax/libs/jquery/2.1.1/jquery.min.js"></script>
    <script>
      $(document).ready(function() {
        $('#clear').on('click',function (e) {
            var $el = $('#sgRNAControl');
            $el.wrap('<form>').closest(
                'form').get(0).reset();
            $el.unwrap();
        });
      });
    </script> 
  </head>
  
  <document.getElementById('sgRNA').value='';
  <body>
  <style>
  body {padding-left: 10px; font-family: Arial !important}
  </style>


  <img src="/images/logo4.png" height="100" width="379"/>
    <br><br>
    %= form_for upload => (enctype => 'multipart/form-data') => begin
      Wildtype/reference sequence<br>
      %= text_area 'wildtype',cols => 80, rows => 5, maxlength => 2000, placeholder => 'Wildtype/reference DNA sequence in plaintext. Must share 5\' and 3\' ends with the edited/desired sequence. Recommend 100 bp flanks around target edit. Max total length 500 bp.'
      <br><br>Edited/desired sequence <br>
      %= text_area 'edited', cols => 80, rows => 5, maxlength => 2000, placeholder => 'Edited/desired DNA sequence in plaintext. Must share 5\' and 3\' ends with the wildtype sequence. Recommend 100 bp flanks around target edit. Max total length 500 bp.'
      <br>
      <br><b>1.</b> Find PE3/PE3b secondary nicking sgRNAs: &nbsp;
     
      %=radio_button 'PE3cb', value => 1, checked => 1
      Yes&nbsp;&nbsp;&nbsp
      %=radio_button 'PE3cb', value => 0
      No<p>

      &nbsp;&nbsp;&nbsp * Min nick distance (pegRNA::PE3):&nbsp;
      %= text_field minNickDist => 40, size => 6
      <br>

      &nbsp;&nbsp;&nbsp * Max nick distance (pegRNA::PE3):&nbsp;
      %= text_field maxNickDist => 150, size => 6

      <br><br><b>2.</b> Select CRISPR Enzyme: &nbsp;&nbsp;
      %= select_field enzyme => ['Cas9-NGG','Cas9-NG','Cas9-SpRY']

      <br> <br>
      <b>3.</b> Number of top candidate sgRNAs to export in results table:
      &nbsp;&nbsp;
      %= text_field maxSgCt => 3, size => 6
      

      <br><br><b>4.</b> (Optional) Use a preselected sgRNA: 
      %=text_field 'c_sgRNA', placeholder => 'Preselected sgRNA sequence', size=>29
      <br><br>

      <b>5.</b> (Optional) Incorporate on-target/off-target predictions. (Cas9-NGG only)<p>
      &nbsp;&nbsp;&nbsp;&nbsp;* Upload <a href="https://portals.broadinstitute.org/gpp/public/analysis-tools/sgrna-design" target="_blank"> Broad sgRNA finder </a> or <a href="https://www.crisprscan.org/?page=sequence" target="_blank"> CRISPRscan </a>results: &nbsp;&nbsp;&nbsp;  
      %= file_field 'sgRNA', id => 'sgRNAControl', value => ''
      <button id = "clear" type="button"> Clear </button>
      <br>&nbsp;&nbsp&nbsp;&nbsp;* Use the wildtype sequence as input and report all possible sgRNAs.<br> <br>
      
      

      %= submit_button 'Design pegRNAs'
      &nbsp;&nbsp;&nbsp;&nbsp;&nbsp &nbsp;&nbsp;&nbsp;&nbsp;&nbsp
      <input type="reset" value="Reset" align="right" />
      
    % end

  <!-- Default Statcounter code for Pegfinder
  http://pegfinder.sidichenlab.org/ -->
  <script type="text/javascript">
  var sc_project=12160803; 
  var sc_invisible=1; 
  var sc_security="a9f5687d"; 
  </script>
  <script type="text/javascript"
  src="https://www.statcounter.com/counter/counter.js"
  async></script>
  <noscript><div class="statcounter"><a title="Web Analytics"
  href="https://statcounter.com/" target="_blank"><img
  class="statcounter"
  src="https://c.statcounter.com/12160803/0/a9f5687d/1/"
  alt="Web Analytics"></a></div></noscript>
  <!-- End of Statcounter Code -->

  </body>

  <footer>
    <br> If you use pegFinder in your research, please consider citing:
    <br> <a href="https://www.nature.com/articles/s41551-020-00622-8" target="_blank"> Chow RD*, Chen JS*, Shen J, and Chen S. Nature Biomedical Engineering 2020. </a> <br>
    <br><p> &copy 2019-2020, Laboratory of Sidi Chen</p>
  </footer>
  
</html>