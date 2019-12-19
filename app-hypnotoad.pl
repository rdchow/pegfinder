#!/usr/bin/env perl
use Mojolicious::Lite;
#plugin AutoReload => {};
require "./lib/sub-needleman-wunsch.pl";
require "./lib/sub-sgRNA-finder.pl";
require "./lib/sub-sgRNA-validator.pl";
require "./lib/sub-PBS-RT-design.pl";

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

  # Process form input
 
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

  #Needleman-Wunsch alignment
  my $seq1 = uc($wt);
  my $seq2 = uc($edit);
  my @nw_out = needleman_wunsch($seq1,$seq2);
  my ($align1,$align2,$minEditPos,$maxEditPos,$trimmingStatus5p,$trimmingStatus3p,$wtdelcounter) = @nw_out;
  my $align1f = $align1; #formatted $align1 to have 50 chars per line
  my $align2f = $align2; #formmated $align2 to have 50 chars per line

  $align1f =~ s/(.{50})\K(?=.)/<br>/g;
  $align2f =~ s/(.{50})\K(?=.)/<br>/g;
  
  #if alignments are gapped, kill the process
  if ($trimmingStatus5p != 0 || $trimmingStatus3p != 0){
     return $c->render(text => 'Wildtype and desired sequence are not aligned on the 5\' or 3\' ends.<br>Please trim off the hanging DNA bases and rerun pegFinder (see alignment below). <br><br>If you used Broad sgRNA results, please also rerun the Broad sgRNA finder using the revised wildtype sequence (if modified). <br><hr><br>Wildtype:<br>'.$align1f.'<br><br>Edited:<br>'.$align2f.'<br> ', status => 200);
  }
  #Process the sgRNA input
  #Can either take a file upload or a chosen sgRNA
  my $sgfile = $c->param('sgRNA');
  my $c_sg = $c->param('c_sgRNA');
  my $pe3Bool;
  if (defined $c->param('PE3cb')){
    $pe3Bool = $c->param('PE3cb');
  }
  else {
    $pe3Bool = 0;
  }

  my ($chosenSG,$chosenCutPos,$chosenOrientation,$chosenDistance,$gcPctg, %sghash, $sghashpt, $data_text,$sgfoundStatus, $maxEditDistance, $chosenNickSGPos,$chosenNickSG,$chosenNickSGDist,$nicksghashpt, %nicksghash,$chosenNickOrientation); # declare all the variables

  my $sgtable = '<table style ="width:40%">';
  my $nicksgtable = '<table style ="width:40%">';
  return $c->redirect_to('form') unless ($sgfile ne ""|| $c_sg ne "");
  
  # Check input size
  return $c->render(text => 'Input is too big. Try reducing sequence length.', status => 200)
    if $c->req->is_limit_exceeded;

  $maxEditDistance = 150;


  #if sgRNA was pre-chosen:
  if (defined $c_sg && $c_sg ne ""){
    #if an sgRNA file was uploaded:
    my $boolfile = 0; #tells us if the uploaded sgRNA file is the expected Broad output
    $data_text = $sgfile->asset->slurp;
    $c_sg =~ s/[\r\n]//g;
    $c_sg = uc ($c_sg);
    $chosenSG = $c_sg;
    #if sgRNA is not 20nt, or contains non DNA characters
    if (length($c_sg) != 20 || $c_sg !~ /\A[ACGT]+\z/i){ 
      return $c->render(text => 'Preselected sgRNA is invalid, please only enter the 20nt spacer sequence.', status => 200);
    }
    else {
      #if sgRNA is not in the wildtype sequence:
      if ($seq1 !~ /$c_sg.GG/ && reverse_complement($seq1) !~ /$c_sg.GG/ ){
        return $c->render(text => 'Preselected sgRNA is not present in wildtype sequence (using NGG PAM). Please use a different sgRNA or use the Broad sgRNA finder and upload the results.', status => 200);
      }
      #if sgRNA is present, then go ahead and calculate the sgRNA metrics
      else {
        my @sgdata = process_chosen_sgRNA($c_sg,$seq1,$minEditPos,$maxEditPos,$maxEditDistance,$wtdelcounter); #returns ($chosenCutPos,$chosenOrientation,$chosenDistance,$gcPctg,$sgfoundStatus);
        #if the chosen sgRNA matches to multiple positions, report back
        return $c->render(text => 'Preselected sgRNA matches to multiple positions in wildtype sequence. Please use a different sgRNA or use the Broad sgRNA finder and upload the results.', status => 200)
          unless ($sgdata[0] ne "non-unique");

        #if chosen sgRNA is otherwise okay, report that back
        ($chosenCutPos,$chosenOrientation,$chosenDistance,$gcPctg,$sgfoundStatus) = @sgdata;
        if ($sgfoundStatus == 0){
           return $c->render(text=>'Preselected sgRNA is incompatible with desired edit: it cuts 3\' (downstream) of the alterations. Please use a different sgRNA or use the Broad sgRNA finder and upload the results.');
        }
        elsif ($sgfoundStatus == 1){
           return $c->render(text=>'Preselected sgRNA cuts >'.$maxEditDistance.'nt from the desired edits, and is predicted to be lower efficiency. <br>Please choose a different sgRNA, or use the Broad sgRNA finder and upload the results to the tool.');
        }
        $sgtable = "Not calculated (preselected sgRNA was used).<br>";

        #evaluate whether we need to find PE3 secondary sgRNAs

        if ($pe3Bool == 1){ # check mark for finding PE3 secondary guides
          if (defined $sgfile && $data_text ne ""){
            #Validate the sgRNA file
            $boolfile = validate_sgRNA($data_text,"100");
          }
          else { # if sgRNA file is invalid:
            $c->render(text => 'Broad sgRNA file is absent. Please use the results file from the Broad sgRNA designer, or rerun without searching for PE3 secondary sgRNAs.', status => 200);
          }

          if ($boolfile == 0){
              $c->render(text => 'sgRNA file is invalid. Please use the results file from the Broad sgRNA designer, or rerun without searching for PE3 secondary sgRNAs.', status => 200);
          }
          #find a nicking sgRNA if the box was ticked on and the sgRNA file is valid
          else {
            my @nickdata = find_nicksgRNA($data_text, $minEditPos,$maxEditPos,$maxEditDistance,$chosenCutPos,$chosenOrientation);
            ($chosenNickSG,$chosenNickSGPos,$chosenNickOrientation,$chosenNickSGDist,$nicksghashpt) = @nickdata;
            if ($chosenNickSG ne "none found"){
              %nicksghash = %$nicksghashpt;

              $nicksgtable .= "<tr><th>Nicking sgRNA_Seq</th><th>OnTargetScore</th><th>CutPosition</th><th>Orientation</th><th>DistanceToPrimary_sgRNACut</th><th>Chosen</th></tr>";
              foreach my $k (sort {$nicksghash{$b}[0]<=>$nicksghash{$a}[0]} keys %nicksghash){
                $nicksgtable .= "<tr><td>$k</tdh><td>$nicksghash{$k}[0]</td><td>$nicksghash{$k}[1]</td><td>$nicksghash{$k}[2]</td><td>$nicksghash{$k}[3]</td>";
                if ($k eq $chosenNickSG){
                  $nicksgtable .= "<td>X</td></tr>";
                }
                else { 
                    $nicksgtable .= "<td></td></tr>";
                }
              }
            }
            else {
              $nicksgtable .= '<tr><br>No nicking sgRNA found</tr>';
            }
            $nicksgtable .= "</table>";
          }
        }
      }
    }
  }

  #if the sgRNA was not pre-chosen:
  else {
    my $boolfile = 0; #tells us if the uploaded sgRNA file is the expected Broad output
    $data_text = $sgfile->asset->slurp;

    #if no sgRNA file was uploaded either:
    if(!defined $data_text || $data_text eq ""){
      $c->render(text => 'No preselected sgRNA or Broad sgRNA finder results were included. Please resubmit the form.', status => 200);
    }

    #if an sgRNA file was uploaded:
    if (defined $sgfile && $data_text ne ""){
      #Validate the sgRNA file
      $boolfile = validate_sgRNA($data_text,"100");
    }
    if ($boolfile eq "0"){
        $c->render(text => 'sgRNA file is invalid or absent. Please use the results file from the Broad sgRNA designer on the wildtype sequence, or enter a preselected sgRNA.', status => 200);
    }

    if ($boolfile == 1){ #if valid sgRNA file
      #Choose an sgRNA
      #print "$data_text\n\n$minEditPos\n$maxEditPos\t$maxEditDistance\n";
      my @sgdata = process_sgRNA($data_text, $minEditPos,$maxEditPos,$maxEditDistance,$wtdelcounter);
      
      #if no sgRNAs were found, report back
      return $c->render(text => 'No candidate sgRNAs found.', status => 200)
        unless ($sgdata[0] ne "none");
      
      #Otherwise, if an sgRNA was found:
      ($sghashpt,$chosenSG,$chosenCutPos,$chosenOrientation,$chosenDistance,$gcPctg,) = @sgdata;
      %sghash = %$sghashpt;
      $sgfoundStatus = 2;
      $sgtable .= "<tr><th>sgRNA_Seq</th><th>OnTargetScore</th><th>CutPosition</th><th>Orientation</th><th>DistanceToEditStart</th><th>Chosen</th></tr>";
      foreach my $k (sort {$sghash{$b}[0]<=>$sghash{$a}[0]} keys %sghash){
        $sgtable .= "<tr><td>$k</tdh><td>$sghash{$k}[0]</td><td>$sghash{$k}[1]</td><td>$sghash{$k}[2]</td><td>$sghash{$k}[3]</td>";
        if ($k eq $chosenSG){
          $sgtable .= "<td>X</td></tr>";
        }
        else { 
            $sgtable .= "<td></td></tr>";
        }
      }
      $sgtable .= "</table>";

      #find a nicking sgRNA if the box was ticked on
      if ($pe3Bool == 1){
        my @nickdata = find_nicksgRNA($data_text, $minEditPos,$maxEditPos,$maxEditDistance,$chosenCutPos,$chosenOrientation);
        ($chosenNickSG,$chosenNickSGPos,$chosenNickOrientation,$chosenNickSGDist,$nicksghashpt) = @nickdata;
        

        if ($chosenNickSG ne "none found"){
          %nicksghash = %$nicksghashpt;
          $nicksgtable .= "<tr><th>Nicking sgRNA_Seq</th><th>OnTargetScore</th><th>CutPosition</th><th>Orientation</th><th>DistanceToPrimary_sgRNACut</th><th>Chosen</th></tr>";
          foreach my $k (sort {$nicksghash{$b}[0]<=>$nicksghash{$a}[0]} keys %nicksghash){
            $nicksgtable .= "<tr><td>$k</tdh><td>$nicksghash{$k}[0]</td><td>$nicksghash{$k}[1]</td><td>$nicksghash{$k}[2]</td><td>$nicksghash{$k}[3]</td>";
            if ($k eq $chosenNickSG){
              $nicksgtable .= "<td>X</td></tr>";
            }
            else { 
                $nicksgtable .= "<td></td></tr>";
            }
          }
        }
        else {
          $nicksgtable .= '<tr><br>No nicking sgRNA found</tr>';
        }
      }
      $nicksgtable .= "</table>";

    }
  }

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
    
    #compile the pegRNA and get the oligos
    my $pegRNA = $chosenSG.$scaffold.$chosenRT.$chosenPBS;
    my $oligotable ='<table style ="width:80%; float = left">';
    $oligotable .= "<tr><th>OligoName</th><th>Sequence</th><th>Description</th></tr>";

    $oligotable .= '<tr><td>sgF</td><td>caccg'.$chosenSG.'gttttaga</td><td>sgRNA, forward</tr>';
    $oligotable .= '<tr><td>sgR</td><td>tagctctaaaac'.reverse_complement($chosenSG).'c</td><td>sgRNA, reverse</td></tr>';
    $oligotable .= '<tr><td>scaffF</td><td>GCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCG</td><td>Scaffold, forward (invariant)</td></tr>';
    $oligotable .= '<tr><td>scaffR</td><td>GCACCGACTCGGTGCCACTTTTTCAAGTTGATAACGGACTAGCCTTATTTTAACTTGCTATTTC</td><td>Scaffold, reverse (invariant)</td></tr>';
    $oligotable .= '<tr><td>extensF</td><td>gtgc'.$extension.'</td><td>3\' extension, forward</td></tr>';
    $oligotable .= '<tr><td>extensR</td><td>aaaaa'.reverse_complement($extension).'</td><td>3\' extension, reverse</td></tr>';

    #if a PE3 nicking sgRNA was chosen:
    if ($pe3Bool == 1 && defined $chosenNickSG){
      if ($chosenNickSG ne "none found"){
        $oligotable .= '<tr><td>nick_sgF</td><td>caccg'.$chosenNickSG.'</td><td>PE3 nick sgRNA, forward</td></tr>';
        $oligotable .= '<tr><td>nick_sgF</td><td>aaac'.reverse_complement($chosenNickSG).'c</td><td>PE3 nick sgRNA, reverse</td></tr>';
      }
    }

    $oligotable .= '</table>';

    #Print the output:
    if (defined $sgfile && ($data_text ne "" || $sgfoundStatus == 2)){
      if ($pe3Bool == 0){
        $c->render(text=>"
        <img src=\"/images/logo4.png\" height=\"100\" width=\"379\" />
        <br><br><hr><br>
        <style>table, th, td {border: 1px solid black;}
        table {border-collapse:collapse;}tr:hover {background-color:#f5f5f5;}
        body {padding-left: 10px; font-family: Arial !important}
        th {background-color: #18a2bf;color: white;}
        th,td{padding-left:10px; padding-right:10px; padding-top: 3px; padding-bottom:3px; text-align: left;}
        </style>
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
        br><hr>
        <p><b><u>Sequence alignment:</b></u></p>
        Wildtype:<br><tt>$align1f</tt><br><br>Edited:<br><tt>$align2f</tt><br>  

        ");
      }
      elsif (defined $chosenNickSG){
        if ($chosenNickSG ne ""){
          $c->render(text=>"
        <img src=\"/images/logo4.png\" height=\"100\" width=\"379\" />
        <br><br><hr><br>
        <style>table, th, td {border: 1px solid black;}
        table {border-collapse:collapse;}tr:hover {background-color:#f5f5f5;}
        body {padding-left: 10px; font-family: Arial}
        th {background-color: #18a2bf;color: white;}
        th,td{padding-left:10px; padding-right:10px; padding-top: 3px; padding-bottom:3px; text-align: left;}
        </style>
        <b><u>Recommended selections for pegRNA design</u></b><br>
          sgRNA: $chosenSG<br>
          RT template ($chosenRTlen nt): $chosenRT<br>
          PBS ($chosenPBSlen nt): $chosenPBS<br>
          Sense 3' extension: $extension<br>
          Full-length pegRNA: $pegRNA<br><br>

          PE3 nicking sgRNA: $chosenNickSG<br><br>
          
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

#app->config(hypnotoad => {proxy => 1}};
app->start;
__DATA__

@@ form.html.ep
<!DOCTYPE html>
<html>
  <head>
  <title>pegFinder: pegRNA designer for Prime Editing</title></head>
  <document.getElementById('sgRNA').value='';
  <body>
  <style>
  body {padding-left: 10px; font-family: Arial !important}
  </style>


  <img src="/images/logo4.png" height="100" width="379"/>
    <br><br>
    %= form_for upload => (enctype => 'multipart/form-data') => begin
      Wildtype/reference sequence<br>
      %= text_area 'wildtype',cols => 80, rows => 5, maxlength => 10000, placeholder => 'Wildtype/reference DNA sequence in plaintext. Must share 5\' and 3\' ends with the edited/desired sequence. Recommend >100 bp flanks around target edit.'
      <br><br>Edited/desired sequence <br>
      %= text_area 'edited', cols => 80, rows => 5, maxlength => 10000, placeholder => 'Edited/desired DNA sequence in plaintext. Must share 5\' and 3\' ends with the wildtype sequence. Recommend >100 bp flanks around target edit.'
      <br>
      <br>Upload <a href="https://portals.broadinstitute.org/gpp/public/analysis-tools/sgrna-design"> Broad sgRNA finder </a> results (tab-delimited *.txt): &nbsp;&nbsp;&nbsp;&nbsp;&nbsp  
      %= file_field 'sgRNA'
      <br>Use the wildtype sequence as input and make sure to report unpicked sequences. <br>
      (<b>Required</b> for designing PE3 secondary nicking sgRNA) <br> <br>
      
      (Optional) Use a preselected sgRNA: 
      %=text_field 'c_sgRNA', placeholder => 'Preselected sgRNA sequence', size=>29
      <br><br>
      Find PE3 secondary nicking sgRNAs: &nbsp;
     
      %=radio_button 'PE3cb', value => 1
      Yes&nbsp;&nbsp;&nbsp
      %=radio_button 'PE3cb', value => 0
      No<br><br>
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
    <p> &copy 2019, Laboratory of Sidi Chen</p>
  </footer>
  
</html>