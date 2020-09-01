use strict;
use warnings;

#############################################################################
# called by annotate_rna; opens and parses output from tRNAscan-SE

sub parse_trnascan{
  my( $infile, $file_ref, $project, $query_taxon, $out_directory, $missed_by_blast_ref, $missed_by_tscan_ref, $missed_ref, $rna_html_files_Href, $all_found_ref, $tscan_html_files_Href, $tscan_html_paths ) = @_;
  my $outfile;
  my $tscan_out = $project . "_tRNAscan.out";
  my $tscan_struct_out = $project . "_tRNAscan.struct";
  my $tscan = 1;
  my %tscan_results;

  my %trna_letters = (
		      Phe => "F",
		      Leu => "L",
		      Ile => "I",
		      Met => "M",
		      Val => "V",
		      Ser => "S",
		      Pro => "P",
		      Thr => "T",
		      Ala => "A",
		      Tyr => "Y",
		      His => "H",
		      Gln => "Q",
		      Asn => "N",
		      Lys => "K",
		      Asp => "D",
		      Glu => "E",
		      Cys => "C",
		      Trp => "W",
		      Arg => "R",
		      Gly => "G",
		      Undet => "UNDETERMINED tRNA"
		     );
  
  # if necessary, run tRNA-scan
  unless( -e "$out_directory/$tscan_out" ){
    print "\n\n\t\'$tscan_out\' not found, running tRNA-scan now (this might take a while) . . .\n";
    system( "tRNAscan/tRNAscan-SE -q -O -o $out_directory/$tscan_out -f $out_directory/$tscan_struct_out $infile" );
  }

  # parse tRNA-scan output
  print "\nParsing tRNA-scan output . . . \n";
  open( TSCAN, "$out_directory/$tscan_out" ) || die "Couldn't open $out_directory/$tscan_out: $!\n"; 
  
  while( <TSCAN> ){
    if( /^($query_taxon)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/ ){
      $5 eq "Sup" and next;
      $tscan_html_files_Href->{$5} = $out_directory . "/" . $5 . "_trnascan.html";
      $tscan_html_paths->{$5}      = "\.\.\/$project" . "_out/$5" . "_trnascan.html";
    }
  }
  close TSCAN;
  
  # array of all possible tRNAs (Phe, Cys, etc.)
  my @trna_3Letter = keys( %trna_letters );
  
  # array of tRNAs found by tRNA-scan (Phe, Cys, etc.)
  my @found_by_tscan = keys( %{$tscan_html_files_Href} );
  
  # array of tRNAs found by BLAST (Phe, Cys, etc.)
  my @found_by_blast = keys( %{$rna_html_files_Href} );
  
  # @missed_by_tscan will contain tRNAs missed by tRNA-scan by calculating difference between 
  # @found_by_tscan and @trna_3Letter
  my %found; #lookup table
  
  #build lookup table
  @found{@found_by_tscan} = ();
		
  foreach( @trna_3Letter ){
    push( @{$missed_by_tscan_ref}, $_ ) unless exists $found{$_};
  }

  # @missed will contain contain AAs missed by both BLAST and tRNA-scan; this is the intersection of @missed_by_blast and @missed_by_tscan
  my @union = my @diff = ();
  my %union = my %isect = ();
  my $e;
  foreach $e ( @{$missed_by_blast_ref} ){ $union{$e} = 1 }
  foreach $e ( @{$missed_by_tscan_ref} ){
    if( $union{$e} ){ $isect{$e} = 1 }
    $union{$e} = 1;
  }
  @{$missed_ref} = keys %isect;
  @union         = keys %union;
  
  # @all_rna will contain AAs found by BLAST and tRNA-scan; this the union of @found_by_blast and @found_by_tscan; 
  my @isect_found = my @diff_found  = ();
  my %union_found = my %isect_found = ();
  my $f;
  foreach $f ( @found_by_blast ){ $union_found{$f} = 1 }
  foreach $f ( @found_by_tscan ){
    if( $union_found{$f} ){ $isect_found{$f} = 1 }
    $union_found{$f} = 1;
  }

  @{$all_found_ref} = keys %union_found;
  @isect_found      = keys %union_found;
  

  # make HTML output for each amino acid, prepending each with HTML header
  foreach( keys(%{$tscan_html_files_Href}) ){
    print_trnascan_output_header( $_, $tscan_html_files_Href->{$_}, $query_taxon, $project );
  }

  # parse tabular tRNAScan-SE output
  open( TSCAN, "$out_directory/$tscan_out" ) || die "Couldn't open $out_directory/$tscan_out: $!\n"; 
  while( <TSCAN> ){
    /\tSup\t/ and next;
    /^\-{8}/ and next;

    my( $species, $trna_num, $start, $end, $aa, $anti, $intron_start, $intron_end, $cove_score ) = split( /\s+/, $_ );
    $trna_num eq "tRNA" and next;
    print "\t$aa\n";
    
    # write to tscan (bottom) panel of HTML summary page
    print_trnascan( $project, $missed_ref, \%trna_letters, $file_ref, $query_taxon, $out_directory, $start, $end, $aa, $anti, $intron_start, $intron_end, $cove_score, $tscan );
  }
  
  close TSCAN;
  
}
#############################################################################
# called by parse_trnascan (above)

sub print_trnascan_output_header{
  my( $aa, $outfile, $query_taxon, $project ) = @_;
  
  open( HTML,">$outfile" ) || die "Couldn't open $outfile: $!\n";
  
  print HTML "<html>\n<head>\n<title>$query_taxon: $aa</title>\n</head>\n";
  print HTML "<style>\n";
  print HTML 
		"p1{\n".
		"\tfont-family: Courier, monospace;\n".
		"\twhite-space: nowrap;\n".
		"\tfont-style: normal;\n".
		"\tfont-variant: normal;\n".
		"\tfont-weight: normal;\n".
		"\tfont-size: 14;\n".
		"\tword-spacing: normal;\n".
		"\tletter-spacing: normal;\n".
		"\ttext-align: left;\n".
		"\tcolor: black ;\n".
		"\talign: left;\n".
		"}\n\n";
  print HTML 
		"p2{\n".
		"\tfont-family: Courier, monospace;\n".
		"\twhite-space: nowrap;\n".
		"\tfont-style: normal;\n".
		"\tfont-variant: normal;\n".
		"\tfont-weight: normal;\n".
		"\tfont-size: 14;\n".
		"\tword-spacing: normal;\n".
		"\tletter-spacing: normal;\n".
		"\ttext-align: left;\n".
		"\tcolor: purple ;\n".
		"\talign: left;\n".
		"}\n\n";
  print HTML 
		"h1{\n".
		"\tfont-family: Courier, monospace;\n".
		"\tfont-style: normal;\n".
		"\tfont-variant: normal;\n".
		"\tfont-weight: normal;\n".
		"\tfont-size: 24;\n".
		"\ttext-align: left;\n".
		"\tcolor: blue;\n".
		"\talign: left;\n".
		"}\n";
  print HTML "</style>\n\n\n\n\n";
  close( HTML );
}
#############################################################################
# called by parse_trnascan (above); does multiple things:
  # 1. extracts nt sequence of tRNA and corresponding upstream and downstream sequences
  # 2. formats and prints tRNA sequence found by tscan for html output

sub print_trnascan{
  my( $project, $missed_ref, $trna_letters_ref, $file_ref, $query_taxon, $out_directory, $start, $end, $aa, $anti, $intron_start, $intron_end, $cove_score, $tscan ) = @_;
  my( $strand, $length, $trna, $up_offset, $down_offset, $upstream, $downstream, $sense, $antisense, $outfile, $coords );
  my $flanking_length = 60; #in addition to gene, extract 60 nt each of up- and downstream sequence
  my $x = 0;
  my $y = 0;
  my $z = 0;
  my $a = 0;

  #for tRNAs on minus strand, extract tRNA and upstream and downstream sequence
  if( ($end - $start) < 0 ){
    $strand = "-";
    $length = $start - $end + 1;
    $start = $end;
    $end = $start + $length - 1;
    $trna = substr($$file_ref, $start-1, $length);
    $up_offset = $end;
    $upstream = substr($$file_ref, $up_offset, $flanking_length);
    $down_offset = ($start-1)-$flanking_length;
    $downstream = substr($$file_ref, $down_offset, $flanking_length);
    $sense  = "$downstream"."$trna"."$upstream";
    $downstream = complement( $downstream );
    $upstream = complement( $upstream );
    $trna = complement( $trna );
    $sense = complement( $sense );
  
  #for trna's on + strand, extract trna sequence, upstream and downstream sequences
  }else{
    $strand = "+";
    $length = $end - $start + 1;
    $trna = uc substr($$file_ref, $start-1, $length);
    $up_offset = ($start-1)-$flanking_length;
    $upstream = substr($$file_ref, $up_offset, $flanking_length);
    
    # this "if" is for genes at the beginning of the genome, e.g., if the genome is oriented
    # to cox1; if the gene starts at position 0 of the genome, the substr function will
    # (properly) extract the last 60 nt of the genome; if, however, the gene starts at
    # position 9, for example, the substr function will extract the last 50 nt of the genome,
    # that's it, it won't grab the 10 nt at the beginning of the genome immediately preceding 
    # the gene
    if( $start < $flanking_length ){
      my $threePrime_up_length = $flanking_length - (length $upstream);
      my $threePrime_up_offset = $start - $threePrime_up_length - 1;
      my $threePrime_upstream = substr( $$file_ref, $threePrime_up_offset, $threePrime_up_length );
      $upstream = "$upstream" . "$threePrime_upstream";
    }
    $down_offset = $end;
    $downstream = substr($$file_ref, $down_offset, $flanking_length);
    $sense = "$upstream"."$trna"."$downstream";
  }
  
  $antisense = $sense;
  $antisense =~ tr/ACGTacgt/TGCAtgca/;
  
  #Layout of query nt printout:
  #|---------------------$a---------------|AA sequence#
  #|--------$x-------|$taxon_name|---$z---|AA sequence#
  
  my $prefix = 60; #number of characters preceding AA sequence
  $y = length( $query_taxon ); #number of characters in taxon name
  $z = 6; #number of characters between taxon name and AA sequence
  $x = ($prefix -$y) - $z; #number of characters preceding taxon name
  $a = $x + $y + $z;	
  
  $query_taxon =~ s/ /&nbsp;/g;
  
  $outfile = "$out_directory/" . $aa . "_trnascan.html";
   
  if( $strand eq "+" ){
    $coords = "$start..$end";
  }else{
    $coords = "$end..$start";	
  }
  
  open( HTML,">>$outfile" ) || die "Couldn't open $outfile: $!\n";
  print_trn_title( $project, $query_taxon, $trna_letters_ref, $aa, $anti, $length, $coords );
  
  #print tRNA coords
  print HTML "<p1>", '&nbsp;' x $a, "<font color=333399>$coords (COVE SCORE: $cove_score)&nbsp;<\/font>", "</p1><br>\n";
  #print "|" below start coord
  print HTML "<p1>", '&nbsp;' x $a, "&#124;", "<\/p1><br>\n";
  #print query taxon
  print HTML "<p2>", '&nbsp;' x $x, "<font COLOR=333399 FACE=COURIER font-weight=700>$query_taxon<\/font>";
  #print tRNA sequence
  print HTML '&nbsp;' x $z, "<font COLOR=333399 FACE=COURIER font-weight=700>", uc $trna, "<\/p2><br><\/font>\n";
  
  tscan_rna_print_nt( $sense, $antisense, $strand, $start, $end, $upstream, $trna, $downstream, $file_ref, $flanking_length, $tscan );
  
  close HTML;
}
#############################################################################
# called by extract_gapped_nt (BLAST OUTPUT) and &print_trnascan (TSCAN OUTPUT); formats and prints 
# nucleotide sequence of concatenated upstream, gene, and downstream sequence below the hit sequences

sub tscan_rna_print_nt {
  my( $sense, $antisense, $strand, $query_start, $query_end, $upstream, $gene, $downstream, $file_ref, $flanking_length, $tscan ) = @_;
  my $i;
  my $position = 12;
  my $position_string = "";
  my $genome = $$file_ref;
  my $HTML_sense;
  my @POSITIONS;
  my $name_color = "333399";
  my $query_color = "333399";
  
  print HTML "<p2>\n";
  
  if( $tscan ){
    #build string of nucleotide coordinates, to be printed below the nucleotide sequences
    for($i = 0; $i < length($sense); $i++){
      if( $strand eq "+" ){
	$position = ($query_start-60) + $i;
	if( $position <= 0 ){
	  $position = ((length "$genome") + $position);
	}
	if( ($i % 12) == 0 ){
	  $position_string .= $position . '&nbsp;' x (12 - length "$position");
	}
      }elsif( $strand eq "-" ){
	$position = ($query_end+60) - $i;               
	if( ($i % 12) == 0 ){
	  $position_string .= $position . '&nbsp;' x (12 - length "$position");
	}
      }
    }
  }
  # build sense strand, formatted for HTML output, with upstream, downstream, and nonsense 
  # strands all gray, and the BLAST hit sequence purple upstream sense strand
  $HTML_sense .= "<font color=gray font-weight=1000>$upstream<\/font>";
  
  #gene
  $HTML_sense .= "<font color=$query_color font-weight=1000>$gene<\/font>";
  
  #downstream sense strand
  $HTML_sense .= "<font color=gray font-weight=1000>$downstream<\/font>";
  
  #print sense strand
  print HTML $HTML_sense;
  #print HTML $HTML_sense;
  print HTML "<\/p2>\n";
  
  #print antisense strand
  print HTML "<p2><font color=gray font-weight=1000>" . $antisense . "<\/font><\/p2>\n";
  
  if( $tscan ){
    #print "^" characters below nt sequences
    print HTML "<p2>";
    for($i = 0; $i < length $sense; $i+=12){
      print HTML "<font color=663399 font-weight=1000>", "&#94;", '&nbsp;' x 11, "<\/font>";
    }
    print HTML"<\/p2>\n";
  }
  #print coordinates below nt sequences
  print HTML "<p2><font color=$query_color font-weight=1000>", $position_string, "<\/font><\/p2>\n";
  print HTML "<br><br><br>\n";
  print HTML "<hr width=10000 size=4>\n\n\n";
}
#############################################################################
# called by print_trnascan
sub print_trn_title{
  my( $project, $query_taxon, $trna_letters, $aa, $anti, $length, $coords ) = @_;
  my $gene_name;
  
  if( $aa eq "Undet" ){
    $gene_name = $trna_letters->{$aa};
    print HTML "<table valign=\"middle\">\n\t<tr nowrap valign=\"middle\">\n";
    print HTML "\t\t<td><font face=\"Times\">", '&nbsp;' x 85, "<\/font><\/td>\n";
    print HTML "\t\t<td><table><form METHOD=\"GET\"><INPUT TYPE=\"button\" VALUE=\"Annotate\" onClick=\"window.open(\'http:\/\/127.0.0.1\/cgi-bin\/trna_make_annotation_form.cgi?gene_name=$gene_name&gene_coords=$coords&product=tRNA-$aa&note=anticodon:$anti&taxon=$query_taxon&project=$project\',\'Sequin_form\',\'width=550,height=600,resizable=yes,menubar=no\')\"><\/form><\/table><\/td>\n";
    print HTML "\t\t<td><font face=\"Courier\" color=\"blue\" size=\"6\">&nbsp;&nbsp;&nbsp;$gene_name&nbsp;($length&nbsp;nt)<\/font><\/td>\n";
    print HTML "\t<\/tr>\n<\/table>\n<br><br>\n";

  }elsif( $aa eq "Cys" ){
    $gene_name = "trn" . $trna_letters->{$aa} . "-" . $anti;
    print HTML "<table valign=\"middle\">\n\t<tr nowrap valign=\"middle\">\n";
    print HTML "\t\t<td><font face=\"Times\">", '&nbsp;' x 85, "<\/font><\/td>\n";
    print HTML "\t\t<td><table><form METHOD=\"GET\"><INPUT TYPE=\"button\" VALUE=\"Annotate\" onClick=\"window.open(\'http:\/\/127.0.0.1\/cgi-bin\/trnC_make_annotation_form.cgi?gene_name=$gene_name&gene_coords=$coords&product=tRNA-$aa&note=anticodon:$anti&taxon=$query_taxon&project=$project\',\'Sequin_form\',\'width=550,height=600,resizable=yes,menubar=no\')\"><\/form><\/table><\/td>\n";
    print HTML "\t\t<td><font face=\"Courier\" color=\"blue\" size=\"6\">&nbsp;&nbsp;&nbsp;$gene_name&nbsp;($length&nbsp;nt)<\/font><\/td>\n";
    print HTML "\t<\/tr>\n<\/table>\n<br><br>\n";

  }else{
    $gene_name = "trn" . $trna_letters->{$aa} . "-" . $anti;
    print HTML "<table valign=\"middle\">\n\t<tr nowrap valign=\"middle\">\n";
    print HTML "\t\t<td><font face=\"Times\">", '&nbsp;' x 85, "<\/font><\/td>\n";
    print HTML "\t\t<td><table><form METHOD=\"GET\"><INPUT TYPE=\"button\" VALUE=\"Annotate\" onClick=\"window.open(\'http:\/\/127.0.0.1\/cgi-bin\/trna_make_annotation_form.cgi?gene_name=$gene_name&gene_coords=$coords&product=tRNA-$aa&note=anticodon:$anti&taxon=$query_taxon&project=$project\',\'Sequin_form\',\'width=550,height=600,resizable=yes,menubar=no\')\"><\/form><\/table><\/td>\n";
    print HTML "\t\t<td><font face=\"Courier\" color=\"blue\" size=\"6\">&nbsp;&nbsp;&nbsp;$gene_name&nbsp;($length&nbsp;nt)<\/font><\/td>\n";
    print HTML "\t<\/tr>\n<\/table>\n<br><br>\n";
  }
}
#############################################################################
  
1;
