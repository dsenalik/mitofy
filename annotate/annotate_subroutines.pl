use strict;
use warnings;

#############################################################################
# Called by parse_matching_sequence (above), formats and prints subject AA
# sequence to HTML output
#############################################################################
sub print_sbjct_AA{
  my ( $sbjct_taxon, $AA_sequence, $sbjct_start, $sbjct_end, $query_start, $query_end, $strand, $percent_match, $gene_length ) = @_;
  my $query_coords;
  my( $x, $y, $z, $a );
  
  #Layout of AA printout:
  #|---------------------$a---------------|AA sequence#
  #|--------$x-------|$taxon_name|---$y---|AA sequence#

  my $prefix = 60; #number of characters preceding AA sequence
  $y = length $sbjct_taxon; #number of characters in taxon name
  $z = 6; #number of characters between taxon name and AA sequence
  $x = ($prefix -$y) - $z; #number of characters preceding taxon name
  $a = $x + $y +$z +1;	
  
  if( $strand eq "+" ){
    $query_coords = "$query_start..$query_end";
  }else{
    $query_coords = "$query_end..$query_start";	
  }
  
  print HTML "<p1>", '&nbsp;' x $a, "$sbjct_start..$sbjct_end&nbsp;//&nbsp;<font color=663399>$query_coords<\/font>", " (", "match length = ", $query_end-$query_start+1, " nt, ", "$percent_match\% identity within match, $strand strand)", "</p1><br>\n";
  print HTML "<p1>", '&nbsp;' x $a, "&#124;", "<\/p1><br>\n";
  
  my @AA = split( //, $AA_sequence );
  $AA_sequence = join( '&nbsp;&nbsp;', @AA );
  
  $sbjct_taxon =~ s/ /&nbsp;/g;
  $sbjct_taxon =~ s/_/&nbsp;/g;
  $sbjct_taxon =~ s/-/&#8211;/g;
  
  print HTML "<p1>", '&nbsp;' x $x, "<font COLOR=663399 FACE=COURIER>$sbjct_taxon&nbsp;<\/font>";
  print HTML '&nbsp;' x $z, "<font COLOR=CC66CC FACE=COURIER>", $AA_sequence, "<\/p1><br>\n";
}
#############################################################################
# Called by parse_matching_sequence (above), formats and prints query AA
# sequence to HTML output
#############################################################################
sub print_query_AA{
  my ( $query_taxon, $AA_sequence, $query_start, $query_end ) = @_;
  my( $x, $y, $z, $a );
  
  #Layout of AA printout:
  #|---------------------$a---------------|AA sequence#
  #|--------$x-------|$taxon_name|---$y---|AA sequence#
  
  my $prefix = 60;             #number of characters preceding AA sequence
  $y = length( $query_taxon ); #number of characters in taxon name
  $z = 6;                      #number of characters between taxon name and AA sequence
  $x = ($prefix -$y) - $z;     #number of characters preceding taxon name
  $a = $x + $y +$z +1;
  
  my @AA = split( //, $AA_sequence );
  $AA_sequence = join( '&nbsp;&nbsp;', @AA );
  
  $query_taxon =~ s/ /&nbsp;/g;
  
  print HTML "<p2>", '&nbsp;' x $x, "<font COLOR=663399 FACE=COURIER>$query_taxon&nbsp;<\/font>";
  print HTML '&nbsp;' x $z, "<font COLOR=663399 FACE=COURIER>", $AA_sequence, "<\/p2><br><\/font>\n";
}
####################################################################################################
# Called by Main, extract nucleotide sequence from query genome for each BLAST hit
####################################################################################################
sub extract_nt_sequence{
  my ( $query_match, $query_start, $query_end, $strand, $file_ref, $rna_bool ) = @_;
  my( $up_offset, $down_offset, $upstream, $downstream, $sense, $antisense );

  my $genome_length = length $$file_ref;
  my $hit_length    = $query_end - $query_start + 1;
  my $hit_offset    = $query_start - 1;
  my $query_nt      = substr( $$file_ref, $hit_offset, $hit_length ); #genomic nt sequence for the BLAST hit
  my $flanking_length = 60; #extract 60 nt each of up- and downstream sequence
  
  # IMPORTANT NOTE ABOUT HOW SUBSTR WORKS: if the BLAST hit starts at position 10 of the genome, and you ask
  # substr to extract the 60 nt upstream of this, it will extract the last 50 nt of the genome only, i.e.,
  # it will not also grab the 10 nt at the beginning of the genome immediately preceding the hit, so I have
  # to grab those separately

  if( $strand eq "+" ){

    #EXTRACT UPSTREAM SEQUENCE
    $up_offset = $query_start - $flanking_length - 1;
    $upstream  = substr( $$file_ref, $up_offset, $flanking_length );

    #check whether BLAST hit is sufficiently close to beginning of genome that extracting the
    #full 60 nt of upstream sequence will go past the beginning of the genome
    if( ($query_start - $flanking_length) < 0 ){
      my $additional_upstream_length  = $flanking_length - (length $upstream);
      my $additional_upstream_offset = $query_start - $additional_upstream_length - 1;
      my $additional_upstream  = substr( $$file_ref, $additional_upstream_offset, $additional_upstream_length );
      $upstream = "$upstream" . "$additional_upstream";
    }
    
    #EXTRACT DOWNSTREAM SEQUENCE
    $down_offset = $query_end;
    $downstream  = substr( $$file_ref, $down_offset, $flanking_length );

    #check whether BLAST hit is sufficiently close to end of genome that extracting the
    #full 60 nt of downstream sequence will go past the end of the genome
    if( ($query_end + $flanking_length - 1) > $genome_length ){
      my $additional_downstream_length = $flanking_length - (length $downstream);
      my $additional_downstream_offset = 0;
      my $additional_downstream  = substr( $$file_ref, $additional_downstream_offset, $additional_downstream_length );
      $downstream = "$downstream" . "$additional_downstream";
    }
    
    $sense = $upstream . $query_nt . $downstream;

  }elsif( $strand eq "-" ){
    #print "qstart $query_start qend $query_end\n";
    #EXTRACT UPSTREAM SEQUENCE (equivalent to 5' UTR for a gene on the reverse strand)
    $up_offset = $query_end;
    $upstream  = substr( $$file_ref, $up_offset, $flanking_length );
    
    #check whether BLAST hit is sufficiently close to end of genome that extracting the
    #full 60 nt of "upstream" (its a "-" strand hit) sequence will go past the end of the genome
    if( ($query_end + $flanking_length - 1) > $genome_length ){
      my $additional_upstream_length = $flanking_length - (length $upstream);
      my $additional_upstream_offset = 0;
      my $additional_upstream  = substr( $$file_ref, $additional_upstream_offset, $additional_upstream_length );
      $upstream = "$upstream" . "$additional_upstream";
    }
    
    #EXTRACT DOWNSTREAM SEQUENCE
    $down_offset = $query_start - $flanking_length - 1;
    $downstream  = substr( $$file_ref, $down_offset, $flanking_length );
    
    #check whether BLAST hit is sufficiently close to beginning of genome that extracting the
    #full 60 nt of upstream sequence will go past the beginning of the genome
    if( ($query_start - $flanking_length) < 0 ){
      my $additional_downstream_length  = $flanking_length - (length $downstream);
      my $additional_downstream_offset = $query_start - $additional_downstream_length - 1;
      my $additional_downstream  = substr( $$file_ref, $additional_downstream_offset, $additional_downstream_length );
      $downstream = "$downstream" . "$additional_downstream";
    }
    
    #have to complement these individually because $gene is already complemented, so can't 
    #concatenate them and complement together
    $query_nt   = complement( $query_nt );
    $upstream   = complement( $upstream );
    $downstream = complement( $downstream );

    $sense = $upstream . $query_nt . $downstream;

  }
  
  $antisense = $sense;
  $antisense =~ tr/53ACGTUMRWSYKVHDBNacgtumrwsykvhdbn/35TGCAAKYWSRMBDHVNtgcaakywsrmbdhvn/;

  if( $rna_bool ){
    rna_print_nt( $sense, $antisense, $strand, $query_start, $query_end, $upstream, $query_nt, $downstream, $file_ref, $flanking_length );

  # does query AA sequence have gaps? if yes, print_nt_with_gaps, so you can insert gap
  # characters into the html-formatted nt sequence printed below the AA sequence
  }elsif( $query_match =~ /-/ ){
    # print $query_match, "\n";
    print_nt_with_gaps( $query_match, $sense, $antisense, $strand, $query_start, $query_end, $upstream, $query_nt, $downstream, $file_ref, $flanking_length );
  }else{
    print_nt( $sense, $antisense, $strand, $query_start, $query_end, $upstream, $query_nt, $downstream, $file_ref, $flanking_length );
  }
}
####################################################################################################
# Called by extract_nt_sequence, formats and prints nucleotide sequence, prints nucleotide sequence
# in-frame, dispalying putative start and stop codons in red and green font, respectively
####################################################################################################
sub print_nt{
  my( $sense, $antisense, $strand, $query_start, $query_end, $upstream, $gene, $downstream, $file_ref, $flanking_length ) = @_;
  my( $i, $position_string, $HTML_sense, @start_stops );
  my $position = 12;
  my $genome   = $$file_ref;

  print HTML "<p2>\n";
  
  #build array of nucleotide coordinates (@caret_coordinates), each one to be printed below a caret (^)
  my @caret_coordinates;
  if( $strand eq "+" ){
    for( $i = 0; $i < (length $sense); $i += 12 ){
      $position = $query_start - 60 + $i;
      if( $position < 0 ){
	$position = (length $genome) + $position; #plus because $position is a negative number at this point
	push( @caret_coordinates, $position );
      }elsif( $position > length $genome ){
	$position = $query_end + $i - (length $genome);
	push( @caret_coordinates, $position );
      }else{
	push( @caret_coordinates, $position );
      }
    }
  }elsif( $strand eq "-" ){
    for( $i = 0; $i < (length $sense); $i += 12 ){
      $position = $query_end + 60 - $i;
      if( $position < 0 ){
	$position = (length $genome) + $position; #plus because $position is a negative number at this point
	push( @caret_coordinates, $position );
      }elsif( $position > length $genome ){
	$position = $position - (length $genome);
	push( @caret_coordinates, $position );
      }else{
	push( @caret_coordinates, $position );
      }
    }
  }
  
  #build string of nucleotide positions, which will be printed below the nucleotide sequences
  foreach my $p ( @caret_coordinates ){
    $position_string .= $p . '&nbsp;' x (12 - length $p);
  }

  # build sense strand, formatted for HTML output, with start codons=green, stop codons=red,
  # upstream, downstream, and nonsense strands all gray, and the nt sequence corresponding to the
  # AA sequence from BLAST output
 
  #print upstream sense strand
  for( $i = 0; $i < ((length $upstream) - 2 ); $i += 3 ){
    my $codon = substr( $upstream, $i, 3 );
    if( $codon =~ /ATG|ACG/ ){
      $HTML_sense .= "<font color=green font-weight=700><a title=\"POSITION!!!\">$codon<\/font><\/a>";
    }elsif( $codon =~ /TAA|TAG|TGA/ ){
      $HTML_sense .= "<font color=red font-weight=700><a title=\"POSITION!!!\">$codon<\/font><\/a>";
    }elsif( $codon =~ /CAA|CAG|CGA/ ){
      $HTML_sense .= "<font color=CCCC font-weight=700><a title=\"POSITION!!!\">$codon<\/font><\/a>";
    }else{
      $HTML_sense .= "<font color=gray font-weight=700>$codon<\/font>";
    }
  }
  #print gene
  for( $i = 0; $i < ((length $gene) - 2 ); $i+=3 ){
    my $codon = substr( $gene, $i, 3 );
    if( $codon =~ /ATG|ACG/ ){
      $HTML_sense .= "<font color=green font-weight=700><a title=\"POSITION!!!\">$codon<\/font><\/a>";
    }elsif( $codon =~ /TAA|TAG|TGA/ ){
      $HTML_sense .= "<font color=red font-weight=700><a title=\"POSITION!!!\">$codon<\/font><\/a>";
    }elsif( $codon =~ /CAA|CAG|CGA/ ){
      $HTML_sense .= "<font color=CCCC font-weight=700><a title=\"POSITION!!!\">$codon<\/font><\/a>";
    }else{
      $HTML_sense .= "<font color=663399 font-weight=700>$codon<\/font>";
    }
  }
  #print downstream sense strand
  for( $i = 0; $i < ((length $downstream) - 2 ); $i+=3 ){
    my $codon = substr( $downstream, $i, 3 );
    if( $codon =~ /ATG|ACG/ ){
      $HTML_sense .= "<font color=green font-weight=700><a title=\"POSITION!!!\">$codon<\/font><\/a>";
    }elsif( $codon =~ /TAA|TAG|TGA/ ){
      $HTML_sense .= "<font color=red font-weight=700><a title=\"POSITION!!!\">$codon<\/font><\/a>";
    }elsif( $codon =~ /CAA|CAG|CGA/ ){
      $HTML_sense .= "<font color=CCCC font-weight=700><a title=\"POSITION!!!\">$codon<\/font><\/a>";
    }else{
      $HTML_sense .= "<font color=gray font-weight=700>$codon<\/font>";
    }
  }

  # this code makes an array of the genome position of all start and stop codons in $sense
  if( $strand eq "+" ){
    for( $i = 0; $i <= (length $sense); $i += 3 ){
      my $codon    = substr( $sense, $i, 3 );
      my $position = $query_start - 60 + $i;
      if( $position < 0 ){
	$position = (length $genome) + $position; #plus because $position is a negative number at this point
	$codon =~ /ATG|ACG|TAA|TAG|TGA|CAA|CAG|CGA/ and push( @start_stops, "$position" );
      }elsif( $position > length $genome ){
	$position = $query_end + $i - (length $genome);
	$codon =~ /ATG|ACG|TAA|TAG|TGA|CAA|CAG|CGA/ and push( @start_stops, "$position" );
      }elsif( $codon =~ /ATG|ACG|TAA|TAG|TGA|CAA|CAG|CGA/ ){
	push( @start_stops, "$position" );
      }else{
	next;
      }
    }
  }elsif( $strand eq "-" ){
    for( $i = 0; $i <= (length $sense); $i += 3 ){
      my $codon    = substr( $sense, $i, 3 );
      my $position = $query_end + 60 - $i;
      if( $position < 0 ){
	$position = (length $genome) + $position; #plus because $position is a negative number at this point
	$codon =~ /ATG|ACG|TAA|TAG|TGA|CAA|CAG|CGA/ and push( @start_stops, $position );
      }elsif( $position > length $genome ){
	$position = $position - (length $genome);
	$codon =~ /ATG|ACG|TAA|TAG|TGA|CAA|CAG|CGA/ and push( @start_stops, $position );
      }elsif( $codon =~ /ATG|ACG|TAA|TAG|TGA|CAA|CAG|CGA/ ){
	push( @start_stops, $position );
      }else{
	next;
      }
    }
  }


  # Right now, $HTML_sense has a placeholder (POSITION!!!) for the mouse-over title showing the 
  # nucleotide positions of start and stop codons; this code replaces that placeholder with the
  # actual position, which is stored in the @start_stops array
  my $counter = 0;
  my $print_HTML_sense;
  for( split( /!!!/, $HTML_sense ) ){
    if( /POSITION/ ){
      $_ =~ s/POSITION/$start_stops[$counter]/g;
    }
    $print_HTML_sense .= $_;
    $counter++;
  }
  
  
  # print sense strand
  print HTML $print_HTML_sense;
  print HTML "<\/p2>\n";
  
  #print antisense strand
  print HTML "<p2><font color=gray font-weight=700>" . $antisense . "<\/font><\/p2>\n";

  #print "^" characters below nt sequences
  print HTML "<p2>";
  for($i = 0; $i < length $sense; $i+=12){
    print HTML "<font color=663399 font-weight=700>", "&#94;", '&nbsp;' x 11, "<\/font>";
  }
  print HTML"<\/p2>\n";
  
  #print coordinates below nt sequences
  print HTML "<p2><font color=663399 font-weight=700>", $position_string, "<\/font><\/p2>\n";
  print HTML "<br><br><br>\n";
  print HTML "<hr width=10000 size=4>";
}
####################################################################################################
#Called by Main, prints HTML header for BLAST summary page for each gene
####################################################################################################
sub print_output_header{
  my ( $project, $query_taxon, $gene, $gene_length, $out_directory, $gene_lengths_ref, $rna_bool ) = @_;
  
  open( HTML,">$out_directory/$gene.html" ) || die "Couldn't open $gene.html: $!\n";
  
  print HTML "<html>\n<head>\n<title>$query_taxon: $gene</title>\n</head>\n";
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

  print HTML "</style>\n";
  print HTML "<body onLoad=\"window.open(\'http:\/\/127.0.0.1\/cgi-bin\/make_annotation_form.cgi?gene_name=$gene&taxon=$query_taxon&project=$project\',\'mywindow\',\'width=550,height=680,resizable=yes,menubar=no,screenX=0,screenY=0\')\">\n";

  # @$gene_lengths_ref = sort {$a <=> $b}( @$gene_lengths_ref);
  # my $gl = join( ',&nbsp;', @$gene_lengths_ref);
  my $gl;
  my $min = min( @$gene_lengths_ref );
  my $max = max( @$gene_lengths_ref );
  if( $min == $max ){
    $gl = $min;
  }else{
    $gl = "$min\-$max";
  }

  print HTML "<table valign=\"middle\">\n\t<tr nowrap valign=\"middle\">\n";
  print HTML "\t\t<td>", '&nbsp;' x 85, "<\/td>\n";
  print HTML "\t\t<td><table><form METHOD=\"GET\"><INPUT TYPE=\"button\" VALUE=\"Annotate\" onClick=\"window.open(\'http:\/\/127.0.0.1\/cgi-bin\/make_annotation_form.cgi?gene_name=$gene&taxon=$query_taxon&project=$project\',\'Sequin_form\',\'width=550,height=680,resizable=yes,menubar=no\')\"><\/form><\/table><\/td>\n";
  if( $rna_bool ){
    print HTML "\t\t<td><font face=\"Courier\" color=\"blue\" size=\"6\">&nbsp;&nbsp;&nbsp;$gene&nbsp;($gl&nbsp;nt)<\/font><\/td>\n";
  }else{
    print HTML "\t\t<td><font face=\"Courier\" color=\"blue\" size=\"6\">&nbsp;&nbsp;&nbsp;$gene&nbsp;($gl&nbsp;AA)<\/font><\/td>\n";
  }
  print HTML "\t<\/tr>\n<\/table>\n<br>\n<hr size=4>";
  
}
####################################################################################################
sub print_output_footer{
  print HTML "</body>\n</html>";
  close( HTML );
}
####################################################################################################
# Called by Main, prints summary html page for the analyzed genome
####################################################################################################
  sub print_gene_summary{
    my( $project, $query_taxon, $seqlen, $html_files, $blast_files, $out_directory, $missed_ref ) = @_;
    my $gene;
    my $output_file = ">". "$out_directory" . "/" . "$project" . "_summary.html";
    my $n_missed_rows = @{$missed_ref};
    
    
    open( HTML_SUM,"$output_file" ) || die "Couldn't open $output_file: $!\n";
    
    #print header
    print HTML_SUM "<HTML>\n";
    print HTML_SUM "<TITLE>$query_taxon BLAST output</TITLE>\n";
    print HTML_SUM "<HEAD>\n";
    print HTML_SUM "<body BGCOLOR=#C2C2C2>\n\n";
    
    #header table
    print HTML_SUM "<table cellspacing=6 cellpadding=3 valign=\"middle\">\n";
    print HTML_SUM "\t<tr valign=\"middle\">\n";
    print HTML_SUM "\t\t<td>", '&nbsp;' x 20, "<\/td>\n";
    print HTML_SUM "\t\t<td><font face=\"Arial\" color=\"black\" size=\"5\">BLAST results for $query_taxon mtDNA<\/font><\/td>\n";
    print HTML_SUM "\t<\/tr>\n";
    print HTML_SUM "<\/table>\n\n<hr>\n\n";
    
    #results table
    print HTML_SUM "<table cellspacing=6 cellpadding=3 valign=\"middle\">\n";
    print HTML_SUM "\t<tr valign=\"middle\">\n";
    print HTML_SUM "\t\t<td>", '&nbsp;' x 20, "<\/td>\n";
    print HTML_SUM "\t\t<th align=\"left\"><font face=\"Arial\">Gene<\/font><\/th>\n";
    print HTML_SUM "\t\t<th align=\"center\"><font face=\"Arial\">Annotate<\/font><\/th>\n";
    print HTML_SUM "\t\t<th align=\"right\"><font face=\"Arial\">Raw BLAST output<\/font><\/th>\n";
    print HTML_SUM "\t\t<td>", '&nbsp;' x 20, "<\/td>\n";
    
    #print Missing genes portion of results table
    print HTML_SUM "\t\t<td valign=\"top\" rowspan=$n_missed_rows>\n";
    
    print HTML_SUM "\t\t\t<table valign=\"top\" valign=\"middle\">\n";
    print HTML_SUM "\t\t\t\t<tr><th align=\"left\"><font face=\"Arial\"><u>No hits<\/u><\/font><\/th><\/tr>\n";
    foreach( sort @{$missed_ref} ){
      print HTML_SUM "\t\t\t<tr><td align=\"left\"><font face=\"Arial\">$_<\/font><\/td><\/tr>\n";
    }
    
    print HTML_SUM "\t\t\t<\/table>\n";
    print HTML_SUM "\t\t<\/td>\t\n";
    print HTML_SUM "\t<\/tr>\n";
    
    my @mt_genes = qw( atp1 atp4 atp6 atp8 atp9 ccmB ccmC ccmFc ccmFn cob cox1 cox2 cox3 matR mttB nad1 nad2 nad3 nad4 nad4L nad5 nad6 nad7 nad9 rpl2 rpl5 rpl10 rpl16 rps1 rps2 rps3 rps4 rps7 rps10 rps11 rps12 rps13 rps14 rps19 sdh3 sdh4 );

    #print results, one row per gene found
    # foreach $gene( sort keys ( %{$html_files} ) ){
    foreach $gene( @mt_genes ){
      # my $blast_output = "$out_directory" . "/" . "$gene" . ".blastx";
      print HTML_SUM "\t<tr valign=\"middle\">\n";
      print HTML_SUM "\t\t<td>", '&nbsp;' x 20, "<\/td>\n";
      #print "gene" row header
      print HTML_SUM "\t\t<td align=\"left\"><font face=\"Arial\" color=\"purple\" size=\"4\">$gene&nbsp;<\/font><\/td>\n";
      #print "Annotate" button
      #print HTML_SUM "\t\t<td><table><tr align=\"left\"><form><INPUT TYPE=\"button\" VALUE=\"Annotate\" onClick=\"window.open(\'$out_directory/$gene.html\',\'BLAST_summary\')\"><\/form><\/table><\/td>\n";
      print HTML_SUM "\t\t<td><table><tr align=\"left\"><form><INPUT TYPE=\"button\" VALUE=\"Annotate\" onClick=\"window.open(\'$html_files->{$gene}\',\'BLAST_summary\')\"><\/form><\/table><\/td>\n";
      #print link to "Raw BLAST output"
      print HTML_SUM "\t\t<td align=\"right\"><a href=\"$blast_files->{$gene}\"><font face=\"Arial\">$gene BLAST output<\/font><\/a><\/td>\n";
      print HTML_SUM "\t<\/tr>\n";	
    }
    print HTML_SUM "<\/table>\n<br>\n<hr>";
    
  }
####################################################################################################
#Called by &annotate_rna

sub make_dna_no_hits{
  my( $gene, $proteinNoHits, $out_directory, $project ) = @_;
  my $no_hits = "No significant BLAST hits for $gene.";
  my $no_hits_page = "$out_directory/" . $gene . "_NoHits.html";
  $proteinNoHits->{$gene} = "\.\.\/$project" . "_out/" . $gene . "_NoHits.html";
  
  open( NO_HITS, ">$no_hits_page" ) || die "Couldn't open $no_hits_page: $!\n";
  print NO_HITS "<html>\n";
  print NO_HITS "<head></head>\n";
  print NO_HITS "\t<table valign=\"middle\">\n\t<tr nowrap valign=\"middle\">\n";
  print NO_HITS "\t\t<td>" . '&nbsp;' x 101;
  print NO_HITS "<\/td>\n";
  print NO_HITS "\t\t<td><table><\/table><\/td><br><br><br><br>\n";
  print NO_HITS "\t\t<td width=520 align=center><font face=\"Courier\" color=\"blue\" size=\"5\">$no_hits<\/font><\/td>\n";
  print NO_HITS "\t<\/tr>\n<\/table>\n";
  close NO_HITS;
  
}
####################################################################################################
sub complement{
  my $sequence = pop;
  $sequence = reverse $sequence;
  $sequence =~ tr/53ACGTUMRWSYKVHDBNacgtumrwsykvhdbn/35TGCAAKYWSRMBDHVNtgcaakywsrmbdhvn/;
  return $sequence;
}
####################################################################################################
sub max {
  my $max = shift(@_);
  foreach my $foo (@_) {
    $foo > $max and $max = $foo;
  }
  return $max;
}
####################################################################################################
sub min {
  my $min = shift(@_);
  foreach my $foo (@_) {
    $foo < $min and $min = $foo;
  }
  return $min;
}
####################################################################################################
# depracated as of 17 June 2012
# sub parse_matching_sequence{
#   my ( $query_match, $sbjct_match, $match, $query_taxon, $query_start, $query_end, $strand, $sbjct_taxon, $sbjct_start, $sbjct_end, $percent_match, $gene_length ) = @_;
#   my @lines = split( /[\n\n]/, $match );

#   foreach my $line ( @lines){
#     if( $line =~ /Query  (\d+)\s+(\S+)  (\d+)/ ){
#       $query_match .= $2;
#     }elsif( $line =~ /Sbjct  (\d+)\s+(\S+)  (\d+)/ ){
#       $sbjct_match .= $2;
#     }
#   }

#   print_sbjct_AA( $sbjct_taxon, uc $sbjct_match, $sbjct_start, $sbjct_end, $query_start, $query_end, $strand, $percent_match, $gene_length );
#   print_query_AA( $query_taxon, uc $query_match, $query_start, $query_end );
  
#   return $query_match;
# }
####################################################################################################

1;
