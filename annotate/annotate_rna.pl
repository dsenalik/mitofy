use strict;
use warnings;

#############################################################################
# main script for annotating RNAs by BLAST; also calls parse_trnascan, which parses
# tRNA-Scan output

sub annotate_rna{

  my( $infile, $project, $query_taxon, $out_directory, $seqlen, $file, $RNA_EMAX, $RNA_PMIN, $RNA_MINLEN ) = @_;
  my $rna; #rna gene name (e.g., Val, rrnS )
  my %rna_html_files;   #key = three-letter AA abbreviation, value = filename of BLAST html results
  my %rna_blast_files;  #key = three-letter AA abbreviation, value = path of raw BLAST output
  my %tscan_html_files; #key = three-letter AA abbreviation, value = filename of tscan html results
  my %tscan_html_paths; #key = three-letter AA code; value = path of tscan HTML output file (e.g., "$project_out/Phe_trnascan.html")
  my %rnaNoHits;        #key = three-letter AA code; value = path of HTML "No significant BLAST hits" page; key and value exist only
                        #for those RNAs with no significant BLAST HITS, as determined by $counter boolean at end of annotate_rna
  my @missed_by_blast; #RNAs not found by BLAST
  my @missed_by_tscan; #tRNAs not found by tRNA-Scan
  my @missed;          #tRNAs missed by both BLAST and tRNA-Scan
  my @all_found;       #tRNAs found by both BLAST and tRNA-Scan
  my %trna_frame_files;
  my $rna_bool = 1;

  my @rna_order = qw(Arg Asn Asn-cp Asp Cys-bacterial Cys-cp Cys-mt Gln Glu Gly His-cp Ile Ile-cp Leu Leu-cp Lys Met-cp Met-f Phe Phe-cp Pro Pro-cp Ser Ser-cp Trp-cp Tyr Val-cp rrn5 rrnL rrnS );

  my %trna = (
	      Phe   => "F",
	      Leu   => "L",
	      Ile   => "I",
	      Met   => "M",
	      Val   => "V",
	      Ser   => "S",
	      Pro   => "P",
	      Thr   => "T",
	      Ala   => "A",
	      Tyr   => "Y",
	      His   => "H",
	      Gln   => "Q",
	      Asn   => "N",
	      Lys   => "K",
	      Asp   => "D",
	      Glu   => "E",
	      Cys   => "C",
	      Trp   => "W",
	      Arg   => "R",
	      Gly   => "G",
	      Undet => "UNDETERMINED tRNA"
	     );

  # open the rna blast db directory to get all rna gene names
  opendir( DBDIR, "blast_dbs/mt_rna" );
  my @RNA = readdir( DBDIR );
  closedir DBDIR;
  
  print "\nBLASTing RNAs . . . \n";

  foreach $rna ( @rna_order ){
    # $rna =~ /^\./ and next; #ignores the "." and ".." elements in the directory

    my $blast_data;   #the entire BLAST output
    my @gene_lengths; #stores lengths subject sequences in BLAST database
    my $sbjct_taxon;  #taxon name of subject sequence from BLAST hit
    my $gene_length;  #total length of sbjct sequence
    my( $gene_name, $taxa_count, $hit, %HITS );
    my $counter	= 0;  #hit counter, treats each non-redundant hit
                      #as an exon, for exporting hit data to the Sequin Information window
    
    #print rna being annotated to stdout
    print "\t$rna\n";

    #run BLAST for current gene
    # open( BLASTRUN, "blastall -p blastn -F F -d blast_dbs/mt_rna/$rna/$rna -i $infile  | " );
    open( BLASTRUN, "blast/blastn -query $infile -db blast_dbs/mt_rna/$rna/$rna  | " );

    while( <BLASTRUN> ){ #write BLAST output to $blast_data variable
      $blast_data .= $_ ;
    }

    #save raw BLAST output
    open( BLASTOUT, ">$out_directory/$rna.blastn" );
    print BLASTOUT $blast_data;
    close BLASTOUT;
      
    #store locations of html output file for each gene
    if( $rna =~ /fm/i ){
      $rna_html_files{Metf}  = "\.\./$project" . "_out/" . "$rna.html";
      $rna_blast_files{Metf} = "\.\./$project" . "_out/" . "$rna.blastn";
    }else{
      $rna_html_files{$rna}   = "\.\./$project" . "_out/" . "$rna.html";
      $rna_blast_files{$rna}  = "\.\./$project" . "_out/" . "$rna.blastn";
    }
    
    #start parsing BLAST output
    if( $blast_data =~ /No hits found/ ){ 
      print MISSED ("No hits found for $rna.\n");
      $rna_html_files{$rna} = not_found( "BLAST", $rna, $out_directory, $project, \%trna );
      push( @missed_by_blast, $rna ) and next;
    }
    
    #split BLAST output by taxon
    for( split( /> /, $blast_data ) ){
      /greedy/ and next; #skip BLAST header
      
      #parse taxon name and gene length
      # if( /^(\S+).*[\n,\t,\s]*Length = (\d+)/ ){
      if( /\[(.*)\]\s*Length=(\d+)/ ){
	$taxa_count++;
	$sbjct_taxon = $1;
	$gene_length = $2;
	$sbjct_taxon =~ s/,//;

	#store non-redundant gene lengths in @gene_lengths
	if( $taxa_count == 1 ){
	  push( @gene_lengths, $gene_length );
	}else{
	  unless( grep( /$gene_length/, @gene_lengths ) ){
	    push( @gene_lengths, $gene_length );
	  }
	}
      # }elsif( /^(\S+)\s*(\S+)[\n,\t,\s]*Length = (\d+)/ ){
      }elsif( /^(\S+).*\s+Length=(\d+)/ ){
	$taxa_count++;
	$sbjct_taxon = $1;
	$gene_length = $2;
	$sbjct_taxon =~ s/,//;

	#store non-redundant gene lengths in @gene_lengths
	if( $taxa_count == 1 ){
	  push( @gene_lengths, $gene_length );
	}else{
	  unless( grep( /$gene_length/, @gene_lengths ) ){
	    push( @gene_lengths, $gene_length );
	  }
	}
      }

      
      #parse gene name
      if( /ref\|.*\|\s(.*)\s\[/ ){
	$gene_name = $1;
      }
      parse_rna_hits_within_taxon( $query_taxon, $sbjct_taxon, $file, $rna, \%HITS, $gene_length );
    }

    my $rna_bool = 1;
    rna_output_header( $project, $query_taxon, $rna, $gene_length, $out_directory, \@gene_lengths, $rna_bool );
    
    #sort unique hits by $sbjct_start, then parse each hit
    foreach $hit ( sort {$HITS{$b}[10] <=> $HITS{$a}[10]} keys ( %HITS ) ){ #sorts by hit length
      my( $sbjct_taxon, $sbjct_start, $sbjct_end, $strand, $percent_match, $sbjct_match, $hit_e_value, $query_start, $query_end, $query_match, $match_length, $gene_length ) = @{$HITS{$hit}};

      #parse matching nt sequence for subject and query
      if( ($percent_match >= $RNA_PMIN) && ($match_length >= $RNA_MINLEN) && ($hit_e_value <= $RNA_EMAX) ){
	$counter++;

	print_sbjct_nt( $rna, $sbjct_taxon, uc $sbjct_match, $sbjct_start, $sbjct_end, $query_start, $query_end, $strand, $percent_match, $gene_length );
	print_query_nt( $rna, $query_taxon, uc $query_match, $query_start, $query_end );
	extract_gapped_nt( $rna, $query_start ,$query_end, $strand, $file, $query_match, $rna_bool );

      }
    }
    
    print_output_footer();
    $counter == 0 and make_rna_no_hits( $rna, \%rnaNoHits, $out_directory, $project, \%trna );
    $counter == 0 and $rna_html_files{$rna} = $rnaNoHits{$rna};
    
  }
  
  parse_trnascan( $infile, $file, $project, $query_taxon, $out_directory, \@missed_by_blast, \@missed_by_tscan, \@missed, \%rna_html_files, \@all_found, \%tscan_html_files, \%tscan_html_paths );

  make_framesets( $project, $out_directory, $query_taxon, \%rna_html_files, \%tscan_html_files, \@all_found, \%trna_frame_files, \%tscan_html_paths, \%rna_blast_files, \%rnaNoHits, \%trna );

  print_rna_summary( $query_taxon, $project, $seqlen, \%rna_html_files, \%rna_blast_files, \%tscan_html_files, $out_directory, \@missed_by_blast, \@missed, \@all_found, \%trna_frame_files, \%tscan_html_paths );

}
#############################################################################
#Called by annotate_rna, parses the multiple BLAST hits for a given taxon

sub parse_rna_hits_within_taxon{
  my( $query_taxon, $sbjct_taxon, $file, $rna, $HITS_ref, $gene_length ) = @_;
  
  #split by match within each gene	
  for( split( / Score/ ) ){
    my( $strand, $query_strand, $sbjct_strand );
    my $hit_e_value;   #e-value of hit
    my $match_length;  #length of match between subject & query
    my $percent_match; #the %identity between subject & query
    my $gaps = 0;     #boolean, if there are gaps in the match, $gaps = 1; necessary for subsequent
                      #printing of the upstream, gene, and downstream sequences corresponding to the subject match
    my $num;    #first coord of query sequence, start or end coord depending on strand
    my $match;  #the entire match, including AA sequences of sbjct and query
    my $query_start;  #start coordinate of hit in query sequence
    my $query_end;    #end coordinate of hit in query sequence
    my $query_coords; #coordinates of hit in "xx..xx" format
    my $query_match;
    my $sbjct_start;  #start coordinate of hit in subject sequence
    my $sbjct_end;    #end coordinate of hit in subject sequence
    my $sbjct_coords; #coordinates of hit in "xx..xx" format
    my $sbjct_match;
    my $upstream;     #nt sequence upstream of the hit, for identifying alternative start codons, etc.
    my $downstreamn;  #nt sequence upstream of the hit, for identifying alternative start codons, etc.
    my $gene_nt_sequence; #nt sequence of gene
    my $start;
    my @hit_data;

    /Length=/ and next; #skip header

    #parse e-value
    if( /Expect[\(\d+\)]* =\s*(\de-\d+)/ ){
      $hit_e_value = $1;
    }elsif( /Expect[\(\d+\)]* =\s*(e-\d+)/ ){
      $hit_e_value =  "1.$1";
    }elsif( /Expect =\s*(\d+\.\d*)/ ){
      $hit_e_value =  $1;
    }elsif( /Expect =\s*(\d+)/ ){
      $hit_e_value =  $1;
    }
    
    # print $hit_e_value, "\n";

    #parse length of match from number of identities
    if( /Identities = (\d+)\/(\d+) \((\d+)\%\)/ ){
      $match_length = $2;
      $percent_match = $3;
    }
    
    # print $match_length, "\n";
    # print $percent_match, "\n";

    # parse forward or reverse strand
    if( /Strand=(\w+)\/(\w+)/ ){
      $query_strand = $1;
      $sbjct_strand = $2;
      
      if ($query_strand eq "Plus" && $sbjct_strand eq "Minus"){
	$strand = "-";
      }
      if ($query_strand eq "Minus" && $sbjct_strand eq "Plus"){
	$strand = "-";
      }
      if ($query_strand eq "Plus" && $sbjct_strand eq "Plus"){
	$strand = "+";
      }
    }
    
    # print $strand, "\n";

    #parse gap information and adjust match length
    if( /Gaps = (\d+)\/(\d+)/ ){
      $gaps = 1;
      $match_length -= $1;
    }

    if( /Query  (\d+)/ ){ #parse match
      $num = $1;
      #record text until end of part, this captures nt sequence of
      #query and subject sequence, which will be separately parsed
      #with &parse_matching_sequence
      $match = "Query  ".$num.$';
      if( $match =~ (/\n\n\nLambda/) ){ # if it's last taxon,
	$match = $`; # get rid of stuff at very end.
      }
    }

    ( $query_start, $query_end, $query_coords, $query_match, $sbjct_start, $sbjct_end, $sbjct_coords, $sbjct_match ) = match_coords( $sbjct_taxon, $query_taxon, $match, $strand );
    
    #store hit data to split and print later on, MIGHT NOT NEED THIS NOW
    @hit_data = ( $sbjct_taxon, $sbjct_start, $sbjct_end, $strand, $percent_match, $sbjct_match, $hit_e_value, $query_start, $query_end, $query_match, $match_length, $gene_length );
    
    #store hit data in %HITS, indexed by $sbjct_start; this effectively
    #filters redundant hits *ACROSS TAXA* and allows hits to later be 
    #sorted by $sbjct_start, so hits can by displayed by the order of
    #the exons. This is all done in Main.
    
    #not sure why I did it this way, but match_coords returns $query_start, $query_end,
    #$sbjct_start, and $sbjct_end as if it all hits are on "+" strand
    
    my $hit;
    #want to sort rRNAs by $sbjct_start
    if( $rna =~ /rrn/ ){
      $hit = "$query_start.$query_end";
      #want to sort tRNAs by $query_start
    }else{
      if( $strand eq "+" ){
	$hit = "$query_start.$query_end";
      }else{
	$hit = "$query_end.$query_start";
      }
    }
    
    #if $HITS{$hit} is already defined, replace it with current redundant hit only if it
    #has a lower e-value
    if( $HITS_ref->{$hit} ){
      if( $hit_e_value < $HITS_ref->{$hit}->[6] ){
	push @{$HITS_ref->{$hit}}, @hit_data;
      }
    }else{
      push @{$HITS_ref->{$hit}}, @hit_data;
    }
    # $HITS_ref->{ $hit } = $hit_data;
    # print $HITS_ref->{ $sbjct_coords } . "\n";
  }
}
#############################################################################
# Called by annotate_rna (above), formats and prints subject nt to HTML output

sub print_sbjct_nt{
  my ( $rna, $sbjct_taxon, $sbjct_sequence, $sbjct_start, $sbjct_end, $query_start, $query_end, $strand, $percent_match, $gene_length ) = @_;
  my( $query_coords, $name_color, $sbjct_color, $query_color, $x, $y, $z, $a );
  
  if( $rna =~ /rrn/ ){
    $name_color  = "336600";
    $query_color = "336600";
    $sbjct_color = "99CC99";
  }else{
    $name_color  = "333399";
    $query_color = "333399";
    $sbjct_color = "9999FF";
  }

  if( $sbjct_taxon =~ (/trn/) ){
    $sbjct_taxon = $`;
  }
  if( $sbjct_taxon =~ (/tRNA/) ){
    $sbjct_taxon = $`;
  }
  
  #Layout of nt printout:
  #|---------------------$a---------------|nt sequence#
  #|--------$x-------|$taxon_name|---$z---|nt sequence#
  
  my $prefix = 60; #number of characters preceding nt sequence
  $y = length $sbjct_taxon; #number of characters in taxon name
  $z = 6; #number of characters between taxon name and nt sequence
  $x = $prefix -$y - $z; #number of characters preceding taxon name
  $a = $x + $y + $z;	
  
  $sbjct_taxon =~ s/ /&nbsp;/g;
  $sbjct_taxon =~ s/_/&nbsp;/g;
  $sbjct_taxon =~ s/-/&#8211;/g;
  
  if( $strand eq "+" ){
    $query_coords = "$query_start..$query_end";
  }else{
    $query_coords = "$query_end..$query_start";	
  }
  
  print HTML "<p1>", '&nbsp;' x $a, "<font color=$sbjct_color>$sbjct_start..$sbjct_end&nbsp;<\/font>//&nbsp;<font color=$name_color>$query_coords<\/font>", "<font color=black> (match length = ", $query_end-$query_start+1, " nt, $percent_match\% identity within match, $strand strand)", "<\/font></p1><br>\n";
  print HTML "<p1>", '&nbsp;' x $a, "&#124;", "<\/p1><br>\n";
  
  
  print HTML "<p1>", '&nbsp;' x $x, "<font COLOR=$name_color FACE=COURIER font-weight=1000>$sbjct_taxon<\/font>";
  print HTML '&nbsp;' x $z, "<font COLOR=$sbjct_color FACE=COURIER font-weight=1000>", uc $sbjct_sequence, "<\/p1><br>\n";
}
#############################################################################
# Called by annotate_rna (above), formats and prints query nt to HTML output

sub print_query_nt{
  my ( $rna, $query_taxon, $query_sequence, $query_start, $query_end ) = @_;
  my( $query_coords, $name_color, $query_color, $x, $y, $z, $a );
  
  if( $rna =~ /rrn/ ){
    $name_color  = "336600";
    $query_color = "336600";
  }else{
    $name_color  = "333399";
    $query_color = "333399";
  }
  
  #Layout of query nt printout:
  #|---------------------$a---------------|RNA sequence#
  #|--------$x-------|$taxon_name|---$z---|RNA sequence#
  
  my $prefix = 60; #number of characters preceding AA sequence
  $y = length $query_taxon ; #number of characters in taxon name
  $z = 6; #number of characters between taxon name and AA sequence
  $x = $prefix - $y - $z; #number of characters preceding taxon name
  
  $query_taxon =~ s/ /&nbsp;/g;
  
  print HTML "<p2>", '&nbsp;' x $x, "<font COLOR=$name_color FACE=COURIER font-weight=1000>$query_taxon<\/font>";
  print HTML '&nbsp;' x $z, "<font COLOR=$query_color FACE=COURIER font-weight=1000>", uc $query_sequence, "<\/p2><br><\/font>\n";
}
#############################################################################
# Called by annotate_rna, extract nucleotide sequence from query genome for each BLAST hit

sub extract_gapped_nt{
  my ( $rna, $query_start ,$query_end, $strand, $file_ref, $query_match, $rna_bool ) = @_;
  my( $up_offset, $down_offset, $upstream, $downstream, $sense, $antisense );

  my $genome_length = length $$file_ref;
  my $gene_length   = $query_end - $query_start + 1;
  my $gene_offset   = $query_start - 1;
  my $gene = $query_match;
  $gene =~ s/-/&#8209/g;
  
  my $flanking_length = 60; #in addition to gene, extract 60 nt each of upstream and downstream sequence
    
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
    
    $sense = "$upstream" . "$gene" . "$downstream";

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
    $upstream   = complement( $upstream );
    $downstream = complement( $downstream );

    $sense = "$upstream" . "$gene" . "$downstream";

  }
  
  $antisense = $sense;
  $antisense =~ tr/53ACGTUMRWSYKVHDBNacgtumrwsykvhdbn/35TGCAAKYWSRMBDHVNtgcaakywsrmbdhvn/;

  rna_print_nt( uc $sense, uc $antisense, $strand, $query_start, $query_end, uc $upstream, uc $gene, uc $downstream, $file_ref, $flanking_length, $rna, $query_match );

}
#############################################################################
# Called by extract_gapped_nt (for BLAST output) and print_trnascan (for TSCAN output); formats and prints
# nucleotide sequence of concatenated upstream, gene, and downstream sequence below the hit sequences

sub rna_print_nt {
  my( $sense, $antisense, $strand, $query_start, $query_end, $upstream, $gene, $downstream, $file_ref, $flanking_length, $rna, $query_match ) = @_;
  my( $i, $HTML_sense, @POSITIONS, $name_color, $query_color );
  my $genome = $$file_ref;
  
  if( $rna =~ /rrn/ ){
    $name_color  = "336600";
    $query_color = "336600";
  }else{
    $name_color  = "333399";
    $query_color = "333399";
  }		
  
  print HTML "<p2>\n";
  
  # build sense strand, formatted for HTML output, with upstream, downstream, and nonsense 
  # strands all gray, and the BLAST hit sequence purple
  #upstream sense strand
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
  
  my( $position, $position_string );
  
  
  ################################IF $query_match HAS NO GAPS################################
  if( $query_match !~ /-+/ ){
    print HTML "<p2>";
    
    #print "^" characters below nt sequences
    for( $i = 0; $i < length $sense; $i += 12 ){
      if( $rna =~ /rrn/ ){
	print HTML "<font color=336600 font-weight=700>", "&#94;", '&nbsp;' x 11, "<\/font>";
      }else{
	print HTML "<font color=333399 font-weight=700>", "&#94;", '&nbsp;' x 11, "<\/font>";
      }	
    }
    print HTML"<\/p2>\n";
    
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
    
    #print coordinates below nt sequences
    if( $rna =~ /rrn/ ){
      print HTML "<p2><font color=336600 font-weight=700>", $position_string, "<\/font><\/p2>\n";
      print HTML "<br><br><br>\n";
      print HTML "<hr width=10000 size=4>\n\n\n";
    }else{
      print HTML "<p2><font color=333399 font-weight=700>", $position_string, "<\/font><\/p2>\n";
      print HTML "<br><br><br>\n";
      print HTML "<hr width=10000 size=4>\n\n\n";
    }
    
  ################################IF $query_match HAS GAPS################################
  
  }elsif( $query_match =~ /-/ ){
    my $gapped_sense = $sense;
    $gapped_sense =~ s/&#8209/-/g;
    
    #create string ($caret_string), which is the string of "^" characters, printed every 12 nt
    #below the nt sequences
    my @gapped_sense_array = split( //, $gapped_sense );
    my( $caret_string, $caret_counter );
    
    for($i = 0; $i < @gapped_sense_array; $i++){
      if( $i == 0 ){
	$caret_string .= "<font color=663399 font-weight=700>&#94;<\/font>";
      }elsif( $gapped_sense_array[$i] =~ /-/ ){
	$caret_string .= "&nbsp;";
      }elsif( $gapped_sense_array[$i] !~ /-/ ){
	$caret_counter++;
	if( ($caret_counter % 12) == 0 ){
	  $caret_string .= "<font color=663399 font-weight=700>&#94;<\/font>";
	}else{
	  $caret_string .= "&nbsp;";
	}
      }
    }
    
    print HTML "<p2>";
    if( $rna =~ /rrn/ ){
      print HTML "<font color=336600 font-weight=700>", $caret_string, "<\/font>";
    }else{
      print HTML "<font color=333399 font-weight=700>", $caret_string, "<\/font>";
    }	
    print HTML"<\/p2>\n";
    
    
    ####the following block of code creates the coordinate string for gapped sequences####
    ####//////////////////////////////////////////////////////////////////////////////####
    my $coordinate_string = $caret_string;
    $coordinate_string =~ s/<font color=663399 font-weight=700>&#94;<\/font>/^/g;

    # "!" is placeholder for "&nbsp;" in the coordiante string below
    $coordinate_string =~ s/&nbsp;/!/g;

    #@spaces contains strings of "!", with length of the string specifiying the number of spaces
    #between carets (^); need to know this, so I can subtract the number of digits in the coordinate
    my @spaces = split( /\^/, $coordinate_string );
    
    #@spacing contains the number of spaces between each caret
    my $spcLen = 0;
    my @spacing;
    for( $i = 1; $i < @spaces; $i++ ){
      $spcLen = (length $spaces[$i]) + 1;
      push @spacing, $spcLen;
    }
    
    my $total_gap_length = 0;
    for( my $i = 0; $i < (length $gapped_sense); $i++ ){
      my $nt = substr( $gapped_sense, $i, 1 );
      $nt =~ /-/ and $total_gap_length++;
    }
      
    #build array of nucleotide coordinates (@caret_coordinates), each one to be printed below a caret (^)
    my @caret_coordinates;
    if( $strand eq "+" ){
      for( $i = 0; $i < (length($gapped_sense) - $total_gap_length); $i += 12 ){
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
      for( $i = 0; $i < (length($gapped_sense) - $total_gap_length); $i += 12 ){
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
    
    #calculate length (i.e., number of characters) of each $coordinate in @caret_coordinates; necessary
    #to determine how many spaces go between each coordinate in HTML-formatted output
    my @coord_size;
    foreach( @caret_coordinates ){
      push( @coord_size, length $_ );
    }
    
    #calculate the number of spaces between coordinates in $coord_string, which depends on the
    #number of digits in the coordinate and the number of spaces between adjacent carets (by
    #default this is 12 spaces, more if there are gaps in the sequence)
    my $nSpaces;
    my @num_spaces;
    for( $i = 0; $i < @coord_size; $i++){
      
      #this is necessary for cases when there are no spaces after the last caret/coordinate
      #pair; in these cases, $spacing[$i] will be undefined
      if( $spacing[$i] && $coord_size[$i]){
	$nSpaces = $spacing[$i] - $coord_size[$i];
      }
      push @num_spaces, $nSpaces;
    }
    
    #build the actual coordinate string, by concatenating the coordinates and spaces
    for( $i = 0; $i < @caret_coordinates; $i++ ){
      $position_string .= $caret_coordinates[$i];
      $position_string .= "&nbsp;" x $num_spaces[$i];
    }
    
    #print coordinates below nt sequences
    if( $rna =~ /rrn/ ){
      print HTML "<p2><font color=336600 font-weight=700>", $position_string, "<\/font><\/p2>\n";
      print HTML "<br><br><br>\n";
      print HTML "<hr width=10000 size=4>";
    }else{
      print HTML "<p2><font color=333399 font-weight=700>", $position_string, "<\/font><\/p2>\n";
      print HTML "<br><br><br>\n";
      print HTML "<hr width=10000 size=4>\n\n\n";
    }	

  }
}
#############################################################################
#Called by parse_rna_hits_within_taxon
sub match_coords{

  my ( $sbjct_taxon, $query_taxon, $match, $strand ) = @_;
  my( $query_match, $sbjct_match, $query_start, $query_end, $query_coords );
  my( $sbjct_start, $sbjct_end, $sbjct_coords );
  my( @q_starts, @q_ends, @s_starts, @s_ends );

  my @lines = split( /[\n\n]/, $match );

  foreach( @lines ){
    if( /Query  (\d+)\s+(\S+)  (\d+)/ ){
      push( @q_starts, $1 );
      push( @q_ends, $3 );
      $query_match .= $2;
    }elsif( /Sbjct  (\d+)\s+(\S+)  (\d+)/ ){
      push( @s_starts, $1 );
      push( @s_ends, $3 );
      $sbjct_match .= $2;
    }
  }

  $query_start  = shift @q_starts;
  $query_end    = pop @q_ends;
  $query_coords = "$query_start..$query_end";
  
  if( $strand eq "+" ){
    $sbjct_start = shift @s_starts;
    $sbjct_end	 = pop @s_ends;
  }else{
    $sbjct_end	 = shift @s_starts;
    $sbjct_start = pop @s_ends;
  }
  $sbjct_coords = "$sbjct_start..$sbjct_end";
  
  if( $strand eq "-" ){
    $query_match = complement( $query_match );
    $sbjct_match = complement( $sbjct_match );	
  }
  
  return( $query_start, $query_end, $query_coords, $query_match, $sbjct_start, $sbjct_end, $sbjct_coords, $sbjct_match );
}
#############################################################################

1;
