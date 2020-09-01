use strict;
use warnings;

#############################################################################
sub nad5_ex3{
  my( $query_taxon, $blast_data, $out_directory ) = @_;
  my $match_length;  #length of match between subject & query
  my $percent_match; #the %identity between subject & query
  my $strand;        #is hit on + or - strand
  my $num;           #first coord of query sequence, start or end coord depending on strand
  my $match;         #the entire match, including AA sequences of sbjct and query
  my $hit_data;      #a concatenated string of hit data, each datum separated by "!!", so the string can easily be split later on
  my $query_start;   #start coordinate of hit in query sequence
  my $query_end;     #end coordinate of hit in query sequence
  my $query_strand;
  my $query_coords;  #coordinates of hit in "xx..xx" format
  my $query_match;
  my $sbjct_taxon;   #start coordinate of hit in subject sequence
  my $sbjct_start;   #start coordinate of hit in subject sequence
  my $sbjct_end;     #end coordinate of hit in subject sequence
  my $sbjct_strand;
  my $sbjct_match;
  my $sbjct_coords;  #coordinates of hit in "xx..xx" format
  my $gene_length;
  my %nadHits;

  # parse taxon name and gene length
  # if( $blast_data =~ /^> (\S+)\s*(\S+)[\n,\t,\s]*Length = (\d+)/ ){
  #   $sbjct_taxon = $1;
  #   $gene_length = $3;
  #   $sbjct_taxon =~ s/,//;
  # }

  $sbjct_taxon = "Beta";
  $gene_length = 22;
  
  for( split( / Score/, $blast_data ) ){
    /Length=/ and next; #skip database sequence header
    
    #parse just the identical, full-length hits
    if( /Identities = (22)\/(22) \((100)\%\)/ ){
      $match_length = $2;
      $percent_match = $3;
      
      # parse out forward or reverse strand
      #if (/Strand = (\w+) \/ (\w+)/){
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
      
      #parse coords and strand of query sequence from hit
      if( /Query  (\d+)/ ){ #parse start and end coords of query sequence
	$num = $1;

	# record text until end of part, this captures query and subject sequence, which will be separately parsed
	# with parse_matching_sequence
	$match = "Query  ".$num.$';
	if( $strand eq "-" ){
	  $query_end = $num;
	  $query_start = $query_end - (3 * $match_length) + 1;
	  $query_coords = "$query_start..$query_end";
	  #print "query coords = " . $query_coords . "\n";
	}else{
	  $query_start = $num;
	  $query_end = $query_start + (3 * $match_length) -1;
	  $query_coords = "$query_start..$query_end";
	}
      if( $match =~ (/\n\n\nLambda/) ){ # if it's last taxon,
	$match = $`; # get rid of stuff at very end.
      }
      }
      
      #parse coordinates of subject sequence from hit
      if( /Sbjct  (\d+)/ ) { #parse sbjct start coord
	$sbjct_start = $1;
      }
      if( /(\d*)\s*\n\n\n/ ){ #parse sbjct end coord
	$sbjct_end = $1;
      }
      $sbjct_coords = "$sbjct_start..$sbjct_end";
      
      ( $query_start, $query_end, $query_coords, $query_match, $sbjct_start, $sbjct_end, $sbjct_coords, $sbjct_match ) = match_coords( $sbjct_taxon, $query_taxon, $match, $strand );
      
      #store hit data to split and print later on, MIGHT NOT NEED THIS NOW
      $hit_data = "Beta!!$sbjct_start!!$sbjct_end!!$sbjct_match!!$strand!!$match!!$query_taxon!!$query_start!!$query_end!!$query_match!!$match_length";
      #print $hit_data . "\n";
      
      #store hit data in %nadHits, indexed by "$sbjct_start.$sbjct_end"
      my $hit = "$sbjct_start.$sbjct_end";
      $nadHits{$hit} = $hit_data;
      
    }else{
      next;
    }
  }
  
  
  #sort unique hits by $sbjct_start, then parse each hit
  my @keys = keys( %nadHits );
  foreach ( sort {$a <=> $b} keys (%nadHits) ){
    my( $sbjct_taxon, $sbjct_start, $sbjct_end, $sbjct_match, $strand, $match, $query_taxon, $query_start, $query_end, $query_match, $match_length ) = split( /!!/, $nadHits{$_} );
    #print $nadHits{$_} . "\n";
    # extract corresponding nt sequence from query genome, to display
    nad5_ex3_print_sbjct_nt( "nad5_exon3", $sbjct_taxon, $sbjct_match, $sbjct_start, $sbjct_end, $query_start, $query_end, $strand, $out_directory );
    nad5_ex3_print_query_nt( "nad5_exon3", $query_taxon, $query_match, $query_start, $query_end, $out_directory );
  }
  
}
#############################################################################
# called by nad5_ex3 (above), formats and prints subject nt sequence to HTML output

sub nad5_ex3_print_sbjct_nt{
  my ( $gene, $sbjct_taxon, $sbjct_sequence, $sbjct_start, $sbjct_end, $query_start, $query_end, $strand, $out_directory ) = @_;
  my $query_coords;
  my $x = 0;
  my $y = 0;
  my $z = 0;
  my $a = 0;
  
  # print $out_directory . "\n";
  # print $sbjct_sequence . "\n";
  
  my $name_color  = "333399";
  my $query_color = "333399";
  my $sbjct_color = "CC66CC";
  
  #Layout of nt printout:
  #|---------------------$a---------------|nt sequence#
  #|--------$x-------|$taxon_name|---$z---|nt sequence#
  
  my $prefix = 60; #number of characters preceding nt sequence
  $y = length( $sbjct_taxon ); #number of characters in taxon name
  $z = 6; #number of characters between taxon name and nt sequence
  $x = ($prefix -$y) - $z; #number of characters preceding taxon name
  $a = $x + $y + $z;	
  my $b = length "nad5 exon3";
  my $c = ($prefix -$b) - $z; #number of characters preceding taxon name
  
  $sbjct_taxon =~ s/ /&nbsp;/g;
  $sbjct_taxon =~ s/_/&nbsp;/g;
  $sbjct_taxon =~ s/-/&#8211;/g;
  
  if( $strand eq "+" ){
    $query_coords = "$query_start..$query_end";
  }else{
    $query_coords = "$query_end..$query_start";	
  }
  
  open( HTML,">>$out_directory/nad5.html" ) || die "Couldn't open nad5.html: $!\n";
  
  print HTML "<p1>", "<font COLOR=$name_color FACE=COURIER font-weight=1000 size=\"4\">&nbsp;&nbsp;&nbsp;&nbsp;nad5&nbsp;exon3:<\/font>", "</p1><br>\n";
  print HTML "<p1>", '&nbsp;' x $a, "<font color=$name_color>$query_coords<\/font>", " (", $sbjct_end-$sbjct_start+1, " nt, $strand strand)", "</p1><br>\n";
  print HTML "<p1>", '&nbsp;' x $a, "&#124;", "<\/p1><br>\n";
  
  
  print HTML "<p1>", '&nbsp;' x $x, "<font COLOR=$name_color FACE=COURIER font-weight=1000>$sbjct_taxon<\/font>";
  print HTML '&nbsp;' x $z, "<font COLOR=$sbjct_color FACE=COURIER font-weight=1000>", uc $sbjct_sequence, "<\/p1><br>\n";

  close HTML;
}
#############################################################################
sub nad5_ex3_print_query_nt{
  my( $gene, $query_taxon, $query_match, $query_start, $query_end, $out_directory ) = @_;
  my $name_color = "333399";
  my $query_color = "663399";
  
  #Layout of query nt printout:
  #|---------------------$a---------------|RNA sequence#
  #|--------$x-------|$taxon_name|---$z---|RNA sequence#
  
  my $prefix = 60; #number of characters preceding AA sequence
  my $y = length( $query_taxon ); #number of characters in taxon name
  my $z = 6; #number of characters between taxon name and AA sequence
  my $x = ($prefix -$y) - $z; #number of characters preceding taxon name
  
  $query_taxon =~ s/ /&nbsp;/g;
  
  open( HTML,">>$out_directory/nad5.html" ) || die "Couldn't open nad5.html: $!\n";

  print HTML "<p2>", '&nbsp;' x $x, "<font COLOR=$name_color FACE=COURIER font-weight=1000>$query_taxon<\/font>";
  print HTML '&nbsp;' x $z, "<font COLOR=$query_color FACE=COURIER font-weight=1000>", uc $query_match, "<\/p2><br><\/font>\n";
  print HTML "<br><br><br>\n";
  print HTML "<hr width=10000 size=4>\n\n\n";
  print_output_footer();
  close HTML;
  
}
#############################################################################

1;
