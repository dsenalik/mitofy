use strict;
use warnings;

####################################################################################################
# Called by extract_nt_sequence, formats and prints nucleotide sequence, prints nucleotide sequence
# in-frame, dispalying putative start and stop codons in red and green font, respectively
####################################################################################################

sub print_nt_with_gaps {
  my( $query_match, $sense, $antisense, $strand, $query_start, $query_end, $upstream, $query_nt, $downstream, $genome, $flanking_length ) = @_;
  my $position = 12; # nt coordinates will be printed every 12 nt below the nt sequence
  my $counter  = 0;
  my @nt_gaps;	   # holds nt gap strings, in order, first to last
  my @aa_segments; # holds AA strings split out from gapped AA sequence; AA strings are stored in order, first to last
  my $nt_length;   # length of each AA segment, multiplied by 3, to get the corresponding nt length
  my @nt_lengths;  # holds nt lengths for each subsequence ("subsequence" refers to the nt sequences in between gaps)
  my $cumulative_nt_offset = 0; # first *nt* position (i.e., offset) of each subsequence; "subsequence" refers to nt sequences in between gaps
  my @nt_offsets = ( 0 ); # holds offset for each subsequence, in order, first to last; start at zero because nt sequence 
                          # will be extracted from $query_nt, which does not include the upstream and downstream sequences; therefore, 
                          # must build the gapped sequence, then add the upstream and downstream sequences
  my $gapped_nt;   # the gapped gene sequence, no upstream or downstream sequence
  my( $caret_counter, $anti_compl, @antisense_gaps, @antisense_seqs, @antisense_compl );
  my( $i, $position_string, $HTML_sense, @start_stops );
  
  print HTML "<p2>\n";

  # find gap characters in AA sequence, create nt gap strings by making the AA gaps 3x longer
  while( $query_match =~ /(-+)/g ){
    my $nt_gap = '-' x (3 * (length $1));
    push( @nt_gaps, $nt_gap );
  }
  
  # pull out the AA sequence, which is split by gaps
  @aa_segments = split( /-+/, $query_match);
  
  # record length of each AA segment, multiply by three to get length of corresponding nt sequence
  for( $i = 0; $i < @aa_segments; $i++ ){
    my $nt_length = 3*(length $aa_segments[$i] );
    push( @nt_lengths, $nt_length );
  }
  
  # calculate the first *nt* position (i.e., the offset) for each subsequence, store in @nt_offsets
  for( $i = 0; $i < (@nt_lengths - 1); $i++ ){
    $cumulative_nt_offset += $nt_lengths[$i];
    push( @nt_offsets, $cumulative_nt_offset );
  }
  
  #build nt sequence, with gaps
  for( $i = 0; $i < @nt_offsets; $i++){
    $gapped_nt .= substr( $query_nt, $nt_offsets[$i], $nt_lengths[$i] );
    $nt_gaps[$i] and $gapped_nt .= $nt_gaps[$i];
  }
  
  my $gap_count = 0;
  my $gapped_sense = $upstream . $gapped_nt . $downstream;
  
  #create string ($caret_string), which is the string of "^" characters, printed every 12 nt
  #below the nt sequences
  my @gapped_sense_array = split(//, $gapped_sense);
  my $caret_string;
  for($i = 0; $i < @gapped_sense_array; $i++){
    if( $i == 0 ){
      $caret_string .= "<font color=663399 font-weight=700>&#94;";
    }elsif( $gapped_sense_array[$i] =~ /-/ ){
      $caret_string .= "&nbsp;";
    }elsif( $gapped_sense_array[$i] !~ /-/ ){
      $caret_counter++;
      if( ($caret_counter % 12) == 0 ){
	$caret_string .= "<font color=663399 font-weight=700>&#94;";
      }else{
	$caret_string .= "&nbsp;";
      }
    }
  }

####the following block of code creates the coordinate string for gapped sequences####
####//////////////////////////////////////////////////////////////////////////////####

  my $coordinate_string = $caret_string;
  $coordinate_string    =~ s/<font color=663399 font-weight=700>&#94;/^/g;

  # "!" is placeholder for "&nbsp;" in the coordiante string below
  $coordinate_string =~ s/&nbsp;/!/g;
  
  #@spaces contains strings of "!", with length of the string specifiying the number of spaces
  #between carets (^); need to know this, so I can subtract the number of digits in the coordinate
  my @spaces = split( /\^/, $coordinate_string );
  
  #@spacing contains the number of spaces between each caret
  my $spcLen;
  my @spacing;
  for( $i = 1; $i < @spaces; $i++ ){
    $spcLen = (length $spaces[$i]) + 1;
    push( @spacing, $spcLen );
  }
  
  #build array of nucleotide coordinates (@caret_coordinates), each one to be printed below a caret (^)
  my @caret_coordinates;
  if( $strand eq "+" ){
    for( $i = 0; $i < (length $sense); $i += 12 ){
      $position = $query_start - 60 + $i;
      if( $position < 0 ){
	$position = (length $$genome) + $position; #plus because $position is a negative number at this point
	push( @caret_coordinates, $position );
      }elsif( $position > length $$genome ){
	$position = $query_end + $i - (length $$genome);
	push( @caret_coordinates, $position );
      }else{
	push( @caret_coordinates, $position );
      }
    }
  }elsif( $strand eq "-" ){
    for( $i = 0; $i < (length $sense); $i += 12 ){
      $position = $query_end + 60 - $i;
      if( $position < 0 ){
	$position = (length $$genome) + $position; #plus because $position is a negative number at this point
	push( @caret_coordinates, $position );
      }elsif( $position > length $$genome ){
	$position = $position - (length $$genome);
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
    # print "spacing: " . $spacing[$i] . "\n";
    # print "coord_size: " . $coord_size[$i] . "\n";
    $nSpaces = $spacing[$i] - $coord_size[$i];
    push( @num_spaces, $nSpaces );
  }
  
  #build the actual coordinate string, by concatenating the coordinates and spaces
  for( $i = 0; $i < @caret_coordinates; $i++ ){
    $position_string .= $caret_coordinates[$i];
    $position_string .= "&nbsp;" x $num_spaces[$i];
  }
####//////////////////////////////////////////////////////////////////////////////####
####/////////////////////////////////////END//////////////////////////////////////####
	
  # build sense strand, formatted for HTML output, with start codons=green, stop codons=red,
  # upstream, downstream, and nonsense strands all gray, and the nt sequence corresponding to the
  # AA sequence from BLAST output
  if( $strand eq "+" || $strand eq "-" ){
    #print upstream sense strand

    for($i = 0; $i < ((length $upstream) - 2 ); $i += 3){
      my $codon = substr($upstream, $i, 3);
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
    for($i = 0; $i < ((length $gapped_nt) - 2 ); $i += 3){
      my $codon = substr($gapped_nt, $i, 3);
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
    for($i = 0; $i < ((length $downstream) - 2 ); $i += 3){
      my $codon = substr($downstream, $i, 3);
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
  }
  
######the following block of code is for formatting the gapped antisense strand######
######/////////////////////////////////////////////////////////////////////////######
  #find gap characters in AA sequence, create nt gap strings by making the AA gaps 3x longer
  while ($gapped_nt =~ /(-+)/g){
    push( @antisense_gaps, $1 );
  }
  
  #pull out the AA sequence, which is split by gaps
  @antisense_seqs = split( /-+/, $gapped_nt);
  
  foreach( @antisense_seqs ){
    $_ =~ tr/53ACGTUMRWSYKVHDBNacgtumrwsykvhdbn/35TGCAAKYWSRMBDHVNtgcaakywsrmbdhvn/;
    push( @antisense_compl, $_);
  }
  
  for( $i = 0; $i < @antisense_compl; $i++ ){
    $anti_compl .= $antisense_compl[$i];
    if($antisense_gaps[$i]){
      $anti_compl .= $antisense_gaps[$i];
    }
  }
  my $gapped_antisense;
  my $up_trans = $upstream;
  my $down_trans = $downstream;
  
  if( $strand eq "+" ){
    $up_trans   =~ tr/53ACGTUMRWSYKVHDBNacgtumrwsykvhdbn/35TGCAAKYWSRMBDHVNtgcaakywsrmbdhvn/;
    $down_trans =~ tr/53ACGTUMRWSYKVHDBNacgtumrwsykvhdbn/35TGCAAKYWSRMBDHVNtgcaakywsrmbdhvn/;
    $gapped_nt  =~ tr/53ACGTUMRWSYKVHDBNacgtumrwsykvhdbn/35TGCAAKYWSRMBDHVNtgcaakywsrmbdhvn/;
    $gapped_antisense = $up_trans . $gapped_nt . $down_trans;
  }else{
    $up_trans   =~ tr/53ACGTUMRWSYKVHDBNacgtumrwsykvhdbn/35TGCAAKYWSRMBDHVNtgcaakywsrmbdhvn/;
    $down_trans =~ tr/53ACGTUMRWSYKVHDBNacgtumrwsykvhdbn/35TGCAAKYWSRMBDHVNtgcaakywsrmbdhvn/;
    $gapped_nt  =~ tr/53ACGTUMRWSYKVHDBNacgtumrwsykvhdbn/35TGCAAKYWSRMBDHVNtgcaakywsrmbdhvn/;
    $gapped_antisense = $down_trans . $gapped_nt . $up_trans;
  }
######/////////////////////////////////////////////////////////////////////////######
######//////////////////////////////////END////////////////////////////////////######


######the following block of code creates mouse-over coordinates for the sense strand######
######///////////////////////////////////////////////////////////////////////////////######

  # print $sense, "\n";
  # this code makes an array of the genome position of all start and stop codons in $sense
  if( $strand eq "+" ){
    for( $i = 0; $i < ((length $sense) - 2); $i += 3 ){
      my $codon    = substr( $sense, $i, 3 );
      my $position = $query_start - $flanking_length + $i;
      #print "$codon $position\n";
      if( $position < 0 ){
	$position = (length $$genome) + $position; #plus because $position is a negative number at this point
	$codon =~ /ATG|ACG|TAA|TAG|TGA|CAA|CAG|CGA/ and push( @start_stops, $position );
      }elsif( $position > length $$genome ){
	$position = $query_end + $i - (length $$genome);
	$codon =~ /ATG|ACG|TAA|TAG|TGA|CAA|CAG|CGA/ and push( @start_stops, $position );
      }elsif( $codon =~ /ATG|ACG|TAA|TAG|TGA|CAA|CAG|CGA/ ){
	#print "   pushing $codon at position $position\n";
	push( @start_stops, $position );
      }else{
	next;
      }
    }
  }elsif( $strand eq "-" ){
    for( $i = 0; $i < ((length $sense) - 2); $i += 3 ){
      my $codon    = substr( $sense, $i, 3 );
      my $position = $query_end + $flanking_length - $i;
      if( $position < 0 ){
	$position = (length $$genome) + $position; #plus because $position is a negative number at this point
	$codon =~ /ATG|ACG|TAA|TAG|TGA|CAA|CAG|CGA/ and push( @start_stops, $position );
      }elsif( $position > length $$genome ){
	$position = $position - (length $$genome);
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
  my $print_HTML_sense;

  # print "$query_start, $query_end\n";

  for( split( /!!!/, $HTML_sense ) ){
    if( /POSITION/ ){
      #print $counter, " ", $start_stops[$counter], "\n";
      $_ =~ s/POSITION/$start_stops[$counter]/g;
    }
    $print_HTML_sense .= $_;
    $counter++;
  }
######///////////////////////////////////////////////////////////////////////////////######
######///////////////////////////////////////////////////////////////////////////////######

############################write all of this out in HTML format############################
######////////////////////////////////////////////////////////////////////////////////######
  # print sense strand
  print HTML $print_HTML_sense;
  print HTML "<\/p2>\n";
  
  #print antisense strand
  print HTML "<p2><font color=gray font-weight=700>" . $gapped_antisense . "<\/font><\/p2>\n";
  
  #print "^" characters below nt sequences
  print HTML "<p2>";
  print HTML $caret_string;
  print HTML"<\/p2>\n";
  
  #print coordinates below nt sequences
  print HTML "<p2><font color=663399 font-weight=700>", $position_string, "<\/font><\/p2>\n";
  print HTML "<br><br><br>\n";
  print HTML "<hr width=10000 size=4>";
}
######////////////////////////////////////////////////////////////////////////////////######
############################################################################################
1;
