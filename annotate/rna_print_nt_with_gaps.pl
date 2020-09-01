use strict;
use warnings;

sub gapped_rna_coordinate_string{
####the following block of code creates the coordinate string for gapped sequences####
####//////////////////////////////////////////////////////////////////////////////####
  my $coordinate_string = $caret_string;
  $coordinate_string =~ s/<font color=663399 font-weight=700>&#94;<\/font>/^/g;

  # "!" is placeholder for "&nbsp;" in the coordiante string below
  $coordinate_string =~ s/&nbsp;/!/g;
  
  # @spaces contains strings of "!", with length of the string specifiying the number of spaces
  # between carets (^); need to know this, so I can subtract the number of digits in the coordinate
  my @spaces = split( /\^/, $coordinate_string );
  
  # @spacing contains the number of spaces between each caret
  my $spcLen = 0;
  my @spacing;

  for( $i = 1; $i < @spaces; $i++ ){
    $spcLen = length($spaces[$i]) + 1;
    push( @spacing, $spcLen);
  }
  
  #build array of nucleotide coordinates (@positions), each one to be printed below a caret (^)
  my @positions;
  for( $i = 0; $i < length $sense; $i += 12 ){
    if( $strand eq "+" ){
      $position = ($query_start-60) + $i;
      if( $position <= 0 ){
	$position = ((length $genome) + $position);
	push( @positions, $position );
      }else{
	push( @positions, $position );
      }
    }elsif( $strand eq "-" ){
      $position = ($query_end+60) - $i;               
      push( @positions, $position );
    }
  }
	
  #calculate length (i.e., number of characters) of each $coordinate in @positions; necessary
  #to determine how many spaces go between each coordinate in HTML-formatted output
  my @coord_size;
  foreach( @positions ){
    push( @coord_size, length( $_) );
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
  for( $i = 0; $i < @positions; $i++ ){
    $position_string .= $positions[$i];
    $position_string .= "&nbsp;" x $num_spaces[$i];
  }
####//////////////////////////////////////////////////////////////////////////////####
####/////////////////////////////////////END//////////////////////////////////////####

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

