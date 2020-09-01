#!/usr/bin/perl

use strict;
use warnings;
use CGI qw( :standard );

# my $title	  = param('title');
my $title	  = "AndyTitle";
my $blast_out = param('blast');
my $tscan_out = param('tscan');

print "Content-type: text/html\r\n\r\n";
open( TEMPLATE,"frameset.html");
frame_content( $title, $blast_out, $tscan_out );

exit;
 
#################################SUBROUTINES#################################
sub frame_content{
  my( $title, $blast_out, $tscan_out ) = @_;
  while( <TEMPLATE> ){
    if( /TITLE/ ){
      $_ =~ s/TITLE/$title/g;
    }if( /BLAST/ ){
      $_ =~ s/BLAST/$blast_out/g;
    }if( /TRNA_SCAN/ ){
      $_ =~ s/TRNA_SCAN/$tscan_out/g;
    }else{
    }
    print;
  }
  close TEMPLATE;
}
#############################################################################
