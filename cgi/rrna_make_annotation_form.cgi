#!/usr/bin/perl

use strict;
use warnings;
use CGI qw( :standard );

my %gene_names = (
		  rrn5  => "5S ribosomal RNA",
		  rrn16 => "16S ribosomal RNA",
		  rrnS  => "small subunit ribosomal RNA",
		  rrnL  => "large subunit ribosomal RNA",
		 );

my $gene = param('gene_name');
my $product = $gene_names{$gene};
my $taxon = param('taxon');
my $project = param( 'project' );

print "Content-type: text/html\r\n\r\n";
open( TEMPLATE,"rrna_template.html");
sequin_popup_content( $gene, $product, $project );
exit;
 
#################################SUBROUTINES#################################
sub sequin_popup_content{
  my( $gene, $product, $project ) = @_;
  while( <TEMPLATE> ){
    if( /GENE_NAME/ ){
      $_ =~ s/GENE_NAME/$gene/g;
    }elsif( /PRODUCT/ ){
      $_ =~ s/PRODUCT/"$product"/g;
    }elsif( /TAXON/ ){
      $_ =~ s/TAXON/$taxon/g;
    }elsif( /PROJECT/ ){
      $_ =~ s/PROJECT/$project/g;
    }
    print;
  }
  close( TEMPLATE );
}
#############################################################################
