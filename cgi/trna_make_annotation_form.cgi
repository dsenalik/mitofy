#!/usr/bin/perl

use strict;
use warnings;
use CGI qw( :standard );

our $gene    = param('gene_name');
our $coords  = param('gene_coords');
our $product = param('product');
our $note    = param('note');
our $taxon   = param('taxon');
our $project = param('project');
our $cplast  = param('cplast');

$gene    or $gene    = "";
$coords  or $coords  = "";
$product or $product = "";
$note    or $note    = "";

print "Content-type: text/html\r\n\r\n";
open( TEMPLATE,"trna_template.html");

while( <TEMPLATE> ){
  if( /GENE_NAME/ ){
    $_ =~ s/GENE_NAME/$gene/g;
  }elsif( /PRODUCT/ ){
    $_ =~ s/PRODUCT/$product/g;
  }elsif( /NOTE/ ){
    $_ =~ s/NOTE/$note/g;
  }elsif( /GENE_COORDS/ ){
    $_ =~ s/GENE_COORDS/$coords/g;
  }elsif( /TAXON/ ){
    $_ =~ s/TAXON/$taxon/g;
  }elsif( /PROJECT/ ){
    $_ =~ s/PROJECT/$project/g;
  }elsif( /NAME=cplast VALUE=yes/ ){
    if( $cplast ){
      if( $cplast eq "yes" ){
	$_ =~ s/VALUE=yes/VALUE=yes CHECKED/g;
	print;
	while( <TEMPLATE> ){
	  $_ =~ s/NAME=cplast VALUE=no CHECKED/NAME=cplast VALUE=no/;
	  last;
	}
      }
    }
  }
  print;
}

close TEMPLATE;

exit;
