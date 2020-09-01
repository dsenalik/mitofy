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
our $type    = param('type');

$gene    or $gene    = "";
$coords  or $coords  = "";
$product or $product = "";
$note    or $note    = "";

print "Content-type: text/html\r\n\r\n";
open( TEMPLATE,"trnC_template.html");

while( <TEMPLATE> ){
  if( /GENE_NAME/ ){
    $_ =~ s/GENE_NAME/trnC/;
  }elsif( /PRODUCT/ ){
    $_ =~ s/PRODUCT/$product/;
  }elsif( /NOTE/ ){
    $_ =~ s/NOTE/$note/;
  }elsif( /GENE_COORDS/ ){
    $_ =~ s/GENE_COORDS/$coords/;
  }elsif( /TAXON/ ){
    $_ =~ s/TAXON/$taxon/;
  }elsif( /PROJECT/ ){
    $_ =~ s/PROJECT/$project/;
  }elsif( $type ){
    if ( $type eq "mito" && /NAME=type VALUE=mito/ ){
      $_ =~ s/VALUE=mito/VALUE=mito CHECKED/;
    }elsif( $type eq "cplast" && /NAME=type VALUE=cplast/ ){
      $_ =~ s/VALUE=cplast/VALUE=cplast CHECKED/;
    }elsif( $type eq "bacterial" && /NAME=type VALUE=bacterial/ ){
      $_ =~ s/VALUE=bacterial/VALUE=bacterial CHECKED/;
    }
  }
  print;
}

close TEMPLATE;

exit;
