#!/usr/bin/perl

use strict;
use warnings;
use CGI qw( :standard );

our $gene    = param('gene_name');
our $taxon   = param('taxon');
our $project = param('project');

my %gene_names = (
        atp1  => "ATPase subunit 1",
        atp4  => "ATPase subunit 4",
        atp6  => "ATPase subunit 6",
        atp8  => "ATPase subunit 8",
        atp9  => "ATPase subunit 9",
        ccmB  => "cytochrome c biogenesis B",
        ccmC  => "cytochrome c biogenesis C",
        ccmFc => "cytochrome c biogenesis FC",
        ccmFn => "cytochrome c biogenesis FN",
        cob   => "apocytochrome b",
        cox1  => "cytochrome c oxidase subunit 1",
        cox2  => "cytochrome c oxidase subunit 2",
        cox3  => "cytochrome c oxidase subunit 3",
        matR  => "maturase",
        matr  => "maturase",
        mttB  => "transport membrane protein",
        nad1  => "NADH dehydrogenase subunit 1",
        nad2  => "NADH dehydrogenase subunit 2",
        nad3  => "NADH dehydrogenase subunit 3",
        nad4  => "NADH dehydrogenase subunit 4",
        nad4L => "NADH dehydrogenase subunit 4L",
        nad5  => "NADH dehydrogenase subunit 5",
        nad6  => "NADH dehydrogenase subunit 6",
        nad7  => "NADH dehydrogenase subunit 7",
        nad9  => "NADH dehydrogenase subunit 9",
        rpl2  => "ribosomal protein L2",
        rpl5  => "ribosomal protein L5",
        rpl10 => "ribosomal protein L10",
        rpl16 => "ribosomal protein L16",
        rps1  => "ribosomal protein S1",
        rps2  => "ribosomal protein S2",
        rps3  => "ribosomal protein S3",
        rps4  => "ribosomal protein S4",
        rps7  => "ribosomal protein S7",
        rps10 => "ribosomal protein S10",
        rps12 => "ribosomal protein S12",
        rps13 => "ribosomal protein S13",
        rps14 => "ribosomal protein S14",
        rps19 => "ribosomal protein S19",
	sdh3  => "succinate dehydrogenase subunit 3",
	sdh4  => "succinate dehydrogenase subunit 4",
);

print "Content-type: text/html\r\n\r\n";
if(( $gene eq "nad1" ) || ($gene eq "nad2" ) || ($gene eq "nad5" )){
  open( TEMPLATE,"trans_template.html");
  sequin_popup_content( \%gene_names );
  close TEMPLATE;
}else{
  open( TEMPLATE,"cis_template.html");
  sequin_popup_content( \%gene_names );
  close TEMPLATE;
}
exit;
 
#################################SUBROUTINES#################################
sub sequin_popup_content{
  my $gene_names = shift;
  while( <TEMPLATE> ){
    if( /GENE_NAME/ ){
      $_ =~ s/GENE_NAME/$gene/g;
    }elsif( /PRODUCT/ ){
      $_ =~ s/PRODUCT/"$gene_names->{$gene}"/g;
    }elsif( /TAXON/ ){
      $_ =~ s/TAXON/$taxon/g;
    }elsif( /PROJECT/ ){
      $_ =~ s/PROJECT/$project/g;
    }else{
    }
    print;
  }
}
#############################################################################
