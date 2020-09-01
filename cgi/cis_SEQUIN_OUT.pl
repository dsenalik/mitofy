#!/usr/bin/perl

use CGI qw( :standard );
use strict;
use warnings;

######################################MAIN###################################
my( $cis_start, $cis_end );
my( $ex1_start, $ex2_start, $ex3_start, $ex4_start, $ex5_start );
my( $ex1_end, $ex2_end, $ex3_end, $ex4_end, $ex5_end );
my @exons;

# get passed CGI parameters
my $gene        = param('gene');
my $taxon       = param('taxon');
my $project     = param('project');
my $feature_key = param('feature_key');
my $product     = param('product');
my $notes       = param('notes');
my $pseudo      = param('pseudo');

my $sequin_outfile = $project . ".tbl";
my $outfile_path   = "cgi_out/"; #DIRECTORY

# for transl_except
my $te_1s  = param('te_1s');
my $te_1e  = param('te_1e');
my $te_1aa = param('te_1aa');
my $te_2s  = param('te_2s');
my $te_2e  = param('te_2e');
my $te_2aa = param('te_2aa');
my $te1 = ""; #transl_except for start codon
my $te2 = ""; #transl_except for stop codon

my $cis_coords = param( 'gene_coords' );
if( $cis_coords ){
  ( $cis_start, $cis_end ) = split( /\D+/, $cis_coords );
}

# for transl_except
if( $te_1s ){
  $te1 = "(pos:" . $te_1s . ".." . $te_1e . ",aa:" . $te_1aa . ")";
}
if( $te_2s ){
  $te2 = "(pos:" . $te_2s . ".." . $te_2e . ",aa:" . $te_2aa . ")";
}

my $exon1 = param( 'ex1' );
if( $exon1 ){
  push( @exons, $exon1 );
  if( $exon1 =~ /(\d+)\D+(\d+)/ ){
    $ex1_start = $1;
    $ex1_end = $2;
  }
}
my $exon2 = param( 'ex2' );
if( $exon2 ){
  push( @exons, $exon2 );
  if( $exon2 =~ /(\d+)\D+(\d+)/ ){
    $ex2_start = $1;
    $ex2_end	= $2;
  }
}
my $exon3 = param( 'ex3' );
if( $exon3 ){
  push( @exons, $exon3 );
  if( $exon3 =~ /(\d+)\D+(\d+)/ ){
    $ex3_start = $1;
    $ex3_end	= $2;
  }
}
my $exon4 = param( 'ex4' );
if( $exon4 ){
  push( @exons, $exon4 );
  if( $exon4 =~ /(\d+)\D+(\d+)/ ){
    $ex4_start = $1;
    $ex4_end	= $2;
  }
}
my $exon5 = param( 'ex5' );
if( $exon5 ){
  push( @exons, $exon5 );
  if( $exon5 =~ /(\d+)\D+(\d+)/ ){
    $ex5_start = $1;
    $ex5_end	= $2;
  }
}

if( -e "$outfile_path$sequin_outfile" ){
  print_SEQUIN( $cis_start, $cis_end, \@exons, $exon1, $exon2, $exon3, $exon4, $exon5, $ex1_start, $ex2_start, $ex3_start, $ex4_start, $ex5_start, $ex1_end, $ex2_end, $ex3_end, $ex4_end, $ex5_end, $gene, $taxon, $sequin_outfile, $outfile_path, $cis_coords, $feature_key, $product, $notes, $te1, $te2, $pseudo );
}else{
  open( SEQUIN, ">$outfile_path$sequin_outfile" ) || print "Can't open $!";
  print SEQUIN ">Feature $taxon\n";
  close SEQUIN;
  print_SEQUIN( $cis_start, $cis_end, \@exons, $exon1, $exon2, $exon3, $exon4, $exon5, $ex1_start, $ex2_start, $ex3_start, $ex4_start, $ex5_start, $ex1_end, $ex2_end, $ex3_end, $ex4_end, $ex5_end, $gene, $taxon, $sequin_outfile, $outfile_path, $cis_coords, $feature_key, $product, $notes, $te1, $te2, $pseudo );
}

#-NEW LINES----------
print "Content-type: text/html\n\n";
print "<HTML><BODY onload=\"window.close()\"></BODY</HTML>\n";
#--------------------


exit;

#################################SUBROUTINES#################################
sub print_SEQUIN{
  my( $cis_start, $cis_end, $exons_ref, $exon1, $exon2, $exon3, $exon4, $exon5, $ex1_start, $ex2_start, $ex3_start, $ex4_start, $ex5_start, $ex1_end, $ex2_end, $ex3_end, $ex4_end, $ex5_end, $gene, $taxon, $sequin_outfile, $outfile_path, $cis_coords, $feature_key, $product, $notes, $te1, $te2, $pseudo ) = @_;
	
  open( SEQUIN, ">>$outfile_path$sequin_outfile" ) || die "Couldn't open $sequin_outfile: $!\n";
  
  # print misc_features
  if( $feature_key eq "misc_feature" ){
    print SEQUIN "$cis_start\t$cis_end\t$feature_key\n";
    if( $gene ){
      print SEQUIN "\t\t\tgene\t$gene\n";
    }
    if( $pseudo eq "yes" ){
      print SEQUIN "\t\t\tpseudo\n";
    }
    if( $notes ){
      print SEQUIN "\t\t\tnote\t$notes\n";
    }
  # print introns
  }elsif( $feature_key eq "intron" ){
    print SEQUIN "$cis_start\t$cis_end\t$feature_key\n";
    print SEQUIN "\t\t\tgene\t$gene\n";
    if( $pseudo eq "yes" ){
      print SEQUIN "\t\t\tpseudo\n";
    }
    if( $notes ){
      print SEQUIN "\t\t\tnote\t$notes\n";
    }
  # print genes
  }else{
    print SEQUIN "$cis_start\t$cis_end\t$feature_key\n";
    print SEQUIN "\t\t\tgene\t$gene\n";
    if( $pseudo eq "yes" ){
      print SEQUIN "\t\t\tpseudo\n";
    }
    if( $notes ){
      print SEQUIN "\t\t\tnote\t$notes\n";
    }
  }
  
  if( @{$exons_ref} == 5 ){
    print SEQUIN "$ex1_start\t$ex1_end\tCDS\n";
    print SEQUIN "$ex2_start\t$ex2_end\n";
    print SEQUIN "$ex3_start\t$ex3_end\n";
    print SEQUIN "$ex4_start\t$ex4_end\n";
    print SEQUIN "$ex5_start\t$ex5_end\n";
    print SEQUIN "\t\t\tproduct\t$product\n";
    if( $te1 ){
      print SEQUIN "\t\t\ttransl_except\t$te1\n";
    }
    if( $te2 ){
      print SEQUIN "\t\t\ttransl_except\t$te2\n";
    }
    print SEQUIN "$ex1_start\t$ex1_end\texon\n";
    print SEQUIN "\t\t\tnumber 1\n";
    print SEQUIN "$ex2_start\t$ex2_end\texon\n";
    print SEQUIN "\t\t\tnumber 2\n";
    print SEQUIN "$ex3_start\t$ex3_end\texon\n";
    print SEQUIN "\t\t\tnumber 3\n";
    print SEQUIN "$ex4_start\t$ex4_end\texon\n";
    print SEQUIN "\t\t\tnumber 4\n";
    print SEQUIN "$ex5_start\t$ex5_end\texon\n";
    print SEQUIN "\t\t\tnumber 5\n";
  }elsif( @{$exons_ref} == 4 ){
    print SEQUIN "$ex1_start\t$ex1_end\tCDS\n";
    print SEQUIN "$ex2_start\t$ex2_end\n";
    print SEQUIN "$ex3_start\t$ex3_end\n";
    print SEQUIN "$ex4_start\t$ex4_end\n";
    print SEQUIN "\t\t\tproduct\t$product\n";
    if( $te1 ){
      print SEQUIN "\t\t\ttransl_except\t$te1\n";
    }
    if( $te2 ){
      print SEQUIN "\t\t\ttransl_except\t$te2\n";
    }
    print SEQUIN "$ex1_start\t$ex1_end\texon\n";
    print SEQUIN "\t\t\tnumber 1\n";
    print SEQUIN "$ex2_start\t$ex2_end\texon\n";
    print SEQUIN "\t\t\tnumber 2\n";
    print SEQUIN "$ex3_start\t$ex3_end\texon\n";
    print SEQUIN "\t\t\tnumber 3\n";
    print SEQUIN "$ex4_start\t$ex4_end\texon\n";
    print SEQUIN "\t\t\tnumber 4\n";
  }elsif( @{$exons_ref} == 3 ){
    print SEQUIN "$ex1_start\t$ex1_end\tCDS\n";
    print SEQUIN "$ex2_start\t$ex2_end\n";
    print SEQUIN "$ex3_start\t$ex3_end\n";
    print SEQUIN "\t\t\tproduct\t$product\n";
    if( $te1 ){
      print SEQUIN "\t\t\ttransl_except\t$te1\n";
    }
    if( $te2 ){
      print SEQUIN "\t\t\ttransl_except\t$te2\n";
    }
    print SEQUIN "$ex1_start\t$ex1_end\texon\n";
    print SEQUIN "\t\t\tnumber 1\n";
    print SEQUIN "$ex2_start\t$ex2_end\texon\n";
    print SEQUIN "\t\t\tnumber 2\n";
    print SEQUIN "$ex3_start\t$ex3_end\texon\n";
    print SEQUIN "\t\t\tnumber 3\n";
  }elsif( @{$exons_ref} == 2 ){
    print SEQUIN "$ex1_start\t$ex1_end\tCDS\n";
    print SEQUIN "$ex2_start\t$ex2_end\n";
    print SEQUIN "\t\t\tproduct\t$product\n";
    if( $te1 ){
      print SEQUIN "\t\t\ttransl_except\t$te1\n";
    }
    if( $te2 ){
      print SEQUIN "\t\t\ttransl_except\t$te2\n";
    }
    print SEQUIN "$ex1_start\t$ex1_end\texon\n";
    print SEQUIN "\t\t\tnumber 1\n";
    print SEQUIN "$ex2_start\t$ex2_end\texon\n";
    print SEQUIN "\t\t\tnumber 2\n";
  }else{
    if( ($pseudo eq "no") && ($feature_key ne "misc_feature") && ($feature_key ne "intron") ){
      print SEQUIN "$cis_start\t$cis_end\tCDS\n";
      print SEQUIN "\t\t\tproduct\t$product\n";
      if( $te1 ){
	print SEQUIN "\t\t\ttransl_except\t$te1\n";
      }
      if( $te2 ){
	print SEQUIN "\t\t\ttransl_except\t$te2\n";
      }
    }
  }
  close SEQUIN;
}
#############################################################################
# this subroutine determines whether a sequence is on the + or - strand
sub strand{
  my( $start, $end ) = @_;
  my( $strand, $tmp );
  
  if( $start < $end ){
    $strand = "+";
  }else{
    $tmp = $end;
    $end = $start;
    $start = $tmp;
    $strand = "-";
  } 
  return( $start, $end, $strand );
}
#############################################################################

