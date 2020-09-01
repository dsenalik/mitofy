#!/usr/bin/perl

use CGI qw( :standard );
use strict;
use warnings;

######################################MAIN###################################
my $trans_start1;
my $trans_end1;
my $trans_start2;
my $trans_end2;
my $trans_start3;
my $trans_end3;
my @exons;
my $ex1_start;
my $ex2_start;
my $ex3_start;;
my $ex4_start;
my $ex5_start;
my $ex1_end;
my $ex2_end;
my $ex3_end;
my $ex4_end;
my $ex5_end;

my $gene        = param('gene');
my $taxon       = param('taxon');
my $project     = param('project' );
my $feature_key = param('feature_key' );
my $product     = param('product' );
my $notes       = param('notes' );
my $pseudo      = param('pseudo' );
my $sequin_outfile = $project . ".tbl";
my $outfile_path   = "cgi_out/"; #DIRECTORY

#for transl_except
my $te_1s  = param('te_1s');
my $te_1e  = param('te_1e');
my $te_1aa = param('te_1aa');
my $te_2s  = param('te_2s');
my $te_2e  = param('te_2e');
my $te_2aa = param('te_2aa');
my $te1 = ""; #transl_except for start codon
my $te2 = ""; #transl_except for stop codon

my $trans_coords1 = param( 'trans_coords1' );
if( $trans_coords1 ){	
  if ( $trans_coords1 =~ /(\d+)\D+(\d+)/ ){
    $trans_start1 = $1;
    $trans_end1 = $2;
  }
}
my $trans_coords2 = param( 'trans_coords2' );
if( $trans_coords2 ){
  if ( $trans_coords2 =~ /(\d+)\D+(\d+)/ ){
    $trans_start2 = $1;
    $trans_end2 = $2;
  }
}
my $trans_coords3 = param( 'trans_coords3' );
if( $trans_coords3 ){
  if ( $trans_coords3 =~ /(\d+)\D+(\d+)/ ){
    $trans_start3 = $1;
    $trans_end3 = $2;
  }
}

#for transl_except
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
    $ex1_end	= $2;
  }
}
my $exon2 = param( 'ex2' );
if( $exon2 ){
  push( @exons, $exon1 );
  if( $exon2 =~ /(\d+)\D+(\d+)/ ){
    $ex2_start = $1;
    $ex2_end	= $2;
  }
}
my $exon3 = param( 'ex3' );
if( $exon3 ){
  push( @exons, $exon1 );
  if( $exon3 =~ /(\d+)\D+(\d+)/ ){
    $ex3_start = $1;
    $ex3_end	= $2;
  }
}
my $exon4 = param( 'ex4' );
if( $exon4 ){
  push( @exons, $exon1 );
  if( $exon4 =~ /(\d+)\D+(\d+)/ ){
    $ex4_start = $1;
    $ex4_end	= $2;
  }
}
my $exon5 = param( 'ex5' );
if( $exon5 ){
  push( @exons, $exon1 );
  if( $exon5 =~ /(\d+)\D+(\d+)/ ){
    $ex5_start = $1;
    $ex5_end	= $2;
  }
}


if( -e "$outfile_path$sequin_outfile" ){
  print_SEQUIN( $trans_start1, $trans_end1, $trans_start2, $trans_end2, $trans_start3, $trans_end3, \@exons, $exon1, $exon2, $exon3, $exon4, $exon5, $ex1_start, $ex2_start, $ex3_start, $ex4_start, $ex5_start, $ex1_end, $ex2_end, $ex3_end, $ex4_end, $ex5_end, $gene, $taxon, $sequin_outfile, $outfile_path, $trans_coords1, $trans_coords2, $trans_coords3, $feature_key, $product, $notes, $te1, $te2, $pseudo );
}else{
  open( SEQUIN, ">$outfile_path$sequin_outfile" );
  print SEQUIN ">Feature $taxon\n";
  close SEQUIN;
  print_SEQUIN( $trans_start1, $trans_end1, $trans_start2, $trans_end2, $trans_start3, $trans_end3, \@exons, $exon1, $exon2, $exon3, $exon4, $exon5, $ex1_start, $ex2_start, $ex3_start, $ex4_start, $ex5_start, $ex1_end, $ex2_end, $ex3_end, $ex4_end, $ex5_end, $gene, $taxon, $sequin_outfile, $outfile_path, $trans_coords1, $trans_coords2, $trans_coords3, $feature_key, $product, $notes, $te1, $te2, $pseudo );
}

#-NEW LINES----------
print "Content-type: text/html\n\n";
print "<HTML><BODY onload=\"window.close()\"></BODY</HTML>\n";
#--------------------

exit;


#################################SUBROUTINES#################################
sub print_SEQUIN{

  my( $trans_start1, $trans_end1, $trans_start2, $trans_end2, $trans_start3, $trans_end3, $exons_ref, $exon1, $exon2, $exon3, $exon4, $exon5, $ex1_start, $ex2_start, $ex3_start, $ex4_start, $ex5_start, $ex1_end, $ex2_end, $ex3_end, $ex4_end, $ex5_end, $gene, $taxon, $sequin_outfile, $outfile_path, $trans_coords1, $trans_coords2, $trans_coords3, $feature_key, $product, $notes, $te1, $te2, $pseudo ) = @_;
	
  open( SEQUIN, ">>$outfile_path$sequin_outfile" );
  
  # print intron
  if( $feature_key eq "intron" ){
    print SEQUIN "$trans_start1\t$trans_end1\t$feature_key\n";
    print SEQUIN "\t\t\tgene\t$gene\n";
    if( $pseudo eq "yes" ){
      print SEQUIN "\t\t\tpseudo\n";
    }
    if( $notes ){
      print SEQUIN "\t\t\tnote\t$notes\n";
    }
  # print gene
  }else{
    print SEQUIN "$trans_start1\t$trans_end1\t$feature_key\n";
    print SEQUIN "$trans_start2\t$trans_end2\n";
    if( $trans_coords3 ){
      print SEQUIN "$trans_start3\t$trans_end3\n";
    }
    print SEQUIN "\t\t\tgene\t$gene\n";
    if( $pseudo eq "yes" ){
      print SEQUIN "\t\t\tpseudo\n";
    }
    print SEQUIN "\t\t\texception\ttrans-splicing\n";
    if( $notes ){
      print SEQUIN "\t\t\tnote\t$notes\n";
    }
    
    if( @{$exons_ref} == 5 ){
      print SEQUIN "$ex1_start\t$ex1_end\tCDS\n";
      print SEQUIN "$ex2_start\t$ex2_end\n";
      print SEQUIN "$ex3_start\t$ex3_end\n";
      print SEQUIN "$ex4_start\t$ex4_end\n";
      print SEQUIN "$ex5_start\t$ex5_end\n";
      print SEQUIN "\t\t\tproduct\t$product\n";
      print SEQUIN "\t\t\texception\ttrans-splicing\n";
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
    }
    elsif( @{$exons_ref} == 4 ){
      print SEQUIN "$ex1_start\t$ex1_end\tCDS\n";
      print SEQUIN "$ex2_start\t$ex2_end\n";
      print SEQUIN "$ex3_start\t$ex3_end\n";
      print SEQUIN "$ex4_start\t$ex4_end\n";
      print SEQUIN "\t\t\tproduct\t$product\n";
      print SEQUIN "\t\t\texception\ttrans-splicing\n";
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
    }
    elsif( @{$exons_ref} == 3 ){
      print SEQUIN "$ex1_start\t$ex1_end\tCDS\n";
      print SEQUIN "$ex2_start\t$ex2_end\n";
      print SEQUIN "$ex3_start\t$ex3_end\n";
      print SEQUIN "\t\t\tproduct\t$product\n";
      print SEQUIN "\t\t\texception\ttrans-splicing\n";
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
    }
    elsif( @{$exons_ref} == 2 ){
      print SEQUIN "$ex1_start\t$ex1_end\tCDS\n";
      print SEQUIN "$ex2_start\t$ex2_end\n";
      print SEQUIN "\t\t\tproduct\t$product\n";
      print SEQUIN "\t\t\texception\ttrans-splicing\n";
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
    }
  }
  
  close SEQUIN;
  
}
#############################################################################
# this subroutine determines whether a sequence is on the + or - strand
sub strand{
  my( $start, $end ) = @_;
  my $strand;
  my $tmp;
  
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

