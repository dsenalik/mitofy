#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

require "annotate_subroutines.pl";
require "annotate_rna.pl";
require "parse_trnascan.pl";
require "rna_html_pages.pl";
require "print_nt_with_gaps.pl";
require "nad5_ex3.pl";

######################################MAIN###################################
our $PROT_EMAX = 1e-3;
our $PROT_PMIN = 60;
our $RNA_EMAX  = 1e-3;
our $RNA_PMIN  = 70;
our $RNA_MINLEN  = 30;

parseArgs();

my( $infile, $project ) = @ARGV;

my $seqlen;	  #number of nucleotides in input genome, from BLAST output
my $query_taxon;  #taxon name of query sequence, from the fasta header
my $file;	  #stores entire input genome as a string
my $gene;
my $out_directory;
my %proteinNoHits;  #key = three-letter AA code; value = path of HTML "No significant BLAST hits" page; key and value exist only
                    #for those RNAs with no significant BLAST HITS, as determined by $counter boolean at end of annotate_rna
my %html_files;     #hash, key=gene, value=path to html output file for that gene; use this in &print_gene_summary
		    #to create html gene summary page
my %blast_files;    #hash, key=gene, value=path to raw blast output file for that gene; use this in &print_gene_summary
		    #to create html gene summary page
my @missed;
my $fasta_count;    #keep track of number of FASTA files in input file

open( INPUTFILE, "$infile" ) || die "Couldn't open $infile: $!\n"; 
while( <INPUTFILE> ){
  /^$/ and next;
  /^#/ and next;
  if( /^>(\S+)/ ){
    if( $fasta_count ){
      $fasta_count++;
      next;
    }else{
      $query_taxon = $1;
      $fasta_count++;
      next;
    }
  }else{
    chomp( $file .= uc $_ ); #store concatenated genome sequence to $file variable
  }
}
close INPUTFILE;

# print $query_taxon . "\n";
# print $file, "\n";
# print length $file, "\n";
# exit;

$seqlen = length $file;

if( $fasta_count > 1 ){
  $query_taxon = $query_taxon . "_" . $fasta_count;
  open( CONCAT_OUT, ">$infile.concat" );
  print CONCAT_OUT ">$query_taxon\n";
  print CONCAT_OUT $file . "\n";
  close BLASTOUT;
  $infile = "$infile.concat";
}

# create output directory
$out_directory = "blast_output/$project"."_out";

unless( -e "$out_directory" ){
  `mkdir $out_directory`;
}

open( MISSED,">$out_directory/missed.txt" ) || die "Couldn't open missed.txt: $!\n";

my @mt_genes = qw( atp1 atp4 atp6 atp8 atp9 ccmB ccmC ccmFc ccmFn cob cox1 cox2 cox3 matR mttB nad1 nad2 nad3 nad4 nad4L nad5 nad5_ex3 nad6 nad7 nad9 rpl2 rpl5 rpl10 rpl16 rps1 rps2 rps3 rps4 rps7 rps10 rps11 rps12 rps13 rps14 rps19 sdh3 sdh4 );

print "\n\nBLASTing genes . . . \n";

foreach $gene ( @mt_genes ){
  $gene =~ /^\./ and next; #ignores the "." and ".." elements in the directory
  # $gene eq "ccmFc" or next;

  my $blast_data  ; #the entire BLAST output
  my @gene_lengths; #stores lengths subject sequences in BLAST database
  my $sbjct_taxon ; #taxon name of subject sequence from BLAST hit
  my $gene_length ; #total length of sbjct sequence
  my $taxa_count  ;
  my $hit         ;
  my %HITS        ;
  my $counter = 0 ; #hit counter, treats each non-redundant hit as an exon, 

  #print gene being annotated to stdout
  print "\t$gene\n" unless $gene eq "nad5_ex3";

  #run BLAST for current gene
  if( $gene eq "nad5_ex3" ){
    # open( BLASTRUN, "blastall -p blastn -d blast_dbs/mt_genes/$gene/$gene -i $infile | " );
    open( BLASTRUN, "blast/blastn -word_size 9 -reward 2 -penalty -3 -gapopen 5 -gapextend 2 -query $infile -db blast_dbs/mt_genes/$gene/$gene  | " );

    while( <BLASTRUN> ){ #write BLAST output to $blast_data variable
      $blast_data .= $_ ;
    }
    close BLASTRUN;
    
    nad5_ex3( $query_taxon, $blast_data, $out_directory );

    #save raw BLAST output
    open( BLASTOUT, ">$out_directory/$gene.blastn" );	#NEW
    print BLASTOUT $blast_data;
    close BLASTOUT and next;

  }else{ #all other protein genes
    # open( BLASTRUN, "blastall -F F -e 50 -p blastx -d blast_dbs/mt_genes/$gene/$gene -i $infile | " );
    open( BLASTRUN, "blast/blastx -query $infile -db blast_dbs/mt_genes/$gene/$gene  | " );
    while( <BLASTRUN> ){ #write BLAST output to $blast_data variable
      $blast_data .= $_ ;
    }
    close BLASTRUN;
       
    #save raw BLAST output
    open( BLASTOUT, ">$out_directory/$gene.blastx" );	#NEW
    print BLASTOUT $blast_data;
    close BLASTOUT;
    
    #store locations of html and raw blast output files for each gene
    $html_files{$gene} = "\.\.\/$project" . "_out/" . "$gene.html";
    $blast_files{$gene} = "\.\.\/$project" . "_out/" . "$gene.blastx";
    
    #start parsing BLAST output
    if( $blast_data =~ /No hits found/ ){ 
      print MISSED "No hits found for $gene.\n";
      $html_files{$gene} = not_found( "BLAST", $gene, $out_directory, $project );
      push @missed, $gene;
      next;
    }
    
    #split BLAST output by taxon
    for( split( />/, $blast_data ) ){
      /Altschul/ and next; #skip BLAST header

      # parse taxon name and gene length, for NCBI-formatted protein sequence headers
      # if( /.*\[(\w*\s*\w*)\]\s*Length = (\d+)/ ){
      # if( /.*\[([\s\S]*)\]\s*Length = (\d+)/ ){
      if( /\[(.*)\]\s*Length=(\d+)/s ){
	$sbjct_taxon = $1;
	$gene_length = $2;
	$sbjct_taxon =~ s/\s+/_/g;
	$sbjct_taxon =~ s/\.//g;
	$taxa_count++;
	
	#store non-redundant gene lengths in @gene_lengths
	if( $taxa_count == 1 ){
	  push( @gene_lengths, $gene_length );
	}else{
	  unless( grep( /$gene_length/, @gene_lengths ) ){
	    push( @gene_lengths, $gene_length );
	  }
	}
	
      #parse taxon name and gene length, for simple FASTA headers (e.g., ">Cucurbita") 
      # }elsif( /^(\S+)\s*Length = (\d+)/ ){
      }elsif( /^ (\S+)\s*Length=(\d+)/ ){
	$sbjct_taxon = $1;
	$gene_length = $2;
	$taxa_count++;
	#store non-redundant gene lengths in @gene_lengths
	if( $taxa_count == 1 ){
	  push( @gene_lengths, $gene_length );
	}else{
	  unless( grep( /$gene_length/, @gene_lengths ) ){
	    push( @gene_lengths, $gene_length );
	  }
	}
      }
      
      # parse all of the hits for this taxon
      parse_hits_within_taxon( $query_taxon, $sbjct_taxon, \$file, \%HITS, $gene, $gene_length );
      
    }

    my $rna_bool = 0;
    
    print_output_header( $project, $query_taxon, $gene, $gene_length, $out_directory, \@gene_lengths, $rna_bool );
    
    # sort unique hits by $query_start, then parse each hit
    foreach $hit ( sort {$HITS{$a}[1] <=> $HITS{$b}[1]} keys ( %HITS ) ){ #sorts by sbjct_begin
      my( $sbjct_taxon, $sbjct_begin, $sbjct_end, $sbjct_match, $strand, $percent_match, $hit_e_value, $query_begin, $query_end, $query_match, $gene_length, $match_length ) = @{$HITS{$hit}};
      
      if( $percent_match >= $PROT_PMIN && $hit_e_value <= $PROT_EMAX ){
	$counter++;
	print_sbjct_AA( $sbjct_taxon, uc $sbjct_match, $sbjct_begin, $sbjct_end, $query_begin, $query_end, $strand, $percent_match, $gene_length );
	print_query_AA( $query_taxon, uc $query_match, $query_begin, $query_end );
	
	# print "\n$sbjct_taxon $sbjct_begin $sbjct_end\n";
	# print $query_match, "\n";
	# print $sbjct_match, "\n";
	
	# extract corresponding nt sequence from query genome, to display below the AA sequences; don't need to worry about gaps in AA sequence
	extract_nt_sequence( $query_match, $query_begin ,$query_end, $strand, \$file );
      }
    }
    print_output_footer();
    make_dna_no_hits( $gene, \%proteinNoHits, $out_directory, $project ) if $counter == 0;
    $html_files{$gene} = $proteinNoHits{$gene} if $counter == 0;;
  }
}

annotate_rna( $infile, $project, $query_taxon, $out_directory, $seqlen, \$file, $RNA_EMAX, $RNA_PMIN, $RNA_MINLEN );
print_gene_summary( $project, $query_taxon, $seqlen, \%html_files, \%blast_files, $out_directory, \@missed );

print "\n\n** Open \'$out_directory/$project\_summary.html\' and \'$out_directory/$project\_rna\_summary.html\' in a web browser (preferably Safari) to see results. **\n\n";

close MISSED;
exit;

#################################SUBROUTINES#################################
# Main (above) splits the BLAST output by taxon and loops over all the taxa.
# For each taxon, parse_hits_within_taxon is called, which loops over the 
# multiple hits within each taxon

sub parse_hits_within_taxon{

  my( $query_taxon, $sbjct_taxon, $file_ref, $HITS_ref, $gene, $gene_length ) = @_;

  #split by match within each gene	
  foreach my $hit( split( /Score/ ) ){
    my $hit_e_value  ; #e-value of hit
    my $match_length ; #length of match between subject & query
    my $percent_match; #the %identity between subject & query
    my $strand	     ; #is hit on + or - strand
    my $query_match  ; #query sequence parsed directly from BLAST hit
    my $query_begin  ; #start coordinate of hit in query sequence
    my $query_end    ; #end coordinate of hit in query sequence
    my $sbjct_match  ; #sbjct sequence parsed directly from BLAST hit
    my $sbjct_begin  ; #start coordinate of hit in subject sequence
    my $sbjct_end    ; #end coordinate of hit in subject sequence
    my $upstream     ; #nt sequence upstream of the hit, for identifying alternative start codons, etc.
    my $downstream   ; #nt sequence upstream of the hit, for identifying alternative start codons, etc.
    my $gene_nt_sequence; #nt sequence of gene
    my $start;
    my $gaps;
    my @hit_data;
    my @sbjct_coords;
    my @query_coords;

    $hit =~ /Length=/ and next; #skip database sequence header

    foreach my $line ( split /\n/, $hit ){
      $line =~ /\S+/ or next;
      $line =~ /Lambda/ and last;

      # parse e-value
      if( $line =~ /Expect[\(\d+\)]* =\s*(\de-\d+)/ ){
	$hit_e_value = $1;
      }elsif( $line =~ /Expect(\(.*\))\s*=\s*(.+)/ ){
	$hit_e_value =  $2;
	if( $2 =~ /^e/ ){
	  $hit_e_value = "1$hit_e_value";
	}
      }elsif( $line =~ /Expect[\(\d+\)]* =\s*(e-\d+)/ ){
	$hit_e_value =  "1$1";
      }elsif( $line =~ /Expect =\s*(\d+\.\d*)/ ){
	$hit_e_value =  $1;
      }elsif( $line =~ /Expect =\s*(\d+)/ ){
	$hit_e_value =  $1;
      }
      
      # parse length of match from number of identities
      if( $line =~ /Identities = (\d+)\/(\d+) \((\d+)\%\),/ ){
	$match_length = $2;
	$percent_match = $3;
      }
      
      # parse forward or reverse strand
      if( $line =~ /Frame = ([+-]+)(\d)/ ){
	$strand = $1;
      }
      
      # get number of gaps
      if( $line =~ /Gaps = (\d+)\/(\d+)/ ){
	$gaps = $1;
      }
      
      #parse query coords
      if( $line =~ /Query  (\d+)\s+(\S+)\s+(\d+)/ ){ #parse begin and end coords of query sequence
	push @query_coords, $1, $3;
	$query_match .= $2;
      }

      #parse sbjct coords
      if( $line =~ /Sbjct  (\d+)\s+(\S+)\s+(\d+)/ ){ #parse begin and end coords of query sequence
	push @sbjct_coords, $1, $3;
	$sbjct_match .= $2;
      }
    }
    
    $sbjct_begin = min( @sbjct_coords );
    $sbjct_end   = max( @sbjct_coords );
    $query_begin = min( @query_coords );
    $query_end   = max( @query_coords );
 
    my $score = $percent_match * (($query_end - $query_begin) + 1);
    
    #store hit data in %HITS, indexed by "$sbjct_begin.$sbjct_end"; this effectively filters redundant hits *ACROSS TAXA* and allows hits to later be 
    #sorted by $sbjct_begin, so hits can by displayed by the order of the exons. Sorting done in Main.
    #if $HITS{"$query_begin.$query_end"} is already defined, replace it with current redundant hit only if it has a better score
    if( $HITS_ref->{"$query_begin.$query_end"} ){
      if( $score > $HITS_ref->{"$query_begin.$query_end"}->[12] ){
	push @{$HITS_ref->{"$query_begin.$query_end"}}, ( $sbjct_taxon, $sbjct_begin, $sbjct_end, $sbjct_match, $strand, $percent_match, $hit_e_value, $query_begin, $query_end, $query_match, $gene_length, $match_length, $score );
      }else{
	next;
      }
    }else{
      push @{$HITS_ref->{"$query_begin.$query_end"}}, ( $sbjct_taxon, $sbjct_begin, $sbjct_end, $sbjct_match, $strand, $percent_match, $hit_e_value, $query_begin, $query_end, $query_match, $gene_length, $match_length, $score );
    }
  }
}
#############################################################################
sub parseArgs{

my $usage = "\nusage:$0 [options] genome.fasta projectID\n" .
"      projectID = unique project name; all output files with have this prefix

options
      --prot_emax - maximum BLAST expect value for protein genes (default: 1e-3)
      --prot_pmin - minimum percent identity for protein genes (default: 60)
      --rna_emax  - maximum BLAST expect value for RNA genes (default: 1e-3)
      --rna_pmin  - minimum percent identity for RNA genes (default: 70)
      --rna_mlen  - minimum length for RNA BLASTN matches (default: 30)
\n";

  my $result = GetOptions
    (
     'prot_emax=s' => \$PROT_EMAX,
     'prot_pmin=s' => \$PROT_PMIN,
     'rna_emax=s'  => \$RNA_EMAX,
     'rna_pmin=s'  => \$RNA_PMIN,
     'rna_mlen=s'  => \$RNA_MINLEN,
    );
  
print $usage and exit unless $ARGV[1];

}
#############################################################################
