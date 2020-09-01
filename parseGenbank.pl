#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

###############################################MAIN#################################################
our( $NTSEQ, $AASEQ, $GFF, $FASTA, $LIST );

parseArgs();

my $file = shift @ARGV;

$file =~ /(.*)\..*/ and my $outname = $1;

open( INFILE, "$file" ) || die "Couldn't open $file: $!\n";

my( $line_number, $in_features, $in_sequence, $genome_size, $dna_sequence, $accession, $features, $binomial, $genus, $definition );
my @gff; #each index = GFF line

while( <INFILE> ){
  if( m|^//| ){
    last;
  }elsif( $in_sequence ){
    $dna_sequence .= $_;
  }elsif( /^ORIGIN/ ){
    $in_sequence = 1;
  }elsif( $in_features ){
    $features .= $_;
  }elsif( /^FEATURES/ ){
    $in_features = 1;
  }elsif( /\s*ORGANISM\s+((.+)\s+(.+))/i ){ #read the organism name from the "ORGANISM" line
    $binomial = $1;
    $genus = $2;
    $binomial =~ s/\s+/_/g;
    $binomial =~ s/_$//g;
  }elsif( /DEFINITION\s*(.+)/ ){
    $definition = $1;
    #print "$definition\n";
  }elsif( /\s*ACCESSION\s+(.+)/i ){ #read the accession number from the "ACCESSION" line
    $accession = $1;
  }
}

$dna_sequence =~ s/[0-9\s]//g;
if( $FASTA ){
  open( OUT, ">$outname.fasta" ) || die "Couldn't open $outname.fasta: $!\n";
  print OUT ">$outname\n";
  print OUT uc $dna_sequence . "\n";
  close OUT;
}

#calculate statistics on the genome sequence
genome_statistics( \$dna_sequence, $binomial );

#parse features
parse_features( $outname, $features, $definition, \$dna_sequence, $binomial, $genus, \@gff, $NTSEQ, $AASEQ, $GFF, $LIST );

close INFILE;

############################################SUBROUTINES#############################################
sub parseArgs{
  
  my $usage = "\nUsage: $0 [options] file.gb
	
options
      --nt   output nucleotide fasta files (default: false)
      --aa   output amino acid fasta files (default: false)
      --gff  output GFF file (default: false)
      --fas  output genome sequence in FASTA format (default: false)
      --list output list of protein and rRNA genes (default: false)
\n

*Note: This script does not capture \"promoter\", \"mRNA\", \"exon\", or \"gene\" annotations.
       It should capture all \"CDS\", \"rRNA\", \"tRNA\", \"intron\", and \"misc_feature\" annotations.\n\n"
;
 
  my $result = GetOptions
    (
     'nt!'    => \$NTSEQ,
     'aa!'    => \$AASEQ,
     'gff!'   => \$GFF,
     'fas!'   => \$FASTA,
     'list!'  => \$LIST,
    );
  
  if( @ARGV == 0 ){
    print $usage and exit;
  }
  
}
####################################################################################################
# Called by parseGenbank.pl; the following subroutine parses out the features table of a GenBank-  #
# formatted file of a complete genome sequence                                                     #
####################################################################################################
sub parse_features{
  my( $outname, $features, $definition, $dna_sequence, $binomial, $genus, $gff_ref, $NTSEQ, $AASEQ, $GFF, $LIST ) = @_;
  my @features = ();
  my $feature_key;
  my %aa;  #key = gene name, value = AA sequence
  my %rna; #key = gene name, value = defined (for outputting gene list)
  my %cds; #key = gene name, value = defined (for outputting gene list)

  # this regular expression matches an entire feature, then saves the entire feature into the 
  # variable $feature; every $feature is then stored in the @features array
  while( $features =~ /^( {5}\S+.*\n(^ {21}\S+.*\n)*)/gm ){
    push( @features, $1 );
  }

  # starting with the array of features for this genome (each element is a single, full feature), 
  # loop over features, extracting the feature name and parsing relevant features
  
  my $i = 0; #increment $i, and append $i to $gene in case there are duplicate copies of the gene
             #this is solely for dealing with the AA files, to avoid overwriting in %aa and
             #outputting duplicate file names
	
  foreach my $feature ( @features ){
    $i++;
    my $cds_coord = "";
    my $note      = "";
    my $gene      = "";
    my $rpt_type  = "";
    
    # extract feature name into $feature_key
    $feature =~ /^ {5}(\S+)/ and $feature_key = $1;
    #print $feature_key . "\n";
      
    # the first capture gets features of type "CDS", "rRNA", "tRNA", and "misc_feature" with "gene" qualifier key; do
    # not want to capture features of type "gene" as this captures the entire gene+intron ranges
    # for intron-containing genes, and creates redundant captures for non-intron-containing genes;
    # it also captures rRNAs as both "gene" and "rRNA" features, and this messes up printing out
    # the GFF file; also skips features of type "exon", which if captured, results in redundant
    # entries (CDS and exon) for the same exon in the GFF file

		
    # skip header
    next if $feature_key eq "source";
    # skip "promoter" features
    next if( $feature =~ /promoter {4,13}((\S+\n)( {21}\S+\n)*)( {21}\/gene=\"(.+)\")/gm );
    # skip single-nucleotide misc_features, i.e., RNA edits
    next if ( $feature =~ /misc_feature {4,13}((\S+\n)( {21}\S+\n)*)/gm ) and $1 !~ /\.\./;
    # skip "mRNA" features
    next if( $feature =~ /mRNA {4,13}((\S+\n)( {21}\S+\n)*)( {21}\/gene=\"(.+)\")/gm );
    # skip "exon" features
    next if( $feature =~ /exon {4,13}((\S+\n)( {21}\S+\n)*)( {21}\/gene=\"(.+)\")/gm );
    # skip "gene" features
    # next if( $feature =~ /gene {4,13}((\S+\n)( {21}\S+\n)*)( {21}\/gene=\"(.+)\")/gm );

    # capture pseudogenes
    # work on following code to include all gene features
    if( $feature =~ /pseudo/ ){
      if( $feature =~ /gene {4,13}((\S+\n)( {21}\S+\n)*)( {21}\/gene=\"(.+)\")/gm ){
	$cds_coord = $1; # extract coordinates
	$gene      = $5; # extract gene name
	$cds_coord =~ s/[\s\n]//g; # eliminate whitespace from $cds_coord
	$feature_key = "pseudo";

	# print $feature, "\n";
	# print $gene . "\n";
	# print $cds_coord . "\n";

	#parse one-line notes
	if( $feature =~ /note="(.+)"/ ){
	  $note = $1;
	  $note =~ s/\n/ /g;
	  $note =~ s/\s+/ /g;
	  $note =~ s/\;/,/g;
	  #print "[[$note]]\n";
	  
	  #parse multi-line notes
	}elsif( $feature =~ /note="(.+)"/s ){
	  $note = $1;
	  $note =~ s/\n/ /g;
	  $note =~ s/\s+/ /g;
	  $note =~ s/\;/,/g;
	  #print "[[$note]]\n";
	}
      }else{
	next;
      }
    
    #capture features of type "CDS", "rRNA", "tRNA", and "misc_feature" with "gene" qualifier key
    }elsif( $feature =~ /$feature_key {4,13}((\S+\n)( {21}\S+\n)*)( {21}\/gene=\"(.+)\")/gm ){
      $cds_coord = $1;		 # extract exon coordinates
      $gene	 = $5;  	 # extract gene name
      $cds_coord =~ s/[\s\n]//g; # eliminate whitespace from $cds_coord
      if( $gene =~ /rrn/i ){
	$rna{$gene} = 1;
      }

       # print $feature, "\n";
       # print $gene . "\n";
       # print $cds_coord . "\n";

      #parse one-line notes
      if( $feature =~ /note="(.+)"/ ){
	$note = $1;
	$note =~ s/\n/ /g;
	$note =~ s/\s+/ /g;
	$note =~ s/\;/,/g;
	#print "[[$note]]\n";

      #parse multi-line notes
      }elsif( $feature =~ /note="(.+)"/s ){
	$note = $1;
	$note =~ s/\n/ /g;
	$note =~ s/\s+/ /g;
	$note =~ s/\;/,/g;
	#print "[[$note]]\n";
      }

      #parse different flavors of codon information for tRNAs
      if( $feature =~ /(codon_recognized=\"([A-Z]{3})\")/ ){
	my $tmp = $1;
	if( $note =~ /\S/ ){
	  $note .= ", $tmp";
	}else{
	  $note = $tmp;
	}
      }elsif( /anticodon:([A-Z]{3})/ ){
	if( $note =~ /\S/ ){
	  $note .= ", $1";
	}else{
	  $note = $1;
	}
      }

#     if( $LIST && $feature !~ /chloroplast/ ){ #for the Vitis mt genome, which has cp genes as CDS features
      if( $LIST && $feature !~ /chloroplast/ ){
	$feature =~ /translation="(.+)"?/s and $cds{$gene} = 1;
      }
      
      #parse amino acid sequence
      if( $AASEQ ){
	if( $feature =~ /translation="(.+)"?/s ){
	  my $pep = $1;
	  $pep =~ s/[\s\n"]//g;
	  if( $aa{$gene} ){
	    $aa{"$gene.$i"} = $pep;
	  }else{
	    $aa{$gene} = $pep;
	  }
	}
      }
      
    #capture features of type "misc_feature" without "gene" as qualifier key
    }elsif( $feature =~ /misc_feature {4,13}((\S+\n)( {21}\S+\n)*)/gm ){
      #print $feature;
      $cds_coord = $1;		 #extract exon coordinates
      $cds_coord =~ s/[\s\n]//g; #eliminate whitespace from $cds_coord
      
      #parse one-line notes
      if( $feature =~ /note="(.+)"/ ){
	$note = $1;
	$note =~ s/\n/ /g;
	$note =~ s/\s+/ /g;
	$note =~ s/\;/,/g;
	#parse multi-line notes
      }elsif( $feature =~ /note="(.+)"/s ){
	$note = $1;
	$note =~ s/\n/ /g;
	$note =~ s/\s+/ /g;
	$note =~ s/\;/,/g;
      }

    #capture repeat_region features, with or without note
    }elsif( $feature =~ /repeat_region {3}((\S+\n)( {21}\S+\n)*)/gm ){
      $cds_coord = $1;			# extract exon coordinates
      $cds_coord =~ s/[\s\n]//g;	# eliminate whitespace from $cds_coord
      
      #parse one-line notes
      if( $feature =~ /note="(.+)"/ ){
	$note = $1;
	$note =~ s/\n/ /g;
	$note =~ s/\s+/ /g;
	$note =~ s/\;/,/g;
	#parse multi-line notes
      }elsif( $feature =~ /note="(.+)"/s ){
	$note = $1;
	$note =~ s/\n/ /g;
	$note =~ s/\s+/ /g;
	$note =~ s/\;/,/g;
      }
      
      #parse repeat type
      if( $feature =~ /rpt_type=(.+)/ ){
	$rpt_type = $1;
	$rpt_type =~ s/\n/ /g;
	$rpt_type =~ s/\s+/ /g;
	$rpt_type =~ s/\;/,/g;
      }
      
    #this is specifically for rRNAs not annotated as genes
    }elsif( $feature =~ /rRNA {12}(\S+\n)( {21}\/product=\"(.+)\")/gm ){
      $cds_coord = $1;	 #extract exon coordinates
      $gene	 = $3;	 #extract gene name
      $cds_coord =~ s/[\s\n]//g; #eliminate whitespace from $cds_coord
      
      #in cases where CDS features do not have a corresponding "gene" feature key, but do have
      #a "product" feature key
    }elsif( $feature =~ /[CDS|tRNA] {12,13}((\S+\n)( {21}\S+\n)*)( {21}\/product=\"(.+)\")/gm ){
      $cds_coord = $1;	  #extract exon coordinates
      $gene	 = $5; 	  #extract gene name
      $gene =~ s/\s+/_/g;
      $cds_coord =~ s/[\s\n]//g;  #eliminate whitespace from $cds_coord

      #parse different flavors of codon information for tRNAs
      if( $feature =~ /(codon.*recognized[=:]\s*\"*([A-Z]{3}))\"/ ){
	if( $note =~ /\S/ ){
	  $note .= ", $1";
	}else{
	  $note = $1;
	}
      }elsif( /anticodon:([A-Z]{3})/ ){
	if( $note =~ /\S/ ){
	  $note .= ", $1";
	}else{
	  $note = $1;
	}
      }
    }else{
      next;
    }
    
    #extract_sequence parses out everything necessary to print out to GFF file and extracts sequences to write out nucleotide FASTA files
    extract_sequence( $outname, $dna_sequence, $definition, $gene, $feature_key, $cds_coord, $genus, $gff_ref, $rpt_type, $NTSEQ, $AASEQ, $GFF, $note, $binomial  );			
  }

  #output AA sequences
  printAA( $binomial, $outname, $genus, \%aa );
  
  #ouput list of gene names
  if( $LIST ){
    my $outfile = "$outname" . ".genes";
    open( GENES, ">$outfile" ) || die "Couldn't open $outfile.\n";
    print GENES "$outname\n";
    foreach( sort keys %cds ){
      print GENES "$_\n";
    }
    foreach( sort keys %rna ){
      print GENES "$_\n";
    }
    close GENES;
  }

#sort and print gff file
  if( $GFF ){
    my $outfile = "$outname" . ".gff"; 
    open( GFFOUT, ">$outfile" ) || die "Couldn't open $outfile.\n";
    
    my $output = [sort {$a->[3] <=> $b->[3]} @$gff_ref];
    
    foreach( @$output ){
      my $line = join( "\t", @$_ );
      print GFFOUT $line . "\n";
    }#foreach $line
    
    close GFFOUT;
  }
  
}

####################################################################################################
sub printAA{
  my( $binomial, $outname, $genus, $aa_ref ) = @_;
  
  foreach( sort keys %{$aa_ref} ){
#   my $outfile = $outname . "_" . $_ . ".fasta";
    my $outfile = $_ . ".fasta";
    chomp( my $dir = `pwd` );
    my $outdir = $dir . "/" . "$outname" . "_aa"; #DIRECTORY
    
    unless( -e $outdir){
      `mkdir $outdir`;
    }
    
    open(OUT, ">$outdir/$outfile") || die "Couldn't open $outfile: $!\n";
    #print OUT ">", $binomial , "_", $_, "\n", $aa_ref->{$_}, "\n";
    print OUT ">", $binomial , "\n", $aa_ref->{$_}, "\n";
    close OUT;
  }
}
####################################################################################################
# The following subroutine calculates the length of the genome and calls the nuc_frequencies       #
# subroutine to calculate nucleotide composition                                                   #
####################################################################################################
sub genome_statistics{
  my( $dna_sequence, $binomial ) = @_;

  $$dna_sequence =~ s/[0-9\s]//g;
  my $genome_size = length $$dna_sequence;

  my $gc = nuc_frequencies( $dna_sequence, $genome_size );

  printf "$binomial  $genome_size nt  %.2f G+C\n", $gc;
}

####################################################################################################
# The following subroutine calculates nucleotide frequencies for the genome it accepts the genome  #
# sequence and the sequence length as inputs from the genome_statistics subroutine		   #
####################################################################################################
sub nuc_frequencies{
	
  my( $dna_sequence, $genome_size ) = @_;
  
  my $count_of_A  = 0;
  my $count_of_C  = 0;
  my $count_of_G  = 0;
  my $count_of_T  = 0;
  my $gc_content  = 0;
  my $errors      = 0;
  
  for( my $position = 0; $position < $genome_size; ++$position ){
    my $nt = substr( $$dna_sequence, $position, 1 );
    if ($nt =~ /A/i){
      ++$count_of_A;
    }elsif ($nt =~ /C/i){
      ++$count_of_C;
    }elsif ($nt =~ /G/i){
      ++$count_of_G;
    }elsif ($nt =~ /T/i){
      ++$count_of_T;
    }else{
      print "The nucleotide \"$nt\" at position $position is ".
	"not recognized.\n\n";
      ++$errors;
    }
  }
  
  my $percent_A   = $count_of_A/$genome_size;
  my $percent_C   = $count_of_C/$genome_size;
  my $percent_G   = $count_of_G/$genome_size;
  my $percent_T   = $count_of_T/$genome_size;
  my $size_verify = $count_of_A + $count_of_C + $count_of_G + $count_of_T + $errors;
  
  $gc_content	= (($count_of_C + $count_of_G) / $genome_size)*100;
  # 	printf ("\tFrequency of A = $count_of_A (%.2f)\n", $percent_A);
  # 	printf ("\tFrequency of C = $count_of_C (%.2f)\n", $percent_C);
  # 	printf ("\tFrequency of G = $count_of_G (%.2f)\n", $percent_G);
  # 	printf ("\tFrequency of T = $count_of_T (%.2f)\n", $percent_T);
  # 	printf ("\tGC content	= %.2f\n", $gc_content);
  
  return $gc_content;
}
####################################################################################################
# Called by parse_features.pl, extracts feature sequences from genome and calls print_sequences 
# and format_gff to output the fasta file and print the feature to a GFF file
####################################################################################################
sub extract_sequence{
  my( $outname, $dna_sequence, $definition, $gene, $feature_key, $cds_coord, $genus, $gff_ref, $rpt_type, $NTSEQ, $AASEQ, $GFF, $note, $binomial ) = @_;
  my( $sequence, $complement, $extracted_sequence, $i, $strand, $start, $end );
  # print $genus;
  # print $gene;
  # print "$gene, extract sequence: " . $cds_coord . "\n\n";
  # print $note . "\n" if defined $note;
  # print $feature_key, "\n";

  my @coords = split( /,/, $cds_coord );

  #if the coordinates begin with "complement", then the entire sequence is
  #the reverse complement, only need get_complement()
  if( $coords[0] =~ /^complement/i ){
    $strand = "-";
    for( $i = 0; $i < @coords; $i++ ){
      ( $start, $end, $extracted_sequence ) = get_sequence( $coords[$i], $dna_sequence );
      $sequence .= $extracted_sequence;
      scalar @coords > 1 and my $exon_number = (scalar @coords) - $i;
      if( $GFF && $feature_key ne "gene" ){
	format_gff( $definition, $feature_key, $gene, $start, $end, $strand, $gff_ref, $rpt_type, $note, $exon_number  );
      }
    }
    $sequence = complement( $sequence );
    #entire set of coordinates not on reverse strand
  }else{
    for( $i = 0; $i < @coords; $i++ ){
      # if coordinate is a reverse complement
      if( $coords[$i] =~ /complement/ ){
	$strand = "-";
	( $start, $end, $complement ) = get_complement( $coords[$i], $dna_sequence );
	$sequence .= $complement;
	scalar @coords > 1 and my $exon_number = $i+1;
	if( $GFF && $feature_key ne "gene" ){
	  format_gff( $definition, $feature_key, $gene, $start, $end, $strand, $gff_ref, $rpt_type, $note, $exon_number  );
	}
      }else{ 
	$strand = "+";
	( $start, $end, $extracted_sequence ) = get_sequence( $coords[$i], $dna_sequence );
	$sequence .= $extracted_sequence;
	scalar @coords > 1 and my $exon_number = $i+1;
	if( $GFF && $feature_key ne "gene" ){
	  format_gff( $definition, $feature_key, $gene, $start, $end, $strand, $gff_ref, $rpt_type, $note, $exon_number  );
	}
      }
    }
  }
  
  if(( $NTSEQ ) && ( $feature_key eq "CDS" || $feature_key eq "tRNA" || $feature_key eq "rRNA" || $feature_key eq "intron" )){
    my $seq_ID = $genus . "_" . $gene;
    print_sequence( $outname, $seq_ID, $gene, $feature_key, $sequence, $genus, $binomial );
  }
  
}
####################################################################################################
sub get_complement {
  my( $coords, $dna_sequence ) = @_;
  my $offset;
  my $start;
  my $end;
  my $complement = "";
  
  # match "number..number"
  if( $coords =~ /(\d+)..[<>]*(\d+)/ ){
    # set $offset to first number
    $start = $1;
    $end = $2;
    $offset = $1 - 1;
    # $length = length of sequence
    my $length = ($2 - $1)+1;
    # extract the sequence
    my $sequence = substr( $$dna_sequence, $offset, $length );
    # convert sequence to its reverse complement
    $complement = reverse( $sequence );
    $complement =~ tr/53ACGTUMRWSYKVHDBNacgtumrwsykvhdbn/35TGCAAKYWSRMBDHVNtgcaakywsrmbdhvn/;
  }
  # return the extracted sequence
  return( $start, $end, $complement );
}
####################################################################################################
sub get_sequence {
  my( $coords, $dna_sequence ) = @_;
  my $offset;
  my $start;
  my $end;
  my $sequence = "";
  
  #match "number..number"
  if( $coords =~ /(\d+)..[<>]*(\d+)/ ){
    # set $offset to first number
    $start = $1;
    $end = $2;
    $offset = $1 - 1;
    
    my $length = ($2 - $1)+1;
    # extract the sequence
    $sequence = substr( $$dna_sequence, $offset, $length );
  }
  # return the extracted sequence
  return( $start, $end, $sequence );
}
####################################################################################################
sub complement{
  my $sequence = pop;
  
  my $complement = reverse( $sequence );
  $complement =~ tr/ACGTNacgtn/TGCANtgcan/;
  
  return( $complement );
}
####################################################################################################
# Called by extract_sequence.pl, formats the feature for GFF output and puts
# the feature into a hash (key="start.$end" coordinate, value=GFF entry), so that GFF
# lines can be sorted numerically by start coordinate of the feature, then printed
# to a GFF file
####################################################################################################
sub format_gff{
  my( $definition, $feature_key, $gene, $start, $end, $strand, $gff_ref, $rpt_type, $note, $exon_number ) = @_;
  my( $type, $color, $attribute_type, @attributes );
  my $length = $end - $start + 1;
  #print $note . "\n" if defined $note;
 
  if( $feature_key eq "tRNA" || $feature_key eq "rRNA" ){
    $attribute_type = "exon";
    $color = "red";
  }elsif( $feature_key eq "CDS" ){
    $attribute_type = "exon";
    $color = "blue";
  }elsif( $feature_key eq "intron" ){
    $attribute_type = "intron";
    $color = "green";
  }elsif( $feature_key eq "misc_feature" ){
    $attribute_type = undef;
    $color = "black";
  }elsif( $feature_key eq "repeat_region" ){
    $attribute_type = undef;
    $color = "grey47";
  }elsif( $feature_key eq "gene" ){
    $color = "pink";
  }elsif( $feature_key eq "promoter" ){
    $color = "yellow";
  }elsif( $feature_key eq "mRNA" ){
    $attribute_type = "mRNA";
    $color = "blue";
  }elsif( $feature_key eq "pseudo" ){
    $attribute_type = "pseudogene";
    $color = "black";
  }else{
    $color = "grey47";
  }
  
  # 	print $definition . "\n";
  # 	print $gene . "\t" . $feature_key . "\t" . $color . "\n";

  my $line = [$definition, "GenBank", $feature_key, $start, $end, $length, $strand, ".", ];
  
  if( $gene ){
    if( $exon_number ){
      push( @attributes, "Gene $gene exon $exon_number" ); #9.1-attribute-gene
    }else{
      push( @attributes, "Gene $gene" ) if $gene;	#9.1-attribute-gene
    }
    
    #if $rpt_type = defined, append to note
    if( $rpt_type ){
      if( $note =~ /\S/ ){
	push( @attributes, "Note $note" );                  #9.2-attribute-note/rpt_type
	push( @attributes, "REPEAT $rpt_type repeat" );     #9.2-attribute-note/rpt_type
      }else{
	push( @attributes, "Note $note" ) if $note =~ /\S/; #9.2-attribute-note
      }
    }else{
      push( @attributes, "Note $note" ) if $note =~ /\S/; #9.2-attribute-note
    }
    
  #else $gene = undef
  }else{
    if( $note =~ /\S/ ){
      push( @attributes, "Gene $note" ) if $note =~ /\S/;	#9.2-attribute-gene (gene = undef && $note = defined)
    }else{
      push( @attributes, "Gene $rpt_type repeat" ) if $rpt_type; #9.2-attribute-note/rpt_type
    }
  }
  #print "$start..$end\n";
  push( @attributes, "COLOR $color" ); #9.2-attribute-color
  push( @attributes, "type $attribute_type" ) if $attribute_type; #9.3-attribute-type
  
  my $at = join( " ; ", @attributes );
  $$line[8] = $at;
  push( @$gff_ref, $line );
  
}
####################################################################################################
sub print_sequence{
  my( $outname, $seq_ID, $gene, $feature_key, $sequence, $genus, $binomial ) = @_;
  my $outfile;
  
  chomp( my $dir = `pwd` );
  my $outdir = $dir . "/" . "$outname" . "_nt"; #DIRECTORY
  #my $outdir = $dir . "/nt"; #DIRECTORY
  unless( -e $outdir ){
    `mkdir $outdir`;
  }
  
  #if( $copy_number > 0 ){
  #$copy_number++;
  #$gene = "$gene.$copy_number";
  #print $gene . "\n";
  #}
  
  if( $feature_key eq "intron" ){
    $outfile = $gene . "_" . $genus . "_" . $feature_key . ".fasta";
    $seq_ID = $seq_ID . "-intron";
    $outfile =~ s/ /_/g;
  }else{
    #print $gene . "\n";
    #$outfile = $gene . "_" . $genus . ".fasta";
    $outfile = $gene . ".fasta";
    $outfile =~ s/ /_/g;

    open(OUT, ">>$outdir/$outfile") || die "Couldn't open $outfile: $!\n";
    $sequence = uc( $sequence );
#   print OUT ">$binomial\n";
#   print OUT ">$binomial $gene\n";
    print OUT ">$genus", ", ", "$gene\n";
#   print OUT ">$genus\n";
#   print OUT ">$gene Silene_latifolia genomic\n";
#   print OUT ">$genus", "_sativus\n";
#   print OUT ">$gene" . "_cuke_cDNA\n";
    print OUT $sequence, "\n";
  }

#   open(OUT, ">>$outdir/$outfile") || die "Couldn't open $outfile: $!\n";
#   $sequence = uc( $sequence );
#   #print OUT ">", $genus , "_", $gene, "\n";
#   print OUT ">", $outname . "\n";
#   print OUT $sequence, "\n";
  
  close OUT;
}
####################################################################################################
