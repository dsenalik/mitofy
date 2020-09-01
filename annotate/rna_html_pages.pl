use strict;
use warnings;
my $mitofybinbase = $ENV{"mitofybinbase"};  # VCRU addition this is where mitofy.pl and other called programs are located
my $webbase = $ENV{"webbase"};  # VCRU addition this is the URL of the mitofy cgi directory (not used here, only shown for complete list)

#######################################################################################################################
# called by annotate_rna

sub rna_output_header{
  my( $project, $query_taxon, $rna, $gene_length, $out_directory, $gene_lengths_ref, $rna_bool ) = @_;
  my %trna = (
	      Phe => "F",
	      Leu => "L",
	      Ile => "I",
	      Met => "M",
	      Val => "V",
	      Ser => "S",
	      Pro => "P",
	      Thr => "T",
	      Ala => "A",
	      Tyr => "Y",
	      His => "H",
	      Gln => "Q",
	      Asn => "N",
	      Lys => "K",
	      Asp => "D",
	      Glu => "E",
	      Cys => "C",
	      Trp => "W",
	      Arg => "R",
	      Gly => "G",
	      Undet => "UNDETERMINED tRNA"
	     );

  my $gene_name;
  
  if( $rna =~ /rrn/ ){
    $gene_name = $rna;
  }elsif( ($rna =~ /fM/i) || ($rna =~ /met-f/i) ){
    $gene_name = "trnfM"
  }elsif( $rna =~ /(\w+)-cp/ ){
    $gene_name = "trn" . $trna{$1} . "-cp";
  }elsif( $rna =~ /(\w+)-mt/ ){
    $gene_name = "trn" . $trna{$1}. "-mt";
  }elsif( $rna =~ /(\w+)-bacterial/ ){
    $gene_name = "trn" . $trna{$1} . "-bacterial";
  }else{
    $gene_name = "trn" . $trna{$rna};
  }
  
  open( HTML,">$out_directory/$rna.html" ) || die "Couldn't open $rna.html: $!\n";

  print HTML "<html>\n<head>\n<title>$query_taxon: $rna</title>\n</head>\n";
  print HTML "<style>\n";
  print HTML "p1{\n";
  print HTML "\tfont-family: Courier, monospace;\n";
  print HTML "\twhite-space: nowrap;\n";
  print HTML "\tfont-style: normal;\n";
  print HTML "\tfont-variant: normal;\n";
  print HTML "\tfont-weight: normal;\n";
  print HTML "\tfont-size: 14;\n";
  print HTML "\tword-spacing: normal;\n";
  print HTML "\tletter-spacing: normal;\n";
  print HTML "\ttext-align: left;\n";
  print HTML "\tcolor: black ;\n";
  print HTML "\talign: left;\n";
  print HTML "}\n\n";
  print HTML "p2{\n";
  print HTML "\tfont-family: Courier, monospace;\n";
  print HTML "\twhite-space: nowrap;\n";
  print HTML "\tfont-style: normal;\n";
  print HTML "\tfont-variant: normal;\n";
  print HTML "\tfont-weight: normal;\n";
  print HTML "\tfont-size: 14;\n";
  print HTML "\tword-spacing: normal;\n";
  print HTML "\tletter-spacing: normal;\n";
  print HTML "\ttext-align: left;\n";
  print HTML "\tcolor: purple ;\n";
  print HTML "\talign: left;\n";
  print HTML "}\n\n";
  print HTML "h1{\n";
  print HTML "\tfont-family: Courier, monospace;\n";
  print HTML "\tfont-style: normal;\n";
  print HTML "\tfont-variant: normal;\n";
  print HTML "\tfont-weight: normal;\n";
  print HTML "\tfont-size: 24;\n";
  print HTML "\ttext-align: left;\n";
  print HTML "\tcolor: blue;\n";
  print HTML "\talign: left;\n";
  print HTML "}\n";
  
  print HTML "</style>\n";
  
  # @$gene_lengths_ref = sort {$a <=> $b}( @$gene_lengths_ref);
  # my $gl  =join( '&nbsp;', @$gene_lengths_ref);

  my $gene_lengths;
  my $min = min( @$gene_lengths_ref );
  my $max = max( @$gene_lengths_ref );
  if( $min == $max ){
    $gene_lengths = $min;
  }else{
    $gene_lengths = "$min\-$max";
  }

  print HTML "<table valign=\"middle\">\n\t<tr nowrap valign=\"middle\">\n";
  print HTML "\t\t<td>", '&nbsp;' x 85, "<\/td>\n";

  # for rRNAs
  if( $rna =~ /rrn/ ){
    # print HTML "<body onLoad=\"window.open(\'${webbase}rrna_make_annotation_form.cgi?gene_name=$gene_name&taxon=$query_taxon&project=$project\',\'mywindow\',\'width=550,height=600,resizable=yes,menubar=no,scrollbars=1,screenX=0,screenY=0\')\">\n"; # VCRU change, add path $webbase
    print HTML "\t\t<td><table><form METHOD=\"GET\"><INPUT TYPE=\"button\" VALUE=\"Annotate\" onClick=\"window.open(\'${webbase}rrna_make_annotation_form.cgi?gene_name=$gene_name&taxon=$query_taxon&project=$project\',\'Sequin_form\',\'width=550,height=600,resizable=yes,menubar=no\')\"><\/form><\/table><\/td>\n";  # VCRU change, add path $webbase
    print HTML "\t\t<td><font face=\"Courier\" color=\"blue\" size=\"6\">&nbsp;&nbsp;&nbsp;$gene_name&nbsp;($gene_lengths&nbsp;nt)<\/font><\/td>\n";
    print HTML "\t<\/tr>\n<\/table>\n<br>\n<hr size=4>";

  # for Met-f
  }elsif( $rna =~ /fM/i ){
    $rna = "Met";
    my $label = "Met-f";
    print HTML "\t\t<td><table><form METHOD=\"GET\"><INPUT TYPE=\"button\" VALUE=\"Annotate\" onClick=\"window.open(\'${webbase}trna_make_annotation_form.cgi?gene_name=$gene_name&product=tRNA-$rna&taxon=$query_taxon&project=$project\',\'Sequin_form\',\'width=550,height=550,resizable=yes,menubar=no\')\"><\/form><\/table><\/td>\n"; # VCRU change, add path $webbase
    print HTML "\t\t<td><font face=\"Courier\" color=\"blue\" size=\"6\">&nbsp;&nbsp;&nbsp;$gene_name&nbsp;($gene_lengths&nbsp;nt)<\/font><\/td>\n";
    print HTML "\t<\/tr>\n<\/table>\n<br>\n<hr size=4>";

  # for mitochondrial Cys
  }elsif( $rna eq "Cys-mt" ){
    $rna = "Cys";
    my $product = "tRNA-Cys";
    # print HTML "<body onLoad=\"window.open(\'${webbase}trnC_make_annotation_form.cgi?gene_name=$gene_name&product=$product&type=mito&taxon=$query_taxon&project=$project\',\'mywindow\',\'width=550,height=600,resizable=yes,menubar=no,scrollbars=1,screenX=0,screenY=0\')\">\n"; # VCRU change, add path $webbase
    print HTML "\t\t<td><table><form METHOD=\"GET\"><INPUT TYPE=\"button\" VALUE=\"Annotate\" onClick=\"window.open(\'${webbase}trnC_make_annotation_form.cgi?gene_name=$gene_name&product=$product&type=mito&taxon=$query_taxon&project=$project\',\'Sequin_form\',\'width=550,height=550,resizable=yes,menubar=no\')\"><\/form><\/table><\/td>\n"; # VCRU change, add path $webbase
    print HTML "\t\t<td><font face=\"Courier\" color=\"blue\" size=\"6\">&nbsp;&nbsp;&nbsp;$gene_name&nbsp;($gene_lengths&nbsp;nt)<\/font><\/td>\n";
    print HTML "\t<\/tr>\n<\/table>\n<br><br>\n";

  # for bacterial Cys
  }elsif( $rna =~ "Cys-bacterial" ){
    $rna = "Cys";
    my $product = "tRNA-Cys";
    # print HTML "<body onLoad=\"window.open(\'${webbase}trnC_make_annotation_form.cgi?gene_name=$gene_name&product=$product&type=bacterial&taxon=$query_taxon&project=$project\',\'mywindow\',\'width=550,height=600,resizable=yes,menubar=no,scrollbars=1,screenX=0,screenY=0\')\">\n"; # VCRU change, add path $webbase
    print HTML "\t\t<td><table><form METHOD=\"GET\"><INPUT TYPE=\"button\" VALUE=\"Annotate\" onClick=\"window.open(\'${webbase}trnC_make_annotation_form.cgi?gene_name=$gene_name&product=$product&type=bacterial&taxon=$query_taxon&project=$project\',\'Sequin_form\',\'width=550,height=550,resizable=yes,menubar=no\')\"><\/form><\/table><\/td>\n"; # VCRU change, add path $webbase
    print HTML "\t\t<td><font face=\"Courier\" color=\"blue\" size=\"6\">&nbsp;&nbsp;&nbsp;$gene_name&nbsp;($gene_lengths&nbsp;nt)<\/font><\/td>\n";
    print HTML "\t<\/tr>\n<\/table>\n<br><br>\n";

  # for chloroplast Cys
  }elsif( $rna =~ /Cys-cp/ ){
    $rna = "Cys";
    my $product = "tRNA-Cys";
    # print HTML "<body onLoad=\"window.open(\'${webbase}trnC_make_annotation_form.cgi?gene_name=$gene_name&product=$product&type=cplast&taxon=$query_taxon&project=$project\',\'mywindow\',\'width=550,height=600,resizable=yes,menubar=no,scrollbars=1,screenX=0,screenY=0\')\">\n"; # VCRU change, add path $webbase
    print HTML "\t\t<td><table><form METHOD=\"GET\"><INPUT TYPE=\"button\" VALUE=\"Annotate\" onClick=\"window.open(\'${webbase}trnC_make_annotation_form.cgi?gene_name=$gene_name&product=$product&type=cplast&taxon=$query_taxon&project=$project\',\'Sequin_form\',\'width=550,height=550,resizable=yes,menubar=no\')\"><\/form><\/table><\/td>\n"; # VCRU change, add path $webbase
    print HTML "\t\t<td><font face=\"Courier\" color=\"blue\" size=\"6\">&nbsp;&nbsp;&nbsp;$gene_name&nbsp;($gene_lengths&nbsp;nt)<\/font><\/td>\n";
    print HTML "\t<\/tr>\n<\/table>\n<br><br>\n";

  # for chloroplast tRNAs
  }elsif( $rna =~ /(\w+)-cp/ ){
    $rna = $1;
    my $product = "tRNA-" . $1;
    # print HTML "<body onLoad=\"window.open(\'${webbase}trna_make_annotation_form.cgi?gene_name=$gene_name&product=$product&note=chloroplast-like&taxon=$query_taxon&project=$project&cplast=yes\',\'mywindow\',\'width=550,height=600,resizable=yes,menubar=no,scrollbars=1,screenX=0,screenY=0\')\">\n"; # VCRU change, add path $webbase
    print HTML "\t\t<td><table><form METHOD=\"GET\"><INPUT TYPE=\"button\" VALUE=\"Annotate\" onClick=\"window.open(\'${webbase}trna_make_annotation_form.cgi?gene_name=$gene_name&product=$product&note=chloroplast-like&taxon=$query_taxon&project=$project&cplast=yes\',\'Sequin_form\',\'width=550,height=550,resizable=yes,menubar=no\')\"><\/form><\/table><\/td>\n"; # VCRU change, add path $webbase
    print HTML "\t\t<td><font face=\"Courier\" color=\"blue\" size=\"6\">&nbsp;&nbsp;&nbsp;$gene_name&nbsp;($gene_lengths&nbsp;nt)<\/font><\/td>\n";
    print HTML "\t<\/tr>\n<\/table>\n<br><br>\n";

  # for non-chloroplast, non-Cys mitochondrial tRNAs
  }else{
    print HTML "\t\t<td><table><form METHOD=\"GET\"><INPUT TYPE=\"button\" VALUE=\"Annotate\" onClick=\"window.open(\'${webbase}trna_make_annotation_form.cgi?gene_name=$gene_name&product=tRNA-$rna&taxon=$query_taxon&project=$project&cplast=no\',\'Sequin_form\',\'width=550,height=550,resizable=yes,menubar=no\')\"><\/form><\/table><\/td>\n"; # VCRU change, add path $webbase
    print HTML "\t\t<td><font face=\"Courier\" color=\"blue\" size=\"6\">&nbsp;&nbsp;&nbsp;$gene_name&nbsp;($gene_lengths&nbsp;nt)<\/font><\/td>\n";
    print HTML "\t<\/tr>\n<\/table>\n<br><br>\n";
  }
}
#######################################################################################################################
# Called by annotate_rna in Main; prints summary html page for the analyzed genome

sub print_rna_summary{
  my( $query_taxon, $project, $seqlen, $rna_html_files, $rna_blast_files, $tscan_html_files_Href, $out_directory, $missed_by_blast_Aref, $missed_Aref, $all_found_Aref, $trna_frames_Href, $tscan_html_paths ) = @_;
  my $gene;
  my $output_file = "$out_directory" . "/" . "$project" . "_rna_summary.html";
  my $n_missed_rows = @{$missed_Aref};
  my $frame;
  
  open( HTML_SUM,">$output_file" ) || die "Couldn't open $output_file: $!\n";
  
  # print header
  print HTML_SUM "<HTML>\n";
  print HTML_SUM "<TITLE>$query_taxon RNA results</TITLE>\n";
  print HTML_SUM "<HEAD>\n";
  print HTML_SUM "<body BGCOLOR=#C2C2C2>\n\n";
  # print HTML_SUM "<body BGCOLOR=#CCCC99>\n\n";

  # header table
  print HTML_SUM "<table cellspacing=6 cellpadding=3 valign=\"middle\">\n";
  print HTML_SUM "\t<tr valign=\"middle\">\n";
  print HTML_SUM "\t\t<td>", '&nbsp;' x 20, "<\/td>\n";
  print HTML_SUM "\t\t<td><font face=\"Arial\" color=\"black\" size=\"5\">RNA BLAST results for $query_taxon mtDNA<\/font><\/td>\n";
  print HTML_SUM "\t<\/tr>\n";
  print HTML_SUM "<\/table>\n\n<hr>\n\n";
  
  # begin results table
  print HTML_SUM "<table cellspacing=6 cellpadding=3 valign=\"middle\">\n";
  print HTML_SUM "\t<tr valign=\"middle\">\n";
  print HTML_SUM "\t\t<td>", '&nbsp;' x 20, "<\/td>\n";
  print HTML_SUM "\t\t<th align=\"left\"><font face=\"Arial\">Gene<\/font><\/th>\n";
  print HTML_SUM "\t\t<th align=\"center\"><font face=\"Arial\">Annotate<\/font><\/th>\n";
  print HTML_SUM "\t\t<th align=\"right\"><font face=\"Arial\">Raw BLAST output<\/font><\/th>\n";
  print HTML_SUM "\t\t<td>", '&nbsp;' x 20, "<\/td>\n";
  
  # missing genes portion of results table
  print HTML_SUM "\t\t<td valign=\"top\" rowspan=$n_missed_rows>\n";
  print HTML_SUM "\t\t\t<table valign=\"top\" valign=\"middle\" align=\"center\">\n";
  print HTML_SUM "\t\t\t\t<tr><th align=\"center\"><font face=\"Arial\">Missed RNAs<\/font><\/th><\/tr>\n";
  foreach( sort @{$missed_Aref} ){
    print HTML_SUM "\t\t\t<tr><td align=\"center\"><font face=\"Arial\">$_<\/font><\/td><\/tr>\n";
  }
  
  print HTML_SUM "\t\t\t<\/table>\n";
  print HTML_SUM "\t\t<\/td>\t\n";
  print HTML_SUM "\t<\/tr>\n";
  
  # print results, one row per gene found
  foreach $gene ( sort @{$all_found_Aref} ){
    $gene eq "Cys" and next;
    if( $gene eq "Met-f"){
      $gene = "trnfM";
      my $label = "Met-f";
      print HTML_SUM "\t<tr valign=\"middle\">\n";
      print HTML_SUM "\t\t<td>", '&nbsp;' x 20, "<\/td>\n";
      #print "gene" row header
      print HTML_SUM "\t\t<td align=\"left\"><font face=\"Arial\" color=\"purple\" size=\"4\">$label&nbsp;<\/font><\/td>\n";
      #print "Annotate" button
      print HTML_SUM "\t\t<td><table><tr align=\"left\"><form><INPUT TYPE=\"button\" VALUE=\"Annotate\" onClick=\"window.open(\'\.\.\/$trna_frames_Href->{Metf}\',\'BLAST_summary\')\"><\/form><\/table><\/td>\n";
      #print link to "Raw BLAST output"
      print HTML_SUM "\t\t<td align=\"right\"><a href=\"$rna_blast_files->{$label}\"><font face=\"Arial\">$label BLAST output<\/font><\/a><\/td>\n";
      print HTML_SUM "\t<\/tr>\n";	
    }elsif( $gene =~ /rrn/ ){
      print HTML_SUM "\t<tr valign=\"middle\">\n";
      print HTML_SUM "\t\t<td>", '&nbsp;' x 20, "<\/td>\n";
      #print "gene" row header
      print HTML_SUM "\t\t<td align=\"left\"><font face=\"Arial\" color=\"purple\" size=\"4\">$gene&nbsp;<\/font><\/td>\n";
      #print "Annotate" button
      print HTML_SUM "\t\t<td><table><tr align=\"left\"><form><INPUT TYPE=\"button\" VALUE=\"Annotate\" onClick=\"window.open(\'$rna_html_files->{$gene}\',\'BLAST_summary\')\"><\/form><\/table><\/td>\n";
      #print link to "Raw BLAST output"
      print HTML_SUM "\t\t<td align=\"right\"><a href=\"$rna_blast_files->{$gene}\"><font face=\"Arial\">$gene BLAST output<\/font><\/a><\/td>\n";
      print HTML_SUM "\t<\/tr>\n";	
    }else{
      print HTML_SUM "\t<tr valign=\"middle\">\n";
      print HTML_SUM "\t\t<td>", '&nbsp;' x 20, "<\/td>\n";
      #print "gene" row header
      print HTML_SUM "\t\t<td align=\"left\"><font face=\"Arial\" color=\"purple\" size=\"4\">$gene&nbsp;<\/font><\/td>\n";
      #print "Annotate" button
      print HTML_SUM "\t\t<td><table><tr align=\"left\"><form><INPUT TYPE=\"button\" VALUE=\"Annotate\" onClick=\"window.open(\'\.\.\/$trna_frames_Href->{$gene}\',\'BLAST_summary\')\"><\/form><\/table><\/td>\n";
      #print link to "Raw BLAST output", only if it exists; e.g., BLAST output might not exist for mitochondrial tRNA genes if there is no BLAST
      #database (e.g., Trp), but it is found by tRNA-scan, so there should be a link to Trp output in the summary page, but there should not be a
      #link to BLAST output, because it doesn't exist
      if( $rna_blast_files->{$gene} ){
	print HTML_SUM "\t\t<td align=\"right\"><a href=\"$rna_blast_files->{$gene}\"><font face=\"Arial\">$gene BLAST output<\/font><\/a><\/td>\n";
      }
      print HTML_SUM "\t<\/tr>\n";	
    }
  }
  print HTML_SUM "<\/table>\n<br>\n<hr>";
  
}
#######################################################################################################################
# called by annotate_rna

sub make_framesets{
  my( $project, $out_directory, $query_taxon, $rna_html_files_Href, $tscan_html_files_Href, $all_found_Aref, $trna_frames_Href, $tscan_html_paths, $rna_blast_files, $rnaNoHits, $trna_AA2letter ) = @_;
  my $AA;
  my $frame_file;
  my $blast = "BLAST";
  my $tscan = "tRNA-scan";

  foreach $AA ( @{$all_found_Aref} ){
    #make top (title) frame for tRNAs only
    $AA =~ /rrn/ and next;

    my $top_frame = make_top_frame( $query_taxon, $AA, $out_directory, $project );

    open( FRAME_TEMPLATE, $mitofybinbase."/frameset.html" ) || die "Couldn't open frameset.html: $!\n"; # VCRU change, add path $mitofybinbase
    $frame_file = "$out_directory/$AA" . "_rna_frame.html";
    
    if( $AA eq "Met-f" ){
      $trna_frames_Href->{"Metf"} = $project . "_out/$AA" . "_rna_frame.html";
      $top_frame =~ s/Metf/Met-f/;
    }else{
      $trna_frames_Href->{$AA} = $project . "_out/$AA" . "_rna_frame.html";
    }
    
    #open output file
    open( FRAME, ">$frame_file" ) || die "Couldn't open $frame_file: $!\n";

    while( <FRAME_TEMPLATE> ){
      if( /TOP/ ){
	$AA =~ /rrn/ and next;
	$_ =~ s/TOP/$top_frame/g;
	print FRAME $_;
      }elsif( /BLAST/ ){
	if( $rnaNoHits->{$AA} ){
	  $_ =~ s/BLAST/$rnaNoHits->{$AA}/g;
	  print FRAME $_;
	}elsif( $rna_html_files_Href->{$AA} ){
	  $_ =~ s/BLAST/$rna_html_files_Href->{$AA}/g;
	  print FRAME $_;
	}else{
	  my $not_found = not_found( $blast, $AA, $out_directory, $project, $trna_AA2letter );
	  $_ =~ s/BLAST/$not_found/g;
	  print FRAME $_;
	}

      }elsif( /TSCAN/ ){
	#tscan output irrelevant for rRNA genes
	$AA =~ /rrn/ and next;

	#if found by tscan (NON "*-cp, -mt, -bacterial" tRNAs)
	if( $tscan_html_files_Href->{$AA} ){
	  $_ =~ s/TSCAN/$tscan_html_paths->{$AA}/;
	  print FRAME $_ and next;

	}elsif( $AA eq "Met-f" && $tscan_html_files_Href->{Met} ){
	  $_ =~ s/TSCAN/$tscan_html_paths->{Met}/;
	  print FRAME $_ and next;

	#if found by tscan (for "-cp" tRNAs)
	}elsif( $AA =~ /(\w+)-cp/ && $tscan_html_files_Href->{$1} ){
	  my $a = $1;
	  if( $tscan_html_files_Href->{$a} ){
	    $_ =~ s/TSCAN/$tscan_html_paths->{$a}/;
	    print FRAME $_ and next;
	  }
	
	#if found by tscan (for Cys-mt)
	}elsif( $AA eq "Cys-mt" && $tscan_html_files_Href->{Cys} ){
	  if( $tscan_html_files_Href->{Cys} ){
	    $_ =~ s/TSCAN/$tscan_html_paths->{Cys}/;
	    print FRAME $_ and next;
	  }

	#if found by tscan (for Cys-bacterial)
	}elsif( $AA eq "Cys-bacterial" && $tscan_html_files_Href->{Cys} ){
	  if( $tscan_html_files_Href->{Cys} ){
	    $_ =~ s/TSCAN/$tscan_html_paths->{Cys}/;
	    print FRAME $_ and next;
	  }

        #NOT found by tscan
	}else{
	  my $not_found = not_found( $tscan, $AA, $out_directory, $project, $trna_AA2letter );
	  $_ =~ s/TSCAN/$not_found/g;
	  print FRAME $_ and next;
	}
      }else{
	print FRAME $_ and next;
      }
    }
    close FRAME_TEMPLATE;
    close FRAME;
  }
}
#######################################################################################################################
sub make_top_frame{
  my( $query_taxon, $AA, $out_directory, $project ) = @_;
  my $top_frame = "$out_directory/$AA" . "_top_frame.html";
  my $label;
  
  if( $AA eq "Met-f" ){
    $AA = "Metf";
  }elsif( $AA eq "Cys-mt" ){
    $AA = "Cys";
  }elsif( $AA eq "Cys-bacterial" ){
    $AA = "Cys";
  }
  
  my %rrna_names = (
		    rrn5  => "5S ribosomal RNA",
		    rrn16 => "16S ribosomal RNA",
		    rrnS  => "18S ribosomal RNA",
		    rrnL  => "26S ribosomal RNA"
		   );

  my %trna_letters = (
		      Phe => "F",
		      Leu => "L",
		      Ile => "I",
		      Met => "M",
		      Val => "V",
		      Ser => "S",
		      Pro => "P",
		      Thr => "T",
		      Ala => "A",
		      Tyr => "Y",
		      His => "H",
		      Gln => "Q",
		      Asn => "N",
		      Lys => "K",
		      Asp => "D",
		      Glu => "E",
		      Cys => "C",
		      Trp => "W",
		      Arg => "R",
		      Gly => "G",
		      Metf => "M",
		      Undet => "UNDETERMINED tRNA"
		     );
  
  my %AA_names = (
		  Phe => "Phenylalanine",
		  Leu => "Leucine",
		  Ile => "Isoleucine",
		  Met => "Methionine",
		  Val => "Valine",
		  Ser => "Serine",
		  Pro => "Proline",
		  Thr => "Threonine",
		  Ala => "Alanine",
		  Tyr => "Tyrosine",
		  His => "Histidine",
		  Gln => "Glutamine",
		  Asn => "Asparagine",
		  Lys => "Lysine",
		  Asp => "Aspartic&nbsp;acid",
		  Glu => "Glutamic&nbsp;acid",
		  Cys => "Cysteine",
		  Trp => "Tryptophan",
		  Arg => "Arginine",
		  Gly => "Glycine",
		  Metf => "Methionine",
		 );

  #for mitochondrial rRNAs
  if( $AA =~ /rrn/ ){
    $label = "$AA" . "&nbsp;" . "(" . "$rrna_names{$AA}" . ")";

  #for mitochondrial Metf
  }elsif( $AA eq "Metf" ){ 
    $label = "trnfM" . "&nbsp;" . "(" . "Methionine" . ",&nbsp;Met-f". ")";

  #for chloroplast tRNAs
  }elsif( $AA =~ /(\w+)-cp/ && ($AA !~ /rrn/) ){
    $label = "trn" . "$trna_letters{$1}" . "&nbsp;" . "(" . "$AA_names{$1}" . ",&nbsp;$AA". ")";
    
  #for "Undetermined" tRNAs found by tRNAscan
  }elsif( $AA eq "Undet" ){
    $label = $trna_letters{$AA};

  #for trnC-mt
  }elsif( $AA eq "trnC-mt" ){
    $label = $trna_letters{$AA};

  #for trnC-bacterial
  }elsif( $AA eq "trnC-bacterial" ){
    $label = $trna_letters{$AA};

  #for non-chloroplast, mitochondrial tRNAs
  }else{
    $label = "trn" . "$trna_letters{$AA}" . "&nbsp;" . "(" . "$AA_names{$AA}" . ",&nbsp;$AA". ")";
  }

  open( TOP_FRAME, ">$top_frame" ) || die "Couldn't open $top_frame: $!\n";
  print TOP_FRAME "<html>\n<head>\n<title>$query_taxon: $AA</title>\n</head>\n";
  print TOP_FRAME "<body BGCOLOR=#999999>\n\n";
  
  print TOP_FRAME "<table valign=\"middle\">\n\t<tr nowrap valign=\"middle\">\n";
  print TOP_FRAME "\t\t<td>", '&nbsp;' x 101, "<\/td>\n";
  print TOP_FRAME "\t\t<td><table><\/table><\/td>\n";
  print TOP_FRAME "\t\t<td><font face=\"Courier\" color=\"blue\" size=\"6\">&nbsp;&nbsp;&nbsp;$label<\/font><\/td>\n";
  print TOP_FRAME "\t<\/tr>\n<\/table>\n<hr size=4>";
  
  close TOP_FRAME;
  
  #change $top_frame for use in make_framesets
  $top_frame = "\.\./$project" . "_out/$AA" . "_top_frame.html";

  return $top_frame;
}
#######################################################################################################################
# called by mitofy.pl and annotate_rna

sub not_found{
  my( $nf, $AA, $out_directory, $project, $trna_AA2letter ) = @_;
  
  if( $trna_AA2letter ){ #called by annotate_rna, so it's a missing RNA gene -- determine RNA type
    if( $AA !~ /rrn/ ){
      if( $AA eq "Metf" || $AA eq "Met-f" ){ 
	$AA = "trnfM";

      #for chloroplast tRNAs
      }elsif( $AA =~ /(\w+)-cp/ ){
	$AA = "trn" . $trna_AA2letter->{$1} . "-cp";

      #for Cys-mt
      }elsif( $AA =~ /(\w+)-mt/ ){
	$AA = "trn" . $trna_AA2letter->{Cys} . "-mt";

      #for Cys-bacterial
      }elsif( $AA =~ /(\w+)-bacterial/ ){
	$AA = "trn" . $trna_AA2letter->{Cys} . "-bacterial";

      #for UNDETERMINED tRNAs
      }elsif( $AA eq "Undet" ){
	$AA = $trna_AA2letter->{$AA};

      #for non-chloroplast, mitochondrial tRNAs
      }else{
	$AA = "trn" . $trna_AA2letter->{$AA};
      }
    }
  }else{ #called by mitofy.pl, so it's a missing protein gene
  }

  my $not_found = "$AA not found by $nf.";
  my $not_found_page = "$out_directory/" . $AA . "_$nf" . "_NotFound.html";
  open( NOT_FOUND, ">$not_found_page" ) || die "Couldn't open $not_found_page: $!\n";

  print NOT_FOUND "<html>\n";
  print NOT_FOUND "<head></head>\n";
  print NOT_FOUND "\t<table valign=\"middle\">\n\t<tr nowrap valign=\"middle\">\n";
  print NOT_FOUND "\t\t<td>" . '&nbsp;' x 101;
  print NOT_FOUND "<\/td>\n";
  print NOT_FOUND "\t\t<td><table><\/table><\/td><br><br><br><br>\n";
  print NOT_FOUND "\t\t<td width=520 align=center><font face=\"Courier\" color=\"blue\" size=\"5\">$not_found<\/font><\/td>\n";
  print NOT_FOUND "\t<\/tr>\n<\/table>\n";

  close NOT_FOUND;
  
  # e.g., for Trp, which is found by tRNA-scan but for which there is no BLAST database and therefore
  # no file of raw blast output; this is the best place to assign $rna_blast_files->{$AA} to an 
  # actual file, the "Not found" page
  $not_found_page = "../$project" . "_out/" . $AA . "_$nf" . "_NotFound.html";
  return( $not_found_page );
}
#######################################################################################################################
# called by annotate_rna

sub make_rna_no_hits{
  my( $rna, $rnaNoHits, $out_directory, $project, $trna_AA2letter ) = @_;
  
  if( $rna !~ /rrn/ ){
    if( $rna eq "Metf" || $rna eq "Met-f" ){ 
      $rna = "trnfM";
    # for chloroplast tRNAs
    }elsif( $rna =~ /(\w+)-cp/ ){
      $rna = "trn" . $trna_AA2letter->{$1} . "-cp";
     # for Cys-mt
    }elsif( $rna eq "Cys-mt" ){
      $rna = "trn" . $trna_AA2letter->{Cys} . "-mt";
     # for Cys-bacterial
    }elsif( $rna eq "Cys-bacterial" ){
      $rna = "trn" . $trna_AA2letter->{Cys} . "-bacterial";
    # for non-chloroplast, non-Cys mitochondrial tRNAs
    }else{
      $rna = "trn" . $trna_AA2letter->{$rna};
    }
  }
  
  my $no_hits = "No significant BLAST hits for $rna.";
  my $no_hits_page = "$out_directory/" . $rna . "_NoHits.html";
  $rnaNoHits->{$rna} = "\.\.\/$project" . "_out/" . $rna . "_NoHits.html";
  
  open( NO_HITS, ">$no_hits_page" ) || die "Couldn't open $no_hits_page: $!\n";
  print NO_HITS "<html>\n";
  print NO_HITS "<head></head>\n";
  print NO_HITS "\t<table valign=\"middle\">\n\t<tr nowrap valign=\"middle\">\n";
  print NO_HITS "\t\t<td>" . '&nbsp;' x 101;
  print NO_HITS "<\/td>\n";
  print NO_HITS "\t\t<td><table><\/table><\/td><br><br><br><br>\n";
  print NO_HITS "\t\t<td width=520 align=center><font face=\"Courier\" color=\"blue\" size=\"5\">$no_hits<\/font><\/td>\n";
  print NO_HITS "\t<\/tr>\n<\/table>\n";
  close NO_HITS;
  
}
#######################################################################################################################

1;
