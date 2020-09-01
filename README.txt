Douglas Senalik dsenalik@wisc.edu June 29, 2013

Included files are my modifications to MITOFY to allow use on a web server.
The MITOFY home page is http://dogma.ccbb.utexas.edu/mitofy/
Every line in MITOFY that I modified is marked with a # VCRU tag.

To use on a linux system, you will have to first install the original MITOFY,
and alsoinstall tRNAscan-SE, the NCBI Blast+ programs, and Apache2.
The versions of tRNAscan-SE and blast included with MITOFY do not work as
they seem to be MacOS versions.
My notes for installing tRNAscan-SE can be found at
http://vcru.wisc.edu/simonlab/bioinformatics/programs/install/trnascan-se.htm
Your Apache2 installation will need a directory with suitable permissions
to execute cgi programs. Install MITOFY into that directory. Replace the
following files with the files from this archive:
  annotate/mitofy.pl
  annotate/annotate_subroutines.pl
  annotate/rna_html_pages.pl
  annotate/parse_trnascan.pl
  annotate/annotate_rna.pl
  cgi/cis_template.html
  cgi/rrna_template.html
  cgi/trans_template.html
  cgi/trna_template.html
  cgi/trnC_template.html
  cgi/cis_SEQUIN_OUT.pl
  cgi/rrna_SEQUIN_OUT.pl
  cgi/trans_SEQUIN_OUT.pl
  cgi/trna_SEQUIN_OUT.pl
  cgi/trnC_SEQUIN_OUT.pl
The main cgi program is mitofy.cgi

You will have to modify the directories at the top of the mitofy.cgi program
to reflect your install locations, and locations of tRNAscan-SE and blast+
