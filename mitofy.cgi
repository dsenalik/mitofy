#!/usr/bin/perl -wT
# This is a cgi wrapper for the Mitofy program (version 2012-) from http://dogma.ccbb.utexas.edu/mitofy/
$| = 1;  # autoflush, disable output buffering
use strict;
use CGI;
use CGI::Carp qw ( fatalsToBrowser );
use POSIX ":sys_wait_h";  # for forking

############################################################
# Configuration variables
############################################################
$CGI::POST_MAX = 1024 * 10000; # 10 megabyte maximum post size

$ENV{PATH}="";  # otherwise we get error: Insecure $ENV{PATH} while running with -T switch

# configurable directories
# these will be referenced by modified verisions of mitofy.pl and the programs that mitofy.pl calls
# customize these to appropriate locations, and changes will then propagate to all called programs
my $tmpbase = "/var/www/html/tmp";      # cgi-accessible temporary file directory
my $tmpbaseurl = "/tmp";                # url of above directory
$ENV{mitofybinbase} = "/var/www/cgi-bin/mitofy/annotate";   # the location of mitofy.pl and other called programs
$ENV{webbase} = "/cgi-bin/mitofy/cgi/"; # the URL of the mitofy cgi directory
$ENV{blastplusdir} = "/usr/bin";        # where the NCBI blast+ binaries are located
$ENV{trnascandir} = "/usr/local/bin";   # where tRNAscan-SE program is located
my $zipbinary = "/usr/bin/zip";         # full path to program for creating archive of all files

my $webtitle = "MITOFY Analysis Server";
my $headerinclude = "";  # leave undefined or null string to omit header
my $footerinclude = "";  # leave undefined or null string to omit footer



############################################################
# global variables
############################################################
(my $prognopath = $0) =~ s|^.*[\/\\]||;  # this program's name with the path removed
(my $thiscgiweb = $0) =~ s|^.*(\/cgi-bin.*)$|$1|;  # url of this cgi as seen from web
my $mitofybin = $ENV{"mitofybinbase"} . "/mitofy.pl";
my $sleepinterval = 20; # seconds, status update interval while waiting for mitofy.pl to finish



############################################################
# read CGI passed values
############################################################
my $query = new CGI;
my $projectid   = $query->param("projectid") // '';  $projectid =~ s/\s/_/g; $projectid = HD_untaint ( $projectid );
my $prot_emax   = $query->param("prot_emax") // '';  $prot_emax = HD_untaint ( $prot_emax );
my $prot_pmin   = $query->param("prot_pmin") // '';  $prot_pmin = HD_untaint ( $prot_pmin );
my $rna_emax    = $query->param("rna_emax") // '';   $rna_emax  = HD_untaint ( $rna_emax );
my $rna_pmin    = $query->param("rna_pmin") // '';   $rna_pmin  = HD_untaint ( $rna_pmin );
my $rna_mlen    = $query->param("rna_mlen") // '';   $rna_mlen  = HD_untaint ( $rna_mlen );
my $data        = $query->param("data") // '';
my $datahandle  = $query->upload("datafile") // '';
my $datafile = '';
if ( $datahandle ) { $datafile   = $query->tmpFileName($datahandle); }
if ( $datafile ) { $datafile = HD_untaint ( $datafile, ' \/' ); }
my $debug       = $query->param("debug");
my $submit      = $query->param("submit");
if ( ! $debug ) { $debug = 0; }  # to avoid warning message
if ( ! $projectid ) { $projectid = time(); }

my $projectdir  = "$tmpbase/${projectid}_out";
$ENV{outpath}   = $projectdir;   # where all mitofy.pl output files will be saved
$ENV{TMPDIR}    = $projectdir;  # no trailing slash - this is used by tRNAscan-SE
my $tmpfile1    = $projectdir . "/$projectid.fasta";  # input sequence file
my $tmpfile2    = $projectdir . "/mitofy.log";  # log of mitofy.pl output
my $tmpfile3    = "$tmpbase/$projectid.zip";  # compressed archive of all output files
my $tmpfile3url = "$tmpbaseurl/$projectid.zip";
my $tmpfile4    = "$tmpbase/$projectid.ziplog.txt";  # log of zip file creation
my $tmpfile4url = "$tmpbaseurl/$projectid.ziplog.txt";



############################################################
# initialize web page
############################################################
print "Content-type: text/html\n\n";
print "<HTML>\n<HEAD>\n<TITLE>$webtitle</TITLE>\n</HEAD>\n<BODY>\n";

include ( $headerinclude );
print "<H2>Public MITOFY Analysis Web Server</H2>\n";
print "<p>The MITOFY home page is <a href=\"http://dogma.ccbb.utexas.edu/mitofy/\">http://dogma.ccbb.utexas.edu/mitofy/</a><br>
If you publish an analysis using MITOFY, please cite the publication listed on the
web site shown above. The algorithms of MITOFY have not been modified from the
authors' version, however the versions of tRNAscan-SE and blast+ are different
than the MacOS versions supplied in the original program.<br>
This server uses the March 22, 2012 version of <a href=\"http://dogma.ccbb.utexas.edu/mitofy/\">MITOFY</a>,
version 1.3.1 of <a href=\"http://lowelab.ucsc.edu/tRNAscan-SE/\">tRNAscan-SE</a>,
and version 2.6.0 of the <a href=\"http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download\">NCBI blast+ programs</a> (updated from version 2.2.28 on Aug. 24, 2017)<br>
This cgi program, and MITOFY files modified to allow use as a web-based cgi can be downloaded as
<strike><a href=\"/simonlab/sdata/software/vcrumitofycgi20130701.tar.gz\">vcrumitofycgi20130701.tar.gz</a></strike>
<a href=\"/simonlab/sdata/software/vcrumitofycgi20170824.tar.gz\">vcrumitofycgi20170824.tar.gz</a><br>
Contact information: <a href=\"mailto:dsenalik\@wisc.edu\">dsenalik\@wisc.edu</a></p>\n";

print "<hr size=\"1\">\n";


############################################################
# either: process data
############################################################
if ( $submit )
  {



############################################################
# get input data stream, store in @infiledata
############################################################
  my @infiledata = ();
  my @outfiledata = ();
  if ( ! $datahandle )
    {
      debugmsg ( "Parsing pasted data" );
      @infiledata = split ( /[\n\r]+/, $data );  # side effect is to delete blank lines, but this is ok
      undef $data;
    }
  else
    {
      debugmsg ( "Parsing temporary input file from uploaded file handle \"$datahandle\"" );
      # some extra filtering to handle any possible combination of \n or \r for line breaks
      my $tmp = join ( '', <$datahandle> );
      @infiledata = split ( /[\n\r]+/, $tmp );  # side effect is to delete blank lines, but this is ok
      close $datahandle;
    } # else
  debugmsg ( scalar @infiledata . " lines of input data" );



############################################################
# save input FASTA data to a temporary file
############################################################
  {
  # create output temporary directory
  mkdir ( $projectdir ) or HD_die ( "Error creating temporary directory \"$projectdir\": $!" );

  # NOTE - mitofy will not accept a FASTA header if it contains any spaces!
  open my $OUTF,">",$tmpfile1 or HD_die ( "Could not create temporary FASTA file \"$tmpfile1\": $!" );
  my $nlines = 0;
  my $nlinessaved = 0;
  my $headerlines = 0;
  my $totalbases = 0;
  foreach my $aline ( @infiledata )
    {
      $nlines++;
      next if ( $aline =~ m/^\s*$/ );  # blank line removal
      if ( $aline =~ m/^\s*>/ )
        {
          unless ( $aline =~ m/^>\S/ )
            { HD_die ( "Invalid FASTA header line $nlines: \"$aline\"\n" ); }
          $aline =~ s/^\s*(>[^\s]*)\s.*$/$1/ ;  # remove first white space and any subsequent text in headers to be compatible with MITOFY
          $headerlines++;
        }
      else  # check for only valid characters
        {
          # check that first non-blank line must be FASTA header
          if ( ! $headerlines )
            { HD_die ( "First line is not a FASTA header starting with \">\"\n" ); }

          $aline =~ s/[\s\-]//g;  # remove all white space and "-" gap characters in sequence
          if ( $aline =~ m/([^AaCcTtGgMmRrWwSsYyKkVvHhDdBbNn]+)/ )
            { HD_die ( "Invalid character \"$1\" in line $nlines, not a nucleotide: \"$aline\"\n" ); }
          $totalbases += length($aline);
        }
      print $OUTF $aline, "\n";
      $nlinessaved++;
    }
  close $OUTF;
  debugmsg ( "Created temporary input FASTA file \"$tmpfile1\" with $headerlines sequences and $totalbases total bases, $nlinessaved lines" );
  unless ( $totalbases )
    { HD_die ( "Error, no nucleotides found in input data, no sequence supplied\n" ); }
  }



############################################################
# run MITOFY as child forked process
############################################################
  my $starttime = time;
  my $pid = fork();
  unless ( defined $pid ) { HD_die ( "unable to fork: $!\n" ); }
  if ( $pid == 0 )
    {
      # child process
      my $cmd = $mitofybin;
      $cmd .= " --prot_emax=\"$prot_emax\"";
      $cmd .= " --prot_pmin=\"$prot_pmin\"";
      $cmd .= " --rna_emax=\"$rna_emax\"";
      $cmd .= " --rna_pmin=\"$rna_pmin\"";
      $cmd .= " --rna_mlen=\"$rna_mlen\"";
      $cmd .= " \"$tmpfile1\"";
      $cmd .= " \"$projectid\"";
      $cmd .= " &> \"$tmpfile2\"";
      debugmsg ( "Running command \"$cmd\"" );
      my $result = system ( $cmd );
      if ( $result ) { HD_die ( "Error $result running command \"$cmd\"" ); }
      exit ( 0 );
    }
  else # parent
    {
  


############################################################
# wait for run to finish
############################################################
      print "<p>Run started, waiting for results (this could take an hour or more for large genomes)</p>\n";
      my $kid = 0;
      while ( $kid == 0 )
        {
          my $etime = time() - $starttime;
          print "Analysis is running, elapsed time ".timestr( $etime )."<br>\n";
          sleep $sleepinterval;
          $kid = waitpid($pid, WNOHANG);
        }



############################################################
# print MITOFY log
############################################################
      print "<hr size=\"1\">\n";
      print "<h3>Mitofy Results:</h3>\n";
      if ( -e $tmpfile2 )
        {
          open my $INF,"<",$tmpfile2 or HD_die ( "Error opening \"$tmpfile2\": $!" );
          while ( my $aline = <$INF> )
            { print $aline; }
          close $INF;
        }
      else
        {
          print "<p><b>Error, no output file was generated!</b></p>\n";
          exit ( 0 );
        }

      print "<p>Annotations, <i>after you complete them</i>, will be saved and available in this file:<br>\n";
      print "<a href=\"/tmp/${projectid}_out/${projectid}.tbl\" target=\"_blank\"><b>${projectid}.tbl</b></a></p>\n";


############################################################
# compress output files into a single archive
############################################################
      {
      my $cmd = "$zipbinary -9 -j -r \"$tmpfile3\" \"$projectdir\" &> \"$tmpfile4\"";  # -j to store files without path
      debugmsg ( "Running command \"$cmd\"" );
      my $result = system ( $cmd );
      if ( $result )
        { print ( "<p>Error $result creating archive of all output files with command \"$cmd\"<br>See <a href=\"$tmpfile4url\">log file</a></p>\n" ); }
      else
        {
          print "<p>A compressed archive of all output files can be downloaded as <a href=\"$tmpfile3url\"><b>$projectid.zip</b></a></p>\n";
          unlink ( $tmpfile4 );
        }
      }


    } # else parent
  } # if ( $submit )



else
{



############################################################
# or: print input form
############################################################
print "<form name=\"scheda\" method=\"post\" action=\"$thiscgiweb\" enctype=\"multipart/form-data\">\n";

print "
<table border=\"1\" cellpadding=\"3\" cellspacing=\"0\"><tr><td bgcolor=#E8FFE8>
&nbsp;Enter FASTA file for analysis by <a href=\"http://dogma.ccbb.utexas.edu/mitofy/\">MITOFY</a>:<br>

<textarea name=\"data\" rows=\"15\" cols=\"82\" wrap=\"off\"></textarea></td></tr>

<tr><td>or upload a file: <input type=\"file\" name=\"datafile\"/ size=\"48\"></td></tr>

";#<tr><td>or use this Linux file: <input type=\"text\" name=\"datalinux\"/ size=\"48\"></td></tr>
print "

</td></tr><tr><td bgcolor=#E8E8FF>
<input type=\"text\" name=\"projectid\"/ size=\"12\" value=\"$projectid\">
Project ID (required)
</td></tr><tr><td bgcolor=#E8FFE8>
<input type=\"text\" name=\"prot_emax\"/ size=\"4\" value=\"1e-3\">
--prot_emax - maximum BLAST expect value for protein genes
</td></tr><tr><td bgcolor=#E8E8FF>
<input type=\"text\" name=\"prot_pmin\"/ size=\"4\" value=\"60\">
--prot_pmin - minimum percent identity for protein genes
</td></tr><tr><td bgcolor=#E8FFE8>
<input type=\"text\" name=\"rna_emax\"/ size=\"4\" value=\"1e-3\">
--rna_emax - maximum BLAST expect value for RNA genes
</td></tr><tr><td bgcolor=#E8E8FF>
<input type=\"text\" name=\"rna_pmin\"/ size=\"4\" value=\"70\">
--rna_pmin - minimum percent identity for RNA genes
</td></tr><tr><td bgcolor=#E8FFE8>
<input type=\"text\" name=\"rna_mlen\"/ size=\"4\" value=\"30\">
--rna_mlen - minimum length for RNA BLASTN matches
</td></tr><tr><td>

<input type=\"submit\" name=\"submit\" value=\"Submit\" />&nbsp;&nbsp;&nbsp;&nbsp;<input type=\"reset\" value=\"Reset\">

</td></tr></table>

";


print "</form>\n";



} # else ( ! $submit )



print "<hr size=\"1\">\n";
print "<p><a href=\"http://www.perl.org/\"><img src=\"/images/perlpowered.png\" border=\"0\"></a></p>\n";
include ( $footerinclude );
print "</BODY>\n</HTML>\n";

exit 0;



#//////////////////////////////////////////////////////////#
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#
#//////////////////////////////////////////////////////////#



############################################################
sub trim { my ( $text ) = @_;
############################################################
  $text =~ s/^[\s]+//;
  $text =~ s/[\s]+$//;
  return $text;
} # sub trim



############################################################
sub debugmsg { my ( $msg, $nonewline ) = @_;
############################################################
# print message only if debuggin mode active
  if ( $debug )
    { 
      print $msg; 
      if ( ! $nonewline ) { print "<br>\n"; }
    }
} # sub debugmsg



############################################################
sub errormsg { my ( $errtext ) = @_;
############################################################
# print HTML paragraph of supplied text in bold + red
  print "<p><font color=\"RED\"><B>$errtext</B></font></p>\n";
} # sub errormsg



###############################################################
sub timestr
###############################################################
  {
    my ( $atime ) = @_;
    my $sec = $atime % 60;
    $atime -= $sec;
    my $min = $atime % 3600;
    $atime -= $min;
    @_ = localtime( $atime );
    return(sprintf("%02d:%02d:%02d", $atime / 3600, $min / 60, $sec ));
  } # sub timestr



############################################################
sub HD_die { my ( $message ) = @_;
############################################################
# one parameter is required, a text string containing an error message to display
# unlike die, this prints to STDOUT
    if ( ! $message ) { $message = "undefined error"; }
    # get script name without path
    ( my $prognopath = $0 ) =~ s/^.*\///;
    print "Error in script \"$prognopath\": " . $message . "\n";
    exit 1;
  } # HD_die



############################################################
sub include { my ( $filename ) = @_;
############################################################
  if ( open my $HDRFILE,"<",$filename ) # no error if not present
    {
      print <$HDRFILE>;
      close $HDRFILE;
    }
} # include



############################################################
sub HD_untaint { my ( $filename, $allowchars ) = @_;
############################################################
# untaint method based on http://perldoc.perl.org/perlsec.html#Laundering-and-Detecting-Tainted-Data
# untaints passed string if contains only safe characters, or returns a
# null string if illegal characters are present
# safe characters are alphabetics, numerics, underscore, a hyphen, an at sign, or a dot.
# to allow other characters, add to $allowchars in a regex compatible form,
# e.g. to allow a slash put '\/' in $allowchars
# and you may wish to add space ' ' too
  if ( ! defined $allowchars ) { $allowchars = ''; } $allowchars .= '';  # this line is to prevent a warning in apache log
  if ( $filename )
    {
      if ( $filename =~ m/^([\-\@\w\.$allowchars]+)$/ )
        { $filename = $1; } # untaints $filename
      else
        { $filename = ''; }
    } # if ( $filename )
  return $filename;
} # sub HD_untaint



# end of file
