# mitofy

A fork of mitofy for linux web interface.

The original mitofy was written by Andy Alverson and can be found at
https://dogma.ccbb.utexas.edu/mitofy/

Citation: Alverson et al. 2010. Insights into the evolution of mitochondrial genome size
from complete sequences of Citrullus lanatus and Cucurbita pepo (Cucurbitaceae).
Mol Biol Evol 27, 1436-1448.

## Purpose

The original version was written for MacOS. The purpose of this project is to
create a linux version, that can be run on a public web server, allowing users
to upload their FASTA sequence file, and analyze without the need to install
software on their own computer.

This program is currently available for public analyses at
https://www.vcru.wisc.edu/cgi-bin/mitofy/mitofy.cgi

## Dependencies

The original mitofy came with MacOS versions of tRNAscan-SE and blast+.
These are not included in this project, and will need to be installed.
Dependencies are:

* tRNAscan-SE - http://lowelab.ucsc.edu/tRNAscan-SE/
* NCBI blast+ - http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download
* zip - should be included with any linux distribution
* apache2 web server

## Installation

1. git clone https://github.com/dsenalik/mitofy.git
1. Edit this section of mitofy/mitofy.cgi to configure paths to dependencies
```# configurable directories
# these will be referenced by modified verisions of mitofy.pl and the programs that mitofy.pl calls
# customize these to appropriate locations, and changes will then propagate to all called programs
my $tmpbase = "/var/www/html/tmp";      # cgi-accessible temporary file directory
my $tmpbaseurl = "/tmp";                # url of above directory
$ENV{mitofybinbase} = "/var/www/cgi-bin/mitofy/annotate";   # the location of mitofy.pl and other called programs
$ENV{webbase} = "/cgi-bin/mitofy/cgi/"; # the URL of the mitofy cgi directory
$ENV{blastplusdir} = "/usr/bin";        # where the NCBI blast+ binaries are located
$ENV{trnascandir} = "/usr/local/bin";   # where tRNAscan-SE program is located
my $zipbinary = "/usr/bin/zip";         # full path to program for creating archive of all files
```
