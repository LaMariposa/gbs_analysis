#!/usr/bin/perl

#make_barcode_file.pl by M. Supple
#created 22 April 2015
#last modified 22 April 2015

#script to generate a barcode file for axe demultiplexing from sample info spread sheets
#usage


use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

my @plate_info;

#die "$usage" unless (@ARGV == 4)

#read in command line arguments
my ($plateName, $plateCSV, $sampleCSV, $barcodePlate)=@ARGV;


#get info on barcode plate from plates.csv
open(PLATE, $plateCSV) || die "can't open file with plate information. $!\n";
while (my $line = <PLATE>)
	{
	  chomp $line;
	  my @fields = split ",", $line;
	  if ($fields[0] eq $plateName){@plate_info=@fields}
	}


#open output file

#print header to output file
print ($plate_info[1]);

#read barcode plate layout into hash

#extract sample info from the samples.csv

#for each sample->get SampleID, PlateCoordinates (barcode1 and barcode2)
#print to output file (barcode1<tab>barcode2<tab>SampleID)

