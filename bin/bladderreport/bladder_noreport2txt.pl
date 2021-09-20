#!/usr/bin/perl -w

use strict;
use Data::Dumper;
use File::Basename;

# script to create a tab separated txt file from html file.

my $html_path = $ARGV[0];
my %data;

my @label=( "Prov-ID",
         "Personnummer",
         "Patientnamn",
         "RNA-kvalitet",
         "Koncentration",
         "Ankomstdatum",
         "Rapportdatum",
         "Indatafil",
         "Expressionsprofilen tyder på att tumören är av subtypskategorin",
         "Undergrupp av \“Urotel-lik\” \(om bedömbart\)",
         "Molekylär grad bedöms som",
         "RIN",
         "DV200",
         "Unikt inpassade läsningar (M)",
         "Unikt inpassade läsningar (%)",
         "Multiinpassade läsningar (%)",
         "Missmatchratio per bas (%)",
         "Kanoniska splitsningar \(n\)",
         "Icke-kanoniska splitsningar \(n\)"
);

my @html_content = `lynx -dump -nolist -width 120 $html_path`; 

for my $line (@html_content){
    chomp( $line );
    if( $line =~/(Patientnamn|Personnummer|Prov-ID|Ankomstdatum|Rapportdatum|Indatafil|Ingångsmaterial|Koncentration|RNA-kvalitet|
                  Expressionsprofilen\styder\spå\satt\stumören\sär\sav\ssubtypskategorin|Undergrupp\sav\s\“Urotel-lik\”\s+\(om\s+bedömbart\)|
                  Molekylär\sgrad\sbedöms\ssom|RIN|DV200|Unikt\sinpassade\släsningar\s\(M\)|Unikt\sinpassade\släsningar\s\(\%\)|
                  Multiinpassade\släsningar\s+\(\%\)|Missmatchratio\sper\sbas\s+\(\%\)|Kanoniska\ssplitsningar\s+\(n\)|
                  Icke-kanoniska\ssplitsningar\s+\(n\))[:\s+]*(.*)/ ){
        $data{ $1 } = $2;
    }
}

$data{ 'RNA-kvalitet' } = "NA";

if( defined($data{ 'Indatafil' }) ){ $data{ 'Indatafil' } = '/srv/data/pontus_e/rsem_bladder_reports/'.basename( $data{ 'Indatafil' } ); }

print "Labb-ID\tPersonnummer\tPatientens Namn\tRNA-kvalité\tKoncentration\tAnkomstdatum\tRapportdatum\tIndatafil\tSubtyp\tUndergrupp\t";
print "Grad\tRIN\tDV200\tUnikt inpassade läsningar \(M\)\tUnikt inpassade läsningar \(\%\)\tMultiinpassade läsningar \(\%\)\tMissmatchratio per bas \(\%\)\tKanoniska splitsningar \(n\)\tIcke-kanoniska splitsningar \(n\)\n"; 
print "$data{ $label[0] }"; 
for my $index( 1..$#label ){  
    if( $data{ $label[$index] } ){
        if( $data{ $label[$index] }=~/\s{2,20}(.*)/ ){ 
            print "\t$1"; 
        }else{ 
            print "\t$data{ $label[$index] }"; 
        }
    }else{ 
        print "\t-"; 
    } 
}  
print "\n"; 
