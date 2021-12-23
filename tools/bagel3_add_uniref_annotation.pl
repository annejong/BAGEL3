#!/usr/bin/env perl

#  BAGEL3 program 
# 	Anne de Jong
#	University of Groningen
#	the Netherlands
#	anne.de.jong@rug.nl
#
#  2012-october
#
#  Add uniref annotation
# 

use strict ;
use warnings ;
use lib "/usr/molgentools/genomics";
use genomics ;

my $sessiondir = ".";
my $queryfile ;
my $outputfile = "uniref_added.txt";
my $unirefdb = "/var/uniref50/uniReference.db"; # file created by uniprot_make_reftable.pl
my $uniref_column = "blast_uniref50";
my $usage = "./bagel3_add_uniref_annotation.pl
				-s Sessiondir [default=current folder]
				-i queryfile
				-is_circular [default=1 for yes]
				-o output file [default=results.txt]\n
		e.g.  ./bagel3_add_uniref_annotation.pl -i 18nov2012/ALL_result.container.txt -o ALL_result.container.annotated.txt 		
" ;

my %conf  = genomics::read_conf('/usr/bagel3/bagel3.conf') ; 
&parseparam() ;


my %uniref = genomics::read_table_to_hash($unirefdb);
print "unireftable loaded\n";
my %table = genomics::read_table_to_hash($queryfile);
my @headers = genomics::get_table_header($queryfile);
print "result table loaded\n" ;

open FILE, ">$sessiondir/$outputfile" or die ("Could not write $sessiondir/$outputfile\n") ;
print FILE join ("\t".@headers)."\n" ; 
foreach my $key (sort { $table{$a}{Organism} cmp $table{$b}{Organism} || $table{$a}{AOI} <=> $table{$b}{AOI} || $table{$a}{Gene_start} <=> $table{$b}{Gene_start} } keys %table ) {
	foreach my $header (@headers) { print FILE $table{$key}{$header}."\t" ; }
	my $uniref_key = $table{$key}{$uniref_column} ;
	print FILE $uniref{$uniref_key}{description}."\n";
}
close FILE ;




sub parseparam {
    my $var ;
    my @arg = @ARGV ;
    while(@arg) {
        $var = shift(@arg) ;
        $sessiondir = shift(@arg) if($var eq '-s') ;
		$queryfile	= shift(@arg) if($var eq '-i') ;
        $outputfile = shift(@arg) if($var eq '-o') ;
    }
    die "No filename found\n$usage" if (!$queryfile) ;
	#remove the last / from sessiondir to make it universal
	$sessiondir =~ s/\/$// ;
}