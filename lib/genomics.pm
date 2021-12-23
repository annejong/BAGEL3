package genomics ;

##############################################################
##      Lib for BAGEL3 program 
## 	Anne de Jong
##	University of Groningen
##	the Netherlands
##	anne.de.jong@rug.nl
##
##############################################################
##
##   
##  2012 April, General functions for BAGEL3
##
##

use strict ;
use warnings ;
use Bio::Seq ;
use Bio::Tools::pICalculator;
use Bio::Tools::SeqPattern ;
use Bio::SearchIO ;
use Bio::SeqIO ;
use Bio::SeqIO::fasta;
use DBI ;


BEGIN {
    use Exporter ();
    @genomics::ISA         = qw(Exporter);
    @genomics::EXPORT      = qw();

}
use vars qw( %codontable) ;

my %codontable ;
my $codontablefile ='/data/bagel3/lib/codon_table.txt';
my $sessiondir = '.';

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------


sub table2html {
	# convert a tab delimted file to a HTML table
	my ($file, $htmlfile) = @_;
	open (TABLE,'<',$file) or die ("Could not find $file\n") ;	
	chomp(my(@lines) = <TABLE>);														# Get the tablefile into a array
	close TABLE ;
	my @elements ;
	my $count = 1 ;

	my @rowcolor = ("CCE5FF","99CCFF");
	my $odd = 0 ;
	open (HTML,'>',$htmlfile) or die ("Could not write $htmlfile\n") ;	
	print HTML "<table>\n" ;
	foreach my $line (@lines) {
		@elements = split /\t/, $line ;
		if ($count == 1) { # print header
			print HTML "<thead bgcolor=3399FF><b>\n\t<tr>" ;
			foreach (@elements) { print HTML "<th>$_</th>"; }
			print HTML "</tr>\n</thead></b><tbody>\n" ;
		} else {
			if ($odd) { $odd=0; } else { $odd=1; }
			print HTML "\t<tr bgcolor=$rowcolor[$odd]>" ;
			foreach (@elements) { print HTML "<td>$_</td>"; }
			print HTML "</tr>\n" ;
		}	
		$count++ ;
	}	
	print HTML "</tbody></table>\n" ;
	close HTML ;
}

sub unique_array {
	# remove duplicates / replicates from array
    return keys %{{ map { $_ => 1 } @_ }};
}

sub print_hash_of_arrays {
	my %table = @_ ;
	foreach my $key ( keys %table )  {
		print "Items in $key are: ";
		foreach ( @{$table{$key}} )  {
			print $_;
		}
	}
}

sub print_hash_of_hashes {
	my %table = @_ ;
	foreach my $key (sort keys %table )  {
		my @line ;
		foreach my $key2 (sort keys %{$table{$key}}) {
			push @line, "$key2=$table{$key}{$key2}";
		}
		print "$key:\t".join("\t",@line)."\n"; ;
	}
}

sub get_second_hashkey {
	my %table = @_ ;
	my @keys = keys %table ;
	my @keys2 = (keys %{$table{$keys[1]}}) ;
	return @keys2 ;
}

sub count_char {
    # returns the number of chars (or regular expression) in the string
	# e.g. for charged residues char = (K|R|H|D|E|C|Y)
	my ($str, $char) = @_;
    my $count = 0 ;
    $count++ while $str =~ /$char/g;
    return $count ;
}

sub calculate_pI {
    # returns the pI of the protein
	my ($seq)  = @_ ;
    my $seqobj = Bio::Seq->new(
				-seq => $seq,
				-id  => 'tocalculate',
				-alphabet => 'protein'
			    ) ; 

    my $calc = Bio::Tools::pICalculator->new(-places => 2);
    $calc->seq( $seqobj );
    my $iep    = $calc->iep ;
    my $charge = $calc->charge_at_pH( 7 );
    return sprintf("%.1f", $iep);;
	#return nearest(.1, $charge), nearest(.1, $iep) ;
}

sub RBS_score {
	# here I expect the small upstream region of the start codon, e.g. 
	my $rbs=shift ;
	my @consensus=split //,"AGGAGG";
	my @seq=split //,$rbs ;
	my $spacing = 5 ; # minimum spacing to the start coding
	my $score = 0 ;
	for my $i (0 .. scalar @seq-6-$spacing) { 
		my $count=0;
		for my $c (0..5) { $count++ if ($seq[$i+$c] eq $consensus[$c]) ; }
		$score=$count if ($count>$score) ;
	}
	return $score ;
} 



sub add_progress {
	my $string=shift;
	open (FILE,">>$sessiondir/sessionprogress") or die ("Could not write to $sessiondir/sessionprogress"); ;
	print FILE $string;
	close(FILE);
}

sub write_log {
	my $filename=shift;
	my $string=shift;
	my $print=shift;
	open (FILE,">>$filename") or die ("Could not write to $filename");
	print FILE "$string\n";
	close(FILE);
	print "$string\n" if ($print eq 'true');
	my $tmp = 1;
	return $tmp ;
}

sub write_string {
	my ($filename, $string) = @_;
	open (FILE,">$filename") or die ("Could not write to $filename");
	print FILE $string;
	close FILE ;
}

sub write_lines {
	my ($filename, @lines) = @_ ;
	open (FILE,">$filename") or die ("Could not write to $filename");
	foreach my $line (@lines) { print FILE $line."\n" ; }
	close FILE ;
}

sub append_lines {
	my ($filename, @lines) = @_ ;
	open (FILE,">>$filename") or die ("Could not write to $filename");
	foreach my $line (@lines) { print FILE $line."\n" ; }
	close FILE ;
}

sub read_lines {
	my $filename=shift;
	open (FILE,"<$filename") or die ("Could not read $filename");
	my @lines = <FILE>;	
	close FILE ;
	return @lines ;
}

sub properfilename {
	my $str = shift ;
	# these chars are not allowed in filenames - \ / | : , ; * " ? < > ( )
	$str =~ s/(\-|\ |\,|\s|\t|\\|\/|\||\:|\;|\*|\"|\?|\<|\>)/\_/g ; 
	$str =~ s/(\(|\))//g ;
	$str =~ s/\.//g ;
	$str =~ s/\.$//g ;
	$str = substr ($str, 1, 200) if (length($str) > 200) ; 
	return $str;
}	

sub add_header {
	# add the string to the top of the file
	my $filename = shift ;
	my $str = shift ;
	my @lines =  read_lines($filename) ;
	unshift (@lines, $str."\n" ) ;
	write_string ($filename, join ("", @lines)) ;
}

sub get_table_header {
	my $filename = shift ;
	my @lines =  read_lines($filename) ;
	chomp @lines ;
	my @tmp = split "\t", $lines[0];
	return @tmp ;
}

sub read_conf {
	my $filename = shift ;
	my %tmp ;
	open(FILE,"$filename") or die("could not find configuration file $filename\n") ;
	my(@lines) = <FILE>;
	chomp @lines ;
	foreach my $line (@lines) {
		$line =~ s/(\s|\t)//g ; # remove spaces and tabs
		if ($line =~ m/(.*)\=(.*)/g ) {	$tmp{$1}=$2; }	
	}
	return %tmp ;
}

sub read_table_to_hash { 
# generic routine to fill a hash from a table file 
# the file should be tab delimited and contain a header. The first columns contains the keys of the hash
	my $filename = shift ;
	my %tmp ;
	open(FILE, "$filename" ) or die("could not find the file $filename\n") ;
	my(@lines) = <FILE>;
	chomp @lines ;
	$lines[0] =~ s/\ //g ;
	my @header = split /\t/, $lines[0] ;
	my $headerline = 1 ;
	foreach my $line (@lines) {
		if ($headerline) {
			$headerline = 0 
		} else {
			my @col = split /\t/, $line ;
			if (scalar @col > 1 and $col[0] ne '') { 
				my $col_count=0;
				foreach my $columname (@header) {
					$tmp{$col[0]}{$columname} = $col[$col_count]; 
					$col_count++;
				}
			}	
		}
	}
	close(FILE) ;
	return %tmp ;
}


sub read_condontable {
	open (FILE,'<', $codontablefile) or die ("$codontablefile does not exists") ;
	my @lines = <FILE>;					# Get the codonfile into a array
	chomp(@lines) ; 
	close FILE ;
	foreach my $line (@lines) {
		my @item = split(/\t/, $line) ;
		$codontable{$item[0]} = $item[1];
	}
}

sub translate {
	my $DNA = shift ; 
	#uc $DNA ; 		# convert to uppercase
	$DNA =~ s/(\s|\n|\r|\x0d)//g; # remove line breaks spaces etc
	my $protein ;
	my $codon ;
	&read_condontable() ;
	for ( my $i=0; $i<(length($DNA)-2); $i+=3) {
		$codon=substr($DNA,$i,3);
		if (defined $codontable{$codon}) { 
			$protein.=$codontable{$codon}
		}	else {
			$protein.='x' ;
		}
	}
	return $protein ;
}

sub inverse_complement {
	my $DNA = shift ;
	$DNA =~ s/(\s|\n|\r|\x0d)//g;
	my $revcomp = reverse($DNA);
	$revcomp =~ tr/[A,C,G,T,a,c,g,t]/[T,G,C,A,t,g,c,a]/;
	return $revcomp;
}
	


## the mandatory one (without it no package!!!)
1

