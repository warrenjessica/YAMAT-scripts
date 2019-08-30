#!/usr/bin/perl

##Script to summarize to total number of fastq reads and number that end in typical tRNA CCA tail
##Written by Jess Warren and Dan Sloan


use warnings;
use strict;

my $usage = "\nUSAGE: perl $0 fastq_file\n\n";

my $file_name = shift or die ($usage);

my @fastq_lines = file_to_array ($file_name);

my $total=0;
my $CCA_total=0;
my %DB_hash;

for(my $i=1; $i < scalar(@fastq_lines); $i+=4){
	$total=$total+1;
	chomp $fastq_lines[$i];
	
	if (substr($fastq_lines[$i], -3) eq 'CCA'){
		$CCA_total=$CCA_total+1;
		
		++$DB_hash {substr($fastq_lines [$i], -4, 1)};
	}
}	

print "File_Name\tTotal_sequences\tSeqs_w_CCA";


foreach(sort keys (%DB_hash)){
	print "\t$_";
}

print "\n";

print "$file_name\t$total\t$CCA_total";

foreach(sort keys (%DB_hash)){
	print "\t$DB_hash{$_}";
}
print "\n";
 

# ##subroutines

sub file_to_array {
	use strict;
	use warnings;

    my($filename) = @_;

    #Initialize variables
    my @filedata = (  );

    unless( open(GET_FILE_DATA, $filename) ) {
        print STDERR "Cannot open file \"$filename\"\n\n";
        exit;
    }

    @filedata = <GET_FILE_DATA>;

    close GET_FILE_DATA;

    return @filedata;
}