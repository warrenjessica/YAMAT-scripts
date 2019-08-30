#!/usr/bin/perl
use strict;
use warnings;
use Bio::SearchIO; 

#This script is used to summarize YAMAT tRNA-seq data by producing a table of all blast hits to reads and a matrix of read counts to the Arabidopsis tRNA DB file using e-value and hit_cov cutoffs. 
##Written by Jess Warren and Dan Sloan

my $usage = "\nUSAGE: $0  INPUTFILESUMMARY   tRNA_DB   E-VALUE   OUTPUT_BASE_NAME HIT_COV_CUTOFF\n\n";

#Summary file needs to be a tab delimited list with samplename	blastfilehits	readsfile


my $input = shift or die ($usage); 
my $DB = shift or die ($usage); 
my $evalue_threshold = shift or die ($usage);
my $output = shift or die ($usage);
my $percent_cov = shift or die ($usage);

my %matrixHoH;  
my %matrixHoH_ties; 


my %hashDB= fasta2hash($DB);
my @DBnames= sort(keys(%hashDB));

my @libraries= file_to_array($input); 


my $FH_SUMM = open_output($output.".summary.txt");

print $FH_SUMM "Library\tQuery_Name\tRead_Count\tSequence\tHit_name\tFilter\tCov_hit_filter\tGenome\tAA\tAnticodon\ttRNA_ID\tEvalue\tBitscore\tPercentID\tSNPs\tGaps\tQuery_Len\tQuery_Cov\tHit_Len\tHit_Cov\tTied\n";

my $FH_Ns = open_output ($output."excluded_Ns.txt");

print $FH_Ns "Library\tExcluded YAMAT Families with Ns\tExcluded YAMAT read count\n";

foreach my $line (@libraries){

	chomp $line; 

	my @library_array=split(/\t/, $line);

	foreach my $trna (@DBnames){ 
		$matrixHoH{$library_array[0]}->{$trna} = 0;
		$matrixHoH_ties{$library_array[0]}->{$trna} = 0;
	}

	my %queryhash = fasta2hash($library_array[2]); 
		
	my $SearchIO_obj = new Bio::SearchIO(-format => 'blast',-file => $library_array[1]);

	my $N_fams = 0;
	my $N_count = 0;

	while( my $result_obj = $SearchIO_obj->next_result ) {
		my $query_name = $result_obj->query_name;
		
		my @split_name = split (/\-/, $query_name); 
		
		my $query_length = $result_obj->query_length;

		
		if ($queryhash{$query_name} =~ /N/){
			$N_count += $split_name[-1];
			++$N_fams;
			next;
		}

		my $hit_bitscore = 0;
		my $hit_evalue; 
		my $hit_length;
		my $hsp_percent;
		my $hsp_length_query;
		my $hsp_length_hit;
		my $SNPs;
		my $gaps;
		my @hit_names;
		
		while(my $hit_obj = $result_obj->next_hit ) {
			my $new_hit_bitscore = $hit_obj->raw_score;
			
			if ($hit_bitscore){
				$hit_bitscore == $new_hit_bitscore or last;
			}
			
			$hit_bitscore = $new_hit_bitscore; 
			push (@hit_names, $hit_obj->name); 
					
			if (scalar (@hit_names) == 1){
				$hit_evalue = $hit_obj->significance; 
				$hit_length = $hit_obj->length;
				my $hsp_obj = $hit_obj->next_hsp;
				$hsp_percent = $hsp_obj->percent_identity;
				$hsp_length_hit = $hsp_obj->length('hit');
				$hsp_length_query = $hsp_obj->length('query');
				$gaps = $hsp_obj->gaps;
				my $hsp_length_aln =  $hsp_obj->length('total');
				my @match_array = $hsp_obj->matches('query');
				$SNPs = $hsp_length_aln - $match_array[0] - $gaps;
			}
		}
				
		
		print $FH_SUMM "$library_array[0]\t$query_name\t$split_name[-1]\t$queryhash{$query_name}\t";

		if (@hit_names){
			
			my @split_hit = split (/\-/, $hit_names[0]);
			
			my $genome = $split_hit[0];
			my $trna_species = $split_hit[1];
			my $trna_ID = $split_hit[2];
			
			my $anticodon = substr($trna_species, -3);
			my $AA = substr($trna_species, 0, length($trna_species) - 3);
			
			print $FH_SUMM "$hit_names[0]\t";
						
			if ($hit_evalue <= $evalue_threshold){
				print $FH_SUMM ".\t";
				
				if (scalar($hsp_length_hit/$hit_length) >= $percent_cov){
					print $FH_SUMM ".\t";
					if (scalar(@hit_names) >= 2){
					
						foreach my $tie_hit (@hit_names){ 
						
							$matrixHoH{$library_array[0]}->{$tie_hit} += $split_name[-1] / scalar(@hit_names);
							$matrixHoH_ties{$library_array[0]}->{$tie_hit} += $split_name[-1] / scalar(@hit_names);
						}
					}else{
						$matrixHoH{$library_array[0]}->{$hit_names[0]} += $split_name[-1];
						}
				}else{
					print $FH_SUMM "cov_fail\t";
					}
			}else{
				print $FH_SUMM "FILTER\t\t";
				
			}
			
			unless ($trna_ID){
				print STDERR "Problem strong: $hit_names[0]\n";
			}
			print $FH_SUMM "$genome\t$AA\t$anticodon\t$trna_ID\t$hit_evalue\t$hit_bitscore\t$hsp_percent\t$SNPs\t$gaps\t$hsp_length_query\t";
			
			print $FH_SUMM $hsp_length_query/$query_length, "\t$hit_length\t", $hsp_length_hit/$hit_length, "\t";
			
			if (scalar(@hit_names) >= 2){
				print $FH_SUMM "TIE:$hit_names[0]";
				for (my $i=1; $i < scalar (@hit_names); ++$i){
					print $FH_SUMM ";$hit_names[$i]";
				}
				print $FH_SUMM "\n";
			}else{
				print $FH_SUMM ".\n";
			}
		
		}else{
			print $FH_SUMM "NO_HIT\tFILTER\n";
		}	
	}
	print "$library_array[0]\t$N_fams\t$N_count\n";
}
	

	

my $FH_MAT = open_output ($output."countmatrix.txt");
my $FH_MAT_TIES = open_output ($output."countmatrix_ties.txt");

print $FH_MAT "Library";
print $FH_MAT_TIES "Library";

my @libs = sort keys %matrixHoH;

my @trnas = sort keys %{$matrixHoH{$libs[0]}};
	
foreach my $column_head (@trnas){
	print $FH_MAT "\t$column_head";
	print $FH_MAT_TIES "\t$column_head";
}
print $FH_MAT "\n";	
print $FH_MAT_TIES "\n";	

foreach (@libs){
	print $FH_MAT $_;
	print $FH_MAT_TIES $_;
	
	foreach my $second_key (@trnas){
		print $FH_MAT "\t", $matrixHoH{$_}->{$second_key};
		print $FH_MAT_TIES "\t", $matrixHoH_ties{$_}->{$second_key};
	}
	
	print $FH_MAT "\n";
	print $FH_MAT_TIES "\n";
	
}	










sub get_fasta_names_and_seqs {
	use strict;
	use warnings;
	my ($inputfilename) = @_;
	my @fasta_names = ();
	my @fasta_seqs= ();

		   
	unless ( open(FILEDATA, $inputfilename) ) {
		print STDERR "Cannot open file \"$inputfilename\"\n\n";
		exit;
	}	

	my @filedata = <FILEDATA>;
	close FILEDATA;
	
	my $seq_count = 0;
	foreach my $line (@filedata){
		if ($line =~ /^>/) {
			if ($line =~ /^>.*[\w]+/){
				my $partialLine = substr ($&, 1);
				push (@fasta_names, $partialLine); 
				push (@fasta_seqs, ""); 
				++$seq_count; 
			}
		}else {
			$fasta_seqs[$seq_count-1] .= $line;
		}
	}
	for (my $i = 0; $i < scalar (@fasta_seqs); ++$i){
		$fasta_seqs[$i] =~s/\s//g;
	}
	
	return (\@fasta_names, \@fasta_seqs);
}


sub fasta2hash {
	my $fasta = shift @_ or die ("\nERROR: No fasta file name provided to fasta2hash\n\n");
	my %fastaHash = arrays2hash (get_fasta_names_and_seqs($fasta));
	return %fastaHash;
}


sub arrays2hash {
	use strict;
	use warnings;

	(my $keyarray, my $valuearray) = @_;
	if (scalar(@$keyarray) != scalar(@$valuearray)) {
		die "Arrays differ in size: Mismatched number of keys and values"; 
	}
	
	my %newhash = ( );
	
	@newhash{ @$keyarray } = @$valuearray;


	return (%newhash);

}


sub file_to_array {
	use strict;
	use warnings;

    my($filename) = @_;

    my @filedata = (  );

    unless( open(GET_FILE_DATA, $filename) ) {
        print STDERR "Cannot open file \"$filename\"\n\n";
        exit;
    }

    @filedata = <GET_FILE_DATA>;

    close GET_FILE_DATA;

    return @filedata;
}

sub open_output {
	use strict;
	use warnings;

    my($filename) = @_;
    my $fh_output;

    unless(open($fh_output, ">$filename")) {
        print "Cannot open file $filename\n";
        exit;
    }
    return $fh_output;
}
