#!/usr/bin/perl

##Script to summarize to YAMAT-seq variants based on prior blast output
##Written by Jess Warren and Dan Sloan

use strict;
use warnings;
use Bio::SearchIO; 
use Getopt::Long;
use IPC::Cmd qw(can_run);
use sloan;

my $usage = 
"\nUsage: perl $0 [options/arguments]
   
   This script takes blast output and the corresponding query and database
   files in fasta format. It parses the top hit (or hits if there are ties),
   realigns the sequences with a global aligner (mafft) and reports SNPs and
   indels.
   
   REQUIRED ARGUMENTS
   
   BLAST output file
         --blast
         Blast output file in standard BLAST text format
   
   Query Fasta
         --query
         Fasta file containing the sequences that were used as queries to
         generate BLAST output

   Database Fasta
         --db
         Fasta file containing the sequences that were used as a database
         to generate BLAST output

   Output Name
         --output
         Output file name
         
         
   OPTIONAL ARGUMENTS
 
   Blast E-value Cutoff
         --evalue [default: 1e-10]      
         Maximum e-value for blast hits

   Blast Hit Coverage Cutoff
         --coverage [default: 0.9]      
         Proportion of length of hit sequence covered by the best BLAST HSP
   
   MAFFT Executable       
         --mafft [default: mafft]
         Path and name of mafft executable if not already in your path
     
   Retain sequences with Ns       
         --retain_Ns [default: off]
         Add this flag if you want to keep YAMAT reads that contain Ns.
               
   EXAMPLE
   
         perl $0 --blast=blast.txt --query=query.fas --db=db.fas --output=my_output       
";


our $BLAST;
our $QUERY;
our $DB;
our $OUTPUT;
our $EVALUE = 1e-10;
our $COVERAGE = 0.9;
our $MAFFT = "mafft";
our $RETAIN_NS;

##Extract parameters from command-line flags

GetOptions(
    'blast=s'  => \$BLAST,
    'query=s'  => \$QUERY,
    'output=s'  => \$OUTPUT,
    'db=s'  => \$DB,    
    'evalue=i'  => \$EVALUE,    
    'coverage=i'  => \$COVERAGE,
    'retain_Ns'  => \$RETAIN_NS,
    'mafft=s'  => \$MAFFT
);

##Print out start time and command line call of Perl script
print "\n" . (localtime) . "\nRunning $0 with the following command:\n", qx/ps -o args $$/, "\n";

##Check dependencies
can_run ($MAFFT) or die ("\n$usage\n\nERROR: Could not find following executable in PATH: $MAFFT. Provide a path and executable name with --mafft.\n\n");

##Check file inputs
$BLAST or die ("\n$usage\n\nERROR: Must provide a BLAST filename with --blast.\n\n");
-e $BLAST or die ("\n$usage\n\nERROR: File specified with --blast does not exist: $BLAST.\n\n");
$QUERY or die ("\n$usage\n\nERROR: Must provide the fasta query input filename with --query.\n\n");
-e $QUERY or die ("\n$usage\n\nERROR: File specified with --query does not exist: $QUERY.\n\n");
$DB or die ("\n$usage\n\nERROR: Must provide the fasta database filename with --db.\n\n");
-e $DB or die ("\n$usage\n\nERROR: File specified with --db does not exist: $DB.\n\n");
$OUTPUT or die ("\n$usage\n\nERROR: Must provide output name with --output.\n\n");

#store fasta files in hashes
my %query_fasta = fasta2hash ($QUERY);
my %db_fasta = fasta2hash ($DB);

#open output files
my $FH1 = open_output("$OUTPUT\.variants.txt");
my $FH2 = open_output("$OUTPUT\.hits.txt");
my $FH3 = open_output("$OUTPUT\.matches.txt");


print $FH1 "Hit Name\tQuery Name\tPosition\tReference\tVariant\n";

#open blast output file and loop through one line at a time
my %hit_HoA;
my %cov_HoAoA;
foreach my $key (sort keys %db_fasta){
	my $tRNAlength = length($db_fasta{$key});
	for (my $i=0;$i< $tRNAlength; ++$i){
		@{$cov_HoAoA{$key}[$i]}=(0,0,0,0,0,0,0,0);
	}
}

my $SearchIO_obj = new Bio::SearchIO(-format => 'blast',-file => $BLAST);  
while( my $result_obj = $SearchIO_obj->next_result ){
	my $query_name = $result_obj->query_name;

	unless ($RETAIN_NS){	
		if ($query_fasta{$query_name} =~ /N/){
			next;
		}
	}

	my @splitname= split(/\-/, $query_name);
	my $readcount = $splitname[-1];
	my $hit_bitscore = 0;
	while(my $hit_obj = $result_obj->next_hit ) {
		#only proceed past first hit if additional hits have exact same bitscore
		my $new_hit_bitscore = $hit_obj->raw_score;       
        if ($hit_bitscore){
        	$hit_bitscore == $new_hit_bitscore or last;
        }
        $hit_bitscore = $new_hit_bitscore;
        my $hit_evalue = $hit_obj->significance;
        
        #skip hits that do not meet the evalue threshold
        $hit_evalue <= $EVALUE or last;
        
        #skip hits that do not meet the coverage threshold
        my $hit_length = $hit_obj->length;
        my $hsp_obj = $hit_obj->next_hsp;
        my $hsp_length = $hsp_obj->length('hit');
        $hsp_length / $hit_length >= $COVERAGE or last;
        
        #skip hits that are on the reverse strand
        $hsp_obj->strand('hit') == -1 and next;
        
        #add query to a list associated with each hit sequence
        my $hit_name = $hit_obj->name;
        push (@{$hit_HoA{$hit_name}}, $query_name);
        
        #extract full length sequences from fasta files and run global alignment with mafft        
        exists ($db_fasta{$hit_name}) or die ("\nERROR: Could not find $hit_name in $DB\.\n\n");
        exists ($query_fasta{$query_name}) or die ("\nERROR: Could not find $query_name in $QUERY\.\n\n");
        
        my $FHT = open_output ("TEMPALIGNMENTINPUT");
        print $FHT ">QUERY\n$query_fasta{$query_name}\n>HIT\n$db_fasta{$hit_name}\n";
        system ("mafft TEMPALIGNMENTINPUT > TEMPALIGNMENTOUTPUT 2>> $OUTPUT\.mafft_logs.txt");
        
        my %aligned_fasta = fasta2hash ("TEMPALIGNMENTOUTPUT");

		#summarize variants
		my $query_string = $aligned_fasta{"QUERY"};
		my $hit_string = $aligned_fasta{"HIT"};

		my $curr_pos = 1;
		my $in_insertion = 0;
		my $insertion_string;
		for (my $i = 0; $i < length ($hit_string); ++$i){ # for loop through each position in the hit string
			if (substr($hit_string, $i, 1) eq '-'){ #if the hit has a gap character there is an insertion in the query
				$insertion_string .= uc(substr($query_string, $i, 1));
				$in_insertion = 1;
				next;
			}elsif($in_insertion){
				print $FH1 "$hit_name\t$query_name\t$curr_pos\tN/A\tInsertion (", $insertion_string ,")\n";
				$in_insertion = 0;
				$insertion_string = "";
				$cov_HoAoA{$hit_name}[$curr_pos -1][7] += $readcount;
			}
			
			if (substr($query_string, $i, 1) eq '-'){ #id the query has a gap character there is a deletion in the query
				print $FH1 "$hit_name\t$query_name\t$curr_pos\t",  uc(substr ($hit_string, $i, 1)),"\tDeletion\n";
				$cov_HoAoA{$hit_name}[$curr_pos -1][1] += $readcount;
				$cov_HoAoA{$hit_name}[$curr_pos -1][6] += $readcount;				
			}elsif (uc(substr($hit_string, $i, 1)) ne uc(substr($query_string, $i, 1))){#if there no gap characters and they still don't match, it's a substition
				print $FH1 "$hit_name\t$query_name\t$curr_pos\t",  uc(substr ($hit_string, $i, 1)),"\t", uc(substr ($query_string, $i, 1)),"\n";
				$cov_HoAoA{$hit_name}[$curr_pos -1][1] += $readcount;
				
				my $character = uc(substr($query_string, $i, 1));
				if ($character eq "A"){
					$cov_HoAoA{$hit_name}[$curr_pos -1][2] += $readcount;
				}elsif($character eq "C"){
					$cov_HoAoA{$hit_name}[$curr_pos -1][3] += $readcount;
				}elsif($character eq "G"){
					$cov_HoAoA{$hit_name}[$curr_pos -1][4] += $readcount;
				}elsif($character eq "T"){
					$cov_HoAoA{$hit_name}[$curr_pos -1][5] += $readcount;
				}elsif($character eq "N"){
					$cov_HoAoA{$hit_name}[$curr_pos -1][8] += $readcount;
				}else {
					die("Error non ACGTN character found: $character/n/n");
				}			
			}else{
				$cov_HoAoA{$hit_name}[$curr_pos -1][0] += $readcount;
			}
			
			++$curr_pos;
			
		}
		if ($in_insertion){
			print $FH1 "$hit_name\t$query_name\t$curr_pos\tN/A\tInsertion (", $insertion_string ,")\n";	
			$cov_HoAoA{$hit_name}[$curr_pos -1][7] += $readcount;
		}
			
		unlink ("TEMPALIGNMENTINPUT");        
		unlink ("TEMPALIGNMENTOUTPUT");        
	}   
} 

#print out the queries that align to each hit
print $FH2 "Hit\tQueries\n";
foreach (sort keys %hit_HoA){
	print $FH2 "$_\t";
	my @query_list = @{$hit_HoA{$_}};
	my $first_query = shift (@query_list);
	print $FH2 $first_query;
	foreach my $next_query (@query_list){
		print $FH2 "\;$next_query";
	}
	print $FH2 "\n";
}

#print out the match data
print $FH3 "tRNA\tPosition\tReference\tMatch\tMismatch\tA\tC\tG\tT\tDel\tIns\tN\n";
foreach my $key (sort keys %cov_HoAoA){
	my $tRNAlength = length($db_fasta{$key});
	
	for (my $i=0;$i< $tRNAlength; ++$i){
		print $FH3 "$key\t", $i+1, "\t", substr($db_fasta{$key}, $i, 1);
		my @countarray =  @{$cov_HoAoA{$key}[$i]};
		foreach(@countarray){
			print $FH3 "\t$_";
		}
		print $FH3 "\n";
	}
}