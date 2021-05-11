#!/usr/bin/perl


use Bio::SeqIO;
use Bio::Tools::Run::StandAloneBlast;
use Bio::SearchIO;
use Getopt::Long;
  
use strict;

my $db_name = 'tbgrna';
  


my $palindrome_bin = 'palindrome';
####### Start Maxicircle Sequences

my $download_dir = '/xraid/bioware/gmod/mblweb-gmod/html/kiss/download';
my $cdna_sequence_file = "all_cdna_pass.fasta";
my $gene_file = "maxicircle_genes.fasta";
my $all_gene_database = "db/all";
 
my $gff_outfile = "gff/cdna_maxicircle.gff";

####### End Maxicircle cdna Sequences



####### Start Minicircle Sequences
my $mini_file = "all_minicircles_pass.fasta";
my $mini_gff_outfile = "gff/mini.gff";
my $repeats_gff_outfile = "gff/repeats.gff";
my $mini_fasta_hit_file = "gff/mini_hits.fasta";
my $grna_fasta_file = "gff/grna.fasta";
my $mini_db = "db/all_minicircles";
my $wu_blastn_bin = "blastn";
#my $wu_blast_options = "-w 2 -matrix goldenoxduk-one -warnings -notes";
my $wu_blast_options = "-hspmax 0 -gspmax 0 -span2 -W 2 -matrix goldenoxduk-hitstatic -warnings -notes";
my $temp_file = "/tmp/temp_blast.out";
my $input_temp_file = "/tmp/temp_blast.fasta";

my $cdna_plus_genes_database = "db/all_genes_plus_cdna";
my $mini_remove_hit_file = "/xraid/bioware/gmod/mblweb-gmod/html/kiss/remove_log.tab";
my $non_hit_hash;

my $min_sequence_size = 300;
my $min_hit_size = 20;
my $max_hit_size = 80;

my $percent_t_transcript_max = 80;
my $percent_c_minicircle_max = 80;

my $min_cont_hits = 15;
my $max_non_hits = 2;
my $min_hit_length = 30;

####### End Minicircle




#### Start Rudo Filter

my $five_linker1 = Bio::Seq->new(-display_id=>'fivelinker1', -seq=>'GCCTCCCTCGCGCCATCAGGATC');
my $five_linker2 = Bio::Seq->new(-display_id=>'fivelinker2', -seq=>'GCCTCCCTCGCGCCATCAGGACC');
my $five_linker  = Bio::Seq->new(-display_id=>'fivelinker2', -seq=>'GCCTCCCTCGCGCCATCAGGA');

my $three_linker = Bio::Seq->new(-display_id=>'threelinker', -seq=>'TTTTTTTTTCTGAGCGGGCTGGCAAGGC');
my $grna_file = "original/grna_rudo.fasta";
my $out_file = "grna_rudo_pass.fasta";

my $out_seqio = Bio::SeqIO->new(-format=>'fasta', -file=>">$out_file");
my $in_seqio = Bio::SeqIO->new(-format=>'fasta', -file=>"$grna_file");


my $bl2seq_factory = Bio::Tools::Run::StandAloneBlast->new('program'=>'blastn', 'F'=>'F');
while(my $seq = $in_seqio->next_seq)
{
	print "---------------------------------------------------\n";
	print $seq->display_name .  " " . $seq->length() . "bp\n";
	print $seq->seq() . "\n";
	# Find if we can get the 5' sequence

	my $start_3;
	my $end_3;
	my $strand_3;

	my $start_5;
	my $end_5;
	my $strand_5;

	my $start;
	my $end;
	my $strand;

	my $subseq;

	my $report = $bl2seq_factory->bl2seq($seq, $three_linker);
	my $result = $report->next_result;

        if(my $hit = $result->next_hit)
	{
		if(my $hsp = $hit->next_hsp)
		{
			# Get this hit/hsp
			# The query is the cDNA	
			print join("\t", 'query3', $hsp->query()->start(), $hsp->query()->end(), $hsp->query()->strand()) . "\n";
			# The hit is the linker
			print join("\t", 'hit3', $hsp->hit()->start(), $hsp->hit()->end(), $hsp->hit()->strand()) . "\n";
			$start_3 = $hsp->query()->start();
			$end_3 = $hsp->query()->end();
			$strand_3 = $hsp->query()->strand();
			# Correct for NN's and ambigious bases
			if($hsp->hit()->strand == 1)
			{
				$start_3 -= 2;
			} else
			{
				$end_3 += 2;
			}
		}
	}
	# Now 
	my $report = $bl2seq_factory->bl2seq($seq, $five_linker);
	my $result = $report->next_result;

        if(my $hit = $result->next_hit)
	{
		if(my $hsp = $hit->next_hsp)
		{
			# Get this hit/hsp
			# The query is the cDNA	
			print join("\t", 'query5', $hsp->query()->start(), $hsp->query()->end(), $hsp->query()->strand()) . "\n";
			# The hit is the linker
			print join("\t", 'hit5', $hsp->hit()->start(), $hsp->hit()->end(), $hsp->hit()->strand()) . "\n";
			$start_5 = $hsp->query()->start();
			$end_5 = $hsp->query()->end();
			$strand_5 = $hsp->query()->strand();
			# Correct for NN's and ambigious bases
			if($hsp->hit()->strand == 1)
			{
				$end_5 += 2;
			} else
			{
				$start_5 -= 2;
			}
		}
	}
	if(defined($start_3) && defined($start_5) )
	{
		if(($start_3 < $start_5) && ($end_3 < $start_5) )
		{
			# 33333333GGGGGGGGGGG55555555555
			$start = $end_3;
			$end = $start_5;
			$strand = -1;
			
		} elsif(($start_5 < $start_3) && ($end_5 < $start_3) )
		{
			# 55555555GGGGGGGGGGGGGG3333333333333333
			$start = $end_5;
			$end = $start_3;
			$strand = 1;
		} else
		{
			# They overlap, so ignore it
			
		}
	}
	if(defined($start) && defined($end) )
	{
		$subseq = $seq->trunc($start, $end);
		$subseq->desc("[$start-$end] " . ($end - $start + 1) . "bp");
		if($subseq->length() >= 20)
		{
			$out_seqio->write_seq($subseq);
		}
		print $subseq->seq() . "\n";
	}


}

exit;

=cut

####### Start GRNA SEQUENCED

my $grna_file = "all_grna_pass.fasta";
my $grna_gff_outfile = "gff/grna.gff";
my $grna_fasta_hit_file = "gff/grna_hits.fasta";
my $grna_fasta_file = "gff/grna_sequenced.fasta";
my $grna_db = "db/all_grna";
my $wu_blastn_bin = "blastn";
#my $wu_blast_options = "-w 2 -matrix goldenoxduk-one -warnings -notes";
my $wu_blast_options = "-hspmax 0 -gspmax 0 -span2 -W 2 -matrix goldenoxduk-hitstatic -warnings -notes";
my $temp_file = "/tmp/temp_blast.out";
my $input_temp_file = "/tmp/temp_blast.fasta";

my $cdna_plus_genes_database = "db/all_genes_plus_cdna";
my $mini_remove_hit_file = "/xraid/bioware/gmod/mblweb-gmod/html/kiss/remove_log.tab";
my $non_hit_hash;

my $min_sequence_size = 20;
my $min_hit_size = 20;
my $max_hit_size = 80;

my $percent_t_transcript_max = 80;
my $percent_c_minicircle_max = 80;

my $min_cont_hits = 15;
my $max_non_hits = 2;
my $min_hit_length = 30;

if($process_grna_sequences)
{

	my $grna =  Bio::SeqIO->new(   '-file'         => $grna_file,
                                '-format'       => "fasta");
  
	open(GFF, ">", $grna_gff_outfile);
	open(FASTA, ">", $grna_fasta_hit_file);
	my $grna_fasta_out = Bio::SeqIO->new(-file=>">$grna_fasta_file", -format=>'fasta');


	# Create cdna plus gene database
	system("cat $gene_file $cdna_sequence_file > $cdna_plus_genes_database");
	system("cd db;xdformat -n all_genes_plus_cdna");
	my $hit_id = 1;


	# Search through all of the non_hits file and create a hash of them

	open(NONHIT, $mini_remove_hit_file);
	while(<NONHIT>)
	{
		chomp($_);
		my ($filedb_name, $hit_name, $hit_start, $hit_stop, $mini_name, $mini_start, $mini_stop) = split("\t", $_);
		if($db_name eq $filedb_name)
		{
			my $hit_full  = join("-", $hit_name, $hit_start, $hit_stop, $mini_name, $mini_start, $mini_stop);
			my $mini_full = join("-", $mini_name, $mini_start, $mini_stop, $hit_name, $hit_start, $hit_stop);
			$non_hit_hash->{$hit_full} = 1;
			$non_hit_hash->{$mini_full} = 1;
		}
	}

	while(my $grna_seq = $grna->next_seq)
	{
		if($grna_seq->length < $min_sequence_size)
		{
			next;
		}

		# Print to the gff file for the minicircles  sequence
	        print GFF join("\t", $grna_seq->display_id, "pcr", "grna", "1", $grna_seq->length, ".", "+", "0", "Sequence " . $grna_seq->display_id) . "\n";
	}

	# NEW
	# Now iterate through the cdna's fully edited/preedited sequences and blast that to the minicircles

	my $grna_seq_db = Bio::SeqIO->new(-file=>$grna_file, -format=>'fasta');

	while(my $grna_seq = $grna_seq_db->next_seq)
	{
		my $tmp_fasta_obj = Bio::SeqIO->new(-file=>">$input_temp_file", -format=>'fasta');
		$tmp_fasta_obj->write_seq($grna_seq);

		# Take the minicircle sequence and blast it to the maxicircles (edited and unedited) and to the set of query cdna files using a wublast gapped blast allowing G-U base pairings
		# blast minicircle to all_gene_database of fully edited and pre-edited sequences
		
		system("$wu_blastn_bin $wu_blast_options -d $cdna_plus_genes_database -o $temp_file -i $input_temp_file");
		my $searchio = new Bio::SearchIO( -format=>'blast', -file=>$temp_file);

		while(my $result = $searchio->next_result)
		{
			while( my $hit = $result->next_hit)
			{
				while( my $hsp = $hit->next_hsp)
				{
					# Find out if we pass the minimum/maximum size range
					my $hit_size = abs($hsp->hit->end - $hsp->hit->start) + 1;
					my ($hit_matches) = $hsp->homology_string =~ tr/[+\|]//;
					my $hit_percent = ($hit_matches / $hit_size) * 100;
					my $hit_start = $hsp->hit->start;
					my $hit_stop = $hsp->hit->end;
					my $query_start = $hsp->query->start;
					my $query_stop = $hsp->query->end;
					# Find out how many gaps and mismatches there are (count blanks)
					my $hs = $hsp->homology_string;
					my ($num_blanks) = $hs =~ tr/\ //;
					
					# Check if we are going to use this one again

					my $hit_full = join("-", $hit->name, $hit_start, $hit_stop, $grna_seq->display_id, $query_start, $query_stop);


					my $source = 'grna-pcr';
					my $mini_source = 'grna-blast';
					if($hit_size < 22)
					{
						$source = 'grna-pcr-putitive';
						$mini_source = 'grna-blast-putitive';
					}
					# If it is in the no-hit hash, it was previously deleted, so mark it as such
					if( defined($non_hit_hash->{$hit_full}) )
					{
						$source = 'grna-pcr-deleted';
						$mini_source = 'grna-blast-deleted';
					}
					
					my $pass = 1;
					# Place any code that would make you throw out this hit below this line

					# Check for poly T in the transcript
					my $num_t_transcript = $hsp->hit_string =~ tr/[Tt]//;
					my $num_c_minicircle = $hsp->query_string =~ tr/[Cc]//;
					my $percent_t_transcript = ($num_t_transcript/$hit_size)*100;
					my $percent_c_minicircle = ($num_c_minicircle/$hit_size)*100;
					if($percent_t_transcript > $percent_t_transcript_max)
					{
						$pass = 0;
					}
					if($percent_c_minicircle > $percent_c_minicircle_max)
					{
						$pass = 0;
					}
					# Check for poly-A/T in the hit
					if($hsp->query_string =~ /(AAAAAAAA)|(TTTTTTTT)/i)
					{
						$pass = 0;
						next;
					}
					
					if(!$pass)
					{
						$source = 'grna-pcr-deleted';
						$mini_source = 'grna-blast-deleted';
					}
	
					my $mini_strand = $hsp->query->strand;
					if($mini_strand < 0)
					{
						$mini_strand = "-";
					} else
					{
						$mini_strand = "+";
					}
					my $cdna_strand = $hsp->hit->strand;
					if($cdna_strand < 0)
					{
						$cdna_strand = "-";
					} else
					{
						$cdna_strand = "+";
					}

					my $hit_name = $hit->name . '_' . $grna_seq->display_id . '.hit-' . $hit_id;

					print 	GFF join("\t", $hit->name, $source, "match", $hit_start, $hit_stop, ".", $cdna_strand, "0", "Target grna:" . $grna_seq->display_id . " $query_start $query_stop ; query_start \"$query_start\" ; " 
						. "query_stop \"$query_stop\" ; hit_string \"" . $hsp->hit_string 
						. '" ; query_string "' . $hsp->query_string . '" ; homology_string "' . $hsp->homology_string . '" ; grna_sequence "' . $grna_seq->display_id . '" ; ') . "\n";
					print 	GFF join("\t", $hit->name, $source, "HSP", $hit_start, $hit_stop, ".", $cdna_strand, "0", "Target grna:" . $grna_seq->display_id . " $query_start $query_stop ; query_start \"$query_start\" ; " 
						. "query_stop \"$query_stop\" ; hit_string \"" . $hsp->hit_string 
						. '" ; query_string "' . $hsp->query_string . '" ; homology_string "' . $hsp->homology_string . '" ; grna_sequence "' . $grna_seq->display_id . '" ; ') . "\n";

#					print FASTA ">" . $hit_name . "\n" . uc($grna_seq->seq) . "\n";
	
					# Now map the reverse of the hit , show the cdna/gene on the minicircle

					my $qstart = $query_start;
					my $qstop = $query_stop;
					my $qdir = "+";
					if($qstart > $qstop)
					{
						my $temp = $qstart;
						$qstart = $qstop;
						$qstop = $temp;
						$qdir = "-";
								
					}

					print 	GFF join("\t", $grna_seq->display_id, $mini_source, "mini_match", $qstart, $qstop, ".", $mini_strand, "0", "minicircle_match " . $hit->name . " ; query_start \"$hit_start\" ; " 
						. "query_stop \"$hit_stop\" ; hit_string \"" . $hsp->query_string 
						. '" ; query_string "' . $hsp->hit_string . '" ; homology_string "' . $hsp->homology_string . '" ; grna_sequence "' . $grna_seq->display_id . '" ; grna_hit_name "' . $hit_name . '" ;') . "\n";
	
					$hit_id++;
					my $q_seq = $hsp->query_string;
					$q_seq =~ s/[\-\ ]//i;
					# Strip .ab1 from the names
					my $mini_name = $grna_seq->display_id;
					my $hit_name = $hit->name;
					$mini_name =~ s/\.ab1//i;
					$hit_name =~ s/\.ab1//i;
					my $q_seq_obj = Bio::Seq->new(-display_id=>join("-", $mini_name, $hit_start, $hit_name, $query_start), -seq=>$q_seq);
					#my $q_seq_obj = Bio::Seq->new(-display_id=>join("-", $grna_seq->display_id, $query_start, $query_stop, $hit->name, $hit_start, $hit_stop) . '--' . $hit_name, -seq=>$grna_seq->trunc($query_start, $query_stop)->seq);
					$grna_fasta_out->write_seq($q_seq_obj);

				}
			}
		}
		

	}


	close(GFF);


}
####### End GRNA SEQUENCED

=cut












