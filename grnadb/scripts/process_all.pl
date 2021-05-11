#!/usr/bin/perl


use Bio::SeqIO;
use Bio::Tools::Run::StandAloneBlast;
use Bio::SearchIO;
use Getopt::Long;
  
use strict;

my $db_name = 'tbgrna';
  
my $process_maxicircle_sequences = 0;
my $process_minicircle_sequences = 0;
my $process_grna_sequences = 0;
my $process_minicircle_inverted_repeats = 0;
my $process_minicircle_features = 0;
my $process_predicted_grna_features = 0;
my $process_blast_databases = 0;
my $process_download = 0;
my $run_statistics = 0;
my $load_database = 0;


my $options = GetOptions(
				"maxicircle", \$process_maxicircle_sequences,
				"minicircle", \$process_minicircle_sequences,
				"grna", \$process_grna_sequences,
				"repeats", \$process_minicircle_inverted_repeats,
				"minicircle_features", \$process_minicircle_features,
				"predicted_grna", \$process_predicted_grna_features,
				"blast", \$process_blast_databases,
				"download", \$process_download,
				"statistics", \$run_statistics,
				"load", \$load_database,
				"db_name=s", \$db_name
		);



my $palindrome_bin = 'palindrome';
####### Start Maxicircle Sequences

my $web_dir = '/srv/www/htdocs/kiss';
my $download_dir = "$web_dir/download";
my $bioware_dir = '/usr/local/bioware';

my $cdna_sequence_file = "all_cdna_pass.fasta";
my $gene_file = "maxicircle_genes.fasta";
my $all_gene_database = "db/all";
 
my $gff_outfile = "gff/cdna_maxicircle.gff";


if($process_maxicircle_sequences)
{ 
 
# Create overlap array of hases
my $overlap;
my $overlap_bloodstream;
my $overlap_procyclic;
my $sequence_fasta_hash; 
 
my $genes =  Bio::SeqIO->new(   '-file'         => $gene_file,
                                '-format'       => "fasta");
 
my $sequences = Bio::SeqIO->new('-file'         => $cdna_sequence_file,
                                '-format'       => "fasta");
 
 
# Create blast databases for each sequence under the db directory
 
my $gene_hash;
open(GFF, ">", $gff_outfile);
 
open(GENEALL, ">", $all_gene_database);
 
while(my $gene = $genes->next_seq)
{
        open(GENE, ">", "db/" . $gene->display_id);
 
        print GENE ">" . $gene->display_id . "\n" . $gene->seq();
        print GENEALL ">" . $gene->display_id . "\n" . $gene->seq() . "\n";
        close(GENE);
        system("cd db;formatdb -t " . $gene->display_id . " -i " . $gene->display_id . " -p F");
        $gene_hash->{$gene->display_id} = $gene;
        print GFF join("\t", $gene->display_id, "pcr", "gene", "1", $gene->length, ".", "+", "0", "Sequence " . $gene->display_id) . "\n";
        # initialize overlap array
        for(0..$gene->length)
        {
                $overlap->{$gene->display_id}[$_] = 0;
        }
 
 
}
 
system("cd db;formatdb -t all -i all -p F");
system("cd db;xdformat -n all");
 
my $un_array;
my $ed_array;
 
my @all_params = (  'program' => 'blastn',
                    'F'=>'F',
                   'database' => "db/all",
                   "E" =>1e-3
);
 
my $all_factory = Bio::Tools::Run::StandAloneBlast->new(@all_params);
my $min_cdna_size = 150;

my $est_hit_num = 1;
while(my $seq = $sequences->next_seq)
{
	$sequence_fasta_hash->{$seq->display_id} = $seq;

        if($seq->length < $min_cdna_size)
        {
                next;
        }
        print GFF join("\t", $seq->display_id, "pcr", "est", "1", $seq->length, ".", "+", "0", "Sequence " . $seq->display_id) . "\n";
        my $gff_type = 'dt_primed';
	if($seq->display_id() =~ /BS_11/)
	{
		$gff_type = 'specific_bloodstream_fraction_11';
	}elsif($seq->display_id() =~ /BS_9/)
	{
		$gff_type = 'specific_bloodstream_fraction_9';
	}elsif($seq->display_id() =~ /PC_11/)
	{
		$gff_type = 'specific_procyclic_fraction_11';
	}elsif($seq->display_id() =~ /PC_9/)
	{
		$gff_type = 'specific_procyclic_fraction_9';	
        }elsif($seq->display_id() =~ /BS/)
        {
                $gff_type = 'specific_bloodstream';
        } elsif($seq->display_id() =~ /PC/)
        {
                $gff_type = 'specific_procyclic';
        } elsif($seq->display_id() =~ /TOCo3001/)
	{
		$gff_type = 'specific_bloodstream';
	}
        # First blast if vs the all db to find out which gene it is closest to
        my $all_report = $all_factory->blastall($seq);
        my $all_result = $all_report->next_result;
 
        my $hit = $all_result->next_hit;
        my $gene_name;
 
        if($hit)
        {
                my $hsp = $hit->next_hsp;
                my $unedited;
                my $edited;
                my $single_gene = 0;
                if($hsp->evalue < 1e-3)
                {
                } else
                {
                        next;
                }
 
                my $gene_hit =  $hit->name;
 
                if($gene_hit =~ /ed$/)
                {
                        ($gene_name) = $gene_hit =~ /(.+)ed$/;
                        $unedited = $gene_name . "un";
                        $edited = $gene_name . "ed";
                } elsif($gene_hit =~ /un$/)
                {
                        ($gene_name) = $gene_hit =~ /(.+)un$/;
                        $unedited = $gene_name . "un";
                        $edited = $gene_name . "ed";
                }else
                {
                        $gene_name = $gene_hit;
                        $unedited = $gene_name;
                        $edited = $gene_name;
                        $single_gene = 1;
                }
 
 
                my $unedit_seq = $gene_hash->{$unedited};
                my $edit_seq = $gene_hash->{$edited};
 
 
                my @un_params = (  'program' => 'blastn',
                                        'F'=>'F',
                                'database' => "db/" . $unedited
                );
 
                my @ed_params = (  'program' => 'blastn',
                                        'F'=>'F',
                                'database' => "db/" . $edited
                );
 
                my $un_factory = Bio::Tools::Run::StandAloneBlast->new(@un_params);
                my $ed_factory = Bio::Tools::Run::StandAloneBlast->new(@ed_params);

                my $max_unedit_size = 0;
                my $max_edit_size = 0;
                my $max_edit_start = 0;
                my $max_edit_stop = 0;
                my $max_unedit_start = 0;
                my $max_unedit_stop = 0;
                my $un_hash;
                my $ed_hash;
 
 
                if($gene_hash->{$unedited})
                {
                        # Find area's
                        my $un_report = $un_factory->blastall($seq);
                        my $un_result = $un_report->next_result;
			my $seq_start = 0;
			my $seq_stop = 0;
			my $hit_string;
			my $query_string;
			my $homology_string;
 
                        while(my $hit = $un_result->next_hit)
                        {
                                while(my $hsp = $hit->next_hsp)
                                {
 
                                        if( ($hsp->hit->end - $hsp->hit->start) > $max_unedit_size)
                                        {
                                                $max_unedit_size = $hsp->hit->end - $hsp->hit->start;
                                                $max_unedit_start = $hsp->hit->start;
                                                $max_unedit_stop = $hsp->hit->end;
                         
						$seq_start = $hsp->query->start;
						$seq_stop = $hsp->query->end;
                                                $hit_string = $hsp->hit_string;
                                                $query_string = $hsp->query_string;
                                                $homology_string = $hsp->homology_string;
			               }
 
                                }
                        }
                        $un_hash->{size} = $max_unedit_size;
                        $un_hash->{start} = $max_unedit_start;
                        $un_hash->{stop} = $max_unedit_stop;
                        $un_hash->{name} = $seq->display_id;
 
			if($seq_start)
			{
                        	print GFF join("\t", $unedit_seq->display_id, "HSP", $gff_type, $un_hash->{start}, $un_hash->{stop}, ".", ".", ".", "Target EST:esthit_" . $est_hit_num . " $seq_start $seq_stop ; est_name \"" . $un_hash->{name} . "\" ; gene_type_match \"pre-edited\" ; query_start \"$seq_start\" ; query_stop \"$seq_stop\" ; hit_string \"$hit_string\" ; query_string \"$query_string\" ; homology_string \"$homology_string\" ; " ) . "\n";
				print GFF join("\t", $un_hash->{name}, "HSP", 'gene_match', $seq_start, $seq_stop, ".", ".", ".", "genehit " . $unedit_seq->display_id . '_hit ; hit_start "' . $un_hash->{start} . '" ; hit_stop "' . $un_hash->{stop} . "\" ; hit_string \"$hit_string\" ; query_string \"$query_string\" ; homology_string \"$homology_string\"" . ' ; gene_type_match "pre-edited" ;') . "\n";
				open(GENE, ">>", 'fasta/' . $unedit_seq->display_id . '.fasta');
				print GENE ">" . $un_hash->{name} . "\n" . $seq->subseq($seq_start, $seq_stop) . "\n";
				close(GENE);
			}
                        for($un_hash->{start}..$un_hash->{stop})
                        {
                                $overlap->{$unedit_seq->display_id}[$_]++;
                        }
                }
 
                if($gene_hash->{$edited} && !($single_gene))
                {
                        my $ed_report = $ed_factory->blastall($seq);
                        my $ed_result = $ed_report->next_result;
			my $seq_start = 0;
			my $seq_stop = 0;
			my $hit_string = '';
			my $query_string = '';
			my $homology_string = '';
 
                        while(my $hit = $ed_result->next_hit)
                        {
                                while(my $hsp = $hit->next_hsp)
                                {
                                                                                                                                                                                                                                                    
               
                                        if( ($hsp->hit->end - $hsp->hit->start) > $max_edit_size)
                                        {
                                                $max_edit_size = $hsp->hit->end - $hsp->hit->start;
                                                $max_edit_start = $hsp->hit->start;
                                                $max_edit_stop = $hsp->hit->end;
						$seq_start = $hsp->query->start;
						$seq_stop = $hsp->query->end;
						$hit_string = $hsp->hit_string;
						$query_string = $hsp->query_string;
						$homology_string = $hsp->homology_string;
                                        }
                                }
                        }
                        $ed_hash->{size} = $max_edit_size;
                        $ed_hash->{start} = $max_edit_start;
                        $ed_hash->{stop} = $max_edit_stop;
                        $ed_hash->{name} = $seq->display_id;
 			if($seq_start)
			{
                        	print GFF join("\t", $edit_seq->display_id, "HSP", $gff_type, $ed_hash->{start}, $ed_hash->{stop}, ".", ".", ".", "Target EST:esthit_" . $est_hit_num . " $seq_start $seq_stop ; est_name \"" . $ed_hash->{name} . "\" ; gene_type_match \"edited\" ; query_start \"$seq_start\" ; query_stop \"$seq_stop\" ; hit_string \"$hit_string\" ; query_string \"$query_string\" ; homology_string \"$homology_string\" ; " ) . "\n";
				print GFF join("\t", $ed_hash->{name}, "HSP", 'gene_match', $seq_start, $seq_stop, ".", ".", ".", "genehit " . $edit_seq->display_id . '_hit ; hit_start "' . $ed_hash->{start} . '" ; hit_stop "' . $ed_hash->{stop} . "\" ; hit_string \"$hit_string\" ; query_string \"$query_string\" ; homology_string \"$homology_string\"" . ' ; gene_type_match "edited" ; ') . "\n";
				open(GENE, ">>", 'fasta/' . $edit_seq->display_id . '.fasta');
				print GENE ">" . $ed_hash->{name} . "\n" . $seq->subseq($seq_start, $seq_stop) . "\n";
				close(GENE);
			}
                        for($ed_hash->{start}..$ed_hash->{stop})
                        {
                                $overlap->{$edit_seq->display_id}->[$_]++;
                        }
                }
 
        } # End If hit
	$est_hit_num++;
} # End while seq
 
 
# Now print out the overlap for the GFF
 
while(my ($gene_name, $overlap_array) = each(%$overlap))
{
        my $count = 0;
        foreach my $overlap_score (@$overlap_array)
        {
                if($count == 0)
                {
                } else
                {
                        print GFF join("\t", $gene_name, "overlap", "overlap", $count, $count, $overlap_score, ".", ".") . "\n";
                }
                $count++;
        }
}
}


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
my $mini_remove_hit_file = "$web_dir/remove_log.tab";
my $non_hit_hash;

my $min_sequence_size = 300;
my $min_hit_size = 20;
my $max_hit_size = 80;

my $percent_t_transcript_max = 80;
my $percent_c_minicircle_max = 80;

my $min_cont_hits = 15;
my $max_non_hits = 2;
my $min_hit_length = 30;

if($process_minicircle_sequences)
{

	my $mini =  Bio::SeqIO->new(   '-file'         => $mini_file,
                                '-format'       => "fasta");
  
	open(GFF, ">", $mini_gff_outfile);
	open(FASTA, ">", $mini_fasta_hit_file);
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

	while(my $mini_seq = $mini->next_seq)
	{
		if($mini_seq->length < $min_sequence_size)
		{
			next;
		}

		# Print to the gff file for the minicircles  sequence
	        print GFF join("\t", $mini_seq->display_id, "pcr", "minicircle", "1", $mini_seq->length, ".", "+", "0", "Sequence " . $mini_seq->display_id) . "\n";
	}

	# NEW

	my $mini_seq_db = Bio::SeqIO->new(-file=>$mini_db, -format=>'fasta');

	while(my $mini_seq = $mini_seq_db->next_seq)
	{
		my $tmp_fasta_obj = Bio::SeqIO->new(-file=>">$input_temp_file", -format=>'fasta');
		$tmp_fasta_obj->write_seq($mini_seq);

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

					my $hit_full = join("-", $hit->name, $hit_start, $hit_stop, $mini_seq->display_id, $query_start, $query_stop);


					my $source = 'pcr';
					my $mini_source = 'blast';
					if($hit_size < 22)
					{
						$source = 'pcr-putitive';
						$mini_source = 'blast-putitive';
					}
					# If it is in the no-hit hash, it was previously deleted, so mark it as such
					if( defined($non_hit_hash->{$hit_full}) )
					{
						$source = 'pcr-deleted';
						$mini_source = 'blast-deleted';
					}
					
					my $pass = 1;
					# Place any code that would make you throw out this hit entierly below this line

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
					if(!$pass)
					{
						$source = 'pcr-deleted';
						$mini_source = 'blast-deleted';
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

					my $hit_name = $hit->name . '_' . $mini_seq->display_id . '.hit-' . $hit_id;

					print 	GFF join("\t", $hit->name, $source, "match", $hit_start, $hit_stop, ".", $cdna_strand, "0", "Target grna:" . $mini_seq->display_id . " $query_start $query_stop ; query_start \"$query_start\" ; " 
						. "query_stop \"$query_stop\" ; hit_string \"" . $hsp->hit_string 
						. '" ; query_string "' . $hsp->query_string . '" ; homology_string "' . $hsp->homology_string . '" ; grna_sequence "' . $mini_seq->display_id . '" ; ') . "\n";
					print 	GFF join("\t", $hit->name, $source, "HSP", $hit_start, $hit_stop, ".", $cdna_strand, "0", "Target grna:" . $mini_seq->display_id . " $query_start $query_stop ; query_start \"$query_start\" ; " 
						. "query_stop \"$query_stop\" ; hit_string \"" . $hsp->hit_string 
						. '" ; query_string "' . $hsp->query_string . '" ; homology_string "' . $hsp->homology_string . '" ; grna_sequence "' . $mini_seq->display_id . '" ; ') . "\n";

#					print FASTA ">" . $hit_name . "\n" . uc($mini_seq->seq) . "\n";
	
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

					print 	GFF join("\t", $mini_seq->display_id, $mini_source, "mini_match", $qstart, $qstop, ".", $mini_strand, "0", "minicircle_match " . $hit->name . " ; query_start \"$hit_start\" ; " 
						. "query_stop \"$hit_stop\" ; hit_string \"" . $hsp->query_string 
						. '" ; query_string "' . $hsp->hit_string . '" ; homology_string "' . $hsp->homology_string . '" ; grna_sequence "' . $mini_seq->display_id . '" ; grna_hit_name "' . $hit_name . '" ;') . "\n";
	
					$hit_id++;
					my $q_seq = $hsp->query_string;
					$q_seq =~ s/[\-\ ]//i;
					# Strip .ab1 from the names
					my $mini_name = $mini_seq->display_id;
					my $hit_name = $hit->name;
					$mini_name =~ s/\.ab1//i;
					$hit_name =~ s/\.ab1//i;
					my $q_seq_obj = Bio::Seq->new(-display_id=>join("-", $mini_name, $hit_start, $hit_name, $query_start), -seq=>$q_seq);
					#my $q_seq_obj = Bio::Seq->new(-display_id=>join("-", $mini_seq->display_id, $query_start, $query_stop, $hit->name, $hit_start, $hit_stop) . '--' . $hit_name, -seq=>$mini_seq->trunc($query_start, $query_stop)->seq);
					$grna_fasta_out->write_seq($q_seq_obj);

				}
			}
		}
		

	}


	close(GFF);


}
####### End Minicircle







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










if($process_minicircle_inverted_repeats)
{
	my $mini =  Bio::SeqIO->new(   '-file'         => $mini_file,
                                '-format'       => "fasta");
  
	open(GFF, ">", $repeats_gff_outfile);
	my $minirepeat_id = 1;
	my $repeat_num_mismatches = 3;
	my $repeat_min_distance = 50;
	my $repeat_max_distance = 150;

	while(my $mini_seq = $mini->next_seq)
	{
		if($mini_seq->length < $min_sequence_size)
		{
			next;
		}

		open(TMP, ">", "temp/seq.fasta");
		print TMP ">" . $mini_seq->display_id . "\n" . $mini_seq->seq() . "\n";
		close(TMP);
		system("$palindrome_bin temp/seq.fasta -minpallen 14 -maxpallen 30 -gaplimit $repeat_max_distance -nummismatches $repeat_num_mismatches -overlap -outfile temp/mini.pal >& /dev/null");
		open(MINIPAL, "<", "temp/mini.pal");

		my $start_proc = 0;
		my $fwd_start = 0;
		my $fwd_stop = 0;
		my $rev_start = 0;
		my $rev_stop = 0;
		my $fwd_seq = '';
		my $rev_seq = '';
		my $match_seq = '';


		while(<MINIPAL>)
		{
			my $line = $_;
			if($start_proc)
			{
				# Is this a fwd line or a reverse line
				if($line =~ /^\d/)
				{
					my ($start, $seq, $stop) = $line =~ /^(\d+)\s+(.+)\s+(\d+)$/;
					$seq =~ s/^\s+//;
					$seq =~ s/\s+$//;

					if($stop > $start)
					{
						# This is a fwd part of the repeat
						$fwd_start = $start;
						$fwd_stop = $stop;
						$fwd_seq = $seq;
					} else
					{
						# This is a reverse part of the repeat
						$rev_start = $stop;
						$rev_stop = $start;
						$rev_seq = $seq;
						
						# Check if we are constrained to the min distance 
						if( ($rev_start - $fwd_stop + 1) >= $repeat_min_distance)
						{
							# Print out the gff 
							print GFF join("\t", $mini_seq->display_id, "palindrome", "inverted_repeats", $fwd_start, $fwd_stop, ".", "+", "0", "mini_repeat \"repeat_$minirepeat_id.fwd\" ; rev_start \"$rev_start\" ; rev_stop \"$rev_stop\" ; fwd_start \"$fwd_start\" ; fwd_stop \"$fwd_stop\" ; fwd_string \"$fwd_seq\" ; rev_string \"$rev_seq\" ; homology_string \"$match_seq\" ;" ) . "\n";
							print GFF join("\t", $mini_seq->display_id, "palindrome", "inverted_repeats", $rev_start, $rev_stop, ".", "-", "0", "mini_repeat \"repeat_$minirepeat_id.rev\" ; rev_start \"$rev_start\" ; rev_stop \"$rev_stop\" ; fwd_start \"$fwd_start\" ; fwd_stop \"$fwd_stop\" ; fwd_string \"$fwd_seq\" ; rev_string \"$rev_seq\" ; homology_string \"$match_seq\" ;" ) . "\n";
						}

						# Clear out the variables
						my $fwd_start = 0;
						my $fwd_stop = 0;
						my $rev_start = 0;
						my $rev_stop = 0;
						my $fwd_seq = '';
						my $rev_seq = '';
						my $match_seq = '';
						$minirepeat_id++;

					} 
				} elsif($line =~ /\|/)
				{
						($match_seq) = $line =~ /^\s+(.+)/;
				}
			} elsif($line =~ /^Palindromes\:/)
			{
				$start_proc = 1;
			}
		}

		close(MINIPAL);

	}
}

my $features_gff = 'gff/features.gff';
my $conserved_region_file = 'conserved_region.fasta';

if($process_minicircle_features)
{
	open(FEATUREGFF, ">", $features_gff);


	# Check for Conserved regions
	my $conserved_seq_io = Bio::SeqIO->new(-file=>$conserved_region_file, -format=>'fasta');
	while(my  $conserved_region_seq = $conserved_seq_io->next_seq())
	{

		my $mini =  Bio::SeqIO->new('-file'=> $mini_file, '-format'=>"fasta");
		my $bl2seq_factory = Bio::Tools::Run::StandAloneBlast->new('program'=>'blastn');

		while(my $mini_seq = $mini->next_seq())
		{
			if($mini_seq->length < $min_sequence_size)
			{
				next;
			}
			my $report = $bl2seq_factory->bl2seq($mini_seq, $conserved_region_seq);

			my $result = $report->next_result;
			while(my $hit = $result->next_hit)
			{
				while(my $hsp = $hit->next_hsp)
				{
					my $start = $hsp->query->start;
					my $stop = $hsp->query->end;
					my $dir = '+';
					if($hsp->hit->start > $hsp->hit->end)
					{
						$dir = "-";
					}
					# Write to the gff file
					my $hit_start = $hsp->hit->start;
					my $hit_stop = $hsp->hit->end;
					print FEATUREGFF join("\t", $mini_seq->display_id, 'bl2seq', 'nc_conserved_region', $start, $stop, '.', $dir, '.', 'Sequence "' . $hit->name . '" ;'  
					. " hit_start \"$hit_start\" ; hit_stop \"$hit_stop\" ; hit_string \"" 
					. $hsp->hit_string . '" ; query_string "' 
					. $hsp->query_string . '" ; homology_string "' . $hsp->homology_string . '" ; '
					) . "\n";
				}
			}
	
		}

		my $cdna = Bio::SeqIO->new('-file'=>$cdna_sequence_file, '-format'=>"fasta");
		while(my $cdna_seq = $cdna->next_seq())
		{
			if($cdna_seq->length < $min_sequence_size)
			{
				next;
			}
			my $report = $bl2seq_factory->bl2seq($cdna_seq, $conserved_region_seq);
	
			my $result = $report->next_result;
			while(my $hit = $result->next_hit)
			{
				while(my $hsp = $hit->next_hsp)
				{
					my $start = $hsp->query->start;
					my $stop = $hsp->query->end;
					my $dir = '+';
					if($hsp->hit->start > $hsp->hit->end)
					{
						$dir = "-";
					}
					# Write to the gff file
					my $hit_start = $hsp->hit->start;
					my $hit_stop = $hsp->hit->end;
					print FEATUREGFF join("\t", $cdna_seq->display_id, 'bl2seq', 'nc_conserved_region', $start, $stop, '.', $dir, '.', 'Sequence "' . $hit->name . '" ;'  
					. " hit_start \"$hit_start\" ; hit_stop \"$hit_stop\" ; hit_string \"" 
					. $hsp->hit_string . '" ; query_string "' 
					. $hsp->query_string . '" ; homology_string "' . $hsp->homology_string . '" ; '
					) . "\n";
				}
			}
		}
	
	}
		

}



my $predicted_grna_gff = 'gff/predicted_grna.gff';
my $predicted_grnas = 'all_grna_predicted.fasta';
my $sequenced_unique_grnas = 'all_grna_sequenced_unique.fasta';
my $unique_db = 'all_grna_nr';

if($process_predicted_grna_features)
{
	# First predict the grnas
	system("../scripts/create_nr_merge_grna.pl --db_name=$db_name  --file=$predicted_grnas --predicted");
	system("../scripts/create_nr__merge_grna.pl --db_name=$db_name  --file=$sequenced_unique_grnas --sequenced");
	open(PREDICTEDGFF, ">", $predicted_grna_gff);
	system("cat $predicted_grnas $sequenced_unique_grnas > db/$unique_db");
	system("cd db;formatdb -p F -i $unique_db");


	# Check for Conserved regions
	my $conserved_seq_io = Bio::SeqIO->new(-file=>'db/' . $unique_db, -format=>'fasta');
	while(my  $conserved_region_seq = $conserved_seq_io->next_seq())
	{

		my $mini =  Bio::SeqIO->new('-file'=> $mini_file, '-format'=>"fasta");
		my $bl2seq_factory = Bio::Tools::Run::StandAloneBlast->new('program'=>'blastn');

		while(my $mini_seq = $mini->next_seq())
		{
			if($mini_seq->length < $min_sequence_size)
			{
				next;
			}
			my $report = $bl2seq_factory->bl2seq($mini_seq, $conserved_region_seq);

			my $result = $report->next_result;
			while(my $hit = $result->next_hit)
			{
				while(my $hsp = $hit->next_hsp)
				{
					if($conserved_region_seq->length != length($hsp->hit_string) || $hsp->homology_string =~ /\ /)
					{
						next;
					}
					my $start = $hsp->query->start;
					my $stop = $hsp->query->end;
					my $dir = '+';
					if($hsp->hit->start > $hsp->hit->end)
					{
						$dir = "-";
					}
					# Write to the gff file
					my $hit_start = $hsp->hit->start;
					my $hit_stop = $hsp->hit->end;
					print PREDICTEDGFF join("\t", $mini_seq->display_id, 'bl2seq', 'predicted_grna', $start, $stop, '.', $dir, '.', 'Sequence "' . $hit->name . '" ;'  
					. " hit_start \"$hit_start\" ; hit_stop \"$hit_stop\" ; hit_string \"" 
					. $hsp->hit_string . '" ; query_string "' 
					. $hsp->query_string . '" ; homology_string "' . $hsp->homology_string . '" ; '
					) . "\n";
				}
			}
	
		}

	
	}
		

}




#my $blast_db_dir = '/xraid/bioware/gmod/data/blastdb';
my $blast_db_dir = "$bioware_dir/grnadb/blastdb'";
my $blastdb_minicircles = $blast_db_dir . "/" . $db_name . "_minicircles";
my $blastdb_grna = $blast_db_dir . "/" . $db_name . "_grna";
my $blastdb_transcript = $blast_db_dir . "/" . $db_name . "_transcript";
my $blastdb_genes = $blast_db_dir . "/" . $db_name . "_genes";
my $blastdb_genes_transcript = $blast_db_dir . "/" . $db_name . "_genes_transcript";

if($process_blast_databases)
{

	system("cp $mini_file $blastdb_minicircles;cd $blast_db_dir;formatdb -i $blastdb_minicircles -p F;xdformat -n $blastdb_minicircles");
	system("cp $grna_fasta_file $blastdb_grna;cd $blast_db_dir;formatdb -i $blastdb_grna -p F;xdformat -n $blastdb_grna");
	system("cp $cdna_sequence_file $blastdb_transcript;cd $blast_db_dir;formatdb -i $blastdb_transcript -p F;xdformat -n $blastdb_transcript");
	system("cp $gene_file $blastdb_genes;cd $blast_db_dir;formatdb -i $blastdb_genes -p F;xdformat -n $blastdb_genes");
	system("cp $cdna_plus_genes_database $blastdb_genes_transcript;cd $blast_db_dir;formatdb -i $blastdb_genes_transcript -p F;xdformat -n $blastdb_genes_transcript");
	
}

if($process_download)
{
	
	unlink("$download_dir/minicircle_sequences.fasta.gz");
	system("cp $mini_file $download_dir/minicircle_sequences.fasta;cd $download_dir;gzip minicircle_sequences.fasta");
	unlink("$download_dir/sequenced_grna.fasta.gz");
	system("cp $grna_file $download_dir/sequenced_grna.fasta;cd $download_dir;gzip sequenced_grna.fasta");
	unlink("$download_dir/sequenced_grna_unique.fasta.gz");
	system("cp $sequenced_unique_grnas $download_dir/sequenced_grna_unique.fasta;cd $download_dir;gzip sequenced_grna_unique.fasta");
	unlink("$download_dir/transcript_sequences.fasta.gz");
	system("cp $cdna_sequence_file $download_dir/transcript_sequences.fasta;cd $download_dir;gzip transcript_sequences.fasta");
	unlink("$download_dir/maxicircle_sequences.fasta.gz");
	system("cp $gene_file $download_dir/maxicircle_sequences.fasta;cd $download_dir;gzip maxicircle_sequences.fasta");
	unlink("$download_dir/predicted_grna.fasta.gz");
	system("cp $predicted_grnas $download_dir/predicted_grna.fasta;cd $download_dir;gzip predicted_grna.fasta");
	system("chmod a+r $download_dir/*");
}

if($run_statistics)
{
	system("$bioware_dir/grnadb/scripts/stats/unique_grna_sequences_per_gene.pl --db_name=$db_name --verbose --outdir=$web_dir/download --webdir=/kiss/download --predicted --sequenced");
	system("$bioware_dir/grnadb/scripts/stats/unique_grna_sequences_per_gene.pl --db_name=$db_name --verbose --outdir=$web_dir/download/gcdna --webdir=/kiss/download/gcdna --sequenced");
	system("$bioware_dir/grnadb/scripts/stats/unique_grna_sequences_per_gene.pl --db_name=$db_name --verbose --outdir=$web_dir/download/predicted --webdir=/kiss/download/predicted --predicted");

	system("$bioware_dir/grnadb/scripts/stats/stats.pl --process_minicircle > $web_dir/download/minicircle_statistics.txt");
	system("$bioware_dir/grnadb/scripts/stats/stats.pl --process_cdna > $web_dir/download/transcript_statistics.txt");
	system("$bioware_dir/grnadb/scripts/stats/stats.pl --process_gcdna > $web_dir/download/gcdna_statistics.txt");

	system("$bioware_dir/grnadb/scripts/grna_find_multipurpose.pl --db_name=$db_name --gene_name=tbrps12 > $web_dir/download/multipurpose_rps12.txt");
	system("$bioware_dir/grnadb/scripts/grna_find_multipurpose.pl --db_name=$db_name --gene_name=tbnd9 > $web_dir/download/multipurpose_nd9.txt");
	system("$bioware_dir/grnadb/scripts/grna_find_multipurpose.pl --db_name=$db_name --gene_name=tbnd8 > $web_dir/download/multipurpose_nd8.txt");
	system("$bioware_dir/grnadb/scripts/grna_find_multipurpose.pl --db_name=$db_name --gene_name=tbnd7 > $web_dir/download/multipurpose_nd7.txt");
	system("$bioware_dir/grnadb/scripts/grna_find_multipurpose.pl --db_name=$db_name --gene_name=tbnd3 > $web_dir/download/multipurpose_nd3.txt");
	system("$bioware_dir/grnadb/scripts/grna_find_multipurpose.pl --db_name=$db_name --gene_name=tba6 > $web_dir/download/multipurpose_a6.txt");
	system("$bioware_dir/grnadb/scripts/grna_find_multipurpose.pl --db_name=$db_name --gene_name=tbmurf2 > $web_dir/download/multipurpose_murf2.txt");
	system("$bioware_dir/grnadb/scripts/grna_find_multipurpose.pl --db_name=$db_name --gene_name=tbcyb > $web_dir/download/multipurpose_cyb.txt");
	system("$bioware_dir/grnadb/scripts/grna_find_multipurpose.pl --db_name=$db_name --gene_name=tbcr4 > $web_dir/download/multipurpose_cr4.txt");
	system("$bioware_dir/grnadb/scripts/grna_find_multipurpose.pl --db_name=$db_name --gene_name=tbcr3 > $web_dir/download/multipurpose_cr3.txt");
	system("$bioware_dir/grnadb/scripts/grna_find_multipurpose.pl --db_name=$db_name --gene_name=tbco3 > $web_dir/download/multipurpose_co3.txt");
	system("$bioware_dir/grnadb/scripts/grna_find_multipurpose.pl --db_name=$db_name --gene_name=tbco2 > $web_dir/download/multipurpose_co2.txt");
	system("$bioware_dir/grnadb/scripts/runall_grna_transcript.pl --db_name=$db_name --output_dir=$web_dir/download/");

}

if($load_database)
{
	system("../scripts/load_gff.pl $db_name");
}


sub check_alignment
{
	my $query_string = shift;
	my $homology_string = shift;
	my $hit_string = shift;

	# Find number of mismatches
	my ($num_mismatches) = $homology_string =~ tr/ //; 
	
	# Find number of gaps
	my ($num_gaps) = $query_string =~ tr/-//;
	($num_gaps) += $hit_string =~ tr/-//;
	
	# Find number of true matches
	my ($true_matches) = $homology_string =~ tr/|//;
	
	# Find number of similarity matches
	my ($sim_matches) = $homology_string =~ tr/+//;

	return ($num_mismatches, $num_gaps, $true_matches, $sim_matches);
	

}

sub check_align_pass
{
	my $min_cont_hits = shift;
	my $max_non_hits = shift;
	my $min_length = shift;

	my $query_string = shift;
	my $homology_string = shift;
	my $hit_string = shift;

	# Check for minimum length requirements
	my $length = length($query_string);
	if($length < $min_cont_hits)
	{
		return 0;
	} elsif($length < $min_length)
	{
		return 0;
	}

	# Check first
	# for each $min_cont_hits window, make sure I have at least 1 window with no spaces
	my $temp_start = 0;
	my $temp_stop = $min_cont_hits - 1;
	my $check_cont_hits = 0;

	while($temp_stop < $length)
	{
		my $check_string = substr($homology_string, $temp_start, $min_cont_hits);
		my ($num_spaces) = $check_string =~ tr/ //;
		$temp_stop++;
		$temp_start++;
		if($num_spaces == 0)
		{
			$temp_stop = $length;
			$check_cont_hits = 1;
		}
	}

	# If we failed, return false
	if(!$check_cont_hits)
	{
		return 0;
	}
	
	# Now check for at most $max_non_hits spaces within a $max_hit_length window

	$temp_start = 0;
	$temp_stop = 0;

	my $check_non_hits = 0;
	
	while($temp_stop < $length)
	{
		my $check_string = substr($homology_string, $temp_start, $min_length);
		my($num_spaces) = $check_string =~ tr/ //;
		$temp_stop++;
		$temp_start++;
		if($num_spaces <= $max_non_hits)
		{
			$check_non_hits = 1;
		}
	}

	if(!$check_non_hits)
	{
		return 0;
	}

	# If we passed everything, return true

	return 1;

}



