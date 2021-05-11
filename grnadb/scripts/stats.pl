#!/usr/bin/perl


use Bio::SeqIO;
use Bio::Tools::Run::StandAloneBlast;
use Bio::SearchIO;
use Bio::Tools::GFF;
use Statistics::Descriptive;

use strict;


my $process_mini_stats = 1;
my $process_cdna_stats = 1;
my $process_maxi_stats = 0;
my $process_grna_hits = 1;
my $process_grna_geneclass = 1;
my $mini_file = "all_minicircles_pass.fasta";
my $cdna_file = "all_cdna_pass.fasta";

my $mini =  Bio::SeqIO->new(   '-file'         => $mini_file,
                                '-format'       => "fasta");

my $cdna =  Bio::SeqIO->new(   '-file'         => $cdna_file,
                                '-format'       => "fasta");

my $mini_size_distribution;
my $num_minicircles = 0;
my $num_used_minicircles = 0;
my @mini_size_array;
my $max_size_minicircle = 300;
my $min_size_mini_dist = 300;

my $cdna_size_distribution;
my $num_cdna = 0;
my $num_used_cdna = 0;
my @cdna_size_array;
my $max_size_cdna = 150;

my @num_array;

if($process_mini_stats)
{
	while(my $seq = $mini->next_seq)
	{
		$num_minicircles++;
		# Trim the sequence and find the largest non-X portion
		my $max_seq = $seq->seq();
		
		if(length($max_seq) >= $max_size_minicircle)
		{
			$num_used_minicircles++;
			push(@mini_size_array, length($max_seq));
			$mini_size_distribution->{length($max_seq)}++;
		}
	}


	# Now find out how much redundancy is in the data
	# Minicircle redundancy

	#iterate through all the sequences again and blast them to the all_minicircles database and find out how many hits per sequence we get


	my $mini_new = Bio::SeqIO->new('-file' => $mini_file, '-format' => "fasta");

	my @blast_params = (  'program' => 'blastn',
	                    'F'=>'F',
	                   'database' => "db/all_minicircles",
	                   "E" =>1e-40
	);

	my $blast_factory = Bio::Tools::Run::StandAloneBlast->new(@blast_params);

	my $times_sum = 0;
	my $num_minis = 0;
	print "SEQUENCE_NAME\tNUMBER_SEEN\n";

	while(my $seq = $mini_new->next_seq)
	{
		my $num_seen = 0;
		if($seq->length() >= $min_size_mini_dist)
		{
	
			$num_minis++;
			my $blast_report = $blast_factory->blastall($seq);
			my $result = $blast_report->next_result;
			while(my $hit = $result->next_hit)
			{
				while(my $hsp = $hit->next_hsp)
				{
					# Check if the length of the hsp is close to the size of the sequence
					my $hit_len = abs($hsp->hit->end - $hsp->hit->start);
					if( ($hit_len * 1.1) > $seq->length())
					{
						if($hsp->percent_identity >= 0.95)
						{
							$times_sum++;
							$num_seen++;
						}
					}
				}
			}
			push(@num_array, $num_seen);
			print $seq->display_id() .  "\t" . $num_seen . "\n";
		}	
		
	}

	my $mini_stat = Statistics::Descriptive::Full->new();

	$mini_stat->add_data(@mini_size_array);

	my $miniseen_stat = Statistics::Descriptive::Full->new();
	$miniseen_stat->add_data(@num_array);

	print "MINICIRCLE Statistics\n\n";
	print "Total number of minicircles of size greater then $max_size_minicircle is " . $mini_stat->count . "\n";
	print "Each minicircle is seen a total of " . $miniseen_stat->mean . " times in the dataset for minicircles greater then " . $min_size_mini_dist . "  bp's (" . $miniseen_stat->count . " sequences)\n";
	print "The std deviation of the distribution of sequences seen is " . $miniseen_stat->standard_deviation . "\n";
	print "The median of the distribution of sequenes seen is " . $miniseen_stat->median . "\n";
	print "The mode of the distribution of sequenes seen is " . $miniseen_stat->mode . "\n";
	print "\n";
	print "The mean minicircle size is " . $mini_stat->mean . "\n";
	print "The median minicircle size is " . $mini_stat->median . "\n";
	print "The mode minicircle size is " . $mini_stat->mode . "\n";
	print "The stdev of the minicircle size is " . $mini_stat->standard_deviation . "\n";
	print "This includes data between the values of " . $mini_stat->min . " and " . $mini_stat->max . "\n";


}

if($process_cdna_stats)
{
	while(my $seq = $cdna->next_seq)
	{
		$num_cdna++;
		# Trim the sequence and find the largest non-X portion
		my $max_size = $seq->length();
		my $max_seq = $seq->seq();
	
		if(length($max_seq) >= $max_size_cdna)
		{
			$num_used_cdna++;
			push(@cdna_size_array, length($max_seq));
			$cdna_size_distribution->{length($max_seq)}++;
		}
	}

	my $cdna_stat = Statistics::Descriptive::Full->new();

	print "\n\n\n";
	$cdna_stat->add_data(@cdna_size_array);
	
	print "CDNA Statistics\n\n";
	print "Total number of cdna's of size greater then $max_size_cdna is " . $cdna_stat->count . "\n";
	print "The mean cdna size is " . $cdna_stat->mean . "\n";
	print "The median cdna size is " . $cdna_stat->median . "\n";
	print "The mode cdna size is " . $cdna_stat->mode . "\n";
	print "The stdev of the cdna size is " . $cdna_stat->standard_deviation . "\n";
	print "This includes data between the values of " . $cdna_stat->min . " and " . $cdna_stat->max . "\n";	

}


# Now parse the gff files that were produced with the process_all.pl script and Determine the % Coverage of Guide RNA's per gene

my $maxicircle_genes_file = "maxicircle_genes.fasta";
my $mini_gff = "gff/mini.gff";
my $cdna_gff = "gff/cdna_maxicircle.gff";


if($process_maxi_stats)
{
	# Obtain each maxicircle seq
	my $maxi_genes = Bio::SeqIO->new('-file' => $maxicircle_genes_file, '-format' =>"fasta");

	my $maxicircle_coverage_hash;

	while(my $seq = $maxi_genes->next_seq())
	{

	  	my $gene_coverage;
		# Initialize the coverage array
		for(1..$seq->length())
		{
			$gene_coverage->[$_] = 0;
		}
	
		$maxicircle_coverage_hash->{$seq->display_id} = $gene_coverage;
	
	}
	
	# Iterage through the mini_hits gff file and pull out all hits to this gene
	my $gff = new Bio::Tools::GFF(	-gff_version 	=> 2,
					-file		=> $mini_gff);
	
	while(my $feature = $gff->next_feature())
	{
		my ($parent)  = $feature->gff_string() =~ /^(\w+)/;
		
		if(defined($maxicircle_coverage_hash->{$parent}) )
		{
			if($feature->primary_tag eq 'match')
			{
				# Get the parent
				for($feature->start..$feature->end)
				{
					$maxicircle_coverage_hash->{$parent}->[$_]++;
				}
			}
		}
			
	}
	
	
	
	while(my ($gene, $gene_array) = each(%$maxicircle_coverage_hash))
	{
		print "\n";
		print $gene . "\n";
		my $size = scalar(@$gene_array) - 1;
		my $covered = 0;
		for(1..$size)
		{
			print join("\t", $_, $gene_array->[$_]) . "\n";
			if($gene_array->[$_] > 0)
			{
				$covered++;
				}
		}
		print "Total Coverage for $gene is " . ($covered/$size * 100) . "\n\n";

	}
}

my $grna_fasta = "gff/mini_hits.fasta";

# This will cluster grna hits together into classes

if($process_grna_hits)
{

	my $grna_seqs = Bio::SeqIO->new( -file=>$grna_fasta, -format=>'fasta');
	my $class_hash;
	my $grna_class = 1;

	# Create a blast database of the grna hits
	system("cp $grna_fasta db/grna");
	system("cd db;formatdb -t grna -i grna -p F");

	my @blast_params = ( 'program' => 'blastn', 'F'=>'F', 'database'=>'db/grna');

	my $grna_blast = Bio::Tools::Run::StandAloneBlast->new(@blast_params);
	

	my $grna_size_hash;
	while(my $seq = $grna_seqs->next_seq())
	{
		$grna_size_hash->{$seq->display_id()} = $seq->length();
	}

	$grna_seqs = Bio::SeqIO->new( -file=>$grna_fasta, -format=>'fasta');
	
	while(my $seq = $grna_seqs->next_seq())
	{
		print "Checking " . $seq->display_id() . "\n";
		# Check if we already classified this grna
		if(defined($class_hash->{$seq->display_id}))
		{
			next;
		}

		# Blast this grna against the database of all grna's
		my $report = $grna_blast->blastall($seq);

		while(my $result = $report->next_result)
		{
			while(my $hit = $result->next_hit())
			{
				while(my $hsp = $hit->next_hsp)
				{
					# Do we hit over the full length with at most 2 mismatches
					my $query_size =  abs($hsp->query->end - $hsp->query->start) + 1;
					my $hit_size = abs($hsp->hit->end - $hsp->hit->start) + 1;
					
					# check homology string for more then 2 spaces
					my ($num_spaces) = $hsp->homology_string =~ tr/ //;
					if( ($num_spaces <= 2) && ($query_size == $hit_size) && ($hit_size == $grna_size_hash->{$hit->name}) )
					{
						# We are good
						$class_hash->{$hit->name} = $grna_class;
						print join("\t", $grna_class, $seq->display_id, $hit->name) . "\n";
						
					}
	
				} # END HSP
			} # END HIT
		} # END RESULT

		$grna_class++;
		
	} # END SEQ
}

if($process_grna_geneclass)
{
	# Now for each maxicircle gene, check how many distinct classes hit it

	my $maxicircle_distinct_grna_hash;
	my $maxicircle_grna_seqs;
	my $maxi_genes = Bio::SeqIO->new('-file' => $maxicircle_genes_file, '-format' =>"fasta");

	while(my $seq = $maxi_genes->next_seq())
	{
		$maxicircle_distinct_grna_hash->{$seq->display_id} = 1;
		my $t;
		$maxicircle_grna_seqs->{$seq->display_id} = $t;

	}
	
	# Iterate through the mini_hits gff file and pull out all hits to this gene
	my $gff = new Bio::Tools::GFF(	-gff_version 	=> 2,
					-file		=> $mini_gff);
	
	while(my $feature = $gff->next_feature())
	{
		my ($parent)  = $feature->gff_string() =~ /^(\w+)/;
		
		if(defined($maxicircle_distinct_grna_hash->{$parent}) )
		{
			if($feature->primary_tag eq 'match')
			{
				my $annotation =  $feature->annotation;
				my $annotations = $annotation->hash_tree;
				my $grna_seq =  $annotations->{'query_string'}->[0];
				my $grna_name = $annotations->{'grna_sequence'}->[0];
				my $grna_hit_name = $annotations->{'Target'}->[0];
				$grna_hit_name =~ s/grna\://;

				# Create
				my $seq_new = Bio::Seq->new( -seq=>$grna_seq, -id=>$grna_hit_name);
				push(@{$maxicircle_grna_seqs->{$parent}}, $seq_new);
			}

		}
	}

	system("mkdir tmp;");
	while(my ($gene, $hit_array) = each(%$maxicircle_grna_seqs) )
	{
		open(GENE, ">", "tmp/" . $gene);
		my $grna_size_hash;
		my $grna_class = 1;
		my $class_hash;
		foreach my $seq (@$hit_array)
		{
			print GENE ">" . $seq->display_id . "\n" . $seq->seq . "\n";
			$grna_size_hash->{$seq->display_id} = $seq->length();
		}
		close(GENE);
		system("cd tmp;formatdb -t $gene -i $gene -p F");

		# Blast this grna against the database of all grna's
		my @blast_params = ( 'program' => 'blastn', 'F'=>'F', 'database'=>'tmp/' . $gene);

		my $grna_blast = Bio::Tools::Run::StandAloneBlast->new(@blast_params);

		foreach my $seq (@$hit_array)
		{
			if(defined($class_hash->{$seq->display_id}))
			{
				next;
			}

			my $report = $grna_blast->blastall($seq);

			while(my $result = $report->next_result)
			{
				while(my $hit = $result->next_hit())
				{
					while(my $hsp = $hit->next_hsp)
					{
						# Do we hit over the full length with at most 2 mismatches
						my $query_size =  abs($hsp->query->end - $hsp->query->start) + 1;
						my $hit_size = abs($hsp->hit->end - $hsp->hit->start) + 1;
						
						# check homology string for more then 2 spaces
						my ($num_spaces) = $hsp->homology_string =~ tr/ //;
						if( ($num_spaces <= 2) && ($query_size == $hit_size) && ($hit_size == $grna_size_hash->{$hit->name}) )
						{
							# We are good
							$class_hash->{$hit->name} = $grna_class;
							print join("\t",$gene,  $grna_class, $seq->display_id, $hit->name) . "\n";
						}
					} # END HSP
				} # END HIT
			} # END REPORT
			$grna_class++;	
		} # END FOREACH SEQ
		$grna_class--;
		print "Gene $gene has $grna_class different Guide RNA Sequences\n";
	}
	
}


