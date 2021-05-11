#!/usr/bin/perl


use Bio::SeqIO;
use Bio::Tools::Run::StandAloneBlast;
use Bio::SearchIO;
use Bio::Tools::GFF;
use Statistics::Descriptive;
use Getopt::Long;
use strict;


my $process_mini_stats = 0;
my $process_cdna_stats = 0;
my $process_gcdna_stats = 0;

my $mini_file = "all_minicircles_pass.fasta";
my $cdna_file = "all_cdna_pass.fasta";
my $gcdna_file = 'all_grna_pass.fasta';

my $db_name;
my $verbose = 0;

my $options = GetOptions(
                                "db_name=s", \$db_name,
                                "verbose!", \$verbose,
				"process_minicircle", \$process_mini_stats,
				"process_cdna", \$process_cdna_stats,
				"process_gcdna", \$process_gcdna_stats,
		);



my $mini =  Bio::SeqIO->new(   '-file'         => $mini_file,
                                '-format'       => "fasta");


my $mini_size_distribution;
my $num_minicircles = 0;
my $num_used_minicircles = 0;
my @mini_size_array;
my $max_size_minicircle = 300;
my $min_size_mini_dist = 300;
my $mini_text;

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
	$mini_text .= "SEQUENCE_NAME\tNUMBER_SEEN\n";

	my %seen_hash;
	while(my $seq = $mini_new->next_seq)
	{
		if(defined($seen_hash{$seq->display_id}))
		{
			next;
		}
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
			$mini_text .= $seq->display_id() .  "\t" . $num_seen . "\n";
		}	
		
	}

	my $mini_stat = Statistics::Descriptive::Full->new();

	$mini_stat->add_data(@mini_size_array);

	my $miniseen_stat = Statistics::Descriptive::Full->new();
	$miniseen_stat->add_data(@num_array);

	print "MINICIRCLE Statistics\n\n";
	print "Total number of minicircles of size greater then $max_size_minicircle is " . $mini_stat->count . "\n";
	print "Each minicircle is seen a total of " . sprintf("%0.1f", $miniseen_stat->mean) . " times in the dataset for minicircles greater then " . $min_size_mini_dist . "  bp's (" . $miniseen_stat->count . " sequences)\n";
	print "The std deviation of the distribution of sequences seen is " .  sprintf("%0.1f", $miniseen_stat->standard_deviation) . "\n";
	print "The median of the distribution of sequenes seen is " .  sprintf("%0.1f", $miniseen_stat->median) . "\n";
	print "The mode of the distribution of sequenes seen is " . $miniseen_stat->mode . "\n";
	print "\n";
	print "The mean minicircle size is " .  sprintf("%0.1f", $mini_stat->mean) . "\n";
	print "The median minicircle size is " .  sprintf("%0.1f", $mini_stat->median) . "\n";
	print "The mode minicircle size is " .  sprintf("%0.1f", $mini_stat->mode) . "\n";
	print "The stdev of the minicircle size is " .  sprintf("%0.1f", $mini_stat->standard_deviation) . "\n";
	print "This includes data between the values of " . $mini_stat->min . " and " . $mini_stat->max . "\n";

	print "\nThe following is the number of similar sequences found per minicircle:\n\n";
	print $mini_text;


}

my $mini =  Bio::SeqIO->new(   '-file'         => $gcdna_file,
                                '-format'       => "fasta");


my $mini_size_distribution;
my $num_minicircles = 0;
my $num_used_minicircles = 0;
my @mini_size_array;
my $max_size_minicircle = 20;
my $min_size_mini_dist = 20;

my $gcdna_text;
my @num_array;

if($process_gcdna_stats)
{
	while(my $seq = $mini->next_seq)
	{
		$num_minicircles++;
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


	my $mini_new = Bio::SeqIO->new('-file' => $gcdna_file, '-format' => "fasta");

	my @blast_params = (  'program' => 'blastn',
	                    'F'=>'F',
	                   'database' => "db/all_grna",
	                   "E" =>1e-40
	);

	my $blast_factory = Bio::Tools::Run::StandAloneBlast->new(@blast_params);

	my $times_sum = 0;
	my $num_minis = 0;
	$gcdna_text .= "SEQUENCE_NAME\tNUMBER_SEEN\n";

	my %seen_hash;
	while(my $seq = $mini_new->next_seq)
	{
		if(defined($seen_hash{$seq->display_id}))
		{
			next;
		}
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
			$gcdna_text .= $seq->display_id() .  "\t" . $num_seen . "\n";
		}	
		
	}

	my $mini_stat = Statistics::Descriptive::Full->new();

	$mini_stat->add_data(@mini_size_array);

	my $miniseen_stat = Statistics::Descriptive::Full->new();
	$miniseen_stat->add_data(@num_array);

	print "gcDNA Statistics\n\n";
	print "Total number of gcDNA of size greater then $max_size_minicircle is " . $mini_stat->count . "\n";
	print "Each gcDNA is seen a total of " . sprintf("%0.1f", $miniseen_stat->mean) . " times in the dataset for gcDNA greater then " . $min_size_mini_dist . "  bp's (" . $miniseen_stat->count . " sequences)\n";
	print "The std deviation of the distribution of sequences seen is " . sprintf("%0.1f", $miniseen_stat->standard_deviation) . "\n";
	print "The median of the distribution of sequenes seen is " . $miniseen_stat->median . "\n";
	print "The mode of the distribution of sequenes seen is " . $miniseen_stat->mode . "\n";
	print "\n";
	print "The mean gcDNA size is " . sprintf("%0.1f", $mini_stat->mean) . "\n";
	print "The median gcDNA size is " . sprintf("%0.1f", $mini_stat->median) . "\n";
	print "The mode gcDNA size is " . sprintf("%0.1f", $mini_stat->mode) . "\n";
	print "The stdev of the gcDNA size is " . sprintf("%0.1f", $mini_stat->standard_deviation) . "\n";
	print "This includes data between the values of " . $mini_stat->min . " and " . $mini_stat->max . "\n";

	print "\nThe following is the number of similar sequences found per sequenced gcDNA:\n\n";
	print $gcdna_text;

}



my $cdna =  Bio::SeqIO->new(   '-file'         => $cdna_file,
                                '-format'       => "fasta");

my $cdna_size_distribution;
my $num_cdna = 0;
my $num_used_cdna = 0;
my @cdna_size_array;
my $max_size_cdna = 150;


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

	$cdna_stat->add_data(@cdna_size_array);
	
	print "cDNA Statistics\n\n";
	print "Total number of cdna's of size greater then $max_size_cdna is " . $cdna_stat->count . "\n";
	print "The mean cdna size is " . sprintf("%0.1f", $cdna_stat->mean) . "\n";
	print "The median cdna size is " . sprintf("%0.1f", $cdna_stat->median) . "\n";
	print "The mode cdna size is " . sprintf("%0.1f", $cdna_stat->mode) . "\n";
	print "The stdev of the cdna size is " . sprintf("%0.1f", $cdna_stat->standard_deviation) . "\n";
	print "This includes data between the values of " . $cdna_stat->min . " and " . $cdna_stat->max . "\n";	

}

