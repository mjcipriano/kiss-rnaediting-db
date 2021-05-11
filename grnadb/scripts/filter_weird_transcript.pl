#!/usr/bin/perl

use strict;


use Bio::Seq;
use Bio::SeqIO;
use Bio::Tools::Run::StandAloneBlast;
use Getopt::Long;

my $infile;
my $outfile;
my $verbose = 0;
my $sequence_db;
my $nomatch = 'forward';
my $nomatch_file = 'nomatch.fasta';
# Percentage of the gene that must be covered by pre-edited and edited transcripts
my $perc_cover = 70;

# The number of bases that must match at the begining and end of the transcript with 
my $bases_near_end = 100;

my $options = GetOptions(
				"infile=s", \$infile,
				"outfile=s", \$outfile,
				"verbose!", \$verbose,
				"sequence-db=s", \$sequence_db,
				"nomatch=s", \$nomatch,
				"nomatch-file=s", \$nomatch_file,
				"perc_cover=s", \$perc_cover
);


if($verbose)
{
	print "
Options:
	infile		$infile
	outfile		$outfile
	sequence-db	$sequence_db
	nomatch		$nomatch
	nomatch-file	$nomatch_file
	verbose		$verbose

";

}


if(!defined($infile) || !defined($outfile) ||  !defined($sequence_db) )
{

	print "
Options not given, please supply the following options

--infile		The fasta formated input file
--outfile		The name of the new file to create with the processed infile's results
--sequence-db		A blast database to compare the orientation of the input sequences to.
--perc_cover		The percentage that the sequence must be covered to be included (default 70).
--nomatch		What to do if no match is found (default forward)
			  forward - Print the sequence in the forward direction
			  reverse - Print the sequence in the reverse direction
			  none    - Do not print the sequence
			  file    - Print out the non-matching sequences to a separate file
			            defined by the nomatch-file option.
--nomatch-file		The name of the file to print to if the nomatch=file option is chosen (default:nomatch.fasta)
--verbose		Use this flag to see debugging info


Example: filter_weird_transcript.pl --infile=myfile.fasta --outfile=result.fasta --sequence_db=nt --perc_cover=60 --verbose
";
exit;
}



my $inseqs = Bio::SeqIO->new( 	-file=>$infile,
				-format=>'fasta');

my $outseqs = Bio::SeqIO->new(	-file=>">$outfile",
				-format=>'fasta');

my $nomatch_seqs;
if($nomatch eq 'file')
{
	$nomatch_seqs = Bio::SeqIO->new( -file=>">$nomatch_file", -format=>'fasta');
}


my $genes =  Bio::SeqIO->new(   '-file'         => $sequence_db,
                                '-format'       => "fasta");
my $gene_hash;

while(my $gene = $genes->next_seq)
{
	$gene_hash->{$gene->display_id} = $gene;
}


while(my $seq = $inseqs->next_seq)
{
	if($verbose)
	{
		print "Searching " . $seq->display_id . "\n";
	}


	if($sequence_db)
	{
		# Open the sequence fasta file

		my $blast_factory = Bio::Tools::Run::StandAloneBlast->new('database'=>$sequence_db, 'program'=>'blastn');
		my $report = $blast_factory->blastall($seq);
		my $result = $report->next_result;

		if(my $hit = $result->next_hit)
		{
			# Find out what this hits to
			my $hsp = $hit->next_hsp;
			my $unedited;
	                my $edited;
	                my $single_gene = 0;
			my $gene_name;
			
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
				$un_hash->{query_start} = $seq_start;
				$un_hash->{query_stop} = $seq_stop;
	                        $un_hash->{name} = $seq->display_id;
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
				$ed_hash->{query_start} = $seq_start;
				$ed_hash->{query_stop} = $seq_stop;
	                        $ed_hash->{name} = $seq->display_id;
			}
			if($verbose)
			{
				print join("\t", $seq->display_id, $gene_name, 'pre-edit', $un_hash->{start} . '..' . $un_hash->{stop}, $un_hash->{query_start} . '..' . $un_hash->{query_stop}, 'edited', $ed_hash->{start} . '..' . $ed_hash->{stop}, $ed_hash->{query_start} . '..' . $ed_hash->{query_stop}) . "\n";
			}

			# Find out how much of the sequence is marked as edited or pre-edited
			my $seq_size = $seq->length();
			# Initialize an array of length size
			my @seq_array;
			for(1..$seq_size)
			{
				$seq_array[$_] = 0;
			}

			# Now move from start to stop and mark it as covered
			for($un_hash->{query_start}..$un_hash->{query_stop})
			{
				$seq_array[$_] += 1;
			}
			if(!$single_gene)
			{
				for($ed_hash->{query_start}..$ed_hash->{query_stop})
				{
					$seq_array[$_] += 2;
				}
			}
			# Now find out what percentage is covered
			my $covered = 0;
			for(1..$seq_size)
			{
				if($seq_array[$_] > 0)
				{
					$covered++;
				}
			}
			my $perc_covered = ($covered/$seq_size)*100;
			if($verbose)
			{
				print $perc_covered . "\n";
			}
			# Now check that the ends have pre-edited and edited
			my $fiveprime_end = 0;
			my $threeprime_end = 0;
			for(1..$bases_near_end)
			{
				if($seq_array[$_] > 0)
				{
					$fiveprime_end = 1;
				}
				
			}
			for($seq_size-$bases_near_end..$seq_size)
			{
				if($seq_array[$_] > 0)
				{
					$threeprime_end = 1;
				}
			}

			if( ($perc_covered > $perc_cover) && $threeprime_end && $fiveprime_end)
			{
				$outseqs->write_seq($seq)
			} else
			{
				if($verbose)
				{
					print "sequence did not pass!\n";
					print "5' End: $fiveprime_end\n3' End: $threeprime_end\nCoverage: $perc_covered\n";
				}
				if($nomatch eq 'file')
				{
					$nomatch_seqs->write_seq($seq);
				}
			}

		} else
		{
			if($verbose)
			{
				print "sequence match not found!\n";
			}
			if($nomatch eq 'file')
			{
				$nomatch_seqs->write_seq($seq);
			}
		}
	}
}

