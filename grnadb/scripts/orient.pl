#!/usr/bin/perl

use strict;


use Bio::Seq;
use Bio::SeqIO;
use Bio::Tools::Run::StandAloneBlast;
use Getopt::Long;

my $infile;
my $outfile;
my $verbose = 0;
my $motif;
my $sequence;
my $nomatch = 'forward';
my $nomatch_file = 'nomatch.fasta';

my $options = GetOptions(
				"infile=s", \$infile,
				"outfile=s", \$outfile,
				"verbose!", \$verbose,
				"motif=s", \$motif,
				"sequence=s", \$sequence,
				"nomatch=s", \$nomatch,
				"nomatch-file=s", \$nomatch_file
);


if($verbose)
{
	print "
Options:
	infile		$infile
	outfile		$outfile
	motif		$motif
	sequence	$sequence
	nomatch		$nomatch
	nomatch-file	$nomatch_file
	verbose		$verbose

";

}


if(!defined($infile) || !defined($outfile) || (!defined($motif) && !defined($sequence)) )
{

	print "
Options not given, please supply the following options

--infile		The fasta formated input file
--outfile		The name of the new file to create with the processed infile's results
--motif			A string representing the motif to find in a specific orientation
--sequence		A sequence file (fasta) which will be blasted to each sequence, and the top hit will determine direction
--nomatch		What to do if no match is found (default forward)
			  forward - Print the sequence in the forward direction
			  reverse - Print the sequence in the reverse direction
			  none    - Do not print the sequence
			  file    - Print out the non-matching sequences to a separate file
			            defined by the nomatch-file option.
--nomatch-file		The name of the file to print to if the nomatch=file option is chosen (default:nomatch.fasta)
--verbose		Use this flag to see debugging info


Example: orient.pl --infile=myfile.fasta --outfile=result.fasta --motif=\"AATAG\" --verbose
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

while(my $seq = $inseqs->next_seq)
{

	if($verbose)
	{
		print "Searching " . $seq->display_id . "\n";
	}


	my $rev_seq = $seq->revcom;

	if($motif)
	{
		# Search for the motif given in this sequence
		if( $seq->seq =~ /$motif/i)
		{
			# The motif was found, continue on to the next sequence and print this one out in the forward orientation
			$outseqs->write_seq($seq);
			if($verbose)
			{
				print " Motif found in forward orientation.\n";
			}
		} elsif($rev_seq->seq =~ /$motif/i)
		{
			# If the motif was not found in the forward orientation, check the reverse
			$outseqs->write_seq($rev_seq);
			if($verbose)
			{
				print " Motif found in the reverse orientation.\n";
			}
		} else
		{
			# The motif was not found on this sequence
			if($nomatch eq 'forward')
			{
				$outseqs->write_seq($seq);
			} elsif($nomatch eq 'reverse')
			{
				$outseqs->write_seq($rev_seq);
			} elsif($nomatch eq 'none')
			{
				# Do nothing
			} elsif($nomatch eq 'file')
			{
				$nomatch_seqs->write_seq($seq);
			}

			if($verbose)
			{
				print " Motif not found!\n";
			}

		}

	}

	if($sequence)
	{
		# Open the sequence fasta file
		my $match_seqio = Bio::SeqIO->new(-file=>$sequence, -format=>'fasta');
	
		my $match_seq = $match_seqio->next_seq or die ("No sequence found in file $sequence.\n");

		my $bl2seq_factory = Bio::Tools::Run::StandAloneBlast->new('program'=>'blastn');
		my $report = $bl2seq_factory->bl2seq($match_seq, $seq);
		my $result = $report->next_result;
		if(my $hit = $result->next_hit)
		{
			if(my $hsp = $hit->next_hsp)
			{
				if($hsp->hit->strand == -1)
				{
					$outseqs->write_seq($rev_seq);
					if($verbose)
					{
						print " Sequence found in the reverse orientation.\n";
						print "  " . $hsp->hit->start . " - " . $hsp->hit->end . "\n";
					}
				} else
				{
					$outseqs->write_seq($seq);
					if($verbose)
					{
						print " Sequence found in forward orientation.\n";
						print "  " . $hsp->hit->start . " - " . $hsp->hit->end . "\n";
					}
				}
			}
		} else
		{
			if($nomatch eq 'forward')
			{
				$outseqs->write_seq($seq);
			} elsif($nomatch eq 'reverse')
			{
				$outseqs->write_seq($rev_seq);
			} elsif($nomatch eq 'none')
			{
				# Do nothing
			} elsif($nomatch eq 'file')
			{
				$nomatch_seqs->write_seq($seq);
			}

			if($verbose)
			{
				print " Sequence not found!\n";
			}
		}
		
	}
}




