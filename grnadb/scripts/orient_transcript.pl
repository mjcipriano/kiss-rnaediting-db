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

my $options = GetOptions(
				"infile=s", \$infile,
				"outfile=s", \$outfile,
				"verbose!", \$verbose,
				"sequence-db=s", \$sequence_db,
				"nomatch=s", \$nomatch,
				"nomatch-file=s", \$nomatch_file
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
--nomatch		What to do if no match is found (default forward)
			  forward - Print the sequence in the forward direction
			  reverse - Print the sequence in the reverse direction
			  none    - Do not print the sequence
			  file    - Print out the non-matching sequences to a separate file
			            defined by the nomatch-file option.
--nomatch-file		The name of the file to print to if the nomatch=file option is chosen (default:nomatch.fasta)
--verbose		Use this flag to see debugging info


Example: orient.pl --infile=myfile.fasta --outfile=result.fasta --sequence_db=nt --verbose
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

	if($sequence_db)
	{
		# Open the sequence fasta file

		my $blast_factory = Bio::Tools::Run::StandAloneBlast->new('database'=>$sequence_db, 'program'=>'blastn');
		my $report = $blast_factory->blastall($seq);
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




