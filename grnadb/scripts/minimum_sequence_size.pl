#!/usr/bin/perl

use Bio::SeqIO;

use strict;

my $infile = $ARGV[0];
my $outfile = $ARGV[1];
my $min_len = $ARGV[2];

my $in = Bio::SeqIO->new(-format=>'fasta', -file=>$infile);

my $out = Bio::SeqIO->new(-format=>'fasta', -file=>">$outfile");

while(my $seq = $in->next_seq)
{
	if($seq->length >= $min_len)
	{
		$out->write_seq($seq);
	}
}
