#!/usr/bin/perl
# This script takes in a sequence file in fasta format and prints to the screen a fasta file with the largest non-X portion of the sequence
use Bio::SeqIO;

use strict;

my $seq_file = $ARGV[0];

my $seqs =  Bio::SeqIO->new(   '-file'         => $seq_file,
                                '-format'       => "fasta");
 
while(my $seq = $seqs->next_seq)
{

        # Trim the sequence and find the largest non-X portion

        my @parts = split("X", $seq->seq());
        my $max_size = 0;
        my $max_seq = "";
        foreach my $seq_part(@parts)
        {
                if(length($seq_part) > length($max_seq))
                {
                        $max_seq = $seq_part;
                }
        }

	print ">" . $seq->display_id . "\n" . $max_seq . "\n";
}

