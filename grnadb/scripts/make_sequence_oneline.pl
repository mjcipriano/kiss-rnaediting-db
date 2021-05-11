#!/usr/bin/perl
# This script takes out all newlines and makes the given sequence into a one line sequence. Fasta is the expected input. Reads from stdin and writes to stdout.

use strict;


my $current_sequence;
my $sequence_name = "";

while(<>)
{
	my $line = $_;
	chomp($line);
	if($line =~ /^$/)
	{
		# Do nothing

	} elsif($line =~ /^\>/)
	{
		# write out last sequence
		if($sequence_name ne "")
		{
			print ">" . $sequence_name . "\n" . $current_sequence . "\n";
			$current_sequence = "";
		}

		($sequence_name) = $line =~ /(\w+)/;
	}else
	{
		$current_sequence .= $line;
	}


}
