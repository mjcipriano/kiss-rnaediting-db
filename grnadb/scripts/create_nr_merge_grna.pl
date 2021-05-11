#!/usr/bin/perl

# Written by Michael J. Cipriano, 2005

use Bio::DB::GFF;
use Bio::SeqFeature::Generic;
use Getopt::Long;
use Bio::Tools::Run::StandAloneBlast;
use File::Temp qw/tempfile tempdir/;
use Bio::SeqIO;

use strict;
my $db_name;
my $verbose = 0;
my $webdir = ".";
my $predicted = 0;
my $sequenced = 0;
my $file = 'nr.fasta';

my $options = GetOptions(
                                "db_name=s", \$db_name,
                                "verbose!", \$verbose,
				"webdir=s", \$webdir,
				"file=s", \$file,
				"predicted", \$predicted,
				"sequenced", \$sequenced
);


my $db = Bio::DB::GFF->new(     -adaptor => 'dbi::mysqlopt',
                                -dsn    => "dbi:mysql:$db_name:gmoddb",
                                -user   => 'gid',
                                -password       => 'gidgid123'
                        );
$db->absolute(1);
my @minicircles = ();

my @search_list;

if($predicted)
{
	push(@search_list, 'mini_match:blast');
	push(@minicircles, $db->features("minicircle:pcr"));
}
if($sequenced)
{
	push(@search_list, 'mini_match:grna-blast');
	push(@minicircles, $db->features("grna:pcr"));
}


my ($all_fh, $all_fn) = tempfile(UNLINK=>1);
my ($nrtemp_fh, $nrtemp_fn) = tempfile(UNLINK=>1);
my $allio = Bio::SeqIO->new(-file=>">$all_fn", -format=>'fasta');

foreach my $minicircle(@minicircles)
{
	my ($temp_fh, $temp_fn) = tempfile(UNLINK=>1);

	if($verbose)
	{
		print $minicircle->name . "\n";
	}
        my $mini_segment = $db->segment($minicircle->name);

	my @grnas = $mini_segment->features(@search_list);
	my @forward_hits = ();
	my @reverse_hits = ();

 	
	foreach my $g (@grnas)
	{
		my $seq = $g->seq();
		if($g->strand > 0)
		{
			push(@forward_hits, $g);
		} else
		{
			push(@reverse_hits, $g);
		}
		if($verbose)
		{
			print " " . join("\t", $g->name, $g->start, $g->end, $g->strand) . "\n";
			#print $mini_segment->subseq($g->start, $g->end)->seq()->seq() . "\n";
			#print $seq . "\n";
		}
	}

	# First do forward
	my @merged_hits = ();
	my %visited_hash;
	foreach my $g (@reverse_hits)
	{
		if(defined($visited_hash{$g}))
		{
			next;
		}
		$visited_hash{$g} = 1;
		foreach my $i (@reverse_hits)
		{
			if(defined($visited_hash{$i}))
			{
				next;
			}
			my $new = check_and_merge($g, $i, 5);
			if(defined($new))
			{
				$visited_hash{$i} = 1;
				$g = $new;
			}
		}
		push(@merged_hits, $g);
	}

	foreach my $g (@forward_hits)
	{
		if(defined($visited_hash{$g}))
		{
			next;
		}
		$visited_hash{$g} = 1;
		foreach my $i (@forward_hits)
		{
			if(defined($visited_hash{$i}))
			{
				next;
			}
			my $new = check_and_merge($g, $i, 5);
			if(defined($new))
			{
				$visited_hash{$i} = 1;
				$g = $new;
			}
		}
		push(@merged_hits, $g);
	}


	foreach my $hit (@merged_hits)
	{
		if($verbose)
		{
			print join("\t", $hit->start, $hit->end, $hit->strand) . "\n";
		}
		if($hit->strand > 0)
		{
			my $hit_sequence_obj =  $mini_segment->subseq($hit->start, $hit->end);
			$allio->write_seq($hit_sequence_obj);
		} else
		{
			my $hit_sequence_obj =  $mini_segment->subseq($hit->start, $hit->end)->revcom;
			$allio->write_seq($hit_sequence_obj);
		}
	}

	my $nrseqs = Bio::SeqIO->new(-file=>$nrtemp_fn, -format=>'fasta');
	while(my $nrseq = $nrseqs->next_seq)
	{
		$allio->write_seq($nrseq);
	}

}

# Now get complete matches and iterate through these and write out the file.

system("nrdb -o $nrtemp_fn -l 10 -d '|' $all_fn");

my $seqfinal = Bio::SeqIO->new(-file=>">$file", -format=>'fasta');
my $num = 1;

my $nrio = Bio::SeqIO->new(-file=>$nrtemp_fn, -format=>'fasta');
while(my $nrseq = $nrio->next_seq)
{
	my $name = 'KISSTBB' . sprintf("%06d", $num);
	if($sequenced)
	{
		$name .= 'S';
	}
	$num++;
	$nrseq->description("");
	$nrseq->display_name($name);
	$seqfinal->write_seq($nrseq);
}



sub check_and_merge
{
	my $a = shift;
	my $b = shift;
	my $overlap = shift;

	my $a_length = $a->length;
	if($a->strand ne $b->strand)
	{
		return undef;
	}
	if($a->start >= $b->start)
	{
		
		if($b->end >= $a->end)
		{
			# A is internal to B
			return $b;
		} elsif($a->start <= $b->end)
		{
			# They overlap with b leading a
			my $overlap_len = $b->end - $a->start +1;
			if($overlap_len >= $overlap)
			{
				my $seqfeat = Bio::SeqFeature::Generic->new(-start=>$b->start, -end=>$a->end, -strand=>$b->strand);
				return $seqfeat;
			}
		}
	} elsif($b->start >= $a->start)
	{
		if($a->end >= $b->end)
		{
			return $a;
		} elsif($b->start <= $a->end)
		{
			# They overlap with b leading a
			my $overlap_len = $a->end - $b->start +1;
			if($overlap_len >= $overlap)
			{
				my $seqfeat = Bio::SeqFeature::Generic->new(-start=>$a->start, -end=>$b->end, -strand=>$a->strand);
				return $seqfeat;
			}
		
		}
	}
	return undef;
}

