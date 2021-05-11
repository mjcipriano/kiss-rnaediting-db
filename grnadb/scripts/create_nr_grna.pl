#!/usr/bin/perl

# Written by Michael J. Cipriano, 2005

use Bio::DB::GFF;
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
	my $io = Bio::SeqIO->new(-file=>">$temp_fn", -format=>'fasta');

	if($verbose)
	{
		print $minicircle->name . "\n";
	}
        my $mini_segment = $db->segment($minicircle->name);

	my @grnas = $mini_segment->features(@search_list);
	foreach my $g (@grnas)
	{
		my $seq = $g->seq();
		if($verbose)
		{
			print " " . join("\t", $g->name, $g->start, $g->end, $g->strand) . "\n";
			print $mini_segment->subseq($g->start, $g->end)->seq()->seq() . "\n";
			print $seq->seq() . "\n";
		}
		$io->write_seq($seq);
	}
	# Now that I have all of these, I want to create a non-redundant set of substrings

	system("patdb -o $nrtemp_fn -c '|' -t n -s 10 -r /dev/null $temp_fn");
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
	$num++;
	$nrseq->description("");
	$nrseq->display_name($name);
	$seqfinal->write_seq($nrseq);
}


