#!/usr/bin/perl


use strict;


use Bio::Tools::Run::StandAloneBlast;
use Bio::SeqIO;
use Getopt::Long;
use File::Temp qw/tempfile tempdir/;


my $db_name;
my $verbose;
my $cdna_file;


my $options = GetOptions(
				"db_name=s", \$db_name,
				"cdna_file=s", \$cdna_file,
				"verbose!", \$verbose
		);

if( !defined($db_name) && !defined($cdna_file) )
{

	print "
--db_name
--cdna_file
--verbose

";
exit;
}

my %result_hash;

my $seqio = Bio::SeqIO->new( -file=>$cdna_file, -format=>'fasta');

my $wu_blastn_bin = "blastn";
# -hpmax = maximum hsp's to report, -gspmax = maximum gapped hsp's to report
# -W 2 = word size, -Q = cost to open a gap, -R = cost to extend a gap
my $wu_blast_options = "-hspmax 0 -gspmax 0 -span2 -W 2 -matrix goldenoxduk-one -warnings -notes";
#my $wu_blast_options = "-nogap -w 2 -matrix goldenoxduk-one -warnings -notes";
my ($seq_fh, $seq_fn) = tempfile();
my ($br_fh, $br_fn) = tempfile();

while(my $seq = $seqio->next_seq)
{

	# write seq to file
	my $seqout = Bio::SeqIO->new(-file=>">$seq_fn", -format=>'fasta');
	$seqout->write_seq($seq);
	# For this cdna, blast it with wublast to the db_name database

	system("$wu_blastn_bin $wu_blast_options -d $db_name -o $br_fn -i $seq_fn");
	my $blast_report = Bio::SearchIO->new(-format=>'blast', -file=>$br_fn);
	while(my $result = $blast_report->next_result)
	{
		while(my $hit = $result->next_hit)
		{
			while(my $hsp = $hit->next_hsp)
			{
				# Check for the longest match (allowing both + and |)
				my $homology_string = $hsp->homology_string;
				if($verbose)
				{
					print 	$seq->display_id .  " [" . $hsp->query->start . "-" . $hsp->query->end . "]" . "\n" 
						. $hit->name . " [" . $hsp->hit->start . "-" . $hsp->hit->end . "]\n"   ;
					print $hsp->query_string . "\n" . $homology_string . "\n" . $hsp->hit_string . "\n";
				}
				my @bhits = $homology_string =~ /[+|]+/g;
				my $max = 0;
				foreach my $bhit (@bhits)
				{
					my $len_hit = length($bhit);
					if($len_hit > $max)
					{
						$max = $len_hit;
					}
				}	
				$result_hash{$max}++;
				if($verbose)
				{
					print "MAX:$max\n";
				}
			}
		}
	}	


}

my $total = 0;
for(1..200)
{
	my $num = $_;
	$total += $result_hash{$num};
	print "$num\t" . $result_hash{$num} . "\t" . $total . "\n";

}


