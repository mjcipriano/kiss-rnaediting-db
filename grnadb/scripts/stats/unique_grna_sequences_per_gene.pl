#!/usr/bin/perl

# Written by Michael J. Cipriano, 2005

use Bio::DB::GFF;
use Getopt::Long;
use Bio::Tools::Run::StandAloneBlast;
use File::Temp qw/tempfile tempdir/;
use Bio::SeqIO;

use strict;
my $db_name;
my $gene_name;
my $verbose = 0;
my $outfile = 'report.txt';
my $outdir = ".";
my $webdir = ".";
my $predicted = 0;
my $sequenced = 0;

my $options = GetOptions(
                                "db_name=s", \$db_name,
                                "gene_name=s", \$gene_name,
                                "verbose!", \$verbose,
                                "outfile=s", \$outfile,
				"outdir=s", \$outdir,
				"webdir=s", \$webdir,
				"predicted", \$predicted,
				"sequenced", \$sequenced
);


my $db = Bio::DB::GFF->new(     -adaptor => 'dbi::mysqlopt',
                                -dsn    => "dbi:mysql:$db_name:gmoddb",
                                -user   => 'gid',
                                -password       => 'gidgid123'
                        );
$db->absolute(1);
my @genes = $db->features("gene:pcr");
my @search_list;

open(REPORT, ">$outdir/$outfile");

if($predicted)
{
	push(@search_list, 'match:pcr');
}
if($sequenced)
{
	push(@search_list, 'match:grna-pcr');
}

my $transcript_to_gene;
my $gene_grnas;

my $grna_coverage_hash;

foreach my $gene(@genes)
{
        my $gene_name = $gene;
        my ($gene_orig) = $gene =~ /\((\w+)\)$/;

	my $grna_seq;
	

        if($gene =~ /un\)$/)
        {
                ($gene_name) = $gene =~ /\((\w+)un\)$/;
        } elsif($gene_name =~ /ed\)$/)
        {
                ($gene_name) = $gene =~ /\((\w+)ed\)$/;
        } else
        {
                ($gene_name) = $gene =~ /\((\w+)\)$/;
        }

        $transcript_to_gene->{$gene_orig} = $gene_name;

        if($verbose)
        {
                print $gene->name . "\n";
        }
        my $gene_segment = $db->segment($gene->name);
#       my @transcripts = $gene_segment->features('specific_bloodstream', 'specific_procyclic', 'dt_primed');
        my @transcripts = $gene_segment->features('specific_bloodstream', 'specific_procyclic', 'dt_primed', 'specific_bloodstream_fraction_11', 'specific_bloodstream_fraction_9', 'specific_procyclic_fraction_11', 'specific_procyclic_fraction_9');

	my @grnas = $gene_segment->features(@search_list);
	foreach my $g (@grnas)
	{
		my ($seq) = $g->attributes('query_string');
		$gene_grnas->{$gene_name}->{$seq} = $g->name;
	}
	# Find the coverage if it is a edited gene
	if($gene->name =~ /ed$/)
	{
		my $len = $gene_segment->length;
		my @overlap_array;
		for(1..$len)
		{
			$overlap_array[$_] = 0;
		}
		foreach my $g(@grnas)
		{
			for($g->start..$g->stop)
			{
				$overlap_array[$_]++;
			}
		}
		# Now count the number of locations > 0;
		my $num_cover = 0;
		for(1..$len)
		{
			if($overlap_array[$_] > 0)
			{
				$num_cover++;
			}
		}
		$grna_coverage_hash->{$gene_name} = ($num_cover/$len) * 100;
	}
	
	
        foreach my $transcript (@transcripts)
        {
                my $est_name = $transcript->attributes('est_name');
		if($verbose)
		{
			print " " . $est_name . "\n";
		}

		my $transcript_segment = $db->segment($est_name);
		my @grnas = $transcript_segment->features(@search_list);
		foreach my $g (@grnas)
		{
			my ($seq) = $g->attributes('query_string');
			$gene_grnas->{$gene_name}->{$seq} = $g->name;
		}
                $transcript_to_gene->{$est_name} = $gene_name;
        }
}

# Now for each gene and grna set, now just count them

print REPORT 
'
gRNA Report

The following are the number of unique gRNAs found for each gene.  For each gene, 2 statistics are generated, Unique and Substring.
The Substring statistic will call 2 gRNAs the same if one is a substring of the other.
The Unique statistic will only call 2 gRNAs the same if the full length of the sequence is the same.

You may download the dataset of each of the genes at the locations listed below.

';

while(my ($gene, $grna_hash) = each(%$gene_grnas))
{
	print REPORT "#############################################\n";
	print REPORT "Gene: $gene\n";
	print REPORT "-----------\n";


	my @count = keys %$grna_hash;
	print REPORT "Unique gRNAs: " . scalar @count . "\n";
	print REPORT "Available at http://gmod.mbl.edu/$webdir/unique_$gene.fasta\n\n";
	
	my ($temp_fh, $temp_fn) = tempfile(UNLINK=>1);
	my $io = Bio::SeqIO->new(-file=>">$temp_fn", -format=>'fasta');
	while(my ($seq_val, $name) = each(%$grna_hash))
	{
		my $seqobj = Bio::Seq->new(-display_id=>$name, -seq=>$seq_val);
		$io->write_seq($seqobj);
	
	}
	
	print join("\t", $gene, scalar @count) . "\n";
	system("patdb -o $outdir/unique_substring_$gene.fasta -c '|' -t n -s 10 -r /dev/null $temp_fn");
	my $sub_count = `grep -c '>' $outdir/unique_substring_$gene.fasta`;
	chomp($sub_count);
	print REPORT "Unique Substring gRNAs: $sub_count\n";
	print REPORT "Available at http://gmod.mbl.edu/$webdir/unique_substring_$gene.fasta\n\n";

	if(defined($grna_coverage_hash->{$gene}))
	{
		print REPORT "Coverage is " . sprintf("%.2f", $grna_coverage_hash->{$gene}) . "%\n";
	} else
	{
		print REPORT "Coverage is not available for non-edited genes\n";
	}
	system("cp $temp_fn $outdir/unique_$gene.fasta");
	close($temp_fh);
	unlink($temp_fn);
}



