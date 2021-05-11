#!/usr/bin/perl

# The purpose of this script is to search for gRNA's that guide for more then one gene.
# Written by Michael J. Cipriano, 2005

use Bio::DB::GFF;
use Getopt::Long;
use Bio::Tools::Run::StandAloneBlast;
use File::Temp qw/tempfile tempdir/;
use Bio::SeqIO;

use strict;
my $db_name;
my $gene_name_opt;
my $verbose = 0;

my $options = GetOptions(
                                "db_name=s", \$db_name,
				"gene_name=s", \$gene_name_opt,
				"verbose!", \$verbose
);


if( !defined($db_name) )
{
	print "

Options are:
--db_name	The name of the Bio::DB::GFF database to query
--gene_name	The name of the gene you are using in your search (all to search all genes)
--verbose	Prints out debugging information

Example: grna_find_multipurpose.pl --db_name=tbgrna --gene_name=Tbco3 --verbose 
";
exit;
}


my $db = Bio::DB::GFF->new(     -adaptor => 'dbi::mysqlopt',
                                -dsn    => "dbi:mysql:$db_name:gmoddb",
                                -user   => 'gid',
                                -password       => 'gidgid123'
                        );


if(!$db)
{
	die(" No database $db_name is found!");
}

if($verbose)
{
	print "CONNECTED: $db_name.\n";
}

$db->absolute(1);




# First create a hash of all of the transcripts and what gene they hit
my $transcript_to_gene;

my @genes = $db->features("gene:pcr");

foreach my $gene(@genes)
{
	my $gene_name = $gene->name;
	
	if($gene->name =~ /un$/)
	{
		($gene_name) = $gene->name =~ /(\w+)un$/;
	} elsif($gene->name =~ /ed$/)
	{
		($gene_name) = $gene->name =~ /(\w+)ed$/;
	} 

	$transcript_to_gene->{$gene->name} = $gene_name;

	if($verbose)
	{
		print $gene->name . "\n";
	}
	my $gene_segment = $db->segment($gene->name);
#	my @transcripts = $gene_segment->features('specific_bloodstream', 'specific_procyclic', 'dt_primed');
	my @transcripts = $gene_segment->features('specific_bloodstream', 'specific_procyclic', 'dt_primed', 'specific_bloodstream_fraction_11', 'specific_bloodstream_fraction_9', 'specific_procyclic_fraction_11', 'specific_procyclic_fraction_9');

	foreach my $transcript (@transcripts)
	{
		my $est_name = $transcript->attributes('est_name');
		$transcript_to_gene->{$est_name} = $gene_name;
	}
}

# Get all minicircle sequences

if($verbose)
{
	print "Searching for minicircle sequences\n";
}

my @minicircles = $db->features("minicircle");

my $printed_hash;
print join("\t", 'a_minicircle_name', 'a_mini_start', 'a_mini_end', 'a_transcript', 'a_transcript_start', 'a_transcript_end', 'a_mini_seq', 'a_gene', 'b_minicircle_name', 'b_mini_start', 'b_mini_end', 'b_transcript', 'b_transcript_start', 'b_transcript_end', 'b_mini_seq', 'b_gene') . "\n";
foreach my $minicircle (@minicircles)
{
	if($verbose)
	{
		print "Minicircle:" . $minicircle->name . "\n";
	}
	# Pull out the minicircle segment
	my $minicircle_segment = $db->segment(-name=>$minicircle->name, -class=>'Sequence');

	# Now check this minicircle for guideRNA's (mini_match)
	my @grnas = $minicircle_segment->features("mini_match:blast");
	foreach my $grna(@grnas)
	{
		
		if($verbose)
		{
			print " grna:" . $grna->name . "\n";
		}
		# Now check what transcript this guide RNA hits and what gene it is.
		my $grna_gene = $transcript_to_gene->{$grna->name};
		if($verbose)
		{
			print "  Gene:$grna_gene\n";
		}

		# Now check for overlapping guideRNA's
		my @overlapping_grnas = $grna->contained_features('mini_match:blast');
		foreach my $overlapping_grna(@overlapping_grnas)
		{

			# Now check this grna's parent transcript
			my $overlapping_grna_gene = $transcript_to_gene->{$overlapping_grna->name};

			# How much do they overlap
			my $overlap = 0;
			if( ($grna->start <= $overlapping_grna->start) && ($grna->end >= $overlapping_grna->end) )
			{
				$overlap = $overlapping_grna->end - $overlapping_grna->start + 1;
			} elsif( ($grna->start <= $overlapping_grna->start) && ($grna->end <= $overlapping_grna->end) )
			{
				$overlap = $grna->end - $overlapping_grna->start + 1;
			} elsif( ($grna->start >= $overlapping_grna->start) && ($grna->end <= $overlapping_grna->end) )
			{
				$overlap = $grna->end - $grna->start + 1;
			} elsif( ($grna->start >= $overlapping_grna->start) && ($grna->end >= $overlapping_grna->end) )
			{
				$overlap = $overlapping_grna->end - $grna->start + 1;
			}
			my $grna_length = $grna->end - $grna->start +1;
			my $perc_overlap = ($overlap / $grna_length) * 100;

			if( $overlapping_grna_gene ne $grna_gene && $perc_overlap > 70 && $overlap >= 20)
			{
				if($verbose)
				{
					print "  overlap:" . $overlapping_grna->name . " ($overlapping_grna_gene) percent:$perc_overlap\n";
				}
				my $a_name = join('-', $minicircle->name, $grna->name, $grna->start, $grna->stop);
				my $b_name = join('-', $minicircle->name,$overlapping_grna->name, $overlapping_grna->start, $overlapping_grna->stop);
				if( !defined($grna_gene) || !defined($overlapping_grna_gene) || $grna_gene eq '' || $overlapping_grna_gene eq '' )
				{ 
					if($verbose)
					{
						print "  undefined gene!\n";
					}
					# Do nothing, this is not a gene we want to print
				} elsif( $printed_hash->{$a_name} eq $b_name || $printed_hash->{$b_name} eq $a_name)
				{
					# Do nothing, we already printed this
				} else
				{
					
					if( defined($gene_name_opt) && ($gene_name_opt =~ /($grna_gene|$overlapping_grna_gene)/i) || !defined($gene_name_opt))
					{
						print join("\t", $minicircle->name, $grna->start, $grna->stop, $grna->name, $grna->attributes('query_start'), $grna->attributes('query_stop'), $grna->seq, $grna_gene, $minicircle->name, $overlapping_grna->start, $overlapping_grna->stop, $overlapping_grna->name, $overlapping_grna->attributes('query_start'), $overlapping_grna->attributes('query_stop'), $overlapping_grna->seq, $overlapping_grna_gene ) . "\n";
						$printed_hash->{$a_name} = $b_name;
						$printed_hash->{$b_name} = $a_name;
					} else
					{
						if($verbose)
						{
							print " SKIP!\n";
						}
					}
				}
			}
			
		}

	}
	if($verbose)
	{
		print "\n\n";
	}

}

