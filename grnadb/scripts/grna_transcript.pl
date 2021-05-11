#!/usr/bin/perl

# The purpose of this script is to search through transcripts for gRNA's that appear on intermediate transcripts of a fully edited gene, but do not appear on the fully edited gene themselves. Emphisis is placed on searching the junction region. 
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
my $outfile;

my $options = GetOptions(
                                "db_name=s", \$db_name,
				"gene_name=s", \$gene_name,
				"verbose!", \$verbose,
				"outfile=s", \$outfile
);


if( !defined($db_name) )
{
	print "
The purpose of this script is to search through transcripts for gRNA's that appear 
on intermediate transcripts of a fully edited gene, but do not appear on the fully 
edited gene themselves. Emphisis is placed on searching the junction region.

Options are:
--db_name	The name of the Bio::DB::GFF database to query
--gene_name	The name of the gene you are using in your search
--verbose	Prints out debugging information
--outfile	Prints the results to a tab dilimited file

Example: grna_transcript.pl --db_name=tbgrnadb --gene_name=Tbco3 --verbose 
";
exit;
}

my $print_out = 0;
if(defined($outfile))
{
	open(OUTFILE, ">", $outfile);
	$print_out = 1;

	# Print the header of the file
	print OUTFILE join("\t", "Transcript name", "pre-edit start", "pre-edit stop", "edit start", "edit stop", "junction start", "junction stop", "gRNA name", "gRNA start", "gRNA stop", "junction type", "blast vs. parent") . "\n";

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

my $preedit_name = $gene_name . "un";
my $edit_name = $gene_name . "ed";

# refhash to hold the transcript names
my $transcripts;

if($verbose)
{
	print "PRE-EDITED NAME: $preedit_name\n";
	print "EDITED NAME    : $edit_name\n";
	
}



# Get all guideRNA's that hit the fully edited and fully pre-edited transcript
my %edited_grna;
my %edited_grna_sequence;
my %preedited_grna;
my %preedited_grna_sequence;

# A guideRNA will be named by
# minicircle_name-minicircle_start-minicircle_stop

my @preedit_segments = $db->segment(-name=>$preedit_name, -class=>'minicircle_match');

if($verbose)
{
	print "Obtaining pre-edited guideRNA features ..\n";
}

foreach my $preedit_segment (@preedit_segments)
{
	my $grna_name = join("-", $preedit_segment->refseq, $preedit_segment->start, $preedit_segment->stop);

	# Set the edit_hash for this guide rna
	$preedited_grna{$grna_name}++;
	my $seq = Bio::Seq->new(-display_id=>$grna_name, -seq=> $preedit_segment->dna);
	$preedited_grna_sequence{$grna_name} = $seq;

	if($verbose)
	{
		print join("\t", $preedit_segment->name, $preedit_segment->start, $preedit_segment->stop) . "\n";
	}

}

if($verbose)
{
	print "Obtained pre-edited guideRNA features!\n";
	print "\n";
}


my @edit_segments = $db->segment(-name=>$edit_name, -class=>'minicircle_match');

if($verbose)
{
	print "Obtaining edited guideRNA features ..\n";
}

foreach my $edit_segment (@edit_segments)
{
	my $grna_name = join("-", $edit_segment->refseq, $edit_segment->start, $edit_segment->stop);

	# Set the edit_hash for this guide rna
	$edited_grna{$grna_name}++;
	my $seq = Bio::Seq->new(-display_id=>$grna_name, -seq=> $edit_segment->dna);
	$edited_grna_sequence{$grna_name} = $seq;

	if($verbose)
	{
		print join("\t", $edit_segment->name, $edit_segment->start, $edit_segment->stop) . "\n";
	}

}
if($verbose)
{
	print "Obtained edited guideRNA features!\n";
	print "\n";
}

# Create a database of them

my $tempdir = tempdir( CLEANUP=>1);
my $temp_seqio = Bio::SeqIO->new(-file=>">$tempdir/hits", -format=>'fasta');
	
while( my($edit_grna_seq_name, $edit_grna_seq) = each(%edited_grna_sequence) )
{
	$temp_seqio->write_seq($edit_grna_seq);
}

system("cd $tempdir;formatdb -p F -i hits");


## Now we have all of the guide rna's that hit the fully edited and pre-edited gene

# Get all of the transcripts that hit this gene (both pre-edited and edited) and place the est sequences in the transcript hash
# pre-edited first

my @genes = $db->get_feature_by_name(-class=>'Sequence', -name=>$preedit_name);
foreach my $gene (@genes)
{
	# Now pull out all transcripts of this gene

	if($verbose)
	{
		print "Getting EST transcript sequences for $preedit_name\n\n";
	}
	my $gene_segment = $db->segment($gene->name);
#	my @transcripts = $gene_segment->features('specific_bloodstream', 'specific_procyclic', 'dt_primed');
#	my @transcripts = $gene_segment->features('specific_bloodstream');
	my @transcripts = $gene_segment->features('specific_bloodstream', 'specific_procyclic', 'dt_primed', 'specific_bloodstream_fraction_11', 'specific_bloodstream_fraction_9', 'specific_procyclic_fraction_11', 'specific_procyclic_fraction_9');

	foreach my $transcript (@transcripts)
	{
		my $est_name = $transcript->attributes('est_name');
		$transcripts->{$transcript->attributes('est_name')} = 1;
		if($verbose)
		{
			print "$est_name FOUND.\n";
		}
	}
}

# now the edited ones

my @genes = $db->get_feature_by_name(-class=>'Sequence', -name=>$edit_name);

foreach my $gene (@genes)
{
	# Now pull out all transcripts of this gene

	if($verbose)
	{
		print "Getting EST transcript sequences for $edit_name\n\n";
	}

	my $gene_segment = $db->segment($gene->name);
#	my @transcripts = $gene_segment->features('specific_bloodstream', 'specific_procyclic', 'dt_primed');
	my @transcripts = $gene_segment->features('specific_bloodstream', 'specific_procyclic', 'dt_primed', 'specific_bloodstream_fraction_11', 'specific_bloodstream_fraction_9', 'specific_procyclic_fraction_11', 'specific_procyclic_fraction_9');
	foreach my $transcript (@transcripts)
	{
		my $est_name = $transcript->attributes('est_name');
		$transcripts->{$transcript->attributes('est_name')} = 1;
		if($verbose)
		{
			print "$est_name FOUND\n";
		}
	}
}

# iterate through all of these transcripts and get their feature objects
my $ed;
my $un;

while(my ($est_name, undef) = each(%$transcripts))
{
	my @features = $db->get_feature_by_attribute(est_name=>$est_name);
	foreach my $feat (@features)
	{
		my $vals;
		$vals->{parent} = $feat->refseq;
		if($feat->refseq =~ /ed$/)
		{
			$ed->{$est_name}->{obj} = $feat;
		} else
		{
			$un->{$est_name}->{obj} = $feat;
		}
	}
}

# iterate through them again

if($verbose)
{
	print "\nSetting Statistics for Transcript Sequences\n\n";
}


while(my ($est_name, undef) = each(%$transcripts))
{

	if($verbose)
	{
		print "$est_name : ";
	}

	# These vars describe the transcript
	my $preedit_start;
	my $preedit_stop;
	my $edit_start;
	my $edit_stop;
	my $seq;
	my $size;

	# These vars tell where on the fully edited/preedited gene this transcript hits
	my $gene_preedit_start;
	my $gene_preedit_stop;
	my $gene_edit_start;
	my $gene_edit_stop;

	my @transcripts = $db->get_feature_by_name(-class=>'Sequence', -name=>$est_name);
	my $size = 1;
	if(defined($transcripts[0]))
	{
		$seq = $transcripts->{seq};
		$size = length($seq);
		if($verbose)
		{
			" size:$size";
		}
	} else
	{
		# If the est name was not found, go the next est;
		if($verbose)
		{
			print "EST $est_name NOT found in the database.\n";
		}
	}


	# Get the locations of the pre-edited and edited start and stop locations 
	if(defined($un->{$est_name}))
	{
		# Get the percent pre-edited portion and start/stop
		$gene_preedit_start = $un->{$est_name}->{obj}->start;
		$gene_preedit_stop = $un->{$est_name}->{obj}->stop;

	}
	if(defined($ed->{$est_name}))
	{
		# Get the percent edited portion and start/stop
		$gene_edit_start = $ed->{$est_name}->{obj}->start;
		$gene_edit_stop = $ed->{$est_name}->{obj}->stop;
	}
	if($verbose)
	{
		print " gene-Pre-edit-range: $gene_preedit_start..$gene_preedit_stop gene-Edit-range: $gene_edit_start..$gene_edit_stop";
	}

	# Try for getting the segment
	my $transcript_segment = $db->segment(-class=>'Sequence', -name=>$est_name);


	# now find the gene_hit features on this (gene_match)

	my @features = $transcript_segment->features('gene_match');

	foreach my $feat (@features)
	{
		if($feat->name =~ /$preedit_name/)
		{
			$preedit_start = $feat->start;
			$preedit_stop = $feat->stop;
		} elsif($feat->name =~ /$edit_name/)
		{
			$edit_start = $feat->start;
			$edit_stop = $feat->stop
		} else
		{
			print "ERROR, wrong type found in gene_match\n";
		}
	}
	if($verbose)
	{
		print " Pre-edit-range: $preedit_start..$preedit_stop Edit-range: $edit_start..$edit_stop";
	}

	# Get all guideRNA's that are located within this region


	# Now determine the "junction" region between the end of the pre=edited and begining of edited
	
	my $junction_start;
	my $junction_stop;

	($junction_start, $junction_stop) = return_junction($preedit_start, $preedit_stop, $edit_start, $edit_stop);

	if($verbose)
	{
		print " junction: $junction_start..$junction_stop";
	}
	
	if($verbose)
	{
		print "\n";
	}

	# Now that we have the junction region, we need to find the guideRNA's that overlap that junction region
	
	my @features = $transcript_segment->features('match:pcr');
	foreach my $grna_feature( @features)
	{
		my $junction_type = 'none';

		my $grna_start = $grna_feature->start;
		my $grna_stop = $grna_feature->stop;
		my $grna_seq = Bio::Seq->new(-display_id=>$grna_feature->name, -seq=>$grna_feature->dna);
		if($verbose)
		{
			print "  " . join(" ", "grna:" . $grna_feature->name, "start:" . $grna_feature->start, "stop:" . $grna_feature->stop);
		}

		# The junction type is overlap-nojunction if there is no junction region and the start is in one part and stop in the other
		if( ($grna_start >= $preedit_start && $grna_start <= $preedit_stop) && ($grna_stop >= $edit_start && $grna_stop <= $edit_stop) )
		{
			$junction_type = 'overlap-nojunction';
		} elsif( ($grna_start >= $edit_start && $grna_start <= $edit_stop) && ($grna_stop >= $preedit_start && $grna_stop <= $preedit_stop) )
		{
			$junction_type = 'overlap-nojunction';
		} elsif( $grna_start >= $preedit_start && $grna_stop <= $preedit_stop )
		{
			$junction_type = 'pre-edit-only';
		} elsif( $grna_start >= $edit_start && $grna_stop <= $edit_stop )
		{
			$junction_type = 'edit-only';
		} else
		{
			$junction_type = 'odd-structure';
		}
		# Find out if this feature is part of the junction region, and the structure of it (lead-in, lead-out, internal, overlap)
#		if($junction_start + $junction_stop > 0)
		# Run for all gRNAs instead of just junction ones
		if(1)
		{


			# lead-in is if the guideRNA overlaps the junction region and the pre-edited portion of the transcript
			# lead-out is if the guideRNA overlaps the junction region and the edited portion of the transcript
			# internal is if the guideRNA is completly internal to the junction region of the transcript
			# overlap is if the guideRNA overlaps the junction region, edited, and pre-edited portions of the transcript
			
			my $preedit_overlap = 0;
			my $edit_overlap = 0;
			my $junction_overlap = 0;
			# Do I overlap the pre-edited portion
			# if grna_start is in between preedit_start and preedit_stop or grna_stop is the same then we overlap
			if( 
			    ( ($grna_start >= $preedit_start) && ($grna_start <= $preedit_stop) ) || 
			    ( ($grna_stop  >= $preedit_start) && ($grna_stop  <= $preedit_stop) )
			  )
			{
				$preedit_overlap = 1;
			}
			# Now check if we overlap the edited portion
			if( 
			    ( ($grna_start >= $edit_start) && ($grna_start <= $edit_stop) ) || 
			    ( ($grna_stop  >= $edit_start) && ($grna_stop  <= $edit_stop) )
			  )
			{
				$edit_overlap = 1;
			}
			# Now the junction region
			if( 
			    ( ($grna_start >= $junction_start) && ($grna_start <= $junction_stop) ) || 
			    ( ($grna_stop  >= $junction_start) && ($grna_stop  <= $junction_stop) )
			  )
			{
				$junction_overlap = 1;
			}

			if($preedit_overlap && !$edit_overlap && $junction_overlap)
			{
				$junction_type = 'lead-in';
			} elsif($preedit_overlap && $edit_overlap && $junction_overlap)
			{
				$junction_type = 'overlap';
			} elsif(!$preedit_overlap && $edit_overlap && $junction_overlap)
			{
				$junction_type = 'lead-out';
			} elsif(!$preedit_overlap && !$edit_overlap && $junction_overlap)
			{
				$junction_type = 'internal';
			}

			if($verbose)
			{
				print " preedit_overlap:$preedit_overlap edit_overlap:$edit_overlap junction_overlap:$junction_overlap";
			}
		}

		# Now check if this guideRNA exists within the fully edited sequence
		# The method of comparison we are going to use is blast against all of the edited grna sequences

		my $blast_factory = Bio::Tools::Run::StandAloneBlast->new('program'=>'blastn', 'database'=>"$tempdir/hits", 'F'=>'F');
		my $report = $blast_factory->blastall($grna_seq);
		my $result = $report->next_result;
		my $hit = $result->next_hit;
		my $blast_text = 'NOHIT';
		if($hit)
		{
			my $hsp = $hit->next_hsp;
			if($hsp)
			{
				$blast_text = "blast match: evalue->" . $hsp->evalue . " perc_identity->" . $hsp->percent_identity;
				if($verbose)
				{
					print " blast match: evalue->" . $hsp->evalue . " perc_identity->" . $hsp->percent_identity;
				}
			}
		} else 
		{ 
			if($verbose)
			{
				print " NOHIT";
			}
		}

		if($verbose)
		{
			print " junction_type:" . $junction_type . "\n";
		}
		if($print_out)
		{
			print OUTFILE join("\t", $est_name, $preedit_start, $preedit_stop, $edit_start, $edit_stop, $junction_start, $junction_stop, $grna_feature->name, $grna_start, $grna_stop, $junction_type, $blast_text) . "\n";
		}
	}

	if($verbose)
	{
		print "\n";
	}
}


sub return_junction
{
	my $preedit_start = shift;
	my $preedit_stop = shift;
	my $edit_start = shift;
	my $edit_stop = shift;

	if($edit_stop <= $preedit_start)
	{
		return ($edit_stop, $preedit_start)
	}elsif($preedit_stop < $edit_start)
	{
		return ($preedit_stop, $edit_start);
	} else
	{
		return (0, 0);
	}

}

sub get_seq_from_blastdb
{
	my $blastdb_name = shift;
	my $seq_name = shift;

	my ($tfh, $tfn) = tempfile();
	system("fastacmd -d $blastdb_name -s $seq_name > $tfn");

	my $seqio = Bio::SeqIO->new(-format=>'fasta', -file=>$tfn);
	my $seq = $seqio->next_seq;
	return $seq;
}


sub align_seq
{
	my $seqone = shift;
	my $seqtwo = shift;

	
}
