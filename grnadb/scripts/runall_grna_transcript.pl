#!/usr/bin/perl

# Written by Michael J. Cipriano, 2005

use Bio::DB::GFF;
use Getopt::Long;
use File::Temp qw/tempfile tempdir/;
use Bio::SeqIO;

use strict;

my $web_dir = '/srv/www/htdocs/kiss';
my $download_dir = "$web_dir/download";
my $bioware_dir = '/usr/local/bioware';


my $db_name;
my $output_dir;
my $verbose = 0;
my $grna_transcript_script = "$bioware_dir/grnadb/scripts/grna_transcript.pl";

my $options = GetOptions(
                                "db_name=s", \$db_name,
				"output_dir=s", \$output_dir,
				"verbose!", \$verbose
);


if( !defined($db_name) && !defined($output_dir) )
{
	print "

Options are:
--db_name	The name of the Bio::DB::GFF database to query
--output_dir	The directory to store all of the results. Many new directories will be created under this directory.
--verbose	Prints out debugging information

Example: runall_grna_unique.pl --db_name=tbgrnadb --prefix=t_bruceibrucei --verbose 
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
my $gene_names;

my @genes = $db->features("gene:pcr");

if(-e $output_dir)
{
	# Good
} else
{
	mkdir $output_dir;
}

if(-e $output_dir)
{
	# Still good
} else
{
	die("Output directory does not exist or could not create new output directory.");
}

foreach my $gene(@genes)
{
	my $gene_name = $gene;
	my ($gene_orig) = $gene =~ /\((\w+)\)$/;
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

	$gene_names->{$gene_name} = 1;

}

while(my ($gene_name,undef) = each(%$gene_names))
{
	if($verbose)
	{
		print $gene_name . "\n";
	}

	my $gene_dir = $output_dir . "/" . $gene_name;
	mkdir($gene_dir);
	system($grna_transcript_script . " --db_name=$db_name --gene_name=$gene_name --outfile=$gene_dir/unique_results.tab");
	system("cd $gene_dir;grep NOHIT unique_results.tab > nohit.tab");
	system("cd $gene_dir;grep NOHIT unique_results.tab | grep internal  > nohit_internal.tab");
	system("cd $gene_dir;grep NOHIT unique_results.tab | grep lead > nohit_lead.tab");
	system("cd $gene_dir;grep -v NOHIT unique_results.tab > hit.tab");
	system("cd $gene_dir;grep NOHIT unique_results.tab | grep -v pre-edit-only | grep 'edit-only' > nohit_editonly.tab");
	system("cd $gene_dir;grep NOHIT unique_results.tab | grep pre-edit-only > nohit_preeditonly.tab");


}


