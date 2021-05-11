package Bio::Graphics::Browser::Plugin::OrfFinder;
# $Id: OrfFinder.pm,v 1.1 2005/07/21 18:37:28 mcipriano Exp $

use strict;
use Bio::Graphics::Browser::Plugin;
use Bio::Graphics::Feature;
use DBI;
use CGI qw(:standard *table);
use Getopt::Long;
use Bio::SimpleAlign;
use Bio::Tools::Run::StandAloneBlast;
use Bio::Tools::BPlite::Sbjct;
use Bio::AlignIO;
use Bio::SeqIO;
use DBI;

use Bio::Root::IO;

use vars '$VERSION','@ISA','$blast_executable';
$VERSION = '0.10';

$ENV{PATH} = '$PATH:/bin:/usr/bin:/usr/local/bin:/var/bio/blast/bin:/var/bio/bin';

my $blast_executable = '';
my $blast_db = '';
my $driver = '';
my $hostname = '';
my $port = '';
my $user = '';
my $password = '';
my $database = '';

@ISA = qw(Bio::Graphics::Browser::Plugin);

sub name { "Orf" }

sub description {
  p("The Orf finder plugin will find current Orfs or do a blast search of previously deleted Orfs.",
    "[NOTE TO SYSADMINS: A local blast server must be set up and blastall must be in the servers path.]").
  p("This plugin was written by Michael Cipriano.");
}

sub type { 'finder' }
sub init {
  my $self = shift;
  my $conf = $self->browser_config;
  $blast_executable = $conf->plugin_setting('blastall_executable');
  $blast_db = $conf->plugin_setting('blast_db');
  $driver =  $conf->plugin_setting('driver');
  $hostname = $conf->plugin_setting('hostname');
  $port = $conf->plugin_setting('port');
  $user = $conf->plugin_setting('user');
  $password = $conf->plugin_setting('password');
  $database = $conf->plugin_setting('database');
}

sub config_defaults {
  my $self = shift;
  return { };
}

# we have no stable configuration
# sub reconfigure { }

sub configure_form {
  my $self = shift;
  my $orfid = param('OrfFinder.orfid'); 
  my $msg  =  $orfid 
              ? font({-color=>'red'},"Invalid orfid: This orfid has never existed.")
	      : '';
  return $msg .
    table(TR({-class=>'searchtitle'},
	     th({-colspan=>2,-align=>'LEFT'},
		'Enter an orfid to search the database for.',
		'The browser will find that orf or identify all genomic regions that are similar',
		'to this sequence.')),
	  TR({-class=>'searchbody'},
	     td('Enter orf ID:'),
	     td(textfield(-name=>'OrfFinder.orfid'))));
}

# find() returns undef unless the SequenceFinder.searchsequence parameter
# is specified and valid.  Returning undef signals the browser to invoke the
# configure_form() method.
# If successful, it returns an array ref of Bio::SeqFeatureI objects.
sub find {
  my $self     = shift;
  my $segments = shift; # current segments - can search inside them or ignore
                        # In this example we do a global search.

  my $db    = $self->database or die "I do not have a database";
  my @results;
  my $orfid = lc param('OrfFinder.orfid');
  my $dbs = $blast_db;
  my $dir = "/tmp";
  my @params = (     'database'    => "$dbs",
                     'program'     => "blastn"
                     );
                                                                                           
  my $factory = Bio::Tools::Run::StandAloneBlast->new(@params);


  # First check to see if this orfid exists in the database
#  my $feature = $db->get_feature_by_name('ORF' => $orfid);
 
  my $max_eval = 1e-100;
#  if($feature)
#  {
#    $max_eval = 0;
#  }


  # Else we have not found this orf, so search the database for the orfid, get the sequence, and blast it against the database

  # Connect to the mysql database

  my $dsn = "DBI:$driver:database=$database;host=$hostname;port=$port";
  my $dbh = DBI->connect($dsn, $user, $password);
  my $drh = DBI->install_driver("mysql");

  my $query = "select sequence from orfs where orfid = " . $dbh->quote($orfid) ;
  my $rh = $dbh->prepare($query);
  $rh->execute();
  my $sequence = '';
  my $row = $rh->fetchrow_hashref;
  if($row)
  {
    $sequence = $row->{sequence};
  } else 
  {
    return;
  }

  my $input = Bio::Seq->new( -id  => 'my_query',
                           -seq =>$sequence);
                                                                                           
  my $blast_report = $factory->blastall($input);
                                                                                           
  my $this_result = $blast_report->next_result;

  while( my $hit = $this_result->next_hit())
  {
    my $this_hsp = $hit->hsp;
    my $homolog_string = $this_hsp->homology_string;
    my $query_string = $this_hsp->query_string;
    my $hit_string = $this_hsp->hit_string;
    my $line_length = 60;
    my $current_portion = 0;
    my $desc = '';
    my $spacer = '<BR>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; ';
    my $qsh = 'Query-> :   ';
    my $rsh = '<BR>Result->:   ';

    while(length($homolog_string)+1 >= $current_portion)
    {
        $desc = $desc . $qsh . substr($query_string, $current_portion, $line_length) . $spacer . substr($homolog_string, $current_portion, $line_length) . $rsh . substr($hit_string, $current_portion, $line_length) . '<BR><BR>';
        $current_portion = $line_length + $current_portion;
    }

    my $hit_name = $hit->name();
    my $current_time = localtime();
    my $name = $hit_name;
    my $type = "blast_search:$current_time" ;
    #my $type = "blast_search";
    my $ref = $hit->name();
    if( ($hit->significance <= $max_eval) && ( abs(length($sequence) - ($this_hsp->hit->end - $this_hsp->hit->start)) <= 20) )
    {
      push @results, Bio::Graphics::Feature->new(-ref   => $ref,
                                                   -type  => $type,
                                                   -name  => $name,
                                                   -start => $this_hsp->hit->start,
                                                   -score => $hit->significance,
                                                   -end   => $this_hsp->hit->end,
                                                   -desc  => $desc);
    }
}
                                                                                           
  return \@results;

#  $self->auto_find($sequence);
}

# auto_find() does the actual work
# It is also called by the main page as a last resort when the user
# types something into the search box that isn't recognized.
sub auto_find {
#  my $self  = shift;
#  my $sequence = shift;

#  my $dbs = 'var/bio/blast/data/giardia';
#  my $dir = "/tmp";
#  my @params = (     'database'    => "$dbs",
#                     'program'     => "blastn"
#                     );
 
#  my $factory = Bio::Tools::Run::StandAloneBlast->new(@params);
#  $sequence =~ s/[ \t\n\f\r]//g;
#  my $input = Bio::Seq->new( -id  => 'my_query',
#                           -seq =>$sequence);

#  my $blast_report = $factory->blastall($input);
                                                                                
#  my $this_result = $blast_report->next_result;
#  my @results;
#  while( my $hit = $this_result->next_hit())
#  {
#    my $this_hsp = $hit->hsp;
#    my $homolog_string = $this_hsp->homology_string;
#    my $query_string = $this_hsp->query_string;
#    my $hit_string = $this_hsp->hit_string;
#    my $desc = "
#Query-> :   $query_string<BR>
#&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; $homolog_string<BR>
#Result->:   $hit_string";
#    my $hit_name = $hit->name();
#    my $current_time = localtime();
#    my $name = $hit_name;
#    my $type = "blast_search:$current_time" ;
#    #my $type = "blast_search";
#    my $ref = $hit->name();
#    push @results, Bio::Graphics::Feature->new(-ref   => $ref,
#						   -type  => $type,
#						   -name  => $name,
#						   -start => $this_hsp->hit->start,
#						   -score => $hit->significance,
#						   -end   => $this_hsp->hit->end,
#						   -desc  => $desc);
#}

#  return \@results;
}



1;
