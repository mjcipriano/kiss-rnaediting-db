package Bio::Graphics::Browser::Plugin::SequenceFinder;
# $Id: SequenceFinder.pm,v 1.1 2005/07/21 18:37:28 mcipriano Exp $
# test plugin
use strict;
use Bio::Graphics::Browser::Plugin;
use Bio::Graphics::Feature;
use DBI;
use CGI qw(:all *table);
use Getopt::Long;
use Bio::SimpleAlign;
use Bio::Tools::Run::StandAloneBlast;
use Bio::Tools::BPlite::Sbjct;
use Bio::AlignIO;
use Bio::SeqIO;
use Bio::Root::IO;

use vars '$VERSION','@ISA';
$VERSION = '0.10';

my $blast_executable = '';
my $blast_db = '';

@ISA = qw(Bio::Graphics::Browser::Plugin);

sub name { "Sequences via Local Blast" }

sub description {
  p("The Blast finder plugin finds sequences and matches them with GMOD.",
    "[NOTE TO SYSADMINS: A local blast server must be set up and blastall must be in the servers path.]").
  p("This plugin was written by Michael Cipriano.");
}

sub type { 'finder' }
sub init {
    my $self = shift;
    my $conf = $self->browser_config;
    $blast_executable = $conf->plugin_setting('blastall_executable');
    $blast_db = $conf->plugin_setting('blast_db');
}

sub config_defaults {
  my $self = shift;

  return { program=> 'blastn',
	   expect=>  '10'};
}

# we have no stable configuration
# sub reconfigure { }

sub configure_form {
  my $self = shift;
  my $sequence = param('SequenceFinder.searchsequence'); 
  my $msg  =  $sequence 
              ? font({-color=>'red'},"Invalid sequence: either too short or not DNA")
	      : '';
  return $msg .
    table(TR({-class=>'searchtitle'},
	     th({-colspan=>2,-align=>'LEFT'},
		'Enter an sequence in raw format.',
		'The browser will identify all genomic regions that are similar',
		'to this sequence.')),
	  TR({-class=>'searchbody'},
	     td('Enter sequence:'),
	     td(textarea(-name=>'SequenceFinder.searchsequence',-cols=>80, -rows=>10))),
	  TR({-class=>'searchbody'},
             td('Program:'),
	     td(popup_menu( -name=>'SequenceFinder.program', -values =>["blastn","tblastn", "tblastx"]))),
          TR({-class=>'searchbody'},
             td('Expect:'),
             td(popup_menu( -name=>'SequenceFinder.expect', -values =>["0.0001", "0.01","1", "10", "100", "1000"], -default=>'10')))
);
}

# find() returns undef unless the SequenceFinder.searchsequence parameter
# is specified and valid.  Returning undef signals the browser to invoke the
# configure_form() method.
# If successful, it returns an array ref of Bio::SeqFeatureI objects.
sub find {
  my $self     = shift;
  my $segments = shift; # current segments - can search inside them or ignore
                        # In this example we do a global search.

  my $sequence = lc param('SequenceFinder.searchsequence');
  my $program = param('SequenceFinder.program');
  my $expect =  param('SequenceFinder.expect');
  my $dbs = $blast_db;
  my $dir = "/tmp";

  if($program eq "")
  {
    $program = 'blastn';
  }
  if($expect eq '')
  {
    $expect = '10';
  }
  my @params = (     'database'    => "$dbs",
                     'program'     => $program,
		     'e'	   => $expect,
		     'F'	   => 'F'
                     );

  my $factory = Bio::Tools::Run::StandAloneBlast->new(@params);

  $sequence =~ s/[ \t\n\f\r]//g;
  $sequence =~ s/[^a-z]//g;

  my $input = Bio::Seq->new( -id  => 'my_query',
                           -seq =>$sequence);
                                                                                           
  my $blast_report = $factory->blastall($input);
                                                                                           
  my $this_result = $blast_report->next_result;
  my @results;
  while( my $hit = $this_result->next_hit())
  {
    while( my $hsp = $hit->next_hsp)
    {
    my $this_hsp = $hsp;
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

    my $ref = $hit->name();
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

}

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
