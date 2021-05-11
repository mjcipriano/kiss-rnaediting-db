# $Id: SixFrameDumper.pm,v 1.1.1.1 2005/06/28 22:10:29 mcipriano Exp $
#
# BioPerl module for Bio::Graphics::Browser::Plugin::SixFrameDumper
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich and Cold Spring Harbor Laboratories 2002
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Graphics::Browser::Plugin::SixFrameDumper - A plugin for dumping sequences in various formats

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

This is a plugin to the Generic Model Organism Database browse used by
Bio::Graphics::Browser to dump out an annotated region in a requested
flatfile format.  Currently the feature formats are 

=head1 FEEDBACK

See the GMOD website for information on bug submission http://www.gmod.org.

=head1 AUTHOR - Jason Stajich

Email jason@bioperl.org

Describe contact details here

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Graphics::Browser::Plugin::SixFrameDumper;
# $Id: SixFrameDumper.pm,v 1.1.1.1 2005/06/28 22:10:29 mcipriano Exp $
# Six Frame Dumper plugin

use strict;
use Bio::Graphics::Browser::Plugin;
use Bio::SeqIO;
use vars qw($VERSION @ISA);
use constant DEBUG => 0;
use CGI qw(:standard);

use Mbl;
use Apache::Session::File;

$VERSION = '0.1';

@ISA = qw(Bio::Graphics::Browser::Plugin);

my $db_name = '';

sub name { "Six Frame Translation" }
sub description {
  p("The Six Frame Translation Dumper prints out the displayed sequences' six frame translation",
    "in fasta format").
  p("This plugin was written by Michael Cipriano");
}

sub init {
	my $self = shift;
	my $conf = $self->browser_config;
	$db_name = $conf->plugin_setting('db_name');
}

sub mime_type {
  my $self = shift;
  my $config = $self->configuration;
  return 'text/plain' if $config->{format} eq 'text';
  return 'text/html'  if $config->{format} eq 'html';
  return wantarray ? ('application/octet-stream','dumped_region') : 'application/octet-stream'
    if $config->{format} eq 'todisk';
  return 'text/plain';
}
                                                                                                                                                                                                                                                  
sub config_defaults {
  my $self = shift;
  return { format           => 'html',
           fileformat       => 'fasta',
           wantsorted       => 0,
       };
}
                                                                                                                                                                                                                                                  
sub reconfigure {
  my $self = shift;
  my $current_config = $self->configuration;
                                                                                                                                                                                                                                                  
  foreach my $param ( $self->config_param() ) {
      $current_config->{$param} = $self->config_param($param);
  }
}

sub configure_form {
  my $self = shift;
  my $current_config = $self->configuration;
  my @choices = TR({-class => 'searchtitle'},
                        th({-align=>'RIGHT',-width=>'25%'},"Output",
                           td(radio_group(-name     => $self->config_name('format'),
                                          -values   => [qw(text html todisk)],
                                          -default  => $current_config->{'format'},
                                          -labels   => {html => 'html/xml',
                                                        'todisk' => 'Save to Disk',
                                                       },
                                          -override => 1))));
  my $browser = $self->browser_config();
  # this to be fixed as more general
                                                                                                                                                                                                                                                  
 
                                                                                                                                                                                                                                                  
  my $html= table(@choices);
  $html;
}





sub dump {
  my $self = shift;
  my $segment = shift;
  
  my $config  = $self->configuration;  
  my $browser = $self->browser_config();


        my $mbl = Mbl::new(undef, $db_name);
        my $dbh = $mbl->dbh();
                                                                                                                                                                                                                                                      
        my %session;
        eval {tie %session, "Apache::Session::File", cookie('SESSION_ID' . '_' . $mbl->organism), { Directory => $mbl->session_tmp_dir }; };
                                                                                                                                                                                                                                                      
        my $login_id = 0;
                                                                                                                                                                                                                                                      
        if(!$session{login_id})
        {
                $login_id = undef;;
        } else
        {
                $login_id = $session{login_id}
        }
	


  my $seq  = new Bio::Seq(-display_id       => $segment->display_id,
			  -seq=>$segment->seq
			  );

#  $seq->add_date(strftime("%d-%b-%Y",localtime));
  $segment->absolute(1);

  my $offset     = $segment->start - 1;
  my $segmentend = $segment->length;

  my $out = new Bio::SeqIO(-format => 'fasta');
  my $mime_type = $self->mime_type;

  my $frame1 = $seq->translate();
  my $frame2 = $seq->trunc(2, $seq->length)->translate();
  my $frame3 = $seq->trunc(3, $seq->length)->translate();
  my $frame4 = $seq->revcom->translate();
  my $frame5 = $seq->revcom->trunc(2, $seq->length)->translate;
  my $frame6 = $seq->revcom->trunc(3, $seq->length)->translate;

  $frame1->display_id("F1");
  $frame2->display_id("F2");
  $frame3->display_id("F3");
  $frame4->display_id("R1");
  $frame5->display_id("R2");
  $frame6->display_id("R3");

  my $sequence = uc($seq->seq);

  my @orfs_found;

  # Find all orfs > threshold in the forward direction
  my $threshold = 100;
  
  my $direction = "+";
  while ($sequence =~ /((ATG)([ACTG][ACTG][ACTG])*?(TGA|TAG|TAA))/g) 
  {   
      
      my $orf_length = length($1);
      my $orf_end = pos($sequence)-1+$offset ;
      my $orf_start = $orf_end - $orf_length +1;
      
      
      
      # Transfer it all to contig coordinates
      my $contig_id = $mbl->get_contig_from_supercontig_coord($seq->display_id, $orf_start);
      $orf_start = $mbl->get_contig_coords_from_supercontig($seq->display_id, $orf_start);
      $orf_end =  $mbl->get_contig_coords_from_supercontig($seq->display_id, $orf_end);
      

      next if $orf_length < $threshold;

      push @orfs_found, [$contig_id, $orf_start, $orf_end, $direction];

   }
   
   # Now do the reverse complement
   $sequence = uc($seq->revcom->seq);
   $direction = "-";
   my $seq_length = $seq->length();
   my $supercontig_end = $seq_length + $offset;
   while($sequence =~/((ATG)([ACTG][ACTG][ACTG])*?(TGA|TAG|TAA))/g)
   {
	   my $orf_length = length($1);
	   my $orf_start = pos($sequence)-1;
	   my $orf_end = $orf_start - $orf_length+1;

	   # Flip it
	   $orf_start = $seq_length - $orf_start - 1;
	   $orf_end = $seq_length - $orf_end - 1;

	   
	   # Make it overall supercontig Coordiates
	   $orf_start = $offset + $orf_start;
	   $orf_end = $offset + $orf_end;
	   
	   my $contig_id = $mbl->get_contig_from_supercontig_coord($seq->display_id, $orf_start);
	   $orf_start = $mbl->get_contig_coords_from_supercontig($seq->display_id, $orf_start);
	   $orf_end =  $mbl->get_contig_coords_from_supercontig($seq->display_id, $orf_end);
	   next if $orf_length < $threshold;

	   push @orfs_found, [$contig_id, $orf_start, $orf_end, $direction];
	   
   }
   
  $mime_type = "html";
  
  if ($mime_type =~ /html/) 
  {
    print "<pre>";  
    print start_html($segment->desc),h1($segment->desc);

    $out->write_seq($seq);
    $out->write_seq($frame1);
    $out->write_seq($frame2);
    $out->write_seq($frame3);
    $out->write_seq($frame4);
    $out->write_seq($frame5);
    $out->write_seq($frame6);
    print "\n<h2>ORFs Found</h2>\n";
    foreach my $orf(@orfs_found)
    {
	    my $orfdir = $orf->[3];
	    if($orfdir eq "+")
	    {
		    $orfdir = '%2B';
	    }
	    if($mbl->check_annotation_admin_rights($login_id))
	    {
		    print $mbl->add_orf_link($orf->[0], $orf->[1], $orf->[2], $orfdir) . "\t";
	    }
	    
	    print join("\t", "contig:" . $orf->[0], $orf->[1], $orf->[2], $orf->[3], "length:" . ($orf->[2] -$orf->[1] + 1) );
	    
	    print "<br>\n";  

    }
    print "</pre";
	
    # Print out the start and stop coordinates
    #print end_html;
  } else {
    $out->write_seq($seq);
    $out->write_seq($frame1);
    $out->write_seq($frame2);
    $out->write_seq($frame3);
    $out->write_seq($frame4);
    $out->write_seq($frame5);
    $out->write_seq($frame6);
  
  }
  
  #undef $out;
}



1;
