package Bio::Graphics::Browser::Plugin::AnnotateOrfs;
# $Id: AnnotateOrfs.pm,v 1.1.1.1 2005/06/28 22:10:29 mcipriano Exp $

use strict;
use Bio::Graphics::Browser::Plugin;
use CGI qw(:standard *table);
use vars '$VERSION','@ISA','$db_name';


=head1 NAME

Bio::Graphics::Browser::Plugin::AnnotateOrfs -- a plugin that executes NCBI's bl2seq on the current view

=head1 SYNOPSIS

 in 0X.organism.conf:
     

=head1 DESCRIPTION

This Gbrowse plugin will take a sequence (entered in the configuration screen)
and BLAST it against the current display, with hits as new sequence features.

You must, of course, have the NCBI Blast suite of programs installed,
you must have configured the plugin to be visible, and you must
set a single plugin parameter in the 0X.organism.conf file:
    [AnnotateOrfs:plugin]
    bl2seq_executable = /path/to/your/bl2seq

=cut

$db_name = "";

$VERSION = '0.02';

@ISA = qw(Bio::Graphics::Browser::Plugin);

my @COLORS = qw(red green blue orange cyan black 
		turquoise brown indigo wheat yellow emerald);

sub name { "Annotate All Orfs" }

sub description {
  p("This plugin will take an input DNA sequence - entered using the 'Configure' button - and will run bl2seq (a BLAST sequence alignment) ",
    "against the assembly sequence in the current view").
  p("This plugin was written by Mark Wilkinson.");
}

sub type { 'annotator' }
sub init {
    my $self = shift;
    my $conf = $self->browser_config;
    $db_name = $conf->plugin_setting('db_name');
}

sub config_defaults {
  my $self = shift;
}

sub reconfigure {
  my $self = shift;
  my $current = $self->configuration;
}

sub configure_form {
  my $self = shift;
  my $current_config = $self->configuration;

}
  

sub annotate {
    my $self = shift;
    my $segment = shift;
    my $ref        = $segment->ref;
    my $abs_start  = $segment->start;
    my $dna        = $segment->seq;
    my $conf = $self->configuration;
     my $feature_list   = Bio::Graphics::FeatureFile->new(-smart_features => 1);
     $feature_list->add_type('Orf_Found'=>{glyph   => 'arrow',
                                    key     => "PossibleOrf",
                                    fgcolor => 'red',
                                    bgcolor => "red",
                                    point   => 0,
                                    linewidth => 2,
                                    orient  => 'Y',
                                    strand_arrow=>'1',
                                    key=>'Possible Orfs'
                                   });




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
      my $contig_orf_start = $mbl->get_contig_coords_from_supercontig($seq->display_id, $orf_start);
      my $contig_orf_end =  $mbl->get_contig_coords_from_supercontig($seq->display_id, $orf_end);
      

      next if $orf_length < $threshold;

      push @orfs_found, [$contig_id, $contig_orf_start, $contig_orf_end, $direction, $seq->display_id, $orf_start, $orf_end,];

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
	   my $contig_orf_start = $mbl->get_contig_coords_from_supercontig($seq->display_id, $orf_start);
	   my $contig_orf_end =  $mbl->get_contig_coords_from_supercontig($seq->display_id, $orf_end);
	   next if $orf_length < $threshold;

	   push @orfs_found, [$contig_id, $contig_orf_start, $contig_orf_end, $direction, $seq->display_id, $orf_start, $orf_end,];
	   
   }



    foreach my $orf(@orfs_found)
    {
	    
	    my $orfdir = $orf->[3];
	    if($orfdir eq "+")
	    {
		    $orfdir = '%2B';
	    }
	    my $link;

	    $link = $mbl->add_orf_link($orf->[0], $orf->[1], $orf->[2], $orfdir, 1) . "\t";
		
	    
	    
	    my $strand;
	    if($orf->[3] eq "+")
	    {
		    $strand = 1;
	    } else
	    {
		    $strand = "-1";
	    }
	    
	    my $orf_length = $orf->[6] - $orf->[5];
	    my $feature = Bio::Graphics::Feature->new(-start=>$orf->[5],-stop=>$orf->[6],-strand=>$strand,-ref=>$orf->[4],-name=>'Orf:' . $orf->[0] . ":" . $orf->[1] . ".." . $orf->[2] . ' ' . $orf_length . "bp", -url=>$link, -glyph=>"anchored_arrow");

	    $feature_list->add_feature($feature,'Orf_Found');
	    

    }



    return $feature_list;
}


1;

