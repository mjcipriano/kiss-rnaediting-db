package Bio::Graphics::Browser::Plugin::ORFAnnotator;
# $Id: ORFAnnotator.pm,v 1.1 2005/07/21 18:37:28 mcipriano Exp $
# test plugin
use strict;
use Bio::Graphics::Browser::Plugin;
use CGI qw(:standard *table);

use vars '$VERSION','@ISA';
$VERSION = '0.21';

@ISA = qw(Bio::Graphics::Browser::Plugin);

my %SITES;


sub name { "ORF quality" }

sub description {
  p("The ORF Annotator will color ORFs based on their quality",
    "on the current view.").
  p("This plugin was written by Michael Cipriano");
}

sub type { 'annotator' }


sub config_defaults {
#  my $self = shift;
  return { };
}

sub reconfigure {
#  my $self = shift;
#  my $current_config = $self->configuration;
#  %$current_config = map {$_=>1} param('ORFAnnotator.variable');
#  $current_config->{on} = param('ORFAnnotator.on');
}



sub configure_form {
  my $self = shift;
  my $current_config = $self->configuration;
}
  

sub annotate {
  my $self = shift;
  my $segment = shift;
  my $config  = $self->configuration;

  my $ref        = $segment->ref;
  my $abs_start  = $segment->start;
  my $abs_end    = $segment->end;
#  my $dna        = $segment->seq;
 
         
  my $feature_list   = Bio::Graphics::FeatureFile->new;
  my $type = 'ORFq';


    $feature_list->add_type('ORF_0'=>{glyph   => 'arrow',
				    key     => "ORF bad",
				    fgcolor => "green",
				    bgcolor => "green",
				    point   => 0,
                                    linewidth => 1,
				    orient  => 'N',
				   });

    $feature_list->add_type('ORF_1'=>{glyph   => 'arrow',
                                    key     => "ORF bad",
                                    fgcolor => "blue",
                                    bgcolor => "blue",
                                    point   => 0,
                                    linewidth => 2,
                                    orient  => 'N',
                                   });

    $feature_list->add_type('ORF_2'=>{glyph   => 'arrow',
                                    key     => "ORF low",
                                    fgcolor => "red",
                                    bgcolor => "red",
                                    point   => 0,
                                    linewidth => 3,
                                    orient  => 'N',
                                   });

    $feature_list->add_type('ORF_3'=>{glyph   => 'arrow',
                                    key     => "ORF middle",
                                    fgcolor => "yellow",
                                    bgcolor => "yellow",
                                    point   => 0,
                                    linewidth => 4,
                                    orient  => 'N',
                                   });

    $feature_list->add_type('ORF_4'=>{glyph   => 'arrow',
                                    key     => "ORF high",
                                    fgcolor => "purple",
                                    bgcolor => "purple",
                                    point   => 0,
				    linewidth => 5,
                                    orient  => 'N',
                                   });


#      my $pos = $abs_start + pos($dna) - length($1) + $offset;

         # pull out all orf features
         my @orfs = $segment->features( -type=>'ORF');
                                                                                           
         for my $o (@orfs) {

	my $attr_score = 1;
      if($o->attributes('CodonPreference') eq 'CP_PASS')
      {
	 $attr_score++;
      }
      if($o->attributes('GeneScan') eq 'GS_PASS')
      {
         $attr_score++;
      }
      if($o->attributes('TestCode') eq 'TC_PASS')
      {
         $attr_score++;
      }
      if($o->attributes('SwissProt') eq '"No significant swissprot hit"')
      {
	 $attr_score--;
      }

      my $line_type = "ORF_$attr_score";
      my $feature = Bio::Graphics::Feature->new(-start=>$o->start,-stop=>$o->stop,-strand=>$o->strand,-ref=>$ref,-name=>$o->name, -url=>"/giardia/showorf.php?orf=". $o->name, -glyph=>"anchored_arrow" );
      $feature_list->add_feature($feature,$line_type);
}


  return $feature_list;
}


1;

