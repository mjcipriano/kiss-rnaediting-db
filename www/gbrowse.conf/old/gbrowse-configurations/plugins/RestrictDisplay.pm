package Bio::Graphics::Browser::Plugin::RestrictDisplay;
# $Id: RestrictDisplay.pm,v 1.1 2005/07/21 18:37:28 mcipriano Exp $
# test plugin
use strict;
use Bio::Graphics::Browser::Plugin;
use CGI qw(table a TR td th p popup_menu radio_group checkbox checkbox_group h1 h2 pre);
use Text::Shellwords;


use vars '$VERSION','@ISA';
$VERSION = '0.21';

@ISA = qw(Bio::Graphics::Browser::Plugin);

my %SITES;


sub name { "Restrict Display" }

sub description {
  p("The Restrict Display Plugin will only show features that have a target word in their name or in any of their attributes",
  p("This plugin was written by Michael Cipriano");
}

sub type { 'annotator' }

#sub init {
#  my $self = shift;
#  my $browser_conf = $self->browser_config;
#  my $configuration = $self->configuration;
#  $configuration->{default_url} = "http://mib.mbl.edu/db/showtag.php";

#  $self->{default_url} = $browser_conf->plugin_setting('default_url');
 
#}

sub config_defaults {
#  my $self = shift;
#  return { default_url  => "http://mib.mbl.edu/db/showtag.php" 
#         };
}

sub reconfigure {
#  my $self = shift;
#  my $current_config = $self->configuration;
#  %$current_config = map {$_=>1} param('ORFAnnotator.variable');
#  $current_config->{on} = param('ORFAnnotator.on');
}



sub configure_form {
  my $self = shift;
  my $sequence = param('RestrictDisplay.searchstring');
  my $msg  =  $sequence
              ? font({-color=>'red'},"Invalid Search String")
              : '';
  return $msg .
    table(TR({-class=>'searchtitle'},
             th({-colspan=>2,-align=>'LEFT'},
                'Enter a search string.',
                'The browser will show all sequences that have that search string in their name or any of their attributes')),
          TR({-class=>'searchbody'},
             td('Enter search string:'),
             td(textarea(-name=>'RestrictDisplay.searchstring',-cols=>80, -rows=>10))));

}
  

sub annotate {
  my $self = shift;
  my $segment = shift;


  my $ref        = $segment->ref;
  my $abs_start  = $segment->start;
  my $abs_end    = $segment->end;
#  my $dna        = $segment->seq;
 
         
  my $feature_list   = Bio::Graphics::FeatureFile->new;



    $feature_list->add_type('search'=>{glyph   => 'arrow',
                                    key     => "sagetag",
                                    fgcolor => 	'black'
                                    bgcolor => "purple",
                                    point   => 0,
				    linewidth => 2,
                                    orient  => 'N',
                                   });



	# Find out how many places this sage tag maps
 	my $db    = $self->database or die "I do not have a database";

        # pull out all orf features

        my @orfs = $segment->features( -type=>'sagetag');
                                                                                           
        for my $o (@orfs) {
        my @feature_list = $db->get_feature_by_name('sagetag' => $o->name);
	my $feature_count = 0;
	foreach my $this_feature (@feature_list)
	{
		if($this_feature->method eq 'sagetag')
		{
			$feature_count++;
		}
		
	}
	my $map_number = $feature_count;

      my $feature = Bio::Graphics::Feature->new(-start=>$o->start,-stop=>$o->stop,-strand=>$o->strand,-ref=>$ref,-name=>$this_name, -url=>$url, -glyph=>"anchored_arrow");
      $feature_list->add_feature($feature,'search');
}


  return $feature_list;
}


1;

