package Bio::Graphics::Browser::Plugin::SageResults;
# $Id: OldSageResults.pm,v 1.1 2005/07/21 18:37:28 mcipriano Exp $
# test plugin
use strict;
use Bio::Graphics::Browser::Plugin;
use CGI qw(:all *table);
use Text::Shellwords;
use Mbl;
use Apache::Session::File;


use vars '$VERSION','@ISA','$sage_url', '$db_name';
$VERSION = '0.21';

@ISA = qw(Bio::Graphics::Browser::Plugin);

my %SITES;


sub name { "Sage Results" }

sub description {

	return  "This track allow a user to visualize SAGE tags that map to a particular region in the genome".br,p, 
"'blue' means the tag is assigned to an ORF in the current view,",br,
"'slategray' means the tag is assigned to an ORF elsewhere in the assembly, and",br,
"'red' means the tag has been classified as Unknown.  <br><br>",p,
"Each tag is given an ID number, tag mapping type",br,
"(PS=primary sense, PA=primary antisense, AS=alternate sense, AA=alternate antisense), ORF ID number, and annotation.",br,p,
"Hovering over a tag brings up a box listing the frequency of a tag appearing throughout the SAGE libraries.",br,p,
  "This plugin was written by Michael Cipriano";
}

sub type { 'annotator' }

$sage_url = '';
$db_name = '';

sub init {
    my $self = shift;
    my $conf = $self->browser_config;
    $sage_url = $conf->plugin_setting('sage_url');
    $db_name = $conf->plugin_setting('db_name');
}

sub config_defaults {
  my $self = shift;

  return { };

}

sub reconfigure 
{
	my $self = shift;
	my $current_config = $self->configuration;
	$current_config->{mincount} = param('SageResults.mincount');
	$current_config->{valtype} = param('SageResults.valtype');
	$current_config->{uniquetranscript} = param('SageResults.uniquetranscript');
	$current_config->{tagtypes} = 1;
	my @tagtypes = param('SageResults.tagtypes');

	$current_config->{"Primary Sense Tag"} = undef;
	$current_config->{"Alternate Sense Tag"} = undef;
	$current_config->{"Primary Antisense Tag"} = undef;
	$current_config->{"Alternate Antisense Tag"} = undef;
	$current_config->{"Unknown"} = undef;

	foreach my $val(@tagtypes)
	{
		$current_config->{$val} = 1;
	}
}



sub configure_form 
{
	my $self = shift;
	my $config  = $self->configuration;

	my $valtype = $config->{valtype};
	my $uniquetranscript = $config->{uniquetranscript};
	my $mincount = $config->{mincount};
	if(!$mincount)
	{
		$mincount = 1;
	}

	my $deftagtypes;
	
	if($config->{"Primary Sense Tag"})
	{
		push(@{$deftagtypes}, "Primary Sense Tag");
	}

        if($config->{"Alternate Sense Tag"})
        {
                push(@{$deftagtypes}, "Alternate Sense Tag");
        }
        if($config->{"Primary Antisense Tag"})
        {
                push(@{$deftagtypes}, "Primary Antisense Tag");
        }
        if($config->{"Alternate Antisense Tag"})
        {
                push(@{$deftagtypes}, "Alternate Antisense Tag");
        }
        if($config->{"Unknown"})
        {
                push(@{$deftagtypes}, "Unknown");
        }

	if(!$config->{tagtypes})
	{
		$deftagtypes = ["Primary Sense Tag", "Alternate Sense Tag", "Primary Antisense Tag", "Alternate Antisense Tag", "Unknown"];
	}
	


	return 
	table(TR({-class=>'searchtitle'},
        	th({-colspan=>2,-align=>'LEFT'},
	                'Enter values to restrict what tags are shown, These values will be saved while browsing.'
	                )),
	        TR({-class=>'searchbody'},
	             td('Minumum Sequence Count:'),
	             td(textfield( -name=>'SageResults.mincount', -default=>$mincount))
		),
	        TR({-class=>'searchbody'},
	             td('Show Values as'),
	             td(popup_menu( -name=>'SageResults.valtype', -values =>["Percent", "Raw Count"], -default=>$valtype))
		),
	        TR({-class=>'searchbody'},
	             td('Show only Unique to Transcript'),
	             td(popup_menu( -name=>'SageResults.uniquetranscript', -values=>['No', 'Yes'], -default=>$uniquetranscript))
		),
		TR({-class=>'searchbody'},
	             td('Show'),
	             td(checkbox_group( -name=>'SageResults.tagtypes', -values =>["Primary Sense Tag", "Alternate Sense Tag", "Primary Antisense Tag", "Alternate Antisense Tag", "Unknown"], -defaults=>$deftagtypes))
		),
	    );
}
  

sub annotate {
	my $self = shift;
	my $segment = shift;
	my $config  = $self->configuration;

	my $ref        = $segment->ref;
	my $abs_start  = $segment->start;
	my $abs_end    = $segment->end;
 
         
	my $feature_list   = Bio::Graphics::FeatureFile->new;
	my $type = 'sagetag';


	if(!$config->{uniquetranscript})
	{
		$config->{uniquetranscript} = "No";
	}
	if(!$config->{mincount})
	{
		$config->{mincount} = 1;
	}
	if(!$config->{valtype})
	{
		$config->{valtype} = "Percent";
	}

	my $uniquetranscript = $config->{uniquetranscript};
	my $mincount = $config->{mincount};
	my $tagtypes = $config->{tagtypes};
	my $valtype = $config->{valtype};
	my @tagtypes = $config->{tagtypes};

	my $mbl = Mbl::new(undef, $db_name);
	my $dbh = $mbl->dbh();


	if(!$config->{"Primary Sense Tag"} && !$config->{"Alternate Sense Tag"} && !$config->{'Primary Antisense Tag'} && !$config->{'Alternate Antisense Tag'} && !$config->{'Unknown'})
	{
		$config->{"Primary Sense Tag"} = 1;
		$config->{"Alternate Sense Tag"} = 1;
		$config->{'Primary Antisense Tag'} = 1;
		$config->{'Alternate Antisense Tag'} = 1;
		$config->{'Unknown'} = 1;
	}
	


	my %session;
        eval {tie %session, "Apache::Session::File", cookie('SESSION_ID' . '_' . $mbl->organism), { Directory => "/var/www/sessions/sessions" }; };
	if($@)
	{
		return undef;
	}

	my $login_id = 0;

	if(!$session{login_id})
	{
		$login_id = undef;;
	} else
	{
		$login_id = $session{login_id}
	}

	$feature_list->add_type('sagetag'=>{glyph   => 'arrow',
                                    key     => "sagetag",
                                    fgcolor => 	sub 
						{
							my $feature = shift;
							my @nums = $feature->name =~ /\d+/g;
							my $num_count = 0;
							foreach my $this_val (@nums)
							{
								$num_count++;
								if( ($num_count == 2) && ($this_val > 1) )
								{
									return 'yellow'
								} elsif ($this_val > 10 && $num_count > 2)
								{
									return 'blue';
								}
							}
							return 'orange';
						},
                                    bgcolor => "purple",
                                    point   => 0,
				    linewidth => 2,
                                    orient  => 'N',
                                   });



	# Find out how many places this sage tag maps
 	my $db    = $self->database or die "I do not have a database";

        # pull out all sagetag features

        my @sagetags = $segment->features( -type=>'sagetag');
                                                                                           
        for my $o (@sagetags) {
	        my @feature_list = $db->get_feature_by_name('sagetag' => $o->name);
		my $mapid = $o->attributes('sagemapid');

		my $showme = 1;

		if($mbl->sage_tag_max_expr($o->name, $login_id) < $mincount)
		{
			$showme = 0;
			next;
		}
		my $orfinfo;
		$orfinfo = $mbl->get_sage_orf_info($o->name);
		my $tagtype;
		if($orfinfo)
		{
			$tagtype = $orfinfo->{tagtype};
		} else
		{
			$tagtype = "Unknown";
		}

		if( $uniquetranscript eq "Yes" )
		{
			if($orfinfo)
			{
				if($orfinfo->{tagmapid} == $mapid)
				{
					# Do Nothing
				} else
				{
					$showme = 0;
					next;
				}
				# Do nothing
			} else
			{
				$showme = 0;
				next;
			}
		}

		if($config->{$tagtype})
		{
			# Do Nothing
		} else
		{
			$showme = 0;
			next;
		}

		$showme = 1;
		if($showme)
		{
			my $this_name = $mbl->get_sage_description_line($o->name, $valtype, $login_id);
		      	my $url = $sage_url  . $o->name;
			my $strand;
			if($o->strand == 1)
			{
				$strand = +1;
			} else
			{
				$strand = -1;
			}
		      	my $feature = Bio::Graphics::Feature->new(-start=>$o->start,-stop=>$o->stop,-strand=>$strand,-ref=>$ref,-name=>$this_name, -url=>$url, -glyph=>"anchored_arrow");
		      	$feature_list->add_feature($feature,'sagetag');
		}
	} # End For loop


	untie(%session);
	return $feature_list;
}


1;

