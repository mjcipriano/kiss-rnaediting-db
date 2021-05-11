package Bio::Graphics::Browser::Plugin::SageGraph;
# $Id: SageGraph.pm,v 1.1.1.1 2005/06/28 22:10:29 mcipriano Exp $
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


sub name { "SAGE Graph" }

sub description {
	return "This track allows a user to visualize SAGE tags that map to a particular region of the genome.
<p>
Blue coloring represents tags assigned to an ORF in the current view, slategray coloring represents tags assigned to an ORF elsewhere in the assembly, and red coloring represents tags that have been classified as Unknown (not resolved to an ORF).
<p>
Each tag is given an ID number, tag mapping type (PS=primary sense, PA=primary antisense, AS=alternate sense, AA=alternate antisense, UK=unknown), and assigned ORF ID number.
<p>
Hovering your mouse over a tag reveals expression profile information.
<p>
This plugin was written by Michael Cipriano";

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
	$current_config->{mincount} = param('SageGraph.mincount');
	$current_config->{valtype} = param('SageGraph.valtype');
	$current_config->{uniquetranscript} = param('SageGraph.uniquetranscript');
	$current_config->{mincount} = param('SageGraph.maxgraph');
	$current_config->{tagtypes} = 1;
	my @tagtypes = param('SageGraph.tagtypes');

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
	my $maxgraph = $config->{maxgraph};

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
	             td(textfield( -name=>'SageGraph.mincount', -default=>$mincount))
		),
	        TR({-class=>'searchbody'},
	             td('Show Values as'),
	             td(popup_menu( -name=>'SageGraph.valtype', -values =>["Percent", "Raw Count"], -default=>$valtype))
		),
	        TR({-class=>'searchbody'},
	             td('Show only Unique to Transcript'),
	             td(popup_menu( -name=>'SageGraph.uniquetranscript', -values=>['No', 'Yes'], -default=>$uniquetranscript))
		),
                TR({-class=>'searchbody'},
                     td('Minumum Sequence Count:'),
                     td(textfield( -name=>'SageGraph.maxgraph', -default=>$maxgraph))
                ),
		TR({-class=>'searchbody'},
	             td('Show'),
	             td(checkbox_group( -name=>'SageGraph.tagtypes', -values =>["Primary Sense Tag", "Alternate Sense Tag", "Primary Antisense Tag", "Alternate Antisense Tag", "Unknown"], -defaults=>$deftagtypes))
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
 
         

	if(!$config->{uniquetranscript})
	{
		$config->{uniquetranscript} = "No";
	}
	if(!$config->{mincount})
	{
		$config->{mincount} = 1;
	}

	if(!$config->{maxgraph})
	{
		$config->{maxgraph} = 3000;
	}
	
	if(!$config->{valtype})
	{
		$config->{valtype} = "Percent";
	}

	my $uniquetranscript = $config->{uniquetranscript};
	my $mincount = $config->{mincount};
	my $maxgraph = $config->{maxgraph};
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
        eval {tie %session, "Apache::Session::File", cookie('SESSION_ID' . '_' . $mbl->organism), { Directory => $mbl->session_tmp_dir }; };
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

	my $feature_list = Bio::Graphics::FeatureFile->new;

	my $libh = $mbl->query('get_library_access');
        $libh->execute($login_id);
	my @libs;

        while(my $row = $libh->fetchrow_hashref)
	{
		push(@libs, $row->{library});
	}


	my $type = 'SAGEGraph';
	my $parent_features_hash;
	$maxgraph = 10;
	foreach my $lib (@libs)
	{
		my $lib_name = $mbl->get_sage_library_short_name($lib);
		my $type = $lib;	
		$feature_list->add_type($type => {glyph => 'xyplot',
	                                    graph_type => 'boxes',
	                                    scale => 'both',
	                                    height => 60,
	                                    key => 'SAGE Library ' . $lib_name,
	                                    fgcolor =>'blue',
	                                    bgcolor => 'blue',
	                                    min_score => 0,
	                                    max_score => $maxgraph,
	                                    });
	
		$parent_features_hash->{$lib} = Bio::Graphics::Feature->new(-start => $abs_start,
	                                                   -end => $abs_end,
	                                                   -ref => $segment->ref,
	                                                   -name => $type,
	                                                   -source => 'UCSC',
	                                                   -type => $type,
	                                                   -configurator => $feature_list);
	}


	# Find out how many places this sage tag maps
 	my $db    = $self->database or die "I do not have a database";

        # pull out all sagetag features

        my @sagetags = $segment->features( -type=>'sagetag');
	my %feature_array_hash;


        for my $o (@sagetags) 
	{
	        my @feature_list = $db->get_feature_by_name('sagetag' => $o->name);
		my $mapid = $o->attributes('sagemapid');

		my $glyphtype = 'sagetag';
		my $showme = 1;
		my $max_exp = $mbl->sage_tag_max_expr($o->name, $login_id);
		if($mbl->sage_tag_max_expr($o->name, $login_id) < $mincount)
		{
			$showme = 0;
			next;
		}
		
		my $orfinfo;
		$orfinfo = $mbl->get_sage_orf_info($o->name);
		my $tagtype;
		my $orf;
		if($orfinfo)
		{
			$tagtype = $orfinfo->{tagtype};
			$orf = $mbl->get_orf_attributes_hash($orfinfo->{orfid});
		} else
		{
			$tagtype = "Unknown";
		}

		# Check if this is assigned to this location only
		if($orfinfo)
		{
			# if it is assigned to an orf and this orf is here
			if($orfinfo->{tagmapid} == $mapid)
			{
				$glyphtype = 'sagetag';
			} else # the orf is somewhere else
			{
				$glyphtype = 'sagetag_other';
			}
		} else
		{
			$glyphtype = 'sagetag_unknown';
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

		if($showme)
		{
			my $scores = $mbl->get_sage_results_hash($o->name, 'percent', $login_id);
			foreach my $lib(@libs)
			{
				my $type = $lib;
				my $score = log(1+($scores->{$lib}*1000));
				if($score > $maxgraph)
				{
					$score = $maxgraph;
				}
				if($score < 1)
				{
					$score = 1;
				}
			      	my $feature = Bio::Graphics::Feature->new(-start=>$o->start,-stop=>$o->stop,-ref=>$ref,-name=>$type, -score=>$score, -subtype=>$type);
				push(@{$feature_array_hash{$lib}}, $feature);
			}
		}
	} # End For loop

	foreach my $lib(@libs)
	{
		my $type = $lib;
		$parent_features_hash->{$lib}->add_segment(@{$feature_array_hash{$lib}});
		$feature_list->add_feature($parent_features_hash->{$lib}, $type);
	}
	untie(%session);

	return $feature_list;


}


1;

