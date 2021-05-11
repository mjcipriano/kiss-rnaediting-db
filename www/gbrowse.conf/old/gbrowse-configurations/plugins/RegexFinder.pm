package Bio::Graphics::Browser::Plugin::RegexFinder;
# $Id: RegexFinder.pm,v 1.1.1.1 2005/06/28 22:10:29 mcipriano Exp $

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


@ISA = qw(Bio::Graphics::Browser::Plugin);

my $nt_fasta_file;

sub name { "Sequences with Regular Expression" }

sub description {
  p("The Regular expression finder plugin finds sequences that match a particular regular expression with GMOD.",
  p("This plugin was written by Michael Cipriano."));
}

sub type { 'finder' }
sub init {
    my $self = shift;
    my $conf = $self->browser_config;
    $nt_fasta_file = $conf->plugin_setting('nt_fasta_file');
}

sub config_defaults {
  my $self = shift;

}

# we have no stable configuration
# sub reconfigure { }

sub configure_form {
  my $self = shift;
  return 
    table(TR({-class=>'searchtitle'},
	     th({-colspan=>2,-align=>'LEFT'},
		'Enter a regular expression below',
		'The browser will identify all genomic regions that match that regular expression',
		)),
	  TR({-class=>'searchbody'},
	     td('Enter regex:'),
	     td(textfield(-name=>'RegexFinder.searchregex'))),
          TR({-class=>'searchbody'},
             td('Search Type(AA translation type not working yet)'),
             td(popup_menu( -name=>'RegexFinder.searchtype', -values =>["Contigs as NT", "Contigs as AA 6 frame translation"]))),

);
}

# find() returns undef unless the SequenceFinder.searchsequence parameter
# is specified and valid.  Returning undef signals the browser to invoke the
# configure_form() method.
# If successful, it returns an array ref of Bio::SeqFeatureI objects.
sub find {
	my $self     = shift;
	my $segments = shift; 
	my @results;
	my $seqstring = '';

	my $regex =  param('RegexFinder.searchregex');
	my $search_type =  param('RegexFinder.searchtype');
	my $config  = $self->configuration;
	my $type = 'generic';

	my $sequences = Bio::SeqIO->new('-file'         => "$nt_fasta_file",
        	                        '-format'       => "fasta"); 

	warn($regex);
	warn($search_type);
	
	while(my $seq = $sequences->next_seq)
	{
		my $name = '';
		my $ref = $seq->id;
		my $desc = "regex match";
		my $start = 0;
		my $end = 0;
		my $contig_length = $seq->length();

		# Check the contigs as nt sequence
			

		if($search_type eq "Contigs as NT")
		{
			my $dir = 1;
			$seqstring = $seq->seq;
			while( $seqstring =~ /($regex)/gi)
			{
			        $desc = $1;
			        $start =  pos($seqstring) - length($1) + 1;
				$end = pos($seqstring) + 1;
	
				push @results, Bio::Graphics::Feature->new(-ref   => $ref,
		                	                                   -type  => $type,
		                        	                           -name  => $name,
									   -strand =>$dir,
	                                		                   -start => $start,
		                                                	   -end   => $end,
			                                                   -desc  => $desc);
			}
			# Check the reverse complement
			$dir = -1;
			$seqstring = $seq->revcom()->seq();
	
			while( $seqstring =~ /($regex)/gi)
	                {
	                        $desc = $1;
				my $start = ( length($seqstring) - pos($seqstring) ) + 1 ;
				my $end = $start + length($1) - 1;
	                                                                                                                                                                                                                                                       
	                        push @results, Bio::Graphics::Feature->new(-ref   => $ref,
	                                                                   -type  => $type,
	                                                                   -name  => $name,
	                                                                   -strand =>$dir,
	                                                                   -start => $start,
	                                                                   -end   => $end,
	                                                                   -desc  => $desc);
        	        }
		} elsif ($search_type eq "Contigs as AA 6 frame translation")
		{
			my $dir = 1;
			my $offset = 0;
			
			# First Frame
			$seqstring = $seq->translate()->seq;
			while( $seqstring =~ /($regex)/gi)
                        {
                                $desc = $1 . " - First Frame";;
                                $start =  ( (pos($seqstring) - length($1) )*3 ) + 1 + $offset ;
                                $end = ( pos($seqstring) * 3) + $offset;
                                                                                                                                                                                                                                                       
                                                                                                                                                                                                                                                       
                                push @results, Bio::Graphics::Feature->new(-ref   => $ref,
                                                                           -type  => $type,
                                                                           -name  => $name,
                                                                           -strand =>$dir,
                                                                           -start => $start,
                                                                           -end   => $end,
                                                                           -desc  => $desc);
                        } # End First Frame
			$offset = 1;
			$seqstring = $seq->trunc(2, $seq->length())->translate()->seq();
			while( $seqstring =~ /($regex)/gi)
                        {
                                $desc = $1 . " - Second Frame";
                                $start =  ( (pos($seqstring) - length($1) )*3 ) + 1 + $offset ;
                                $end = ( pos($seqstring) * 3) + $offset;
                                                                                                                                                                                                                                                       
                                                                                                                                                                                                                                                       
                                push @results, Bio::Graphics::Feature->new(-ref   => $ref,
                                                                           -type  => $type,
                                                                           -name  => $name,
                                                                           -strand =>$dir,
                                                                           -start => $start,
                                                                           -end   => $end,
                                                                           -desc  => $desc);
                        } # End Second Frame

                        $offset = 2;
                        $seqstring = $seq->trunc(3, $seq->length())->translate()->seq();

                        while( $seqstring =~ /($regex)/gi)
                        {
                                $desc = $1 . " - Third Frame";
                                $start =  ( (pos($seqstring) - length($1) )*3 ) + 1 + $offset ;
                                $end = ( pos($seqstring) * 3) + $offset;
                                                                                                                                                                                                                                                       
                                                                                                                                                                                                                                                       
                                push @results, Bio::Graphics::Feature->new(-ref   => $ref,
                                                                           -type  => $type,
                                                                           -name  => $name,
                                                                           -strand =>$dir,
                                                                           -start => $start,
                                                                           -end   => $end,
                                                                           -desc  => $desc);
                        } # End Third Frame



			# Now reverse complement
			$dir = -1;
			my $offset_extra;


                        # First Reverse Frame
                        $seqstring = $seq->revcom->translate()->seq;

			$offset_extra = $contig_length % 3;
                        while( $seqstring =~ /($regex)/gi)
                        {
                                $desc = $1 . " - Rev First Frame";
                                $start =  ( (pos($seqstring) - length($1) )*3 ) + 1 + $offset ;
                                $end = ( pos($seqstring) * 3) + $offset;

				my $tmpstart = $start;
				$start = length($seqstring)*3 - $end + $offset_extra;
				$end = length($seqstring)*3 - $tmpstart + $offset_extra;
                                $name = $offset_extra . " " . $offset;



                                push @results, Bio::Graphics::Feature->new(-ref   => $ref,
                                                                           -type  => $type,
                                                                           -name  => $name,
                                                                           -strand =>$dir,
                                                                           -start => $start,
                                                                           -end   => $end,
                                                                           -desc  => $desc);
                        } # End First Frame

			# Start Second Reverse Frame
                        $offset = 1;
			$offset_extra = ($contig_length - 1) % 3 + 2;
                        $seqstring = $seq->revcom->trunc(2, $seq->length())->translate()->seq();
                        while( $seqstring =~ /($regex)/gi)
                        {
                                $desc = $1 . " - Rev";
                                $start =  ( (pos($seqstring) - length($1) )*3 ) + 1 + $offset ;
                                $end = ( pos($seqstring) * 3) + $offset;

                                my $tmpstart = $start;
                                $start = length($seqstring)*3 - $end + $offset_extra;
                                $end = length($seqstring)*3 - $tmpstart + $offset_extra;
                                $name = $offset_extra . " " . $offset;


                                push @results, Bio::Graphics::Feature->new(-ref   => $ref,
                                                                           -type  => $type,
                                                                           -name  => $name,
                                                                           -strand =>$dir,
                                                                           -start => $start,
                                                                           -end   => $end,
                                                                           -desc  => $desc);
                        } # End Second Frame

			# Start Third Reverse Frame
                        $offset = 2;
			$offset_extra = ($contig_length - 2) % 3 + 3;
                        $seqstring = $seq->revcom->trunc(3, $seq->length())->translate->seq();
                        while( $seqstring =~ /($regex)/gi)
                        {
                                $desc = $1 . " - Rev";
                                $start =  ( (pos($seqstring) - length($1) )*3 ) + 1 + $offset ;
                                $end = ( pos($seqstring) * 3) + $offset;

                                my $tmpstart = $start;
                                $start = length($seqstring)*3 - $end + $offset_extra;
                                $end = length($seqstring)*3 - $tmpstart + $offset_extra;
				$name = $offset_extra . " " . $offset;



                                push @results, Bio::Graphics::Feature->new(-ref   => $ref,
                                                                           -type  => $type,
                                                                           -name  => $name,
                                                                           -strand =>$dir,
                                                                           -start => $start,
                                                                           -end   => $end,
                                                                           -desc  => $desc);
                        } # End Third Frame


		}

	}

                                                                                           
	return \@results;

}

sub auto_find {

}


1;
