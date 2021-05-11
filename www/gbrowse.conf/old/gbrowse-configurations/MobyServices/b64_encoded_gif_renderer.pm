package MobyServices::b64_encoded_gif_renderer;
use strict;
use XML::DOM;
our @ISA = qw(Exporter);
use File::Temp qw/ tempfile /;
#our @EXPORT = qw(render type);
our @EXPORT_OK = qw(render type);
use MIME::Base64;

sub type {
    return "b64_encoded_gif";
}

sub render {
    my ($DOM, $htmldir,$imgdir) = @_;
    my $content;
    foreach my $subnode($DOM->getChildNodes){
        next unless  (($subnode->getNodeType == TEXT_NODE) || ($subnode->getNodeType == CDATA_SECTION_NODE));
        $content .=$subnode->toString;
    }
    $content =~ s/^\s+//; $content =~ s/\s+$//;
    my $bindata = decode_base64($content);
    my ($fh, $filename) = tempfile( DIR => "$htmldir/$imgdir/", SUFFIX=> ".gif" );
    $filename =~ /(\w+\.jpg)$/;
    $filename = $1;
    binmode($fh);
    print $fh $bindata;
    close $fh;
    return ("<img src='$imgdir/$filename'>", 0);  # the 0 indicates that we have only rendered the top-level XML of this object
}

1;


=head1 NAME

b64_encoded_gif_renderer.pm - a renderer (HTML) for b64_encoded_gif type MOBY Objects

=head1 AUTHOR

Please report all bugs to Mark Wilkinson (markw at illuminae.com)

=head1 SYNOPSIS

just put the renderer in your gbrowse.conf/MobyServices folder
and it will work.

=head1 DESCRIPTION

This renderer returns HTML that fits between the
<td> tags in a web-page to display the content
of a b64_encoded_gif (or ontological child of) object.

=head1 METHODS

The module has two methods:

=over

=item type

this returns a scalar indicating the MOBY Object
Class that this renderer is designed to handle.  Objects
of this type, or objects that inherit from this type,
will be passed to this renderer.

=item render

This is called with three pieces of data which may or may not
be useful to your script:
    
=over

=item 1)  $data - an XML::DOM object
representing the deserialized MOBY object


=item 2)  $htmldir - the full path to the directory serving your html

 e.g. /usr/local/apache/htdocs/
 
 (this is the HTDOCS parameter you specified
  when you installed Gbrowse)

=item 3)  $imgdir - the additional path information to a writable directory for images

 e.g. /gbrowse/tmp
 
 (this is the folder specified in the tmpimages parameter
  in your organism.conf file)


=back

=back

=head1 RETURNS

The subroutine should return two pieces of data:

=over

=item 1)  An HTML representation of the Object

 this will appear between <td></td> tags in the webpage

=item 2)  A boolean indicating whether the renderer has parsed all sub-objects, or just the top level object

 '1' indicates that the renderer has fully parsed the Object XML
 '0' indicates that you need Gbrowse to follow child objects
     and render them independently

=back

=cut