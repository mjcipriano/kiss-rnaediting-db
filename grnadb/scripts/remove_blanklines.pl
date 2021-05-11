#!/usr/bin/perl

use strict;


while(<>)
{
	if($_ =~ m/\w/)
	{
		print $_;
	} else
	{
		# Do nothing
	}
}
