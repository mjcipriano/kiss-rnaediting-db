#!/usr/bin/perl


foreach (@ARGV)
{
	print $_ . "\n\n";
	my $filename = $_;
	my $new_filename = "temp/$filename";
	open(FILE, $filename);
	open(NEWFILE, ">", $new_filename);
	while(<FILE>)
	{
		my $line = $_;
		if($line =~ /PATH/)
		{
			# Do nothing
		} else
		{
			print NEWFILE $line;
		}
	
	}
}
