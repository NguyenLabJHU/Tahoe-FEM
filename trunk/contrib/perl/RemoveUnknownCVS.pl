#!/usr/bin/perl
#
# run "cvs -q -n update" and prompt user to delete unrecognized files
#
open (IN, "cvs -q -n update |");

while ($line = <IN>)
{
	chomp($line);
	if ($line =~ s/^(\?\s)(.*)/$2/)
	{
		print STDOUT "remove \"$2\" (y|n|q): ";
		$ans = <STDIN>;
		if ($ans =~ /y/) { 
			system("rm $line"); 
		}
		elsif ($ans =~ /q/) { 
			exit; 
		}
	}
}
