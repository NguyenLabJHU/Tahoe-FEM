#!/usr/bin/perl
# $Id: CVSLogPick.pl,v 1.2 2003-06-09 06:23:01 paklein Exp $
#
# Print the non-empty log messages to stdout
#
# To generate a cvs log for all changes committed after a given date D:
#
#    cvs log -N -b -d">Feb 4, 2003" > & cvs.log
#
if (scalar(@ARGV) == 0) {
	print STDOUT "\n usage: CVSLogPick.pl [cvs log file]\n\n";
	exit;	
}

$cvs_log = $ARGV[0];
open(IN, $cvs_log) || die "could not open file $cvs_log";

$working_file = "none";
$echo = 0;
while ($line = <IN>) {

	# echo
	if ($echo == 1) {

		print STDOUT $line;

		# end echo
		if ($line =~ /=====/) { $echo = 0; }
	}
	elsif ($line =~ /Working\ file/) # looking for working file
	{
		$working_file = $line;
	}

	# look for selected revisions
	if ($line =~ /selected\ revisions:\ [1-9]/) {
		$echo = 1;
		print STDOUT $working_file;
		
		# discard 2 lines
		$line = <IN>;
		$line = <IN>;		
	}	
}

close(IN);
