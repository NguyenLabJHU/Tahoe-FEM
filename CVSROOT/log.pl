#!/usr/bin/perl
#
#
# Usage:  log.pl [-u user] [[-m mailto] ...] [-s] -f logfile 'dirname file ...'
#
#	-u user		- $USER passed from loginfo
#	-m mailto	- for each user to receive cvs log reports
#			(multiple -m's permitted)
#	-s		- to prevent "cvs status -v" messages
#	-f logfile	- for the logfile to append to (mandatory,
#			but only one logfile can be specified).

# here is what the output looks like:
#
#    From: woods@kuma.domain.top
#    Subject: CVS update: testmodule
#
#    Date: Wednesday November 23, 1994 @ 14:15
#    Author: woods
#
#    Update of /local/src-CVS/testmodule
#    In directory kuma:/home/kuma/woods/work.d/testmodule
#    
#    Modified Files:
#    	test3 
#    Added Files:
#    	test6 
#    Removed Files:
#    	test4 
#    Log Message:
#    - wow, what a test
#
# (and for each file the "cvs status -v" output is appended unless -s is used)
#
#    ==================================================================
#    File: test3           	Status: Up-to-date
#    
#       Working revision:	1.41	Wed Nov 23 14:15:59 1994
#       Repository revision:	1.41	/local/src-CVS/cvs/testmodule/test3,v
#       Sticky Options:	-ko
#    
#       Existing Tags:
#    	local-v2                 	(revision: 1.7)
#    	local-v1                 	(revision: 1.1.1.2)
#    	CVS-1_4A2                	(revision: 1.1.1.2)
#    	local-v0                 	(revision: 1.2)
#    	CVS-1_4A1                	(revision: 1.1.1.1)
#    	CVS                      	(branch: 1.1.1)

$cvsroot = $ENV{'CVSROOT'};

# turn off setgid
#
$) = $(;

$dostatus = 1;

# strings to replace repository name with URL to CVSWeb link:
# searches for $repo_name and replaces it with $cvsweb_URL:
$repo_short_name = "Sourceforge";
$repo_name = "/cvsroot/tahoe";
$cvsweb_URL = "http://tahoe.cvs.sourceforge.net/tahoe";

# parse command line arguments
#
while (@ARGV) {
        $arg = shift @ARGV;

	if ($arg eq '-m') {
                $users = "$users " . shift @ARGV;
	} elsif ($arg eq '-u') {
		$login = shift @ARGV;
	} elsif ($arg eq '-f') {
		($logfile) && die "Too many '-f' args";
		$logfile = shift @ARGV;
	} elsif ($arg eq '-s') {
		$dostatus = 0;
	} else {
		($donefiles) && die "Too many arguments!\n";
		$donefiles = 1;
		@files = split(/ /, $arg);
	}
}

# the first argument is the module location relative to $CVSROOT
#
$modulepath = shift @files;

# get a login name for the guy doing the commit....
#
if ($login eq '') {
	$login = getlogin || (getpwuid($<))[0] || "paklein";
      }
#
# bounce mail off sourceforge server forwarding
#
$login = $login."\@users.sourceforge.net";
#
#$mailcmd = "| Mail -s 'Tahoe CVS update: $modulepath'";
#
$mailcmd = "| mail -s 'Tahoe CVS update ($repo_short_name): $modulepath'";

# Initialise some date and time arrays
#
@mos = (January,February,March,April,May,June,July,August,September,
        October,November,December);
@days = (Sunday,Monday,Tuesday,Wednesday,Thursday,Friday,Saturday);

($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime;
$year += 1900;

# open log file for appending
#
#open(OUT, ">>" . $logfile) || die "Could not open(" . $logfile . "): $!\n";

# send mail, if there's anyone to send to!
#
if ($users) {
	$mailcmd = "$mailcmd $users";
	open(MAIL, $mailcmd) || die "Could not Exec($mailcmd): $!\n";
}

# print out the log Header
# 
#print OUT "\n";
#print OUT "****************************************\n";
#print OUT "  Date: $days[$wday] $mos[$mon] $mday, $year @ $hour:" . sprintf("%02d", $min) . "\n";
#print OUT "Author: $login\n\n";

if (MAIL) {
#	print MAIL "\n";
	print MAIL '$Id: log.pl,v 1.8 2006-05-27 01:05:57 paklein Exp $' . "\n";
	print MAIL "===================================================================\n";
	print MAIL "  Date: $days[$wday] $mos[$mon] $mday, $year @ $hour:" . sprintf("%02d", $min) . "\n";
	print MAIL "Author: $login\n\n";
}

# print the stuff from logmsg that comes in on stdin to the logfile
#
open(IN, "-");
$mod_dir = " ";
$file_info = 0;
# 0 = don't write
# 1 = write next time
# 2 = write this time

while (<IN>) {

	# get directory
	$next_line = $_;
	if ($next_line =~ s/(Update\sof\s)(.+)/$2/) {
		
		chomp($next_line);
		$next_line =~ s/$repo_name/$cvsweb_URL/;
		#$next_line =~ s/\/usr\/local\/cvsrep/https:\/\/tahoe\.ca\.sandia\.gov\/cgi-bin\/cvsweb\.cgi/;
		$mod_dir = $next_line;
	}

	# file_info switch
	if ($_ =~ /Modified\sFiles/) { $file_info = 1; }
	elsif ($_ =~ /Added\sFiles/) { $file_info = 1; }
	elsif ($_ =~ /Removed\sFiles/) { $file_info = 1; }
	elsif ($_ =~ /Vendor\sTag\:/) { $file_info = 0; }
	elsif ($_ =~ /Tag\:/)        { $file_info = 1; }
	elsif ($_ =~ /Log\sMessage/) { $file_info = 0; }

	# write file and cvsweb link
	if ($file_info == 2) {
		$the_line = $_;
		chomp($the_line);		
		@the_split = split(/\s/, $the_line);
		foreach $one_file (@the_split) {
			
			if ($one_file =~ /[a-zA-Z\.0-9_]+/) {
				if (MAIL) {
					print MAIL "    $one_file: <$mod_dir/$one_file>\n";
				}
			}
		}
	}

	#print OUT $_;
	if (MAIL && $file_info != 2) {
		print MAIL $_;
	}

	# rotate flag
	if ($file_info == 1) { $file_info = 2; }
}
close(IN);

#print OUT "\n";
if (MAIL) {
print MAIL "\n";
}

# after log information, do an 'cvs -Qq status -v' on each file in the arguments.
#
if ($dostatus != 0) {
	while (@files) {
		$file = shift @files;
		if ($file eq "-") {
			#print OUT "[input file was '-']\n";
			if (MAIL) {
				print MAIL "[input file was '-']\n";
			}
			last;
		}
		$pid = open(RCS, "-|");
		if ( !defined $pid )
		{
			die "fork failed: $!";
		}
		if ($pid == 0)
		{
			exec 'cvs', '-nQq', 'status', '-v', $file;
			die "cvs exec failed: $!";
		}
		while (<RCS>) {
			#print OUT;
			if (MAIL) {
				print MAIL;
			}
		}
		close(RCS);
	}
}

#close(OUT);
#die "Write to $logfile failed" if $?;

close(MAIL);
die "Pipe to $mailcmd failed" if $?;

## must exit cleanly
##
exit 0;
