#!/usr/bin/perl
# $Id: extract_KBC_info.pl,v 1.1 2003-10-04 00:41:21 paklein Exp $
#
# This script scan the given files for the position and force information
# written by KBC controllers and echos the results to a .force file. The files
# are assumed to be tahoe .out files produced from a 3D calculation.
#

if (scalar(@ARGV) < 1) {
	print STDOUT "\tusage: extra_KBC_info.pl out_file1 [out_file2 out_file3 ...]\n";
	exit;
}

# find the next occurrence of $clue in the given $stream
# and return the line containing it, or FAILED if it is not
# found
sub ffw {
	my $stream = shift(@_);
	my $clue   = shift(@_);
	
	my $line = 0;
	my $found = 0;
	while ($found == 0 && defined($line = <$stream>)) {
		chomp($line);
		if ($line =~ /$clue/) {
			$found = 1;
		}	
	}

	# not found	
	if ($found == 0) { $line = "FAILED"; }

	return $line;
}

foreach $file (@ARGV) {
	print STDOUT "processing file: $file\n";

	# open streams
	open (IN, $file) || die "Error: could not open file \"$file\"\n";
	open (OUT, ">".$file.".force") || die "Error: could not create file \"$file.force\"\n";

	# fast forward to results section
	$line = ffw(IN, "T i m e");
	if ($line =~ /FAILED/) { 
		die "Error: file \"$file\" contains no results\n"; 
	} else {
		print STDOUT "found results in file \"$file\"\n";
	}

	# echo results
	$line = ffw(IN, "Time =");
	while ($line !~ /FAILED/) {

		# write the time	
		@line = split(/\s/, $line);
		print OUT $line[3];
		
		# write position
		ffw(IN, "Position . . .");
		if (defined($value = <IN>)) { 
			chomp($value);
			print OUT $value; 
		}
		if (defined($value = <IN>)) { 
			chomp($value);
			print OUT $value; 
		}
		if (defined($value = <IN>)) { 
			chomp($value);
			print OUT $value; 
		}

		# write contact force
		ffw(IN, "Contact force");
		if (defined($value = <IN>)) { 
			chomp($value);
			print OUT $value; 
		}
		if (defined($value = <IN>)) { 
			chomp($value);
			print OUT $value; 
		}
		if (defined($value = <IN>)) { 
			chomp($value);
			print OUT $value; 
		}
	
		# next
		print OUT "\n";
		$line = ffw(IN, "Time =");
	}

	# clean up
	close(IN);
	close(OUT);
}

