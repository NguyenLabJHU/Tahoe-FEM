#!/usr/bin/perl -w
#
# $Id: MirrorTree.pl,v 1.3 2004-11-09 21:27:15 paklein Exp $
#
# Mirror all contents of a directory tree using symbolic links.
#

if (scalar(@ARGV) < 2) {
	print STDOUT "\tusage: MirrorTree.pl\n\t\t[absolute path to root of source tree]\n\t\t[root of destination]\n\t\t[ln command (optional string)]\n";
	exit;
}

# source and destination directories
$source = $ARGV[0];
$dest = $ARGV[1];

# the 'ln' command
$ln = 'ln -sf';
if (scalar(@ARGV) > 2) { $ln = $ARGV[2]; }

use Cwd;

sub mirror_sub
{
	# resolve arguments
	my $source = shift(@_);
	my $dest   = shift(@_);
	my $ln     = shift(@_);

	# echo arguments 
	print STDOUT "     source directory: $source\n";
	print STDOUT "destination directory: $dest\n";
	print STDOUT "           ln command: $ln\n";

	# check source directory
	if (-d $source) {
		print STDOUT "found source directory $source\n";
		
		# path must be from root
		if ($source !~ /^\//) {
			print STDOUT "source directory $source must be an absolute path\n";
			exit;
		}	
	} else {
		print STDOUT "could not find source directory $source\n";
		exit;
	}

	# check destination directory
	if (-d $dest) {
		print STDOUT "found destination directory $dest\n";
	} else {
		print STDOUT "creating destination directory $dest\n";
		system ("mkdir -p $dest");
	}

	# get source directory contents - need to collect all files at top because
	# subdirectories are processes recursively
	opendir(CURR, $source) or die "ERROR: could not open directory: $source\n";
	my @all_files = grep {! /^\./ } readdir CURR;
	closedir CURR;

	# process current directory
	foreach $file (@all_files) {

		my $source_path = "$source/$file";
		my $dest_path = "$dest/$file";

		# directory
		if (-d $source_path) {
	
			# filter directories
			if ($file !~ /^CVS$/) {

				print STDOUT "processing directory $source_path\n";
	
				# create destination if needed
				if (-d $dest_path) {
					print STDOUT "found destination directory $dest_path\n";
				} else {
					print STDOUT "creating destination directory $dest_path\n";
					system ("mkdir -p $dest_path");
				}
		
				# process subdirectories
				mirror_sub($source_path, $dest_path, $ln);
			}
		}
		# plain file
		elsif (-f $source_path) {

			# filter files
			if ($file !~ /\.o$/) {

				# create symbolic link
				print STDOUT "$ln $source_path $dest_path\n";
				system ("$ln $source_path $dest_path");
			}	
		}
		else {
			print STDOUT "skipping file $file\n";
		}
	}
}

# call it
mirror_sub($source, $dest, $ln);
