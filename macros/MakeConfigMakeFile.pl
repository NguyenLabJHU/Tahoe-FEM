#!/usr/bin/perl -w
#
# $Id: MakeConfigMakeFile.pl,v 1.5 2003-08-07 18:08:14 paklein Exp $
#
# Takes a configuration file and the output generated by MakeConfigHeaderFile.pl
# from that file and generates a make-format file which collects the directories
# which must be traversed to build the enabled options into the code.
#
# The format of the file fed to this script is:
#
# [root of output file names]
# 
# [root of preprocessor symbol for option 1]
# [ENABLE/DISABLE for option 1]
# [one line description for option 1]
# [source directories for option 1]
#
# [root of preprocessor symbol for option 2]
# [ENABLE/DISABLE for option 2]
# [one line description for option 2]
# [source directories for option 2]
#
# ...and so one with three lines of information for each option
#

if (scalar(@ARGV) == 0) {
	print STDOUT "\tusage: MakeConfigFile.pl [config file] [(optional) output directory]\n";
	exit;
}

$config_file = $ARGV[0];
$destination = ".";
if (scalar(@ARGV) > 1) { $destination = $ARGV[1]; }

print "reading configuration file: $config_file\n";

# skip leading comments
open(IN, $config_file) || die "could not open config file: $config_file\n";
$out_root = <IN>;
chomp($out_root);
while ($out_root =~ /^#/ || $out_root !~ /[a-zA-Z0-9]/ ) {
	$out_root = <IN>;
	chomp($out_root);
}

# look for corresponding header file created by MakeConfigHeaderFile.pl
print "output file root: $out_root\n";
$config_header_file = "$destination" . "/" . $out_root . ".h";
if (-e $config_header_file) {
	print "found header file: $config_header_file\n";
}
else {
	print "could not find header file: $config_header_file\n";
	print "run MakeConfigHeaderFile.pl $config_file\n";
	print "stop.\n";
	exit;
}
open(HEADER, "$config_header_file") || die "could not open file: $config_header_file\n";

# output make file
$config_make_file = "$destination" . "/" . $out_root . ".make";
if (-e $config_make_file) {
	print "output make file already exists: $config_make_file\n";
	print "stop.\n";
	exit;
}
else {
	print "creating output file: $config_make_file\n";
}
open(OUT, ">$config_make_file") || die "could not open output file: $config_make_file\n";

# print header
print OUT <<FIN;
# \$Id\$
# This file was generated by MakeConfigMakeFile.pl from $config_file
FIN
open (DATE, "date |") || die "could not get the date\n";
while (<DATE>) {
	print OUT "# created: ".$_;
}
close(DATE);

print OUT <<FIN;
#
# \\file $out_root.make
# Configuration of optional components within Tahoe.
# Sections of the code are included or excluded in the build of Tahoe depending in 
# this flags in this file and in the file $out_root.h. Each option has
# a #define definition in $out_root.h and a corresponding directory definition
# in this file. The two items must be set consistently to enable or
# disable materials models. To enable an option:
# -# in this file, uncomment the macro.
# -# in $out_root.h, uncomment the #define statement
#
# The naming convention for the definitions in this file and the macros in
# $out_root.h are as follows. For the option [OPTION]:
# -# the macro below defining the corresponding source directory will be DIRECTORY_[OPTION]
# -# the symbol in $out_root.h file will be [OPTION]
FIN

# scan config file for options
@option_list = ();
@directory_list = ();
$scan_line = 1;
$opt_root = "";
$is_active = 0;
while ($line = <IN>) {
	chomp($line);
	
	# non-comment, non-blank line
	if ($line !~ /^#/ && $line =~ /[a-zA-Z0-9]/ ) {
	
		# read name
		if ($scan_line == 1) {
			$opt_root = $line;
			push @option_list, $opt_root;
			print OUT "\n# \\def DIRECTORY_$opt_root\n";
			print "processing option: $opt_root: ";
			$scan_line++;
		}
		# read ENABLE/DISABLE
		elsif ($scan_line == 2) {
			$scan_line++;
		}
		# read description
		elsif ($scan_line == 3) {
		
			# add period
			if ($line !~ /\.$/) { $line = $line."."; }

			print OUT "# $line\n";
			print OUT "# This option must be set in conjunction with #define $opt_root\n";
			print OUT "# in $out_root.h. */\n";

			# scan header file for symbol
			$symbol = $opt_root;
			$found = 0;
			while ($found == 0 && ($line = <HEADER>)) {
			
				chomp($line);
				if ($line =~ /#.+$symbol/) {

					# not active
					if ($line =~ /#undef/ || $line !~ /^#define/) {
						print "INACTIVE\n";
						$is_active = 0;
					}
					else {					
						print "ACTIVE\n";
						$is_active = 1;
					}
					$found = 1;
				}
			}
			
			$scan_line++;
		}
		# read directories
		else {		
			if ($is_active == 0) {
				print OUT "#"; # comment line out
			}
			$dir_macro = "DIRECTORY_".$opt_root; 
			print OUT "$dir_macro = $line\n";
			push @directory_list, $dir_macro;
			$scan_line = 1;
		}
	}
}

# end of user-editable section
print OUT <<FIN;

#############################################################
############## no configuration options below ###############
# Unless you are adding another option, you should not need #
# to edit the contents of below.                            #
#############################################################

FIN

# create directory list macro name
$dir_macro = $out_root."Directories";
$dir_macro =~ s/([a-z])([A-Z])/$1\_$2/g; # replace inter-caps with "_"
$dir_macro = uc($dir_macro);

# print directory list macro
print OUT "$dir_macro = ";
for ($i = 0; $i < scalar(@directory_list); $i++) {
	$macro = $directory_list[$i];
	print OUT "\\\n\t\$($macro) ";
}
print OUT "\n";
