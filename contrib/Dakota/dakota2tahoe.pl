#!/usr/bin/perl -w

# should have 3 command line arguments
if (scalar(@ARGV) != 3) {
	print STDOUT "\tusage: perl dakota2tahoe.pl [DAKOTA parameter file] [XML template file] [output file name]\n";
	exit;
}

# command line arguments
($dakota, $source, $outfile) = @ARGV;

open (DAKOTA, "$dakota") || die "could not open $dakota";

# read number of parameters
defined($line = <DAKOTA>)|| die "file error";
chomp($line);
$num_vars = $line;
$num_vars =~ s/[\s]+([0-9]+)[\s]+.*/$1/;
print STDOUT "found $num_vars variables\n";

# insert parameters into a hash
for ($i = 0; $i < $num_vars; $i++) {
	defined($line = <DAKOTA>) || die "file error";
	chomp($line);
	$line =~ s/^\s+//; # drop leading white space
	($val, $var) = split(/[\s]+/, $line);
	$varlist{$var} = $val;
}
close(DAKOTA);

# print them out
foreach $key (keys %varlist) {
	print STDOUT "$key = $varlist{$key}\n";
}

print "\n";

# echo file making variable subsitutions
open(IN, "$source") || die "could not open input file $source";
open(OUT, ">$outfile") || die "could not open output file $outfile";
while (defined($line = <IN>)) {
	chomp($line);
	# variable substitution
	foreach $key (keys %varlist) {
		$line =~ s/["']$key["']/"$varlist{$key}"/g;
	}
	# substitution & arithmatic (one expression per line)
 	if ($line =~ /{/) {
 		($begin,$end) = split /{/, $line;
 		($exp,$end)   = split /}/, $end;
 		print "exp : $exp - > ";
		foreach $key (keys %varlist) {
                	$exp =~ s/$key/$varlist{$key}/g;
        	}
 		$nexp = eval($exp);
 		print " $nexp\n";
 		$line = "$begin"."$nexp"."$end";
 	}
	
	
	# echo
	print OUT "$line\n";
}
close(OUT);
close(IN);
