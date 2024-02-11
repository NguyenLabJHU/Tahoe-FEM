#!/usr/bin/perl -w
if (scalar(@ARGV) == 2) {
  ($in, $out) = @ARGV;
} else { die; }

#@runs = ("rate350","creep100","creep500"); # all
#@runs = ("rate350"); # for stiffness
#@runs = ("creep100","creep500"); # for viscosity
@runs = ("creep100"); # for viscosity

$errors = 0.0;

foreach $run (@runs) {
	@items= split(/\./,$in); $tag = $items[1];
	$stem = "$run.$tag";
	$stdout = "$stem.stdout";
	if (-e "$stdout") { unlink "$stdout"; }
	print "running $run:$tag ";
# replace variables
	system("dpp $in $run.xml.tmpl $stem.xml parameter.defaults > $stdout");
 	system("tahoe -f $stem.xml >> $stem.stdout");
	system("extract_1D $stem.io1.exo $stem.io1.exo $stem.uvst_Y >> $stdout");
#	$cmd = "plot \"$run.dat\" w l, \"$stem.uvst_Y\" w l; pause 3";
#	system("echo \' $cmd \' | gnuplot - > /dev/null");
	system("least_square $run.dat $stem.uvst_Y $stem.error >> $stdout");
#	open(ERR,"cat $run.error|");
	open(ERR,"$stem.error");
	$line = <ERR>; $line =~ s/^\s+//; @items = split(/\s+/,$line); close ERR;
	$error = $items[0];
	print " error : $error \n";
	$errors += $error;
}

read_parameters();
#print "alpha 1 ", $varlist{"ALPHA1"}, "\n";
#print "alpha 2 ", $varlist{"ALPHA2"}, "\n";
#print "alpha 3 ", $varlist{"ALPHA3"}, "\n";
$constraint = $varlist{"ALPHA1"} + $varlist{"ALPHA2"} + $varlist{"ALPHA3"};

open(OUT,">$out");
print OUT "$errors ERROR\n";
print OUT "$constraint CONSTRAINT\n";
close OUT;

#$error_line = "$error error";
#system("echo $error_line > $out");

#---------------------------------------------------------------------------
sub read_parameters {
	open(IN,"$in");
  # read number of parameters
  defined($line = <IN>)|| die "file error";
  chomp($line);
  $num_vars = $line;
  $num_vars =~ s/[\s]*([0-9]+)[\s]+.*/$1/;
  print STDOUT "found $num_vars variables\n";

  # insert parameters into a hash
  for ($i = 0; $i < $num_vars; $i++) {
    defined($line = <IN>) || die "file error";
    chomp($line);
    $line =~ s/^\s+//; # drop leading white space
    ($val, $var) = split(/[\s]+/, $line);
    $varlist{$var} = $val;
  }
	close IN;
}

