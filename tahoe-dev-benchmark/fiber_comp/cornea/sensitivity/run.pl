#!/usr/bin/perl -w
if (scalar(@ARGV) == 2) {
  ($in, $out) = @ARGV;
} else { die; }

@runs = ("rate350","creep100","creep500"); # all
#@runs = ("rate350"); # for stiffness
#@runs = ("creep100","creep500"); # for viscosity

$wgs{"rate350"} = 4.0;
$wgs{"creep100"} = 1.0;
$wgs{"creep500"} = 1.0;

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
	$wg = $wgs{$run};
	$errors += $error*$wg;
}
open(OUT,">$out");
print OUT "$errors ERROR\n";
close OUT;

#$error_line = "$error error";
#system("echo $error_line > $out");
