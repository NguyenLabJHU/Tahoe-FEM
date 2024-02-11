#!/usr/bin/perl -w

$in = "parameter.defaults";

#@runs = ("creep100","creep350","creep500"); 
@runs = ("creep500"); 
open(GPL,">creep.gpl");
print GPL "set logscale xy\nset xrange [0.5:*]\nplot \\\n";
@rss = ("rs0","rs1","rs2");
foreach $run (@runs) {
	$stem = $run;
	$stdout = "$stem.stdout";
	if (-e "$stdout") { unlink "$stdout"; }
	print "running $run: ";
  # run and extract data
	$rsp = "";
	foreach $rs (@rss) {
		system("dpp $in $run.$rs.xml.tmpl $stem.$rs.xml parameter.defaults $rsp >> $stdout");
		$rsp = $rs;
 		system("tahoe -f $stem.$rs.xml >> $stem.stdout");
  	system("extract_1D $stem.$rs.io1.exo $stem.$rs.io1.exo $stem.$rs.uvst_Y >> $stdout");
	}
  system("cat $stem.rs0.uvst_Y  $stem.rs1.uvst_Y $stem.rs2.uvst_Y > $stem.uvst_Y");
	print "... done\n";
	print GPL "\"$run.dat\" w p,\\\n";
	print GPL "\"$stem.uvst_Y\" w l,\\\n";
};
print GPL "pause -1\n";
print GPL "set term png color\nset output \"creep.png\"\nrepl\n";
print GPL "set term post eps color\nset output \"creep.eps\"\nrepl\n";
close GPL;
exit;

@runs = ("rate3_5","rate35","rate350"); 
#@runs = ("rate350"); 
foreach $run (@runs) {
	$stem = $run;
	$stdout = "$stem.stdout";
	if (-e "$stdout") { unlink "$stdout"; }
	print "running $run: ";
  # run and extract data
	system("dpp $in $run.xml.tmpl $stem.xml parameter.defaults > $stdout");
 	system("tahoe -f $stem.xml >> $stem.stdout");
 	system("extract_1D $stem.io1.exo $stem.io1.exo $stem.uvst_Y >> $stdout");
	print "... done\n";
	open(GPL,">$run.gpl");
	print GPL "set logscale xy\nset xrange [0.5:*]\nplot \\\n";
	print GPL "\"$run.dat\" w p,\\\n";
	print GPL "\"$stem.uvst_Y\" w l,\\\n";
	print GPL "pause -1\n";
	print GPL "set term png color\nset output \"$stem.png\"\nrepl\n";
	print GPL "set term post eps color\nset output \"$stem.eps\"\nrepl\n";
	close GPL;
}
