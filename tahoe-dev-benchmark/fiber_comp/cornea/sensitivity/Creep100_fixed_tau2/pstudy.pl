#!/usr/bin/perl -w
$run = "creep100"; # for viscosity

$errors = 0.0;

#@as = (0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16);
#@as = (0.08,0.1,0.12,0.14,0.16,0.18,0.20,0.22,0.24);
@as = (0.12,0.14,0.16,0.18,0.20,0.22,0.24,0.26);
#@as = (0.161);
open(GPL,">pstudy.gpl");
print GPL "set logscale xy\nset xrange [0.5:*]\nplot \\\n";
foreach $a (@as) {
	$tag = "a=$a";
	$stem = "$run.$tag";
	$stdout = "$stem.stdout";
	if (-e "$stdout") { unlink "$stdout"; }
	print "running $run:$tag ";
  # run and extract data
	$in = "in.$tag";
	open(OUT,">$in");
	$b= 3.0*(0.1610-2.0/3.0*$a);
	print OUT "2 variables\n$a ALPHABAR\n$b ALPHA2\n";
  close OUT;
	system("dpp $in $run.xml.tmpl $stem.xml parameter.defaults > $stdout");
 	system("tahoe -f $stem.xml >> $stem.stdout");
	system("extract_1D $stem.io1.exo $stem.io1.exo $stem.rs0.uvst_Y >> $stdout");
	system("dpp $in $run.rs1.xml.tmpl $stem.rs1.xml parameter.defaults $tag >> $stdout");
 	system("tahoe -f $stem.rs1.xml >> $stem.stdout");
	system("extract_1D $stem.rs1.io1.exo $stem.rs1.io1.exo $stem.rs1.uvst_Y >> $stdout");
	$tag = "$tag.rs1";
	system("dpp $in $run.rs2.xml.tmpl $stem.rs2.xml parameter.defaults $tag >> $stdout");
 	system("tahoe -f $stem.rs2.xml >> $stem.stdout");
	system("extract_1D $stem.rs2.io1.exo $stem.rs2.io1.exo $stem.rs2.uvst_Y >> $stdout");
	system("cat $stem.rs0.uvst_Y  $stem.rs1.uvst_Y $stem.rs2.uvst_Y > $stem.uvst_Y");
	print GPL "\"$stem.uvst_Y\" w l,\\\n";

  # unshifted error
	system("least_square $run.dat $stem.uvst_Y $stem.error >> $stdout");
	open(ERR,"$stem.error");
	$line = <ERR>; $line =~ s/^\s+//; @items = split(/\s+/,$line); close ERR;
	$error = $items[0];
	print " error: $error";
  close(ERR);

	# shift and error
	$shift = shift_data();#"$run.dat","$stem.uvst_Y");
	system("least_square $run.dat $stem.shift $stem.shift.error >> $stdout");
	open(ERR,"$stem.shift.error");
	$line = <ERR>; $line =~ s/^\s+//; @items = split(/\s+/,$line); close ERR;
	$shift_error = $items[0];
	print " shift: $shift, error: $shift_error\n";
  close(ERR);

	#$cmd = "set logscale xy; plot \"$run.dat\" w l, \"$stem.uvst_Y\" w l, \"$stem.shift\" w l; pause 10";

	#plot
#	unless (defined ($pid = fork)) { die "cannot fork: $!"; }
#  unless ($pid) {
#		exec("echo \' $cmd \' | gnuplot - > /dev/null");
#	}

	#report error
	$errors += 0.0*$error + 100.0*$shift_error + 0.0*abs($shift);
}
print GPL "\"creep100.dat\" w p\npause -1\n";
print GPL "set term png color\nset output \"pstudy.png\"\nrepl\n";
close GPL;

#read_parameters();
#print "alpha 1 ", $varlist{"ALPHA1"}, "\n";
#print "alpha 2 ", $varlist{"ALPHA2"}, "\n";
#print "alpha 3 ", $varlist{"ALPHA3"}, "\n";
#$constraint = 0.25*(3.0*$varlist{"ALPHA"} + $varlist{"ALPHA1"});

#$error_line = "$error error";
#system("echo $error_line > $out");
#	$run = $runs[0];
#	$cmd = "plot \"$run.dat\" w l, \"$stem.uvst_Y\" w l, \"$stem.shift\" w l; pause 3";
#	exec("echo \' $cmd \' | gnuplot - > /dev/null");

#---------------------------------------------------------------------------
#sub read_parameters {
#	open(IN,"$in");
#  # read number of parameters
#  defined($line = <IN>)|| die "file error";
#  chomp($line);
#  $num_vars = $line;
#  $num_vars =~ s/[\s]*([0-9]+)[\s]+.*/$1/;
#  print STDOUT "found $num_vars variables\n";
#
#  # insert parameters into a hash
#  for ($i = 0; $i < $num_vars; $i++) {
#    defined($line = <IN>) || die "file error";
#    chomp($line);
#    $line =~ s/^\s+//; # drop leading white space
#    ($val, $var) = split(/[\s]+/, $line);
#    $varlist{$var} = $val;
#  }
#	close IN;
#}
#---------------------------------------------------------------------------
sub shift_data {
	$ref="$run.dat";
	$dat="$stem.uvst_Y";
	$new="$stem.shift";
	open(IN,"$ref"); 
	$line = <IN>;
	$line =~ s/^\s+//; 
	($tref, $valref) = split(/[\s]+/, $line);
	close IN;
	open(IN,"$dat"); 
	open(OUT,">$new"); 
	$shift = 0.0;
	while(defined($line=<IN>)) {
		$line =~ s/^\s+//; 
		($t, $val) = split(/[\s]+/, $line);
		if (abs($t - $tref) < 0.0001) { $shift = $val - $valref; }
		$val = $val - $shift;
		print OUT "$t $val\n";
	}
	close OUT;
	close IN;

	return $shift;
}

