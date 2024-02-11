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
  # run and extract data
	system("dpp $in $run.xml.tmpl $stem.xml parameter.defaults > $stdout");
 	system("tahoe -f $stem.xml >> $stem.stdout");
	system("extract_1D $stem.io1.exo $stem.io1.exo $stem.uvst_Y >> $stdout");

  # unshifted error
	system("least_square $run.dat $stem.uvst_Y $stem.error >> $stdout");
	open(ERR,"$stem.error");
	$line = <ERR>; $line =~ s/^\s+//; @items = split(/\s+/,$line); close ERR;
	$error = $items[0];
	print " error : $error ";
  close(ERR);

	# shift and error
	shift_data();#"$run.dat","$stem.uvst_Y");
	system("least_square $run.dat $stem.shift $stem.shift.error >> $stdout");
	open(ERR,"$stem.shift.error");
	$line = <ERR>; $line =~ s/^\s+//; @items = split(/\s+/,$line); close ERR;
	$shift_error = $items[0];
	print " error : $shift_error \n";
  close(ERR);

	$cmd = "set logscale x; plot \"$run.dat\" w l, \"$stem.uvst_Y\" w l, \"$stem.shift\" w l; pause 10";

	#plot
	unless (defined ($pid = fork)) { die "cannot fork: $!"; }
  unless ($pid) {
		exec("echo \' $cmd \' | gnuplot - > /dev/null");
	}

	#report error
	$errors += 0.0*$error + 10.0*$shift_error;
}

#read_parameters();
#print "alpha 1 ", $varlist{"ALPHA1"}, "\n";
#print "alpha 2 ", $varlist{"ALPHA2"}, "\n";
#print "alpha 3 ", $varlist{"ALPHA3"}, "\n";
#$constraint = 0.25*(3.0*$varlist{"ALPHA"} + $varlist{"ALPHA1"});

open(OUT,">$out");
print OUT "$errors ERROR\n";
#print OUT "$constraint CONSTRAINT\n";
print "error: $errors \n";
#print "constraint: $constraint \n";
close OUT;

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
		if (abs($t - $tref) < 0.0001) { $shift = $val - $valref; print "SHIFT: $shift\n"; }
		$val = $val - $shift;
		print OUT "$t $val\n";
	}
	close OUT;
	close IN;

}

