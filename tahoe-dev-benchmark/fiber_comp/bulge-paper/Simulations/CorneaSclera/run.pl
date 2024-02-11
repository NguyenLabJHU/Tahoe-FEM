#!/usr/bin/perl 
use File::Copy;

$run_sim = 1; 
$run_plot = 1;
if (defined($ARGV[0])) { 
  if ($ARGV[0] =~ /plot/) { $run_sim = 0;}
}

@runs = ("cornea_iso", "cornea_aniso", "sclera_iso", "sclera_aniso" );

# run simulation
if ($run_sim) {
  @files = ("globe.exo","globe.xml.tmpl","disp_curve.gpl");
  foreach $run (@runs) {
    print "starting $run ... ";
    mkdir $run;
    chdir $run;
    foreach $file (@files) {copy("../$file","$file") || die "can't copy $file";}
    $sclera = 0; if ($run =~ /sclera/) {$sclera = 1;}
    $iso    = 1; if ($run =~ /aniso/)  {$iso    = 0;}
    system("ipp iso=$iso sclera=$sclera globe.xml.tmpl -o $run.xml");
    system("tahoe -f $run.xml > tahoe_log");
    chdir "..";
    print "done\n";
  }
}


# post processing
if ($run_plot) {
  foreach $run (@runs) {
    system("cd $run; gnuplot disp_curve.gpl > gnuplot_log; mv u.png $run_u.png;  cd ..");
  }
  system("gnuplot comp_disp_curve.gpl");
  system("exodiff cornea_iso/cornea_iso.io0.exo cornea_aniso/cornea_aniso.io0.exo cornea_diff.io0.exo >& /dev/null");
  system("exodiff sclera_iso/sclera_iso.io0.exo sclera_iso/sclera_iso.io0.exo sclera_iso.io0.exo >& /dev/null");
}
