#!/usr/bin/perl -w
# author : reese jones (rjones@sandia.gov)
if (!(defined($ARGV[0]))) {
  print "usage : getSurfaceMesh <geometry_file> <sset_id> <displacement_file>\n";
  exit(1);
}
$geometryfile = $ARGV[0];
$sset_id = $ARGV[1];
$displacementfile = $ARGV[2];

$basename = (split(/\./,$displacementfile))[0];
$nset_id = $sset_id;
$logfile = $basename."_log";

$verbose = 1;

@faces = extract_mesh();
@coords = get_displacements();
@normals = calc_normals();

#print "faces ", scalar(@faces),"\n";
#print "times ", scalar(@coords),"\n";

calc_refractive_power();

#write_vtu(@faces,@coords);
write_vtu();


exit(0);
#==============================================================================
sub calc_refractive_power {
  $np = 1.3375; # composite refractive power of the cornea
  my $xtol = 0.000000001;
  my @conn = @faces;
  my @pts  = @coords;
  my @nms  = @normals;
  $ntimes = scalar(@pts);
  #print "ntimes $ntimes\n";
  my @Xs = @{$pts[0]};
  my $nnodes = scalar(@Xs);
  #print "nnodes $nnodes\n";

  # find apex
  my $apexNode = -1;
  my $found = 0;
  for (my $i = 0; $i < $nnodes; $i++) {
    if ((abs($Xs[$i][0]) < $xtol) && (abs($Xs[$i][1]) < $xtol)) {
      $apexNode = $i; $found++; 
    }
  }
#print "found $found\n";
  if (!($found)) { print ">>> can not find apex\n"; die;}
  my $apexFace = -1;
  $found = 0;
  for (my $i = 0; $i < $nfaces; $i++) {
    for (my $j = 0; $j < 4; $j++) {
      if ($faces[$i][$j] == $apexNode) { 
        $apexFace = $i; $found++;
        $jj = $j-1; if ($jj < 0) {$jj = 3;}
        $apexNeighborLeft = $faces[$i][$jj];
        $jj = $j+1; if ($jj < 0) {$jj = 0;}
        $apexNeighborRight = $faces[$i][$jj];
      }
    }
  }
  if (!($found)) { print ">>> can not find apex face\n"; die;}
#print "found $found\n";
  
  for (my $index = 0; $index < $ntimes; $index++) {
    my @xs = @{$pts[$index]};
    my @ns = @{$nms[$index]};
    my $zApex = $xs[$apexNode][2];
    @ps = ();
    for (my $i = 0; $i < $nnodes; $i++) {
      if ($i == $apexNode) {
        $rp = 0.0; # dummy for now
      }
      else {
        my $thetaI = acos($ns[$i][2]);
        my $thetaT = asin(1.0/$np*sin($thetaI));
        my $z = $xs[$i][2]-$zApex;
        my $r = sqrt( $xs[$i][0]* $xs[$i][0]+ $xs[$i][1]* $xs[$i][1] );
        $rp = $np/($z + $r/(tan($thetaI - $thetaT)));
        $rp *= 1000; # convert from 1/mm --> 1/m  i.e. diopters
      }
#print "$i $thetaI $thetaT $rp\n";
#$rp = $xs[$i][0]/tan($thetaI - $thetaT);
      push(@ps, $rp);
    }
    # handle apex node
    $ps[$apexNode] = 0.5 * ( $ps[$apexNeighborLeft] + $ps[$apexNeighborRight]);
    push (@Rps, [ @ps ] );
  }
}

#==============================================================================
sub write_vtu {
  #my @conn = shift(@_);
  #my @pts  = shift(@_);
  my @conn = @faces;
  my @pts  = @coords;
  my @nms  = @normals;
  $ntimes = scalar(@pts);
  #print "> number of times: $ntimes\n";
  my @Xs = @{$pts[0]};
  my $nnodes = scalar(@Xs);
  print "> number of nodes: $nnodes\n";
  $nfaces = scalar(@conn);
  print "> number of faces: $nfaces\n";
  for (my $index = 0; $index < $ntimes; $index++) {
    my $time = $ts[$index];
    print "> writing time $time\n";
    my @xs = @{$pts[$index]};
    my @ns = @{$nms[$index]};
    my @ps = @{$Rps[$index]};
    open(OUT,">$basename.$index.vtu");
    print OUT "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">\n";
    print OUT "<UnstructuredGrid>\n";
    print OUT "<Piece NumberOfCells=\"",$nfaces,
                  "\" NumberOfPoints=\"",$nnodes,"\">\n";
    print OUT "<Points>\n";
    print OUT "<DataArray NumberOfComponents=\"3\" format=\"ascii\" type=\"Float32\">\n";
    for (my $i = 0; $i < $nnodes; $i++) {
      print OUT $Xs[$i][0], " ", $Xs[$i][1], " ",$Xs[$i][2], "\n";
    }
    print OUT "</DataArray>\n</Points>\n";
    print OUT "<Cells>\n";
    print OUT "<DataArray Name=\"connectivity\" format=\"ascii\" type=\"Int32\">\n";
    for (my $i = 0; $i < $nfaces; $i++) {
      print OUT $conn[$i][0]," ",$conn[$i][1]," ",$conn[$i][2]," ",$conn[$i][3],"\n";
    }
    print OUT "</DataArray>\n";
    print OUT "<DataArray Name=\"offsets\" format=\"ascii\" type=\"Int32\">\n";
    for (my $i = 0; $i < $nfaces; $i++) {
      print OUT ($i+1)*4,"\n";
    }
    print OUT "</DataArray>\n";
    print OUT "<DataArray Name=\"types\" format=\"ascii\" type=\"UInt8\">\n";
    for (my $i = 0; $i < $nfaces; $i++) {
      print OUT "9\n";
    }
    print OUT "</DataArray>\n";
    print OUT "</Cells>\n";
    print OUT "<PointData>\n";
    print OUT "<DataArray Name=\"D_\"  NumberOfComponents=\"3\" format=\"ascii\" type=\"Float32\">\n";
    for (my $i = 0; $i < $nnodes; $i++) {
      print OUT $xs[$i][0]-$Xs[$i][0], 
           " ", $xs[$i][1]-$Xs[$i][1], 
           " ", $xs[$i][2]-$Xs[$i][2], "\n";
    }
    print OUT "</DataArray>\n";
    print OUT "<DataArray Name=\"N_\"  NumberOfComponents=\"3\" format=\"ascii\" type=\"Float32\">\n";
    for (my $i = 0; $i < $nnodes; $i++) {
      print OUT $ns[$i][0],
           " ", $ns[$i][1], 
           " ", $ns[$i][2], "\n";
    }
    print OUT "</DataArray>\n";
    print OUT "<DataArray Name=\"Rp\"  format=\"ascii\" type=\"Float32\">\n";
    for (my $i = 0; $i < $nnodes; $i++) {
      print OUT $ps[$i],"\n";
    }
    print OUT "</DataArray>\n";
    print OUT "</PointData>\n";
    print OUT "</Piece>\n</UnstructuredGrid>\n</VTKFile>\n";
    close OUT;
  }
}

#==============================================================================
sub calc_normals{
  my $xtol = 0.000000001;
  my @conn = @faces;
  my @pts  = @coords;
  $nfaces = scalar(@conn);
  #print "nfaces $nfaces\n";
  $ntimes = scalar(@pts);
  #print "ntimes $ntimes\n";
  my @Xs = @{$pts[0]};
  my $nnodes = scalar(@Xs);
  #print "nnodes $nnodes\n";
  for (my $index = 0; $index < $ntimes; $index++) {
    my @xs = @{$pts[$index]};
    @ns = ();
    @zeros = (0.0,0.0,0.0);
    for (my $i = 0; $i < $nnodes; $i++) {
      push (@ns, [ @zeros ] );  
    }
    for (my $i = 0; $i < $nfaces; $i++) {
      @vs = ();
      for (my $j = 0; $j < 4; $j++) {
        $n1 = $j; $n2 = $j+1; if ($n2 == 4) {$n2 = 0;}
        $node1 = $conn[$i][$n1];
        $node2 = $conn[$i][$n2];
        @v = ($xs[$node2][0]-$xs[$node1][0],
              $xs[$node2][1]-$xs[$node1][1],
              $xs[$node2][2]-$xs[$node1][2]);
        push (@vs,[ @v ]);
      }
      for (my $j = 0; $j < 4; $j++) {
        $n1 = $j; $n2 = $j+1; if ($n2 == 4) {$n2 = 0;}
        @nm = ($vs[$n1][2]*$vs[$n2][1] - $vs[$n1][1]*$vs[$n2][2],
               $vs[$n1][0]*$vs[$n2][2] - $vs[$n1][2]*$vs[$n2][0],
               $vs[$n1][1]*$vs[$n2][0] - $vs[$n1][0]*$vs[$n2][1]);
        $node = $conn[$i][$n2];
        # correct for symmetry
        my $x = $xs[$node][0];
        my $y = $xs[$node][1];
        if (abs($x) < $xtol) { $nm[0] = 0.0;}
        if (abs($y) < $xtol) { $nm[1] = 0.0;}
        # add normalized normals (as opposed to jac weighted)
        $mag = sqrt($nm[0]*$nm[0] +$nm[1]*$nm[1]+$nm[2]*$nm[2]);
        $ns[$node][0] += $nm[0]/$mag;
        $ns[$node][1] += $nm[1]/$mag;
        $ns[$node][2] += $nm[2]/$mag;
      }
    }
    for (my $i = 0; $i < $nnodes; $i++) {
      $mag = sqrt($ns[$i][0]*$ns[$i][0]
                 +$ns[$i][1]*$ns[$i][1]
                 +$ns[$i][2]*$ns[$i][2]);
      $ns[$i][0] = -$ns[$i][0]/$mag;
      $ns[$i][1] = -$ns[$i][1]/$mag;
      $ns[$i][2] = -$ns[$i][2]/$mag;
#print $ns[$i][0]," ",$ns[$i][1]," ",$ns[$i][2], " $mag\n";
    }
    push (@normals, [ @ns ]);
  }
  return @normals;
}

#==============================================================================
sub get_displacements{
my $nnodes = 0;
my $offset = 0; 
my $tindex = 0;
system("exotxt $displacementfile $displacementfile.txt > /dev/null");
open(IN,"$displacementfile.txt") 
  || die "cannot open displacement file $displacementfile for reading : $!";

while (defined($line = <IN>) ) {
  if ($line =~ /initial/)   {
    $_ = <IN>;
    $_ = <IN> ; chomp ; s/^\s+//; s/\s+$//; @list = split(/ +/);
    $nnodes = $list[0];
  }
  if ($line =~ /Coordinates/)   {
    for (my $i = 0; $i < $nnodes; $i++) {
      $_ = <IN> ; chomp ; s/^\s+//; s/\s+$//; @list = split(/ +/);
      $X = $list[0]; $Y = $list[1]; $Z = $list[2];
      push (@Xs,$X); push (@Ys,$Y); push (@Zs,$Z);
    }
  }
  if ($line =~ /Time/)   {
    $_ = <IN> ; chomp ; s/^\s+//; s/\s+$//; @list = split(/ +/);
    $time = $list[0];
    push (@ts,$time);
  }
  if ($line =~ /Nodal/) {
    @xyzs = ();
    for (my $i = 0; $i < $nnodes; $i++) {
      $_ = <IN> ; chomp ; s/^\s+//; s/\s+$//; @list = split(/ +/);
      $_ = <IN>; # junk
      $u = $list[0]; $v = $list[1]; $w = $list[2];
      push (@us,$u); push (@vs,$v); push (@ws,$w);
      $x = $Xs[$i] + $u; $y= $Ys[$i] + $v; $z= $Zs[$i] + $w;
      @xyz = ($x, $y, $z);
      push @xyzs, [ @xyz ];
    }
    # output to file
    my $time = $ts[$tindex];
    #print "> writing time $time\n";
    open(OUT,">$basename.points.$tindex.dat");
    print OUT "#time : $time\n";
    for (my $j = 0; $j < $nnodes; $j++) {
      $k = $i +$offset;
      print OUT $Xs[$j]," ",$Ys[$j]," ",$Zs[$j],
            " ",$us[$k]," ",$vs[$k]," ",$ws[$k], "\n";
    }
    close OUT; $tindex++;

    $offset += $nnodes;
    push(@coords,[ @xyzs ]);
  }

}
close(IN);

return @coords;
}

#==============================================================================
sub extract_mesh{
	my ($def,$not_found,$ns,$nn);

  $nsetfile = "$basename.ns$nset_id";
  $ssetfile = "$basename.ss$sset_id";
  $meshfile = $geometryfile;
	# generate mesh file if needed
	if (!(-e $meshfile)) {
    $journalfile = ((split(/\./,$meshfile))[0]).".jou";
    if (-e $journalfile) {
		  if ($verbose) {print "> creating mesh file...";}
		  system("cubit -nogui -nographics -nojournal $basename.jou > $logfile");
		  if ($verbose) {print "done\n";}
    }
    else { print "ERROR : no cubit journal file\n"; die; }
	} 

	$txtfile = "$meshfile.txt";
	system("exotxt $meshfile $txtfile > /dev/null");
	open(IN,"$txtfile");

	# read number of nodes
	while (($def = defined($line =<IN>))
		&& !($line =~ /Database initial variables/)) { };
	<IN>;
	@items = get_items($line=<IN>);
	#$nnodes = $items[0];
	#if ($verbose) {print "total number of nodes: $nnodes\n";}

	# read element connectivities
	while (($def = defined($line =<IN>))
		&& !($line =~ /Element block/)) { };
	$ie = 0;
	$continue = 1;
	while ($continue) {
		@items = get_items($line=<IN>);
		$nelems = $items[1];
		<IN>; <IN>;
		for ($i=0; $i < $nelems; $i++) {
			@items = get_items($line=<IN>);
			$ie++;
			$conn[$ie] = [$items[0],$items[1],$items[2],$items[3],
                    $items[4],$items[5],$items[6],$items[7]];
		}
		if (!(defined($line = <IN>)) 
			|| !($line =~ /Element block/)) { $continue = 0;}
	}
	#if ($verbose) {print "read $ie elements\n";}

	# find nset
	$not_found = 1;
	while($not_found) {
		while (($def = defined($line =<IN>)) && !($line =~ /Nodal point set/)) { };
    @items = get_items($line=<IN>);
    if ($items[0] =~ /$nset_id/) { $not_found = 0;}
    $nn = $items[1];
	}
  if ($not_found) { print ">>> nodeset $nset_id not found\n"; die; }

	# write nset (global numbering)
	if ($def) {
		for (my $i=0; $i < $nn; $i++) { 
			@items = get_items($line = <IN>);
 			$nset[$i] = $items[0]; 
		}
		#if ($verbose) {print "> writing $nn nodes...";}
		open(OUT,">$nsetfile");
  	@nset_sorted = sort { $a <=> $b } @nset;
		$i = 1;
		foreach $node (@nset_sorted) {
			print OUT "$node\n";
			$rev_map{$node} = $i++;
		}
 		close OUT;
		#if ($verbose) {print "done\n";}
	} else { print "ERROR : no sideset $sset_id found\n"; die; }

	# find sset
	while (($def = defined($line =<IN>)) 
		&& !($line =~ /Number of Side Sets/)) { };
	$not_found = 1;
	while ($not_found) {
		while (($def = defined($line =<IN>)) && !($line =~ /\#sides/)) { };
		@items = get_items($line);
		if ($items[0] =~ /$sset_id/) { 
			$not_found = 0;
			$ns = $items[1];
		}
	}
  if ($not_found) { print ">>> sideset $sset_id not found\n"; die; }

	# local face ordering
  @side2nodes =
		([0,1,5,4],
		 [1,2,6,5],
		 [2,3,7,6],
		 [3,0,4,7],
		 [3,2,1,0],
		 [4,5,6,7]);

	# write sset (local numbering)
	if ($def) {
		#if ($verbose) {print "> writing $ns facets...";}
 		open(OUT,">$ssetfile");
		<IN>;
		$count = 0;
		$i=0;
		@items = get_items($line = <IN>);
		while ($count < $ns) {
			if (($count != 0) && (($count % 6) == 0)) {
				@items = get_items($line = <IN>);
				$i = 0;
			}
			$elem = $items[2*$i];
			$side = $items[2*$i+1]-1;
			for (my $j=0; $j < 4; $j++) { 
				$face[$j] = $conn[$elem][$side2nodes[$side][$j]];
			}
#			print OUT "$count :: element: $elem, side: $side\n";
			for (my $j=0; $j < 4; $j++) { 
				$face[$j] = $rev_map{$face[$j]} -1;
			}
 			print OUT $face[0]," ",
 			          $face[1]," ",
 			          $face[2]," ",
 			          $face[3],"\n";
      push (@faces, [ @face ]);
			$i++;
			$count++;
		}
		close OUT;
		#if ($verbose) {print "done\n";}
	} else { print "ERROR : no sideset $sset_id found\n"; die; }

  return @faces;
}
#==============================================================================
sub get_items {
  my ($line);
  $line = $_[0];
  chomp($line);
  $line =~ s/^\s+//;
  $line =~ s/\s+$//;
  @items = split(/\s+/,$line);
  return @items;
}
#==============================================================================
sub asin {
  if (abs($_[0]) > 1) { die "input out of range in asin()\n"; }
  atan2($_[0],sqrt(1-$_[0]**2));
}
sub acos {
  if (abs($_[0]) > 1) { die "input out of range in acos()\n"; }
  atan2(sqrt(1-$_[0]**2),$_[0]);
}
sub tan { sin($_[0])/cos($_[0]) }
