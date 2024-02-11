clc

clear
close all
Node401_cls=load('Beam_classic_401.txt');
Node402_cls=load('Beam_classic_402.txt');
Node405_cls=load('Beam_classic_405.txt');
Node408_cls=load('Beam_classic_408.txt');
Node411_cls=load('Beam_classic_411.txt');
Node414_cls=load('Beam_classic_414.txt');
Node417_cls=load('Beam_classic_417.txt');
Node420_cls=load('Beam_classic_420.txt');
Node423_cls=load('Beam_classic_423.txt');
Node426_cls=load('Beam_classic_426.txt');
Node429_cls=load('Beam_classic_429.txt');
Node432_cls=load('Beam_classic_432.txt');
Node435_cls=load('Beam_classic_435.txt');
Node438_cls=load('Beam_classic_438.txt');
Node441_cls=load('Beam_classic_441.txt');
Node444_cls=load('Beam_classic_444.txt');
Node447_cls=load('Beam_classic_447.txt');
Node450_cls=load('Beam_classic_450.txt');
Node453_cls=load('Beam_classic_453.txt');
Node456_cls=load('Beam_classic_456.txt');
Node459_cls=load('Beam_classic_459.txt');
Node462_cls=load('Beam_classic_462.txt');
Node465_cls=load('Beam_classic_465.txt');
Node468_cls=load('Beam_classic_468.txt');
Node471_cls=load('Beam_classic_471.txt');
Node474_cls=load('Beam_classic_474.txt');
Node477_cls=load('Beam_classic_477.txt');
Node480_cls=load('Beam_classic_480.txt');
Node483_cls=load('Beam_classic_483.txt');
Node486_cls=load('Beam_classic_486.txt');
Node489_cls=load('Beam_classic_489.txt');
Node492_cls=load('Beam_classic_492.txt');
Node495_cls=load('Beam_classic_495.txt');


NodeQ27_6_com=load('Beam_Q27P8_comb_0006.txt');
NodeQ27_7_com=load('Beam_Q27P8_comb_0007.txt');
NodeQ27_18_com=load('Beam_Q27P8_comb_0018.txt');
NodeQ27_2312_com=load('Beam_Q27P8_comb_2312.txt');
NodeQ27_2313_com=load('Beam_Q27P8_comb_2313.txt');
NodeQ27_2336_com=load('Beam_Q27P8_comb_2336.txt');
NodeQ27_2338_com=load('Beam_Q27P8_comb_2338.txt');
NodeQ27_2356_com=load('Beam_Q27P8_comb_2356.txt');
NodeQ27_2358_com=load('Beam_Q27P8_comb_2358.txt');
NodeQ27_2376_com=load('Beam_Q27P8_comb_2376.txt');
NodeQ27_2378_com=load('Beam_Q27P8_comb_2378.txt');
NodeQ27_2396_com=load('Beam_Q27P8_comb_2396.txt');
NodeQ27_2389_com=load('Beam_Q27P8_comb_2398.txt');
NodeQ27_2416_com=load('Beam_Q27P8_comb_2416.txt');
NodeQ27_2418_com=load('Beam_Q27P8_comb_2418.txt');
NodeQ27_2436_com=load('Beam_Q27P8_comb_2436.txt');
NodeQ27_2438_com=load('Beam_Q27P8_comb_2438.txt');
NodeQ27_2456_com=load('Beam_Q27P8_comb_2456.txt');
NodeQ27_2458_com=load('Beam_Q27P8_comb_2458.txt');
NodeQ27_2476_com=load('Beam_Q27P8_comb_2476.txt');
NodeQ27_2478_com=load('Beam_Q27P8_comb_2478.txt');
NodeQ27_2496_com=load('Beam_Q27P8_comb_2496.txt');
NodeQ27_2498_com=load('Beam_Q27P8_comb_2498.txt');
NodeQ27_2516_com=load('Beam_Q27P8_comb_2516.txt');
NodeQ27_2518_com=load('Beam_Q27P8_comb_2518.txt');
NodeQ27_2536_com=load('Beam_Q27P8_comb_2536.txt');
NodeQ27_2538_com=load('Beam_Q27P8_comb_2538.txt');
NodeQ27_2556_com=load('Beam_Q27P8_comb_2556.txt');
NodeQ27_2558_com=load('Beam_Q27P8_comb_2558.txt');
NodeQ27_2576_com=load('Beam_Q27P8_comb_2576.txt');
NodeQ27_2578_com=load('Beam_Q27P8_comb_2578.txt');
NodeQ27_2596_com=load('Beam_Q27P8_comb_2596.txt');
NodeQ27_2598_com=load('Beam_Q27P8_comb_2598.txt');
NodeQ27_2616_com=load('Beam_Q27P8_comb_2616.txt');
NodeQ27_2618_com=load('Beam_Q27P8_comb_2618.txt');
NodeQ27_2636_com=load('Beam_Q27P8_comb_2636.txt');
NodeQ27_2638_com=load('Beam_Q27P8_comb_2638.txt');
NodeQ27_2656_com=load('Beam_Q27P8_comb_2656.txt');
NodeQ27_2658_com=load('Beam_Q27P8_comb_2658.txt');
NodeQ27_2676_com=load('Beam_Q27P8_comb_2676.txt');
NodeQ27_2678_com=load('Beam_Q27P8_comb_2678.txt');
NodeQ27_2696_com=load('Beam_Q27P8_comb_2696.txt');
NodeQ27_2698_com=load('Beam_Q27P8_comb_2698.txt');
NodeQ27_2716_com=load('Beam_Q27P8_comb_2716.txt');
NodeQ27_2718_com=load('Beam_Q27P8_comb_2718.txt');
NodeQ27_2736_com=load('Beam_Q27P8_comb_2736.txt');
NodeQ27_2738_com=load('Beam_Q27P8_comb_2738.txt');
NodeQ27_2756_com=load('Beam_Q27P8_comb_2756.txt');
NodeQ27_2758_com=load('Beam_Q27P8_comb_2758.txt');
NodeQ27_2776_com=load('Beam_Q27P8_comb_2776.txt');
NodeQ27_2778_com=load('Beam_Q27P8_comb_2778.txt');
NodeQ27_2796_com=load('Beam_Q27P8_comb_2796.txt');
NodeQ27_2798_com=load('Beam_Q27P8_comb_2798.txt');
NodeQ27_2816_com=load('Beam_Q27P8_comb_2816.txt');
NodeQ27_2818_com=load('Beam_Q27P8_comb_2818.txt');
NodeQ27_2836_com=load('Beam_Q27P8_comb_2836.txt');
NodeQ27_2838_com=load('Beam_Q27P8_comb_2838.txt');
NodeQ27_2856_com=load('Beam_Q27P8_comb_2856.txt');
NodeQ27_2858_com=load('Beam_Q27P8_comb_2858.txt');
NodeQ27_2876_com=load('Beam_Q27P8_comb_2876.txt');
NodeQ27_2878_com=load('Beam_Q27P8_comb_2878.txt');
NodeQ27_2896_com=load('Beam_Q27P8_comb_2896.txt');
NodeQ27_2898_com=load('Beam_Q27P8_comb_2898.txt');
NodeQ27_2914_com=load('Beam_Q27P8_comb_2914.txt');



Node401_com=load('Beam_comb_401.txt');
Node402_com=load('Beam_comb_402.txt');
Node405_com=load('Beam_comb_405.txt');
Node408_com=load('Beam_comb_408.txt');
Node411_com=load('Beam_comb_411.txt');
Node414_com=load('Beam_comb_414.txt');
Node417_com=load('Beam_comb_417.txt');
Node420_com=load('Beam_comb_420.txt');
Node423_com=load('Beam_comb_423.txt');
Node426_com=load('Beam_comb_426.txt');
Node429_com=load('Beam_comb_429.txt');
Node432_com=load('Beam_comb_432.txt');
Node435_com=load('Beam_comb_435.txt');
Node438_com=load('Beam_comb_438.txt');
Node441_com=load('Beam_comb_441.txt');
Node444_com=load('Beam_comb_444.txt');
Node447_com=load('Beam_comb_447.txt');
Node450_com=load('Beam_comb_450.txt');
Node453_com=load('Beam_comb_453.txt');
Node456_com=load('Beam_comb_456.txt');
Node459_com=load('Beam_comb_459.txt');
Node462_com=load('Beam_comb_462.txt');
Node465_com=load('Beam_comb_465.txt');
Node468_com=load('Beam_comb_468.txt');
Node471_com=load('Beam_comb_471.txt');
Node474_com=load('Beam_comb_474.txt');
Node477_com=load('Beam_comb_477.txt');
Node480_com=load('Beam_comb_480.txt');
Node483_com=load('Beam_comb_483.txt');
Node486_com=load('Beam_comb_486.txt');
Node489_com=load('Beam_comb_489.txt');
Node492_com=load('Beam_comb_492.txt');
Node495_com=load('Beam_comb_495.txt');

Node401_rot=load('Beam_rotation_401.txt');
Node402_rot=load('Beam_rotation_402.txt');
Node405_rot=load('Beam_rotation_405.txt');
Node408_rot=load('Beam_rotation_408.txt');
Node411_rot=load('Beam_rotation_411.txt');
Node414_rot=load('Beam_rotation_414.txt');
Node417_rot=load('Beam_rotation_417.txt');
Node420_rot=load('Beam_rotation_420.txt');
Node423_rot=load('Beam_rotation_423.txt');
Node426_rot=load('Beam_rotation_426.txt');
Node429_rot=load('Beam_rotation_429.txt');
Node432_rot=load('Beam_rotation_432.txt');
Node435_rot=load('Beam_rotation_435.txt');
Node438_rot=load('Beam_rotation_438.txt');
Node441_rot=load('Beam_rotation_441.txt');
Node444_rot=load('Beam_rotation_444.txt');
Node447_rot=load('Beam_rotation_447.txt');
Node450_rot=load('Beam_rotation_450.txt');
Node453_rot=load('Beam_rotation_453.txt');
Node456_rot=load('Beam_rotation_456.txt');
Node459_rot=load('Beam_rotation_459.txt');
Node462_rot=load('Beam_rotation_462.txt');
Node465_rot=load('Beam_rotation_465.txt');
Node468_rot=load('Beam_rotation_468.txt');
Node471_rot=load('Beam_rotation_471.txt');
Node474_rot=load('Beam_rotation_474.txt');
Node477_rot=load('Beam_rotation_477.txt');
Node480_rot=load('Beam_rotation_480.txt');
Node483_rot=load('Beam_rotation_483.txt');
Node486_rot=load('Beam_rotation_486.txt');
Node489_rot=load('Beam_rotation_489.txt');
Node492_rot=load('Beam_rotation_492.txt');
Node495_rot=load('Beam_rotation_495.txt');


Node401_str=load('Beam_stretch_401.txt');
Node402_str=load('Beam_stretch_402.txt');
Node405_str=load('Beam_stretch_405.txt');
Node408_str=load('Beam_stretch_408.txt');
Node411_str=load('Beam_stretch_411.txt');
Node414_str=load('Beam_stretch_414.txt');
Node417_str=load('Beam_stretch_417.txt');
Node420_str=load('Beam_stretch_420.txt');
Node423_str=load('Beam_stretch_423.txt');
Node426_str=load('Beam_stretch_426.txt');
Node429_str=load('Beam_stretch_429.txt');
Node432_str=load('Beam_stretch_432.txt');
Node435_str=load('Beam_stretch_435.txt');
Node438_str=load('Beam_stretch_438.txt');
Node441_str=load('Beam_stretch_441.txt');
Node444_str=load('Beam_stretch_444.txt');
Node447_str=load('Beam_stretch_447.txt');
Node450_str=load('Beam_stretch_450.txt');
Node453_str=load('Beam_stretch_453.txt');
Node456_str=load('Beam_stretch_456.txt');
Node459_str=load('Beam_stretch_459.txt');
Node462_str=load('Beam_stretch_462.txt');
Node465_str=load('Beam_stretch_465.txt');
Node468_str=load('Beam_stretch_468.txt');
Node471_str=load('Beam_stretch_471.txt');
Node474_str=load('Beam_stretch_474.txt');
Node477_str=load('Beam_stretch_477.txt');
Node480_str=load('Beam_stretch_480.txt');
Node483_str=load('Beam_stretch_483.txt');
Node486_str=load('Beam_stretch_486.txt');
Node489_str=load('Beam_stretch_489.txt');
Node492_str=load('Beam_stretch_492.txt');
Node495_str=load('Beam_stretch_495.txt');


MicroPolar = load('Beam.txt');
Lateral_deflection_classic = load('Beam_timoshenko_classic.txt'); 
Rotation_Beam = load('Rotation_Beam.txt'); 
Mxy = load('Mxy.txt'); 


x = 0:1/32:1;

E = 20e9;
b = 0.05;
h = 0.2;

I = (1/12)*b*h^3;

Coeff = 3*E*I/1000;

Result_phi11_com = [x(1,1) Node401_com(1,2);x(1,2) Node402_com(1,2);x(1,3) Node405_com(1,2);x(1,4) Node408_com(1,2);x(1,5) Node411_com(1,2);x(1,6) Node414_com(1,2);
    x(1,7) Node417_com(1,2);x(1,8) Node420_com(1,2);x(1,9) Node423_com(1,2);x(1,10) Node426_com(1,2);x(1,11) Node429_com(1,2);x(1,12) Node432_com(1,2);
    x(1,13) Node435_com(1,2);x(1,14) Node438_com(1,2);x(1,15) Node441_com(1,2);x(1,16) Node444_com(1,2);x(1,17) Node447_com(1,2);x(1,18) Node450_com(1,2);
    x(1,19) Node453_com(1,2);x(1,20) Node456_com(1,2);x(1,21) Node459_com(1,2);x(1,22) Node462_com(1,2);x(1,23) Node465_com(1,2);x(1,24) Node468_com(1,2);
    x(1,25) Node471_com(1,2);x(1,26) Node474_com(1,2);x(1,27) Node477_com(1,2);x(1,28) Node480_com(1,2);x(1,29) Node483_com(1,2);x(1,30) Node486_com(1,2);
    x(1,31) Node489_com(1,2);x(1,32) Node492_com(1,2);x(1,33) Node495_com(1,2)];


Result_phi11_rot = [x(1,1) Node401_rot(1,2);x(1,2) Node402_rot(1,2);x(1,3) Node405_rot(1,2);x(1,4) Node408_rot(1,2);x(1,5) Node411_rot(1,2);x(1,6) Node414_rot(1,2);
    x(1,7) Node417_rot(1,2);x(1,8) Node420_rot(1,2);x(1,9) Node423_rot(1,2);x(1,10) Node426_rot(1,2);x(1,11) Node429_rot(1,2);x(1,12) Node432_rot(1,2);
    x(1,13) Node435_rot(1,2);x(1,14) Node438_rot(1,2);x(1,15) Node441_rot(1,2);x(1,16) Node444_rot(1,2);x(1,17) Node447_rot(1,2);x(1,18) Node450_rot(1,2);
    x(1,19) Node453_rot(1,2);x(1,20) Node456_rot(1,2);x(1,21) Node459_rot(1,2);x(1,22) Node462_rot(1,2);x(1,23) Node465_rot(1,2);x(1,24) Node468_rot(1,2);
    x(1,25) Node471_rot(1,2);x(1,26) Node474_rot(1,2);x(1,27) Node477_rot(1,2);x(1,28) Node480_rot(1,2);x(1,29) Node483_rot(1,2);x(1,30) Node486_rot(1,2);
    x(1,31) Node489_rot(1,2);x(1,32) Node492_rot(1,2);x(1,33) Node495_rot(1,2)];


Result_phi11_str = [x(1,1) Node401_str(1,2);x(1,2) Node402_str(1,2);x(1,3) Node405_str(1,2);x(1,4) Node408_str(1,2);x(1,5) Node411_str(1,2);x(1,6) Node414_str(1,2);
    x(1,7) Node417_str(1,2);x(1,8) Node420_str(1,2);x(1,9) Node423_str(1,2);x(1,10) Node426_str(1,2);x(1,11) Node429_str(1,2);x(1,12) Node432_str(1,2);
    x(1,13) Node435_str(1,2);x(1,14) Node438_str(1,2);x(1,15) Node441_str(1,2);x(1,16) Node444_str(1,2);x(1,17) Node447_str(1,2);x(1,18) Node450_str(1,2);
    x(1,19) Node453_str(1,2);x(1,20) Node456_str(1,2);x(1,21) Node459_str(1,2);x(1,22) Node462_str(1,2);x(1,23) Node465_str(1,2);x(1,24) Node468_str(1,2);
    x(1,25) Node471_str(1,2);x(1,26) Node474_str(1,2);x(1,27) Node477_str(1,2);x(1,28) Node480_str(1,2);x(1,29) Node483_str(1,2);x(1,30) Node486_str(1,2);
    x(1,31) Node489_str(1,2);x(1,32) Node492_str(1,2);x(1,33) Node495_str(1,2)];


Result_phi31_com = [x(1,1) Node401_com(1,3);x(1,2) Node402_com(1,3);x(1,3) Node405_com(1,3);x(1,4) Node408_com(1,3);x(1,5) Node411_com(1,3);x(1,6) Node414_com(1,3);
    x(1,7) Node417_com(1,3);x(1,8) Node420_com(1,3);x(1,9) Node423_com(1,3);x(1,10) Node426_com(1,3);x(1,11) Node429_com(1,3);x(1,12) Node432_com(1,3);
    x(1,13) Node435_com(1,3);x(1,14) Node438_com(1,3);x(1,15) Node441_com(1,3);x(1,16) Node444_com(1,3);x(1,17) Node447_com(1,3);x(1,18) Node450_com(1,3);
    x(1,19) Node453_com(1,3);x(1,20) Node456_com(1,3);x(1,21) Node459_com(1,3);x(1,22) Node462_com(1,3);x(1,23) Node465_com(1,3);x(1,24) Node468_com(1,3);
    x(1,25) Node471_com(1,3);x(1,26) Node474_com(1,3);x(1,27) Node477_com(1,3);x(1,28) Node480_com(1,3);x(1,29) Node483_com(1,3);x(1,30) Node486_com(1,3);
    x(1,31) Node489_com(1,3);x(1,32) Node492_com(1,3);x(1,33) Node495_com(1,3)];


Result_phi31_rot = [x(1,1) Node401_rot(1,3);x(1,2) Node402_rot(1,3);x(1,3) Node405_rot(1,3);x(1,4) Node408_rot(1,3);x(1,5) Node411_rot(1,3);x(1,6) Node414_rot(1,3);
    x(1,7) Node417_rot(1,3);x(1,8) Node420_rot(1,3);x(1,9) Node423_rot(1,3);x(1,10) Node426_rot(1,3);x(1,11) Node429_rot(1,3);x(1,12) Node432_rot(1,3);
    x(1,13) Node435_rot(1,3);x(1,14) Node438_rot(1,3);x(1,15) Node441_rot(1,3);x(1,16) Node444_rot(1,3);x(1,17) Node447_rot(1,3);x(1,18) Node450_rot(1,3);
    x(1,19) Node453_rot(1,3);x(1,20) Node456_rot(1,3);x(1,21) Node459_rot(1,3);x(1,22) Node462_rot(1,3);x(1,23) Node465_rot(1,3);x(1,24) Node468_rot(1,3);
    x(1,25) Node471_rot(1,3);x(1,26) Node474_rot(1,3);x(1,27) Node477_rot(1,3);x(1,28) Node480_rot(1,3);x(1,29) Node483_rot(1,3);x(1,30) Node486_rot(1,3);
    x(1,31) Node489_rot(1,3);x(1,32) Node492_rot(1,3);x(1,33) Node495_rot(1,3)];


Result_phi31_str = [x(1,1) Node401_str(1,3);x(1,2) Node402_str(1,3);x(1,3) Node405_str(1,3);x(1,4) Node408_str(1,3);x(1,5) Node411_str(1,3);x(1,6) Node414_str(1,3);
    x(1,7) Node417_str(1,3);x(1,8) Node420_str(1,3);x(1,9) Node423_str(1,3);x(1,10) Node426_str(1,3);x(1,11) Node429_str(1,3);x(1,12) Node432_str(1,3);
    x(1,13) Node435_str(1,3);x(1,14) Node438_str(1,3);x(1,15) Node441_str(1,3);x(1,16) Node444_str(1,3);x(1,17) Node447_str(1,3);x(1,18) Node450_str(1,3);
    x(1,19) Node453_str(1,3);x(1,20) Node456_str(1,3);x(1,21) Node459_str(1,3);x(1,22) Node462_str(1,3);x(1,23) Node465_str(1,3);x(1,24) Node468_str(1,3);
    x(1,25) Node471_str(1,3);x(1,26) Node474_str(1,3);x(1,27) Node477_str(1,3);x(1,28) Node480_str(1,3);x(1,29) Node483_str(1,3);x(1,30) Node486_str(1,3);
    x(1,31) Node489_str(1,3);x(1,32) Node492_str(1,3);x(1,33) Node495_str(1,3)];


Result_phi31_Q27_com = [x(1,1) NodeQ27_2312_com(1,3);x(1,2) NodeQ27_2313_com(1,3);x(1,3) NodeQ27_2336_com(1,3);x(1,4) NodeQ27_2356_com(1,3);...
    x(1,5) NodeQ27_2376_com(1,3);x(1,6) NodeQ27_2396_com(1,3);x(1,7) NodeQ27_2416_com(1,3);x(1,8) NodeQ27_2436_com(1,3);...
    x(1,9) NodeQ27_2456_com(1,3);x(1,10) NodeQ27_2476_com(1,3);x(1,11) NodeQ27_2496_com(1,3);x(1,12) NodeQ27_2516_com(1,3);...
    x(1,13) NodeQ27_2536_com(1,3);x(1,14) NodeQ27_2556_com(1,3);x(1,15) NodeQ27_2576_com(1,3);x(1,16) NodeQ27_2596_com(1,3);...
    x(1,17) NodeQ27_2616_com(1,3);x(1,18) NodeQ27_2636_com(1,3);x(1,19) NodeQ27_2656_com(1,3);x(1,20) NodeQ27_2676_com(1,3);...
    x(1,21) NodeQ27_2696_com(1,3);x(1,22) NodeQ27_2716_com(1,3);x(1,23) NodeQ27_2736_com(1,3);x(1,24) NodeQ27_2756_com(1,3);...
    x(1,25) NodeQ27_2776_com(1,3);x(1,26) NodeQ27_2796_com(1,3);x(1,27) NodeQ27_2816_com(1,3);x(1,28) NodeQ27_2836_com(1,3);...
    x(1,29) NodeQ27_2856_com(1,3);x(1,30) NodeQ27_2876_com(1,3);x(1,31) NodeQ27_2896_com(1,3);x(1,32) NodeQ27_6_com(1,3);...
    x(1,33) NodeQ27_7_com(1,3)];


Result_phi22_com = [x(1,1) Node401_com(1,4);x(1,2) Node402_com(1,4);x(1,3) Node405_com(1,4);x(1,4) Node408_com(1,4);x(1,5) Node411_com(1,4);x(1,6) Node414_com(1,4);
    x(1,7) Node417_com(1,4);x(1,8) Node420_com(1,4);x(1,9) Node423_com(1,4);x(1,10) Node426_com(1,4);x(1,11) Node429_com(1,4);x(1,12) Node432_com(1,4);
    x(1,13) Node435_com(1,4);x(1,14) Node438_com(1,4);x(1,15) Node441_com(1,4);x(1,16) Node444_com(1,4);x(1,17) Node447_com(1,4);x(1,18) Node450_com(1,4);
    x(1,19) Node453_com(1,4);x(1,20) Node456_com(1,4);x(1,21) Node459_com(1,4);x(1,22) Node462_com(1,4);x(1,23) Node465_com(1,4);x(1,24) Node468_com(1,4);
    x(1,25) Node471_com(1,4);x(1,26) Node474_com(1,4);x(1,27) Node477_com(1,4);x(1,28) Node480_com(1,4);x(1,29) Node483_com(1,4);x(1,30) Node486_com(1,4);
    x(1,31) Node489_com(1,4);x(1,32) Node492_com(1,4);x(1,33) Node495_com(1,4)];


Result_phi22_rot = [x(1,1) Node401_rot(1,4);x(1,2) Node402_rot(1,4);x(1,3) Node405_rot(1,4);x(1,4) Node408_rot(1,4);x(1,5) Node411_rot(1,4);x(1,6) Node414_rot(1,4);
    x(1,7) Node417_rot(1,4);x(1,8) Node420_rot(1,4);x(1,9) Node423_rot(1,4);x(1,10) Node426_rot(1,4);x(1,11) Node429_rot(1,4);x(1,12) Node432_rot(1,4);
    x(1,13) Node435_rot(1,4);x(1,14) Node438_rot(1,4);x(1,15) Node441_rot(1,4);x(1,16) Node444_rot(1,4);x(1,17) Node447_rot(1,4);x(1,18) Node450_rot(1,4);
    x(1,19) Node453_rot(1,4);x(1,20) Node456_rot(1,4);x(1,21) Node459_rot(1,4);x(1,22) Node462_rot(1,4);x(1,23) Node465_rot(1,4);x(1,24) Node468_rot(1,4);
    x(1,25) Node471_rot(1,4);x(1,26) Node474_rot(1,4);x(1,27) Node477_rot(1,4);x(1,28) Node480_rot(1,4);x(1,29) Node483_rot(1,4);x(1,30) Node486_rot(1,4);
    x(1,31) Node489_rot(1,4);x(1,32) Node492_rot(1,4);x(1,33) Node495_rot(1,4)];


Result_phi22_str = [x(1,1) Node401_str(1,4);x(1,2) Node402_str(1,4);x(1,3) Node405_str(1,4);x(1,4) Node408_str(1,4);x(1,5) Node411_str(1,4);x(1,6) Node414_str(1,4);
    x(1,7) Node417_str(1,4);x(1,8) Node420_str(1,4);x(1,9) Node423_str(1,4);x(1,10) Node426_str(1,4);x(1,11) Node429_str(1,4);x(1,12) Node432_str(1,4);
    x(1,13) Node435_str(1,4);x(1,14) Node438_str(1,4);x(1,15) Node441_str(1,4);x(1,16) Node444_str(1,4);x(1,17) Node447_str(1,4);x(1,18) Node450_str(1,4);
    x(1,19) Node453_str(1,4);x(1,20) Node456_str(1,4);x(1,21) Node459_str(1,4);x(1,22) Node462_str(1,4);x(1,23) Node465_str(1,4);x(1,24) Node468_str(1,4);
    x(1,25) Node471_str(1,4);x(1,26) Node474_str(1,4);x(1,27) Node477_str(1,4);x(1,28) Node480_str(1,4);x(1,29) Node483_str(1,4);x(1,30) Node486_str(1,4);
    x(1,31) Node489_str(1,4);x(1,32) Node492_str(1,4);x(1,33) Node495_str(1,4)];


Result_phi13_Q27_com = [x(1,1) NodeQ27_2312_com(1,5);x(1,2) NodeQ27_2313_com(1,5);x(1,3) NodeQ27_2336_com(1,5);x(1,4) NodeQ27_2356_com(1,5);...
    x(1,5) NodeQ27_2376_com(1,5);x(1,6) NodeQ27_2396_com(1,5);x(1,7) NodeQ27_2416_com(1,5);x(1,8) NodeQ27_2436_com(1,5);...
    x(1,9) NodeQ27_2456_com(1,5);x(1,10) NodeQ27_2476_com(1,5);x(1,11) NodeQ27_2496_com(1,5);x(1,12) NodeQ27_2516_com(1,5);...
    x(1,13) NodeQ27_2536_com(1,5);x(1,14) NodeQ27_2556_com(1,5);x(1,15) NodeQ27_2576_com(1,5);x(1,16) NodeQ27_2596_com(1,5);...
    x(1,17) NodeQ27_2616_com(1,5);x(1,18) NodeQ27_2636_com(1,5);x(1,19) NodeQ27_2656_com(1,5);x(1,20) NodeQ27_2676_com(1,5);...
    x(1,21) NodeQ27_2696_com(1,5);x(1,22) NodeQ27_2716_com(1,5);x(1,23) NodeQ27_2736_com(1,5);x(1,24) NodeQ27_2756_com(1,5);...
    x(1,25) NodeQ27_2776_com(1,5);x(1,26) NodeQ27_2796_com(1,5);x(1,27) NodeQ27_2816_com(1,5);x(1,28) NodeQ27_2836_com(1,5);...
    x(1,29) NodeQ27_2856_com(1,5);x(1,30) NodeQ27_2876_com(1,5);x(1,31) NodeQ27_2896_com(1,5);x(1,32) NodeQ27_6_com(1,5);...
    x(1,33) NodeQ27_7_com(1,5)];


Result_phi13_com = [x(1,1) Node401_com(1,5);x(1,2) Node402_com(1,5);x(1,3) Node405_com(1,5);x(1,4) Node408_com(1,5);x(1,5) Node411_com(1,5);x(1,6) Node414_com(1,5);
    x(1,7) Node417_com(1,5);x(1,8) Node420_com(1,5);x(1,9) Node423_com(1,5);x(1,10) Node426_com(1,5);x(1,11) Node429_com(1,5);x(1,12) Node432_com(1,5);
    x(1,13) Node435_com(1,5);x(1,14) Node438_com(1,5);x(1,15) Node441_com(1,5);x(1,16) Node444_com(1,5);x(1,17) Node447_com(1,5);x(1,18) Node450_com(1,5);
    x(1,19) Node453_com(1,5);x(1,20) Node456_com(1,5);x(1,21) Node459_com(1,5);x(1,22) Node462_com(1,5);x(1,23) Node465_com(1,5);x(1,24) Node468_com(1,5);
    x(1,25) Node471_com(1,5);x(1,26) Node474_com(1,5);x(1,27) Node477_com(1,5);x(1,28) Node480_com(1,5);x(1,29) Node483_com(1,5);x(1,30) Node486_com(1,5);
    x(1,31) Node489_com(1,5);x(1,32) Node492_com(1,5);x(1,33) Node495_com(1,5)];


Result_phi13_rot = [x(1,1) Node401_rot(1,5);x(1,2) Node402_rot(1,5);x(1,3) Node405_rot(1,5);x(1,4) Node408_rot(1,5);x(1,5) Node411_rot(1,5);x(1,6) Node414_rot(1,5);
    x(1,7) Node417_rot(1,5);x(1,8) Node420_rot(1,5);x(1,9) Node423_rot(1,5);x(1,10) Node426_rot(1,5);x(1,11) Node429_rot(1,5);x(1,12) Node432_rot(1,5);
    x(1,13) Node435_rot(1,5);x(1,14) Node438_rot(1,5);x(1,15) Node441_rot(1,5);x(1,16) Node444_rot(1,5);x(1,17) Node447_rot(1,5);x(1,18) Node450_rot(1,5);
    x(1,19) Node453_rot(1,5);x(1,20) Node456_rot(1,5);x(1,21) Node459_rot(1,5);x(1,22) Node462_rot(1,5);x(1,23) Node465_rot(1,5);x(1,24) Node468_rot(1,5);
    x(1,25) Node471_rot(1,5);x(1,26) Node474_rot(1,5);x(1,27) Node477_rot(1,5);x(1,28) Node480_rot(1,5);x(1,29) Node483_rot(1,5);x(1,30) Node486_rot(1,5);
    x(1,31) Node489_rot(1,5);x(1,32) Node492_rot(1,5);x(1,33) Node495_rot(1,5)];


Result_phi13_str = [x(1,1) Node401_str(1,5);x(1,2) Node402_str(1,5);x(1,3) Node405_str(1,5);x(1,4) Node408_str(1,5);x(1,5) Node411_str(1,5);x(1,6) Node414_str(1,5);
    x(1,7) Node417_str(1,5);x(1,8) Node420_str(1,5);x(1,9) Node423_str(1,5);x(1,10) Node426_str(1,5);x(1,11) Node429_str(1,5);x(1,12) Node432_str(1,5);
    x(1,13) Node435_str(1,5);x(1,14) Node438_str(1,5);x(1,15) Node441_str(1,5);x(1,16) Node444_str(1,5);x(1,17) Node447_str(1,5);x(1,18) Node450_str(1,5);
    x(1,19) Node453_str(1,5);x(1,20) Node456_str(1,5);x(1,21) Node459_str(1,5);x(1,22) Node462_str(1,5);x(1,23) Node465_str(1,5);x(1,24) Node468_str(1,5);
    x(1,25) Node471_str(1,5);x(1,26) Node474_str(1,5);x(1,27) Node477_str(1,5);x(1,28) Node480_str(1,5);x(1,29) Node483_str(1,5);x(1,30) Node486_str(1,5);
    x(1,31) Node489_str(1,5);x(1,32) Node492_str(1,5);x(1,33) Node495_str(1,5)];

Result_phi33_com = [x(1,1) Node401_com(1,6);x(1,2) Node402_com(1,6);x(1,3) Node405_com(1,6);x(1,4) Node408_com(1,6);x(1,5) Node411_com(1,6);x(1,6) Node414_com(1,6);
    x(1,7) Node417_com(1,6);x(1,8) Node420_com(1,6);x(1,9) Node423_com(1,6);x(1,10) Node426_com(1,6);x(1,11) Node429_com(1,6);x(1,12) Node432_com(1,6);
    x(1,13) Node435_com(1,6);x(1,14) Node438_com(1,6);x(1,15) Node441_com(1,6);x(1,16) Node444_com(1,6);x(1,17) Node447_com(1,6);x(1,18) Node450_com(1,6);
    x(1,19) Node453_com(1,6);x(1,20) Node456_com(1,6);x(1,21) Node459_com(1,6);x(1,22) Node462_com(1,6);x(1,23) Node465_com(1,6);x(1,24) Node468_com(1,6);
    x(1,25) Node471_com(1,6);x(1,26) Node474_com(1,6);x(1,27) Node477_com(1,6);x(1,28) Node480_com(1,6);x(1,29) Node483_com(1,6);x(1,30) Node486_com(1,6);
    x(1,31) Node489_com(1,6);x(1,32) Node492_com(1,6);x(1,33) Node495_com(1,6)];


Result_phi33_rot = [x(1,1) Node401_rot(1,6);x(1,2) Node402_rot(1,6);x(1,3) Node405_rot(1,6);x(1,4) Node408_rot(1,6);x(1,5) Node411_rot(1,6);x(1,6) Node414_rot(1,6);
    x(1,7) Node417_rot(1,6);x(1,8) Node420_rot(1,6);x(1,9) Node423_rot(1,6);x(1,10) Node426_rot(1,6);x(1,11) Node429_rot(1,6);x(1,12) Node432_rot(1,6);
    x(1,13) Node435_rot(1,6);x(1,14) Node438_rot(1,6);x(1,15) Node441_rot(1,6);x(1,16) Node444_rot(1,6);x(1,17) Node447_rot(1,6);x(1,18) Node450_rot(1,6);
    x(1,19) Node453_rot(1,6);x(1,20) Node456_rot(1,6);x(1,21) Node459_rot(1,6);x(1,22) Node462_rot(1,6);x(1,23) Node465_rot(1,6);x(1,24) Node468_rot(1,6);
    x(1,25) Node471_rot(1,6);x(1,26) Node474_rot(1,6);x(1,27) Node477_rot(1,6);x(1,28) Node480_rot(1,6);x(1,29) Node483_rot(1,6);x(1,30) Node486_rot(1,6);
    x(1,31) Node489_rot(1,6);x(1,32) Node492_rot(1,6);x(1,33) Node495_rot(1,6)];


Result_phi33_str = [x(1,1) Node401_str(1,6);x(1,2) Node402_str(1,6);x(1,3) Node405_str(1,6);x(1,4) Node408_str(1,6);x(1,5) Node411_str(1,6);x(1,6) Node414_str(1,6);
    x(1,7) Node417_str(1,6);x(1,8) Node420_str(1,6);x(1,9) Node423_str(1,6);x(1,10) Node426_str(1,6);x(1,11) Node429_str(1,6);x(1,12) Node432_str(1,6);
    x(1,13) Node435_str(1,6);x(1,14) Node438_str(1,6);x(1,15) Node441_str(1,6);x(1,16) Node444_str(1,6);x(1,17) Node447_str(1,6);x(1,18) Node450_str(1,6);
    x(1,19) Node453_str(1,6);x(1,20) Node456_str(1,6);x(1,21) Node459_str(1,6);x(1,22) Node462_str(1,6);x(1,23) Node465_str(1,6);x(1,24) Node468_str(1,6);
    x(1,25) Node471_str(1,6);x(1,26) Node474_str(1,6);x(1,27) Node477_str(1,6);x(1,28) Node480_str(1,6);x(1,29) Node483_str(1,6);x(1,30) Node486_str(1,6);
    x(1,31) Node489_str(1,6);x(1,32) Node492_str(1,6);x(1,33) Node495_str(1,6)];

Result_dx_com = [x(1,1) Node401_com(1,7);x(1,2) Node402_com(1,7);x(1,3) Node405_com(1,7);x(1,4) Node408_com(1,7);x(1,5) Node411_com(1,7);x(1,6) Node414_com(1,7);
    x(1,7) Node417_com(1,7);x(1,8) Node420_com(1,7);x(1,9) Node423_com(1,7);x(1,10) Node426_com(1,7);x(1,11) Node429_com(1,7);x(1,12) Node432_com(1,7);
    x(1,13) Node435_com(1,7);x(1,14) Node438_com(1,7);x(1,15) Node441_com(1,7);x(1,16) Node444_com(1,7);x(1,17) Node447_com(1,7);x(1,18) Node450_com(1,7);
    x(1,19) Node453_com(1,7);x(1,20) Node456_com(1,7);x(1,21) Node459_com(1,7);x(1,22) Node462_com(1,7);x(1,23) Node465_com(1,7);x(1,24) Node468_com(1,7);
    x(1,25) Node471_com(1,7);x(1,26) Node474_com(1,7);x(1,27) Node477_com(1,7);x(1,28) Node480_com(1,7);x(1,29) Node483_com(1,7);x(1,30) Node486_com(1,7);
    x(1,31) Node489_com(1,7);x(1,32) Node492_com(1,7);x(1,33) Node495_com(1,7)];


Result_dx_rot = [x(1,1) Node401_rot(1,7);x(1,2) Node402_rot(1,7);x(1,3) Node405_rot(1,7);x(1,4) Node408_rot(1,7);x(1,5) Node411_rot(1,7);x(1,6) Node414_rot(1,7);
    x(1,7) Node417_rot(1,7);x(1,8) Node420_rot(1,7);x(1,9) Node423_rot(1,7);x(1,10) Node426_rot(1,7);x(1,11) Node429_rot(1,7);x(1,12) Node432_rot(1,7);
    x(1,13) Node435_rot(1,7);x(1,14) Node438_rot(1,7);x(1,15) Node441_rot(1,7);x(1,16) Node444_rot(1,7);x(1,17) Node447_rot(1,7);x(1,18) Node450_rot(1,7);
    x(1,19) Node453_rot(1,7);x(1,20) Node456_rot(1,7);x(1,21) Node459_rot(1,7);x(1,22) Node462_rot(1,7);x(1,23) Node465_rot(1,7);x(1,24) Node468_rot(1,7);
    x(1,25) Node471_rot(1,7);x(1,26) Node474_rot(1,7);x(1,27) Node477_rot(1,7);x(1,28) Node480_rot(1,7);x(1,29) Node483_rot(1,7);x(1,30) Node486_rot(1,7);
    x(1,31) Node489_rot(1,7);x(1,32) Node492_rot(1,7);x(1,33) Node495_rot(1,7)];


Result_dx_str = [x(1,1) Node401_str(1,7);x(1,2) Node402_str(1,7);x(1,3) Node405_str(1,7);x(1,4) Node408_str(1,7);x(1,5) Node411_str(1,7);x(1,6) Node414_str(1,7);
    x(1,7) Node417_str(1,7);x(1,8) Node420_str(1,7);x(1,9) Node423_str(1,7);x(1,10) Node426_str(1,7);x(1,11) Node429_str(1,7);x(1,12) Node432_str(1,7);
    x(1,13) Node435_str(1,7);x(1,14) Node438_str(1,7);x(1,15) Node441_str(1,7);x(1,16) Node444_str(1,7);x(1,17) Node447_str(1,7);x(1,18) Node450_str(1,7);
    x(1,19) Node453_str(1,7);x(1,20) Node456_str(1,7);x(1,21) Node459_str(1,7);x(1,22) Node462_str(1,7);x(1,23) Node465_str(1,7);x(1,24) Node468_str(1,7);
    x(1,25) Node471_str(1,7);x(1,26) Node474_str(1,7);x(1,27) Node477_str(1,7);x(1,28) Node480_str(1,7);x(1,29) Node483_str(1,7);x(1,30) Node486_str(1,7);
    x(1,31) Node489_str(1,7);x(1,32) Node492_str(1,7);x(1,33) Node495_str(1,7)];

Result_dx_cls = [x(1,1) Node401_cls(1,7);x(1,2) Node402_cls(1,7);x(1,3) Node405_cls(1,7);x(1,4) Node408_cls(1,7);x(1,5) Node411_cls(1,7);x(1,6) Node414_cls(1,7);
    x(1,7) Node417_cls(1,7);x(1,8) Node420_cls(1,7);x(1,9) Node423_cls(1,7);x(1,10) Node426_cls(1,7);x(1,11) Node429_cls(1,7);x(1,12) Node432_cls(1,7);
    x(1,13) Node435_cls(1,7);x(1,14) Node438_cls(1,7);x(1,15) Node441_cls(1,7);x(1,16) Node444_cls(1,7);x(1,17) Node447_cls(1,7);x(1,18) Node450_cls(1,7);
    x(1,19) Node453_cls(1,7);x(1,20) Node456_cls(1,7);x(1,21) Node459_cls(1,7);x(1,22) Node462_cls(1,7);x(1,23) Node465_cls(1,7);x(1,24) Node468_cls(1,7);
    x(1,25) Node471_cls(1,7);x(1,26) Node474_cls(1,7);x(1,27) Node477_cls(1,7);x(1,28) Node480_cls(1,7);x(1,29) Node483_cls(1,7);x(1,30) Node486_cls(1,7);
    x(1,31) Node489_cls(1,7);x(1,32) Node492_cls(1,7);x(1,33) Node495_cls(1,7)];


Result_dy_com = [x(1,1) Node401_com(1,8);x(1,2) Node402_com(1,8);x(1,3) Node405_com(1,8);x(1,4) Node408_com(1,8);x(1,5) Node411_com(1,8);x(1,6) Node414_com(1,8);
    x(1,7) Node417_com(1,8);x(1,8) Node420_com(1,8);x(1,9) Node423_com(1,8);x(1,10) Node426_com(1,8);x(1,11) Node429_com(1,8);x(1,12) Node432_com(1,8);
    x(1,13) Node435_com(1,8);x(1,14) Node438_com(1,8);x(1,15) Node441_com(1,8);x(1,16) Node444_com(1,8);x(1,17) Node447_com(1,8);x(1,18) Node450_com(1,8);
    x(1,19) Node453_com(1,8);x(1,20) Node456_com(1,8);x(1,21) Node459_com(1,8);x(1,22) Node462_com(1,8);x(1,23) Node465_com(1,8);x(1,24) Node468_com(1,8);
    x(1,25) Node471_com(1,8);x(1,26) Node474_com(1,8);x(1,27) Node477_com(1,8);x(1,28) Node480_com(1,8);x(1,29) Node483_com(1,8);x(1,30) Node486_com(1,8);
    x(1,31) Node489_com(1,8);x(1,32) Node492_com(1,8);x(1,33) Node495_com(1,8)];


Result_dy_rot = [x(1,1) Node401_rot(1,8);x(1,2) Node402_rot(1,8);x(1,3) Node405_rot(1,8);x(1,4) Node408_rot(1,8);x(1,5) Node411_rot(1,8);x(1,6) Node414_rot(1,8);
    x(1,7) Node417_rot(1,8);x(1,8) Node420_rot(1,8);x(1,9) Node423_rot(1,8);x(1,10) Node426_rot(1,8);x(1,11) Node429_rot(1,8);x(1,12) Node432_rot(1,8);
    x(1,13) Node435_rot(1,8);x(1,14) Node438_rot(1,8);x(1,15) Node441_rot(1,8);x(1,16) Node444_rot(1,8);x(1,17) Node447_rot(1,8);x(1,18) Node450_rot(1,8);
    x(1,19) Node453_rot(1,8);x(1,20) Node456_rot(1,8);x(1,21) Node459_rot(1,8);x(1,22) Node462_rot(1,8);x(1,23) Node465_rot(1,8);x(1,24) Node468_rot(1,8);
    x(1,25) Node471_rot(1,8);x(1,26) Node474_rot(1,8);x(1,27) Node477_rot(1,8);x(1,28) Node480_rot(1,8);x(1,29) Node483_rot(1,8);x(1,30) Node486_rot(1,8);
    x(1,31) Node489_rot(1,8);x(1,32) Node492_rot(1,8);x(1,33) Node495_rot(1,8)];


Result_dy_str = [x(1,1) Node401_str(1,8);x(1,2) Node402_str(1,8);x(1,3) Node405_str(1,8);x(1,4) Node408_str(1,8);x(1,5) Node411_str(1,8);x(1,6) Node414_str(1,8);
    x(1,7) Node417_str(1,8);x(1,8) Node420_str(1,8);x(1,9) Node423_str(1,8);x(1,10) Node426_str(1,8);x(1,11) Node429_str(1,8);x(1,12) Node432_str(1,8);
    x(1,13) Node435_str(1,8);x(1,14) Node438_str(1,8);x(1,15) Node441_str(1,8);x(1,16) Node444_str(1,8);x(1,17) Node447_str(1,8);x(1,18) Node450_str(1,8);
    x(1,19) Node453_str(1,8);x(1,20) Node456_str(1,8);x(1,21) Node459_str(1,8);x(1,22) Node462_str(1,8);x(1,23) Node465_str(1,8);x(1,24) Node468_str(1,8);
    x(1,25) Node471_str(1,8);x(1,26) Node474_str(1,8);x(1,27) Node477_str(1,8);x(1,28) Node480_str(1,8);x(1,29) Node483_str(1,8);x(1,30) Node486_str(1,8);
    x(1,31) Node489_str(1,8);x(1,32) Node492_str(1,8);x(1,33) Node495_str(1,8)];

Result_dy_cls = [x(1,1) Node401_cls(1,8);x(1,2) Node402_cls(1,8);x(1,3) Node405_cls(1,8);x(1,4) Node408_cls(1,8);x(1,5) Node411_cls(1,8);x(1,6) Node414_cls(1,8);
    x(1,7) Node417_cls(1,8);x(1,8) Node420_cls(1,8);x(1,9) Node423_cls(1,8);x(1,10) Node426_cls(1,8);x(1,11) Node429_cls(1,8);x(1,12) Node432_cls(1,8);
    x(1,13) Node435_cls(1,8);x(1,14) Node438_cls(1,8);x(1,15) Node441_cls(1,8);x(1,16) Node444_cls(1,8);x(1,17) Node447_cls(1,8);x(1,18) Node450_cls(1,8);
    x(1,19) Node453_cls(1,8);x(1,20) Node456_cls(1,8);x(1,21) Node459_cls(1,8);x(1,22) Node462_cls(1,8);x(1,23) Node465_cls(1,8);x(1,24) Node468_cls(1,8);
    x(1,25) Node471_cls(1,8);x(1,26) Node474_cls(1,8);x(1,27) Node477_cls(1,8);x(1,28) Node480_cls(1,8);x(1,29) Node483_cls(1,8);x(1,30) Node486_cls(1,8);
    x(1,31) Node489_cls(1,8);x(1,32) Node492_cls(1,8);x(1,33) Node495_cls(1,8)];

Result_dz_cls = [x(1,1) Node401_cls(1,9);x(1,2) Node402_cls(1,9);x(1,3) Node405_cls(1,9);x(1,4) Node408_cls(1,9);x(1,5) Node411_cls(1,9);x(1,6) Node414_cls(1,9);
    x(1,7) Node417_cls(1,9);x(1,8) Node420_cls(1,9);x(1,9) Node423_cls(1,9);x(1,10) Node426_cls(1,9);x(1,11) Node429_cls(1,9);x(1,12) Node432_cls(1,9);
    x(1,13) Node435_cls(1,9);x(1,14) Node438_cls(1,9);x(1,15) Node441_cls(1,9);x(1,16) Node444_cls(1,9);x(1,17) Node447_cls(1,9);x(1,18) Node450_cls(1,9);
    x(1,19) Node453_cls(1,9);x(1,20) Node456_cls(1,9);x(1,21) Node459_cls(1,9);x(1,22) Node462_cls(1,9);x(1,23) Node465_cls(1,9);x(1,24) Node468_cls(1,9);
    x(1,25) Node471_cls(1,9);x(1,26) Node474_cls(1,9);x(1,27) Node477_cls(1,9);x(1,28) Node480_cls(1,9);x(1,29) Node483_cls(1,9);x(1,30) Node486_cls(1,9);
    x(1,31) Node489_cls(1,9);x(1,32) Node492_cls(1,9);x(1,33) Node495_cls(1,9)];

Result_dz_com = [x(1,1) Node401_com(1,9);x(1,2) Node402_com(1,9);x(1,3) Node405_com(1,9);x(1,4) Node408_com(1,9);x(1,5) Node411_com(1,9);x(1,6) Node414_com(1,9);
    x(1,7) Node417_com(1,9);x(1,8) Node420_com(1,9);x(1,9) Node423_com(1,9);x(1,10) Node426_com(1,9);x(1,11) Node429_com(1,9);x(1,12) Node432_com(1,9);
    x(1,13) Node435_com(1,9);x(1,14) Node438_com(1,9);x(1,15) Node441_com(1,9);x(1,16) Node444_com(1,9);x(1,17) Node447_com(1,9);x(1,18) Node450_com(1,9);
    x(1,19) Node453_com(1,9);x(1,20) Node456_com(1,9);x(1,21) Node459_com(1,9);x(1,22) Node462_com(1,9);x(1,23) Node465_com(1,9);x(1,24) Node468_com(1,9);
    x(1,25) Node471_com(1,9);x(1,26) Node474_com(1,9);x(1,27) Node477_com(1,9);x(1,28) Node480_com(1,9);x(1,29) Node483_com(1,9);x(1,30) Node486_com(1,9);
    x(1,31) Node489_com(1,9);x(1,32) Node492_com(1,9);x(1,33) Node495_com(1,9)];


Result_dz_rot = [x(1,1) Node401_rot(1,9);x(1,2) Node402_rot(1,9);x(1,3) Node405_rot(1,9);x(1,4) Node408_rot(1,9);x(1,5) Node411_rot(1,9);x(1,6) Node414_rot(1,9);
    x(1,7) Node417_rot(1,9);x(1,8) Node420_rot(1,9);x(1,9) Node423_rot(1,9);x(1,10) Node426_rot(1,9);x(1,11) Node429_rot(1,9);x(1,12) Node432_rot(1,9);
    x(1,13) Node435_rot(1,9);x(1,14) Node438_rot(1,9);x(1,15) Node441_rot(1,9);x(1,16) Node444_rot(1,9);x(1,17) Node447_rot(1,9);x(1,18) Node450_rot(1,9);
    x(1,19) Node453_rot(1,9);x(1,20) Node456_rot(1,9);x(1,21) Node459_rot(1,9);x(1,22) Node462_rot(1,9);x(1,23) Node465_rot(1,9);x(1,24) Node468_rot(1,9);
    x(1,25) Node471_rot(1,9);x(1,26) Node474_rot(1,9);x(1,27) Node477_rot(1,9);x(1,28) Node480_rot(1,9);x(1,29) Node483_rot(1,9);x(1,30) Node486_rot(1,9);
    x(1,31) Node489_rot(1,9);x(1,32) Node492_rot(1,9);x(1,33) Node495_rot(1,9)];


Result_dz_str = [x(1,1) Node401_str(1,9);x(1,2) Node402_str(1,9);x(1,3) Node405_str(1,9);x(1,4) Node408_str(1,9);x(1,5) Node411_str(1,9);x(1,6) Node414_str(1,9);
    x(1,7) Node417_str(1,9);x(1,8) Node420_str(1,9);x(1,9) Node423_str(1,9);x(1,10) Node426_str(1,9);x(1,11) Node429_str(1,9);x(1,12) Node432_str(1,9);
    x(1,13) Node435_str(1,9);x(1,14) Node438_str(1,9);x(1,15) Node441_str(1,9);x(1,16) Node444_str(1,9);x(1,17) Node447_str(1,9);x(1,18) Node450_str(1,9);
    x(1,19) Node453_str(1,9);x(1,20) Node456_str(1,9);x(1,21) Node459_str(1,9);x(1,22) Node462_str(1,9);x(1,23) Node465_str(1,9);x(1,24) Node468_str(1,9);
    x(1,25) Node471_str(1,9);x(1,26) Node474_str(1,9);x(1,27) Node477_str(1,9);x(1,28) Node480_str(1,9);x(1,29) Node483_str(1,9);x(1,30) Node486_str(1,9);
    x(1,31) Node489_str(1,9);x(1,32) Node492_str(1,9);x(1,33) Node495_str(1,9)];

Result_dz_Q27_com = [x(1,1) NodeQ27_2312_com(1,9);x(1,2) NodeQ27_2313_com(1,9);x(1,3) NodeQ27_2336_com(1,9);x(1,4) NodeQ27_2356_com(1,9);...
    x(1,5) NodeQ27_2376_com(1,9);x(1,6) NodeQ27_2396_com(1,9);x(1,7) NodeQ27_2416_com(1,9);x(1,8) NodeQ27_2436_com(1,9);...
    x(1,9) NodeQ27_2456_com(1,9);x(1,10) NodeQ27_2476_com(1,9);x(1,11) NodeQ27_2496_com(1,9);x(1,12) NodeQ27_2516_com(1,9);...
    x(1,13) NodeQ27_2536_com(1,9);x(1,14) NodeQ27_2556_com(1,9);x(1,15) NodeQ27_2576_com(1,9);x(1,16) NodeQ27_2596_com(1,9);...
    x(1,17) NodeQ27_2616_com(1,9);x(1,18) NodeQ27_2636_com(1,9);x(1,19) NodeQ27_2656_com(1,9);x(1,20) NodeQ27_2676_com(1,9);...
    x(1,21) NodeQ27_2696_com(1,9);x(1,22) NodeQ27_2716_com(1,9);x(1,23) NodeQ27_2736_com(1,9);x(1,24) NodeQ27_2756_com(1,9);...
    x(1,25) NodeQ27_2776_com(1,9);x(1,26) NodeQ27_2796_com(1,9);x(1,27) NodeQ27_2816_com(1,9);x(1,28) NodeQ27_2836_com(1,9);...
    x(1,29) NodeQ27_2856_com(1,9);x(1,30) NodeQ27_2876_com(1,9);x(1,31) NodeQ27_2896_com(1,9);x(1,32) NodeQ27_6_com(1,9);...
    x(1,33) NodeQ27_7_com(1,9)];

Result_S11_com = [x(1,1) Node401_com(1,10);x(1,2) Node402_com(1,10);x(1,3) Node405_com(1,10);x(1,4) Node408_com(1,10);x(1,5) Node411_com(1,10);x(1,6) Node414_com(1,10);
    x(1,7) Node417_com(1,10);x(1,8) Node420_com(1,10);x(1,9) Node423_com(1,10);x(1,10) Node426_com(1,10);x(1,11) Node429_com(1,10);x(1,12) Node432_com(1,10);
    x(1,13) Node435_com(1,10);x(1,14) Node438_com(1,10);x(1,15) Node441_com(1,10);x(1,16) Node444_com(1,10);x(1,17) Node447_com(1,10);x(1,18) Node450_com(1,10);
    x(1,19) Node453_com(1,10);x(1,20) Node456_com(1,10);x(1,21) Node459_com(1,10);x(1,22) Node462_com(1,10);x(1,23) Node465_com(1,10);x(1,24) Node468_com(1,10);
    x(1,25) Node471_com(1,10);x(1,26) Node474_com(1,10);x(1,27) Node477_com(1,10);x(1,28) Node480_com(1,10);x(1,29) Node483_com(1,10);x(1,30) Node486_com(1,10);
    x(1,31) Node489_com(1,10);x(1,32) Node492_com(1,10);x(1,33) Node495_com(1,10)];


Result_S11_rot = [x(1,1) Node401_rot(1,10);x(1,2) Node402_rot(1,10);x(1,3) Node405_rot(1,10);x(1,4) Node408_rot(1,10);x(1,5) Node411_rot(1,10);x(1,6) Node414_rot(1,10);
    x(1,7) Node417_rot(1,10);x(1,8) Node420_rot(1,10);x(1,9) Node423_rot(1,10);x(1,10) Node426_rot(1,10);x(1,11) Node429_rot(1,10);x(1,12) Node432_rot(1,10);
    x(1,13) Node435_rot(1,10);x(1,14) Node438_rot(1,10);x(1,15) Node441_rot(1,10);x(1,16) Node444_rot(1,10);x(1,17) Node447_rot(1,10);x(1,18) Node450_rot(1,10);
    x(1,19) Node453_rot(1,10);x(1,20) Node456_rot(1,10);x(1,21) Node459_rot(1,10);x(1,22) Node462_rot(1,10);x(1,23) Node465_rot(1,10);x(1,24) Node468_rot(1,10);
    x(1,25) Node471_rot(1,10);x(1,26) Node474_rot(1,10);x(1,27) Node477_rot(1,10);x(1,28) Node480_rot(1,10);x(1,29) Node483_rot(1,10);x(1,30) Node486_rot(1,10);
    x(1,31) Node489_rot(1,10);x(1,32) Node492_rot(1,10);x(1,33) Node495_rot(1,10)];


Result_S11_str = [x(1,1) Node401_str(1,10);x(1,2) Node402_str(1,10);x(1,3) Node405_str(1,10);x(1,4) Node408_str(1,10);x(1,5) Node411_str(1,10);x(1,6) Node414_str(1,10);
    x(1,7) Node417_str(1,10);x(1,8) Node420_str(1,10);x(1,9) Node423_str(1,10);x(1,10) Node426_str(1,10);x(1,11) Node429_str(1,10);x(1,12) Node432_str(1,10);
    x(1,13) Node435_str(1,10);x(1,14) Node438_str(1,10);x(1,15) Node441_str(1,10);x(1,16) Node444_str(1,10);x(1,17) Node447_str(1,10);x(1,18) Node450_str(1,10);
    x(1,19) Node453_str(1,10);x(1,20) Node456_str(1,10);x(1,21) Node459_str(1,10);x(1,22) Node462_str(1,10);x(1,23) Node465_str(1,10);x(1,24) Node468_str(1,10);
    x(1,25) Node471_str(1,10);x(1,26) Node474_str(1,10);x(1,27) Node477_str(1,10);x(1,28) Node480_str(1,10);x(1,29) Node483_str(1,10);x(1,30) Node486_str(1,10);
    x(1,31) Node489_str(1,10);x(1,32) Node492_str(1,10);x(1,33) Node495_str(1,10)];

Result_S11_cls = [x(1,1) Node401_cls(1,10);x(1,2) Node402_cls(1,10);x(1,3) Node405_cls(1,10);x(1,4) Node408_cls(1,10);x(1,5) Node411_cls(1,10);x(1,6) Node414_cls(1,10);
    x(1,7) Node417_cls(1,10);x(1,8) Node420_cls(1,10);x(1,9) Node423_cls(1,10);x(1,10) Node426_cls(1,10);x(1,11) Node429_cls(1,10);x(1,12) Node432_cls(1,10);
    x(1,13) Node435_cls(1,10);x(1,14) Node438_cls(1,10);x(1,15) Node441_cls(1,10);x(1,16) Node444_cls(1,10);x(1,17) Node447_cls(1,10);x(1,18) Node450_cls(1,10);
    x(1,19) Node453_cls(1,10);x(1,20) Node456_cls(1,10);x(1,21) Node459_cls(1,10);x(1,22) Node462_cls(1,10);x(1,23) Node465_cls(1,10);x(1,24) Node468_cls(1,10);
    x(1,25) Node471_cls(1,10);x(1,26) Node474_cls(1,10);x(1,27) Node477_cls(1,10);x(1,28) Node480_cls(1,10);x(1,29) Node483_cls(1,10);x(1,30) Node486_cls(1,10);
    x(1,31) Node489_cls(1,10);x(1,32) Node492_cls(1,10);x(1,33) Node495_cls(1,10)];


Result_S33_com = [x(1,1) Node401_com(1,11);x(1,2) Node402_com(1,11);x(1,3) Node405_com(1,11);x(1,4) Node408_com(1,11);x(1,5) Node411_com(1,11);x(1,6) Node414_com(1,11);
    x(1,7) Node417_com(1,11);x(1,8) Node420_com(1,11);x(1,9) Node423_com(1,11);x(1,10) Node426_com(1,11);x(1,11) Node429_com(1,11);x(1,12) Node432_com(1,11);
    x(1,13) Node435_com(1,11);x(1,14) Node438_com(1,11);x(1,15) Node441_com(1,11);x(1,16) Node444_com(1,11);x(1,17) Node447_com(1,11);x(1,18) Node450_com(1,11);
    x(1,19) Node453_com(1,11);x(1,20) Node456_com(1,11);x(1,21) Node459_com(1,11);x(1,22) Node462_com(1,11);x(1,23) Node465_com(1,11);x(1,24) Node468_com(1,11);
    x(1,25) Node471_com(1,11);x(1,26) Node474_com(1,11);x(1,27) Node477_com(1,11);x(1,28) Node480_com(1,11);x(1,29) Node483_com(1,11);x(1,30) Node486_com(1,11);
    x(1,31) Node489_com(1,11);x(1,32) Node492_com(1,11);x(1,33) Node495_com(1,11)];


Result_S33_rot = [x(1,1) Node401_rot(1,11);x(1,2) Node402_rot(1,11);x(1,3) Node405_rot(1,11);x(1,4) Node408_rot(1,11);x(1,5) Node411_rot(1,11);x(1,6) Node414_rot(1,11);
    x(1,7) Node417_rot(1,11);x(1,8) Node420_rot(1,11);x(1,9) Node423_rot(1,11);x(1,10) Node426_rot(1,11);x(1,11) Node429_rot(1,11);x(1,12) Node432_rot(1,11);
    x(1,13) Node435_rot(1,11);x(1,14) Node438_rot(1,11);x(1,15) Node441_rot(1,11);x(1,16) Node444_rot(1,11);x(1,17) Node447_rot(1,11);x(1,18) Node450_rot(1,11);
    x(1,19) Node453_rot(1,11);x(1,20) Node456_rot(1,11);x(1,21) Node459_rot(1,11);x(1,22) Node462_rot(1,11);x(1,23) Node465_rot(1,11);x(1,24) Node468_rot(1,11);
    x(1,25) Node471_rot(1,11);x(1,26) Node474_rot(1,11);x(1,27) Node477_rot(1,11);x(1,28) Node480_rot(1,11);x(1,29) Node483_rot(1,11);x(1,30) Node486_rot(1,11);
    x(1,31) Node489_rot(1,11);x(1,32) Node492_rot(1,11);x(1,33) Node495_rot(1,11)];


Result_S33_str = [x(1,1) Node401_str(1,11);x(1,2) Node402_str(1,11);x(1,3) Node405_str(1,11);x(1,4) Node408_str(1,11);x(1,5) Node411_str(1,11);x(1,6) Node414_str(1,11);
    x(1,7) Node417_str(1,11);x(1,8) Node420_str(1,11);x(1,9) Node423_str(1,11);x(1,10) Node426_str(1,11);x(1,11) Node429_str(1,11);x(1,12) Node432_str(1,11);
    x(1,13) Node435_str(1,11);x(1,14) Node438_str(1,11);x(1,15) Node441_str(1,11);x(1,16) Node444_str(1,11);x(1,17) Node447_str(1,11);x(1,18) Node450_str(1,11);
    x(1,19) Node453_str(1,11);x(1,20) Node456_str(1,11);x(1,21) Node459_str(1,11);x(1,22) Node462_str(1,11);x(1,23) Node465_str(1,11);x(1,24) Node468_str(1,11);
    x(1,25) Node471_str(1,11);x(1,26) Node474_str(1,11);x(1,27) Node477_str(1,11);x(1,28) Node480_str(1,11);x(1,29) Node483_str(1,11);x(1,30) Node486_str(1,11);
    x(1,31) Node489_str(1,11);x(1,32) Node492_str(1,11);x(1,33) Node495_str(1,11)];


Result_S33_cls = [x(1,1) Node401_cls(1,11);x(1,2) Node402_cls(1,11);x(1,3) Node405_cls(1,11);x(1,4) Node408_cls(1,11);x(1,5) Node411_cls(1,11);x(1,6) Node414_cls(1,11);
    x(1,7) Node417_cls(1,11);x(1,8) Node420_cls(1,11);x(1,9) Node423_cls(1,11);x(1,10) Node426_cls(1,11);x(1,11) Node429_cls(1,11);x(1,12) Node432_cls(1,11);
    x(1,13) Node435_cls(1,11);x(1,14) Node438_cls(1,11);x(1,15) Node441_cls(1,11);x(1,16) Node444_cls(1,11);x(1,17) Node447_cls(1,11);x(1,18) Node450_cls(1,11);
    x(1,19) Node453_cls(1,11);x(1,20) Node456_cls(1,11);x(1,21) Node459_cls(1,11);x(1,22) Node462_cls(1,11);x(1,23) Node465_cls(1,11);x(1,24) Node468_cls(1,11);
    x(1,25) Node471_cls(1,11);x(1,26) Node474_cls(1,11);x(1,27) Node477_cls(1,11);x(1,28) Node480_cls(1,11);x(1,29) Node483_cls(1,11);x(1,30) Node486_cls(1,11);
    x(1,31) Node489_cls(1,11);x(1,32) Node492_cls(1,11);x(1,33) Node495_cls(1,11)];


Result_Sig33_com = [x(1,1) Node401_com(1,12);x(1,2) Node402_com(1,12);x(1,3) Node405_com(1,12);x(1,4) Node408_com(1,12);x(1,5) Node411_com(1,12);x(1,6) Node414_com(1,12);
    x(1,7) Node417_com(1,12);x(1,8) Node420_com(1,12);x(1,9) Node423_com(1,12);x(1,10) Node426_com(1,12);x(1,11) Node429_com(1,12);x(1,12) Node432_com(1,12);
    x(1,13) Node435_com(1,12);x(1,14) Node438_com(1,12);x(1,15) Node441_com(1,12);x(1,16) Node444_com(1,12);x(1,17) Node447_com(1,12);x(1,18) Node450_com(1,12);
    x(1,19) Node453_com(1,12);x(1,20) Node456_com(1,12);x(1,21) Node459_com(1,12);x(1,22) Node462_com(1,12);x(1,23) Node465_com(1,12);x(1,24) Node468_com(1,12);
    x(1,25) Node471_com(1,12);x(1,26) Node474_com(1,12);x(1,27) Node477_com(1,12);x(1,28) Node480_com(1,12);x(1,29) Node483_com(1,12);x(1,30) Node486_com(1,12);
    x(1,31) Node489_com(1,12);x(1,32) Node492_com(1,12);x(1,33) Node495_com(1,12)];


Result_Sig33_rot = [x(1,1) Node401_rot(1,12);x(1,2) Node402_rot(1,12);x(1,3) Node405_rot(1,12);x(1,4) Node408_rot(1,12);x(1,5) Node411_rot(1,12);x(1,6) Node414_rot(1,12);
    x(1,7) Node417_rot(1,12);x(1,8) Node420_rot(1,12);x(1,9) Node423_rot(1,12);x(1,10) Node426_rot(1,12);x(1,11) Node429_rot(1,12);x(1,12) Node432_rot(1,12);
    x(1,13) Node435_rot(1,12);x(1,14) Node438_rot(1,12);x(1,15) Node441_rot(1,12);x(1,16) Node444_rot(1,12);x(1,17) Node447_rot(1,12);x(1,18) Node450_rot(1,12);
    x(1,19) Node453_rot(1,12);x(1,20) Node456_rot(1,12);x(1,21) Node459_rot(1,12);x(1,22) Node462_rot(1,12);x(1,23) Node465_rot(1,12);x(1,24) Node468_rot(1,12);
    x(1,25) Node471_rot(1,12);x(1,26) Node474_rot(1,12);x(1,27) Node477_rot(1,12);x(1,28) Node480_rot(1,12);x(1,29) Node483_rot(1,12);x(1,30) Node486_rot(1,12);
    x(1,31) Node489_rot(1,12);x(1,32) Node492_rot(1,12);x(1,33) Node495_rot(1,12)];


Result_Sig33_str = [x(1,1) Node401_str(1,12);x(1,2) Node402_str(1,12);x(1,3) Node405_str(1,12);x(1,4) Node408_str(1,12);x(1,5) Node411_str(1,12);x(1,6) Node414_str(1,12);
    x(1,7) Node417_str(1,12);x(1,8) Node420_str(1,12);x(1,9) Node423_str(1,12);x(1,10) Node426_str(1,12);x(1,11) Node429_str(1,12);x(1,12) Node432_str(1,12);
    x(1,13) Node435_str(1,12);x(1,14) Node438_str(1,12);x(1,15) Node441_str(1,12);x(1,16) Node444_str(1,12);x(1,17) Node447_str(1,12);x(1,18) Node450_str(1,12);
    x(1,19) Node453_str(1,12);x(1,20) Node456_str(1,12);x(1,21) Node459_str(1,12);x(1,22) Node462_str(1,12);x(1,23) Node465_str(1,12);x(1,24) Node468_str(1,12);
    x(1,25) Node471_str(1,12);x(1,26) Node474_str(1,12);x(1,27) Node477_str(1,12);x(1,28) Node480_str(1,12);x(1,29) Node483_str(1,12);x(1,30) Node486_str(1,12);
    x(1,31) Node489_str(1,12);x(1,32) Node492_str(1,12);x(1,33) Node495_str(1,12)];


Result_M113_com = [x(1,1) Node401_com(1,13);x(1,2) Node402_com(1,13);x(1,3) Node405_com(1,13);x(1,4) Node408_com(1,13);x(1,5) Node411_com(1,13);x(1,6) Node414_com(1,13);
    x(1,7) Node417_com(1,13);x(1,8) Node420_com(1,13);x(1,9) Node423_com(1,13);x(1,10) Node426_com(1,13);x(1,11) Node429_com(1,13);x(1,12) Node432_com(1,13);
    x(1,13) Node435_com(1,13);x(1,14) Node438_com(1,13);x(1,15) Node441_com(1,13);x(1,16) Node444_com(1,13);x(1,17) Node447_com(1,13);x(1,18) Node450_com(1,13);
    x(1,19) Node453_com(1,13);x(1,20) Node456_com(1,13);x(1,21) Node459_com(1,13);x(1,22) Node462_com(1,13);x(1,23) Node465_com(1,13);x(1,24) Node468_com(1,13);
    x(1,25) Node471_com(1,13);x(1,26) Node474_com(1,13);x(1,27) Node477_com(1,13);x(1,28) Node480_com(1,13);x(1,29) Node483_com(1,13);x(1,30) Node486_com(1,13);
    x(1,31) Node489_com(1,13);x(1,32) Node492_com(1,13);x(1,33) Node495_com(1,13)];


Result_M113_rot = [x(1,1) Node401_rot(1,13);x(1,2) Node402_rot(1,13);x(1,3) Node405_rot(1,13);x(1,4) Node408_rot(1,13);x(1,5) Node411_rot(1,13);x(1,6) Node414_rot(1,13);
    x(1,7) Node417_rot(1,13);x(1,8) Node420_rot(1,13);x(1,9) Node423_rot(1,13);x(1,10) Node426_rot(1,13);x(1,11) Node429_rot(1,13);x(1,12) Node432_rot(1,13);
    x(1,13) Node435_rot(1,13);x(1,14) Node438_rot(1,13);x(1,15) Node441_rot(1,13);x(1,16) Node444_rot(1,13);x(1,17) Node447_rot(1,13);x(1,18) Node450_rot(1,13);
    x(1,19) Node453_rot(1,13);x(1,20) Node456_rot(1,13);x(1,21) Node459_rot(1,13);x(1,22) Node462_rot(1,13);x(1,23) Node465_rot(1,13);x(1,24) Node468_rot(1,13);
    x(1,25) Node471_rot(1,13);x(1,26) Node474_rot(1,13);x(1,27) Node477_rot(1,13);x(1,28) Node480_rot(1,13);x(1,29) Node483_rot(1,13);x(1,30) Node486_rot(1,13);
    x(1,31) Node489_rot(1,13);x(1,32) Node492_rot(1,13);x(1,33) Node495_rot(1,13)];


Result_M113_str = [x(1,1) Node401_str(1,13);x(1,2) Node402_str(1,13);x(1,3) Node405_str(1,13);x(1,4) Node408_str(1,13);x(1,5) Node411_str(1,13);x(1,6) Node414_str(1,13);
    x(1,7) Node417_str(1,13);x(1,8) Node420_str(1,13);x(1,9) Node423_str(1,13);x(1,10) Node426_str(1,13);x(1,11) Node429_str(1,13);x(1,12) Node432_str(1,13);
    x(1,13) Node435_str(1,13);x(1,14) Node438_str(1,13);x(1,15) Node441_str(1,13);x(1,16) Node444_str(1,13);x(1,17) Node447_str(1,13);x(1,18) Node450_str(1,13);
    x(1,19) Node453_str(1,13);x(1,20) Node456_str(1,13);x(1,21) Node459_str(1,13);x(1,22) Node462_str(1,13);x(1,23) Node465_str(1,13);x(1,24) Node468_str(1,13);
    x(1,25) Node471_str(1,13);x(1,26) Node474_str(1,13);x(1,27) Node477_str(1,13);x(1,28) Node480_str(1,13);x(1,29) Node483_str(1,13);x(1,30) Node486_str(1,13);
    x(1,31) Node489_str(1,13);x(1,32) Node492_str(1,13);x(1,33) Node495_str(1,13)];


Result_M131_com = [x(1,1) Node401_com(1,14);x(1,2) Node402_com(1,14);x(1,3) Node405_com(1,14);x(1,4) Node408_com(1,14);x(1,5) Node411_com(1,14);x(1,6) Node414_com(1,14);
    x(1,7) Node417_com(1,14);x(1,8) Node420_com(1,14);x(1,9) Node423_com(1,14);x(1,10) Node426_com(1,14);x(1,11) Node429_com(1,14);x(1,12) Node432_com(1,14);
    x(1,13) Node435_com(1,14);x(1,14) Node438_com(1,14);x(1,15) Node441_com(1,14);x(1,16) Node444_com(1,14);x(1,17) Node447_com(1,14);x(1,18) Node450_com(1,14);
    x(1,19) Node453_com(1,14);x(1,20) Node456_com(1,14);x(1,21) Node459_com(1,14);x(1,22) Node462_com(1,14);x(1,23) Node465_com(1,14);x(1,24) Node468_com(1,14);
    x(1,25) Node471_com(1,14);x(1,26) Node474_com(1,14);x(1,27) Node477_com(1,14);x(1,28) Node480_com(1,14);x(1,29) Node483_com(1,14);x(1,30) Node486_com(1,14);
    x(1,31) Node489_com(1,14);x(1,32) Node492_com(1,14);x(1,33) Node495_com(1,14)];


Result_M131_rot = [x(1,1) Node401_rot(1,14);x(1,2) Node402_rot(1,14);x(1,3) Node405_rot(1,14);x(1,4) Node408_rot(1,14);x(1,5) Node411_rot(1,14);x(1,6) Node414_rot(1,14);
    x(1,7) Node417_rot(1,14);x(1,8) Node420_rot(1,14);x(1,9) Node423_rot(1,14);x(1,10) Node426_rot(1,14);x(1,11) Node429_rot(1,14);x(1,12) Node432_rot(1,14);
    x(1,13) Node435_rot(1,14);x(1,14) Node438_rot(1,14);x(1,15) Node441_rot(1,14);x(1,16) Node444_rot(1,14);x(1,17) Node447_rot(1,14);x(1,18) Node450_rot(1,14);
    x(1,19) Node453_rot(1,14);x(1,20) Node456_rot(1,14);x(1,21) Node459_rot(1,14);x(1,22) Node462_rot(1,14);x(1,23) Node465_rot(1,14);x(1,24) Node468_rot(1,14);
    x(1,25) Node471_rot(1,14);x(1,26) Node474_rot(1,14);x(1,27) Node477_rot(1,14);x(1,28) Node480_rot(1,14);x(1,29) Node483_rot(1,14);x(1,30) Node486_rot(1,14);
    x(1,31) Node489_rot(1,14);x(1,32) Node492_rot(1,14);x(1,33) Node495_rot(1,14)];


Result_M131_str = [x(1,1) Node401_str(1,14);x(1,2) Node402_str(1,14);x(1,3) Node405_str(1,14);x(1,4) Node408_str(1,14);x(1,5) Node411_str(1,14);x(1,6) Node414_str(1,14);
    x(1,7) Node417_str(1,14);x(1,8) Node420_str(1,14);x(1,9) Node423_str(1,14);x(1,10) Node426_str(1,14);x(1,11) Node429_str(1,14);x(1,12) Node432_str(1,14);
    x(1,13) Node435_str(1,14);x(1,14) Node438_str(1,14);x(1,15) Node441_str(1,14);x(1,16) Node444_str(1,14);x(1,17) Node447_str(1,14);x(1,18) Node450_str(1,14);
    x(1,19) Node453_str(1,14);x(1,20) Node456_str(1,14);x(1,21) Node459_str(1,14);x(1,22) Node462_str(1,14);x(1,23) Node465_str(1,14);x(1,24) Node468_str(1,14);
    x(1,25) Node471_str(1,14);x(1,26) Node474_str(1,14);x(1,27) Node477_str(1,14);x(1,28) Node480_str(1,14);x(1,29) Node483_str(1,14);x(1,30) Node486_str(1,14);
    x(1,31) Node489_str(1,14);x(1,32) Node492_str(1,14);x(1,33) Node495_str(1,14)];



box on
Figure1 = figure(1);
hold on
plot(Result_phi11_str(:,1),Result_phi11_str(:,2),'--o',Result_phi11_com(:,1),Result_phi11_com(:,2),'--or','LineWidth',2)
set(Figure1,'defaulttextinterpreter','latex');
ylabel('$\Phi_{11}$','fontsize',18)
xlabel('$X_1$(m)','fontsize',18)
lgnd=legend({'Case B Micromorphic','Case A Micromorphic'},'Interpreter','Latex');
set(gca,'FontName','mwa_cmr10','FontSize',14)
set(lgnd, 'XColor', 'w', 'YColor', 'w', 'Color', 'none');


box on
Figure2 = figure(2);
hold on
plot(Result_phi22_com(:,1),Result_phi22_com(:,2),'--or',Result_phi22_str(:,1),Result_phi22_str(:,2),'--o','LineWidth',2)
set(Figure2,'defaulttextinterpreter','latex');
ylabel('$\Phi_{22}$','fontsize',18)
xlabel('$X_1$(m)','fontsize',18)
lgnd=legend({'Case A Micromorphic','Case B Micromorphic'},'Interpreter','Latex');
set(gca,'FontName','mwa_cmr10','FontSize',14)
set(lgnd, 'XColor', 'w', 'YColor', 'w', 'Color', 'none');



box on
Figure3 = figure(3);
hold on
plot(Result_phi33_com(:,1),Result_phi33_com(:,2),'--or',Result_phi33_str(:,1),Result_phi33_str(:,2),'--o','LineWidth',2)
set(Figure3,'defaulttextinterpreter','latex');
ylabel('$\Phi_{33}$','fontsize',18)
xlabel('$X_1$(m)','fontsize',18)
lgnd=legend({'Case A Micromorphic','Case B Micromorphic'},'Interpreter','Latex');
set(gca,'FontName','mwa_cmr10','FontSize',14)
set(lgnd, 'XColor', 'w', 'YColor', 'w', 'Color', 'none');


box on
Figure4 = figure(4);
hold on
plot(Result_phi13_rot(:,1),0.5*(Result_phi13_rot(:,2)-Result_phi31_rot(:,2)),'--og',Rotation_Beam(:,1),Rotation_Beam(:,2),'--o',...
                        Result_phi13_com(:,1),0.5*(Result_phi13_com(:,2)-Result_phi31_com(:,2)),'--or','LineWidth',2)
set(Figure4,'defaulttextinterpreter','latex');
ylabel('$\Phi_2^{rot} = \frac{1}{2}(\Phi_{13}-\Phi_{31})$','fontsize',18)
xlabel('$X_1$(m)','fontsize',18)
lgnd=legend({'Case C Micromorphic','Micropolar Ramezani et al. Eur J Mech A-Solid-2009','Case A Micromorphic'},'Interpreter','Latex');
set(gca,'FontName','mwa_cmr10','FontSize',14)
set(lgnd, 'XColor', 'w', 'YColor', 'w', 'Color', 'none');


box on
Figure5 = figure(5);
hold on
plot(Result_dx_str(:,1),Result_dx_str(:,2),'--o',Result_dx_rot(:,1),Result_dx_rot(:,2),'--ok'...
    ,Result_dx_cls(:,1),Result_dx_cls(:,2),'--og',Result_dx_com(:,1),Result_dx_com(:,2),'--or','LineWidth',2)
set(Figure5,'defaulttextinterpreter','latex');
ylabel('$u_1$','fontsize',18)
xlabel('$X_1$(m)','fontsize',18)
lgnd=legend({'Case B Micromorphic','Case C Micromorphic',...
    'Classical Continuum','Case A Micromorphic'},'Interpreter','Latex');
set(gca,'FontName','mwa_cmr10','FontSize',14)
set(lgnd, 'XColor', 'w', 'YColor', 'w', 'Color', 'none');


box on
Figure6 = figure(6);
hold on
plot(Result_dy_com(:,1),-Result_dy_com(:,2),'--or',Result_dy_cls(:,1),-Result_dy_cls(:,2),'--og'...
    ,Result_dy_rot(:,1),-Result_dy_rot(:,2),'--ok',Result_dy_str(:,1),-Result_dy_str(:,2),'--o','LineWidth',2)
set(Figure6,'defaulttextinterpreter','latex');
ylabel('$u_2$','fontsize',18)    
xlabel('$X_1$(m)','fontsize',18)
lgnd=legend({'Case A Micromorphic','Classical Continuum','Case C Micromorphic',...
    'Case B Micromorphic',},'Interpreter','Latex');
set(gca,'FontName','mwa_cmr10','FontSize',14)
set(lgnd, 'XColor', 'w', 'YColor', 'w', 'Color', 'none');

box on
Figure7 = figure(7);
hold on

CC = plot(Result_dz_com(:,1),Coeff*Result_dz_com(:,2),'--or',Lateral_deflection_classic(:,1),Lateral_deflection_classic(:,2),'c',Result_dz_cls(:,1),Coeff*Result_dz_cls(:,2),'--og',...
    MicroPolar(:,1),MicroPolar(:,2),Result_dz_rot(:,1),Coeff*Result_dz_rot(:,2),'--ok',Result_dz_str(:,1),Coeff*Result_dz_str(:,2),'--ob', 'LineWidth',2);

set(Figure7,'defaulttextinterpreter','latex');
ylabel('$3EIu_3/(FL^3)$','fontsize',18)
xlabel('$X_1$(m)','fontsize',18)
lgnd=legend({'Case A Micromorphic','Classical Beam Model Timoshenko and Goodier','Classical Continuum','Micropolar Ramezani et al. Eur J Mech A-Solid-2009',...
    'Case C Micromorphic','Case B Micromorphic'},'Interpreter','Latex');
set(CC(4),'Color', [1, 0.5, 0], 'LineStyle', '--','Marker','o')
set(gca,'FontName','mwa_cmr10','FontSize',14)
set(lgnd, 'XColor', 'w', 'YColor', 'w', 'Color', 'none');



box on
Figure8 = figure(8);
hold on
plot(Mxy(:,1),Mxy(:,2),'--ob',Result_M113_com(:,1),Result_M113_com(:,2)-Result_M131_com(:,2),'--om'...
    ,Result_M113_str(:,1),Result_M113_str(:,2)-Result_M131_str(:,2),'--or'...
    ,Result_M113_rot(:,1),Result_M113_rot(:,2)-Result_M131_rot(:,2),'--og','LineWidth',2)
set(Figure8,'defaulttextinterpreter','latex');
ylabel('Couple stress','fontsize',18)
xlabel('$X_1$(m)','fontsize',18)
lgnd=legend({'Micropolar Ramezani et al. Eur J Mech A-Solid-2009','$\Phi_{11}$, $\Phi_{22}$, $\Phi_{33}$, $\Phi_{13}$, $\Phi_{31}$ Micromorphic',...
    '$\Phi_{11}$, $\Phi_{22}$, $\Phi_{33}$ Micromorphic','$\Phi_{13}$, $\Phi_{31}$ Micromorphic'},'Interpreter','Latex');
set(gca,'FontName','mwa_cmr10','FontSize',14)
set(lgnd, 'XColor', 'w', 'YColor', 'w', 'Color', 'none');


box on
Figure9 = figure(9);
hold on
plot(Result_dz_Q27_com(:,1),-Coeff*Result_dz_Q27_com(:,2),'--o',Result_dz_com(:,1),Coeff*Result_dz_com(:,2),'--or','LineWidth',2)
set(Figure9,'defaulttextinterpreter','latex');
ylabel('$3EIu_3/(FL^3)$','fontsize',18)
xlabel('$X_1$(m)','fontsize',18)
lgnd=legend({'Case A Micromorphic Q27P8','Case A Micromorphic Q8P8'},'Interpreter','Latex');
set(gca,'FontName','mwa_cmr10','FontSize',14)
set(lgnd, 'XColor', 'w', 'YColor', 'w', 'Color', 'none');


box on
Figure10 = figure(10);
box on
hold on
plot(Result_phi13_com(:,1),0.5*(Result_phi13_com(:,2)-Result_phi31_com(:,2)),'--or'...
    ,Result_phi13_Q27_com(:,1),-0.5*(Result_phi13_Q27_com(:,2)-Result_phi31_Q27_com(:,2)),'--o','LineWidth',2)
set(Figure10,'defaulttextinterpreter','latex');
ylabel('$\Phi_2^{rot} = \frac{1}{2}(\Phi_{13}-\Phi_{31})$','fontsize',18)
xlabel('$X_1$(m)','fontsize',18)
lgnd=legend({'Case A Micromorphic Q8P8'...
    ,'Case A Micromorphic Q27P8'},'Interpreter','Latex');
set(gca,'FontName','mwa_cmr10','FontSize',14)
set(lgnd, 'XColor', 'w', 'YColor', 'w', 'Color', 'none');
