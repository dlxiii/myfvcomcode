./grrider << EOF
initial_zeta.dat
1
-0.5700E+05  0.5700E+05 -0.5700E+05  0.5700E+05
50 50
initial_zeta.grd
1
24
0
EOF
./grrider << EOF
analysis_day1_day2_zeta.dat
2
-0.5700E+05  0.5700E+05 -0.5700E+05  0.5700E+05
50 50
analysis_day1_day2_zeta.grd
1
24
0
EOF
./grrider << EOF
true_initial_day1_day2_zeta.dat
3
-0.5700E+05  0.5700E+05 -0.5700E+05  0.5700E+05
50 50
true_initial_day1_day2_zeta.grd
1
24
0
EOF
