set title 'Le potentiel exacte de O2'
set autoscale
set xlabel 'r'
set ylabel 'E(eV)'
m = "./exo_2_potentiel_exacte.dat"
set grid
plot m using 1:2 with linespoints