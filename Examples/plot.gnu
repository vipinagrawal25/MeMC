se term x11 0
se log y
plot "< tail -54000 d_sph/mc_log | awk '{print $3 }' | gsl-histogram 1230 1430 64" u (0.5*($1+$2)):3 w lp t "", 'hist_start.dat' u 1:2 w p pt 6 ps 2 t ""


se term x11 1

plot "< tail -44000 out/mc_log | awk '{print $3 }' | gsl-histogram 450 750 64" u (0.5*($1+$2)):3 w lp t " ", 'hist_memc.dat' u 1:2 w p pt 6 ps 2 t ""
pause -1
