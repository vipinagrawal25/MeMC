se term x11 0
se log y
plot "< tail -54000 d_sph/mc_log | awk '{print $3 }' | gsl-histogram 1230 1430 64" u
(0.5*($1+$2)):3 w lp t "local from start", 'hist_start.dat' u 1:2 w p t "hist old"


se term x11 1

plot "< tail -44000 out/mc_log | awk '{print $3 }' | gsl-histogram 450 750 64" u
(0.5*($1+$2)):3 w lp t "local from memc", 'hist_memc.dat' u 1:2 w p t "hist old"

pause -1
