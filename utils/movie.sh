#!/bin/bash
#dir=/home/vikash/opt/mnt2/dardel/shear/solid-fluid2/fluid_nores/00009/viz_all_
dir=/home/vikash/Code/MeMC/utils/avrgSnaps/fluid_nores_wse_
odir="avrgSnaps/fluid_wse_begining_"
for i in {1..399}
do 
    # inf = printf("%s%04d", dir, i)
    printf -v inf '%s%05d.vtk' "$dir" "$i"
    printf -v outf '%s%04d' "$odir" "$i"
    /opt/visit/bin/visit -nowin -cli -s pseudocolor.py $inf $outf
    # magick mogrify -fuzz 4% -define trim:percent-background=100% -trim +repage -format png $outf.png
    # echo $inf
done
