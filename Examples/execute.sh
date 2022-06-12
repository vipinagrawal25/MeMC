echo "Create a surface for N=512"
../bin/exe_start 512 sph d_sph/ 60000
python ../utils/gen_memc_conf.py d_sph/snap_0300.h5
sleep 2
../bin/exe_memc para_file.in out/
