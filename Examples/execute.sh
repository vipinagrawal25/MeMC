echo "Create a surface for N=1024"
../bin/exe_start 1024 sph d_sph/
python ../utils/gen_memc_conf.py d_sph/part_pos0003.bin
sleep 2
../bin/exe_memc para_file.in out/
