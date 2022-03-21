echo "Create a surface for N=512"
../bin/exe_start 512 sph d_sph/
python ../utils/gen_memc_conf.py d_sph/part_pos0300.bin
sleep 2
../bin/exe_memc para_file.in out/
