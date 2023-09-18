import f90nml
import os
import shutil as sh

nml = f90nml.read("infiles/case2.nml")

z_c2 = [0, 5, 10, 15, 20, 25, 30]

z_c2 = [int(i) for i in range(0, 101, 5)]

for z in z_c2:
    nml["nml_px"]["z_injection"][3] = 100 - z
    nml["nml_px"]["z_injection"][4] = z
    nml.write("tmp_infile.nml", force=True)
    os.system("fpm run --profile release -- --infile tmp_infile.nml")
    sh.copytree("output", f"run_outputs/output_{z:03d}")