#!/bin/bash -l

#SBATCH --job-name="olfactory-bulb"
#SBATCH --partition=normal
#SBATCH --constraint=mc
#SBATCH --mail-type=ALL
#SBATCH --account=ich002

module swap PrgEnv-cray PrgEnv-intel
module load daint-mc cray-python/3.9.4.1 
module use /apps/hbp/ich002/hbp-spack-deployments/softwares/23-02-2022/modules/tcl/cray-cnl7-haswell
module load neuron/8.0.2

export PATH=/store/hbp/ich035/kumbhar/neuron_home/install3/bin:$PATH
export PYTHONPATH=/store/hbp/ich035/kumbhar/neuron_home/install3/lib/python/:$PYTHONPATH
export LANG=C
export LC_ALL=en_US


mkdir -p ./ob_sim_all
cp -r /apps/hbp/ich002/cnr-software-utils/olfactory-bulb/olfactory-bulb-3d/ ./ob_sim_all
cd ./ob_sim_all/olfactory-bulb-3d/

git checkout 2to3

cd sim

rm -rf x86x64

nrnivmodl .

python3 /apps/hbp/ich002/cnr-software-utils/olfactory-bulb/olfactory-bulb-utils/create_sim_launcher.py --glomlist $1 --outputfolder $2 --odor $3 --sniffintvl $4

srun ./x86_64/special -python -mpi bulbsimlauncher.py

cp /apps/hbp/ich002/cnr-software-utils/olfactory-bulb/olfactory-bulb-utils/sim_dict_to_json.py .
cp /apps/hbp/ich002/cnr-software-utils/olfactory-bulb/olfactory-bulb-utils/bulbdef_llb.py .
cp /apps/hbp/ich002/cnr-software-utils/olfactory-bulb/olfactory-bulb-utils/bulbdict_llb.py .
cp /apps/hbp/ich002/cnr-software-utils/olfactory-bulb/olfactory-bulb-utils/params_llb.py .

python3 sim_dict_to_json.py --outputfolder sim_out --inputfolder .

cp sim_out/simgloms.json ../../..
cp sim_out/simcells.json ../../..
cp sim_out/connections.json ../../..
cd ../../..

zip -r ob_sim_all.zip ob_sim_all

