import os
import json
import argparse

def main():
    parser = argparse.ArgumentParser(description = '''The create_sim_launcher.py script  ''' \
        '''create the launch file for running the olfactory bulb simulation.''')
    parser.add_argument("--glomlist", type=str, required=False, default="[5, 7, 32, 37, 78]", 
        help="list of the glomeruli to be simulated -> e.g.: [5, 7, 32] - default: '[5, 7, 32, 37, 78]'")
    parser.add_argument("--outputfolder", type=str, required=False, default="./output", 
        help="path to the folder where the file will be created (if not present, it will be created) -> e.g.: ./my_results - default: './'" )
    parser.add_argument("--odor", type=str, required=False, default="Onion", 
        help="odor to be supplied as stimulus -> e.g.: Apple - default: 'Onion'" )
    parser.add_argument("--sniffintvl", type=str, required=False, default="500", 
        help="odor to be supplied as stimulus -> e.g.: 600 - default: '500'" )

    args = parser.parse_args()

    glomlist = args.glomlist
    outputfolder = args.outputfolder
    odor = args.odor
    sniffintvl = args.sniffintvl

    if not os.path.exists(outputfolder):
        os.makedirs(outputfolder)

    with open("bulbsimlauncher.py", "w") as f:
        f.write("import params\n")
        f.write("import runsim\n")
        f.write("from math import sqrt\n")
        f.write("from neuron import h\n")
        f.write("import odors\n\n")
        f.write("params.sniff_invl_min = params.sniff_invl_max = " + sniffintvl + "\n")
        f.write("params.training_exc = params.training_inh = True\n")
        f.write("h('sigslope_AmpaNmda=5')\n")
        f.write("h('sigslope_FastInhib=5')\n")
        f.write("h('sigexp_AmpaNmda=4')\n\n")
        f.write("params.odor_sequence = [ ('" + odor + "', 50, 1000, 1e+9) ]\n")
        f.write("runsim.build_part_model(" + glomlist + ", [])\n")
        f.write("runsim.run()\n")

if __name__ == "__main__":
    main()