import pprint
import argparse
import os
import json 
import shutil
from collections import OrderedDict

glom_cells_dict = {}

def map_glom_to_cells():
    
    glom_id = list(range(0, 127))
    mitr_x_glom = 5
    tuft_x_glom = 10
    num_mitr = 635

    for i in glom_id:
        glom_cells_dict[str(i)] = []
        start_mitr = i * mitr_x_glom
        end_mitr = i * mitr_x_glom + mitr_x_glom
        start_tft = num_mitr + i * tuft_x_glom
        end_tft = num_mitr + i * tuft_x_glom + tuft_x_glom
        mitr = list(range(start_mitr, end_mitr))
        tuft = list(range(start_tft, end_tft))
        for m in mitr:
            glom_cells_dict[str(i)].append(m)
        for t in tuft:
            glom_cells_dict[str(i)].append(t)


def generateGlomPositions(gloms):
    '''
    Reference:
    {
        geom        : {diam:1,nseg:1}
        glom 1      : {position:{x:_,y:_,z:_}}
        glom 2      : {position:{x:_,y:_,z:_}}
        glom 3      : {position:{x:_,y:_,z:_}}
        ...
        glom 127    : {position:{x:_,y:_,z:_}}
    }


    {
        "geom": {"diam": 100,"nseg": 1},
        "0"         : {"pos": {"x": 207.283,"y": 221.71,"z": 866.201}},
        "1"         : {"pos": {"x": 563.378,"y": 1699.8,"z": 1033.58}},
        "2"         : {"pos": {"x": 759.676,"y": 1839.35,"z": 863.402}},
        "3"         : {"pos": {"x": 1002.92,"y": 1799.8,"z": 605.69}},

    '''
    # load glomeruli positions (extracted from bulbdef.py)
    glom_radius=50
    glom_coord=[]
    with open('realgloms.txt', 'r') as f:
        l = f.readline()
        while l:
            tk = l.split()
            p = ()
            for _tk in tk:
                p += (float(_tk), )
            glom_coord.append(p)
            l = f.readline()


    # creating sections with equivalent names in the NetPyNE format to store its properties
    netpyne_cell={}
    netpyne_cell =  {'geom':{'diam':2*glom_radius,'nseg':1}}
    for i in gloms:
        print('generating Glomeruli: ',i, '\tw/ coord-> x:',glom_coord[i][0],'y:',glom_coord[i][1],'z:',glom_coord[i][2],)
        glom_dict =  {   i:  
                                {   'pos':  {   'x':glom_coord[i][0],
                                                'y':glom_coord[i][1],
                                                'z':glom_coord[i][2],
                                            },
                                }
                        }
        netpyne_cell.update(glom_dict)

    with open('./netpyne_cells/netpyne_glomeruli.json', 'w', encoding='utf-8') as f:
        # json.dump(netpyne_cell, f, ensure_ascii=False)
        json.dump(netpyne_cell, f, ensure_ascii=False, indent=4)

def generateMitraCells(gloms,mitrals):
    '''
    Reference:
    {
    "secs": {
        "dend_0": {"geom": {"nseg": 4,"pt3d": [[x0,y0,z0,d0],[x1,y1,z1,d1],[x2,y2,z2,d2],[x3,y3,z3,d3]]}},
        "dend_1": {"geom": {"nseg": 4,"pt3d": [[x0,y0,z0,d0],[x1,y1,z1,d1],[x2,y2,z2,d2],[x3,y3,z3,d3]]}},
        ...
        "soma_0": {"geom": {"nseg": 4,"pt3d": [[x0,y0,z0,d0],[x1,y1,z1,d1],[x2,y2,z2,d2],[x3,y3,z3,d3]]}},
        "apic_1": {"geom": {"nseg": 4,"pt3d": [[x0,y0,z0,d0],[x1,y1,z1,d1],[x2,y2,z2,d2],[x3,y3,z3,d3]]}},
        "tuft_1": {"geom": {"nseg": 4,"pt3d": [[x0,y0,z0,d0],[x1,y1,z1,d1],[x2,y2,z2,d2],[x3,y3,z3,d3]]}},
        ...
    '''
    # parameters extracted from 'params.py'
    # mitral GIDs
    gid_mitral_begin = 0
    Nmitral_per_glom = 5 # mitral per glomerolus
    Nmitral = len(gloms) * Nmitral_per_glom

    #[jv comment: creating GIDs for the middle tufted cells]
    # middle tufted
    gid_mtufted_begin = gid_mitral_begin+Nmitral
    Nmtufted_per_glom = 2*Nmitral_per_glom # twice than mitral!
    Nmtufted = len(gloms) * Nmtufted_per_glom

    # gids
    gids = set()
    for glomid in gloms:
        gids.update(list(range(glomid * Nmitral_per_glom, (glomid+1) * Nmitral_per_glom)) + \
                    list(range(glomid * Nmtufted_per_glom + gid_mtufted_begin, (glomid+1) * Nmtufted_per_glom + gid_mtufted_begin)))
    gids.update(mitrals)

    gid_list=list(gids)
    gid_list.sort()
    print('sim gids: ',gid_list)

    import cellreader
    for cell_gid in gid_list:
        # print('Cell GID: \t',cell_gid)
        cr=cellreader.CellReader('mitral_cell_library.car') # Joao library - same as mccells.car, but 
        # cr=cellreader.CellReader('mccells.car') # default mitral cells library
        cell = cr.readcell(cell_gid)
        
        # Adds a cell.name property in the cell object to facilitate NetPyNE conversion
        for sec_type in cell.__dict__.keys():
            for sec_num in range(len(cell.__dict__[sec_type])):
                cell.__dict__[sec_type][sec_num].name=sec_type+'_'+str(sec_num)

        hist = [0]*int(1040/40)
        import misc
        center = misc.centroid(cell.soma[0].points)
        print('\nCell GID: \t',cell_gid,'\t N of comp | \t soma:', len(cell.soma),'\t dend:', len(cell.dend),'\t apic:', len(cell.apic),'\t tuft:', len(cell.tuft))
        
        printSecInfo=False
        if printSecInfo:
            for ind_dnd,dnd in enumerate(cell.dend):
                # cell.dend[ind_dnd].name='dend_'+str(ind_dnd)
                for ind_p,p in enumerate(dnd.points):
                    d = misc.distance(p, center)
                    # hist[int(d/40)] += 1.
                    print('dend ',ind_dnd,'\t point ',ind_p,'\t pos (x,y,z):', p,'\t soma dist:',d)
            
            for ind_api,api in enumerate(cell.apic):
                for ind_p,p in enumerate(api.points):
                    d = misc.distance(p, center)
                    print('apic ',ind_api,'\t comp ',ind_p,'\t pos (x,y,z):', p,'\t soma dist:',d)
            
            for ind_tuf,tuf in enumerate(cell.tuft):
                for ind_p,p in enumerate(tuf.points):
                    d = misc.distance(p, center)
                    print('tuft ',ind_tuf,'\t comp ',ind_p,'\t pos (x,y,z):', p,'\t soma dist:',d)

        # initializing a dictionary in NetPyNE format to store the cell properties
        netpyne_cell = {'secs': {},'conds':{}}
        # retrieving the compartments in dictionary format
        for sec_type in cell.__dict__.keys():
            # retrieving each compartment for a specific segment type (e.g. soma:(soma_0, soma_1), dend)
            for sec_num in range(len(cell.__dict__[sec_type])):
                # creating a name attribute for each compartment to make it easier to link with NetPyNE formating (e.g. cell.<segment>[index].name)
                sec_name = sec_type+'_'+str(sec_num)
                cell.__dict__[sec_type][sec_num].name=sec_name


                pt3d_points=[]
                [pt3d_points.append([p1,p2,p3,p4]) for p1,p2,p3,p4 in cell.__dict__[sec_type][sec_num].points]

                # # creating sections with equivalent names in the NetPyNE format to store its properties
                netpyne_cell['secs'][sec_name] =    { 
                                                    'geom': {
                                                            # 'L':1,
                                                            # 'diam': 1,
                                                            'nseg': int(len(pt3d_points)-1), 
                                                            'pt3d':pt3d_points
                                                            },
                                                    }
        import json
        with open('./netpyne_cells/netpyne_mitral_cell_'+str(cell_gid)+'.json', 'w', encoding='utf-8') as f:
            json.dump(netpyne_cell, f, ensure_ascii=False, indent=4)

def generateGranulePositions(gloms):
    '''
    Reference:
    {
        geom        : {diam:1,nseg:1}
        granule gid : {position:{x:_,y:_,z:_}}
        granule gid : {position:{x:_,y:_,z:_}}
        granule gid : {position:{x:_,y:_,z:_}}
        ...
        granule gid : {position:{x:_,y:_,z:_}}
    }

    {
        "geom": {"diam": 1, "nseg": 1}, 
        "1905": {"pos": {"x": -646.0, "y": 1173.0, "z": -17.0}}, 
        "1906": {"pos": {"x": -646.0, "y": 1173.0, "z": 0.0}},
        ...
    }

    '''
    ggid2pos={}
    pos2ggid={}
    # load glomeruli positions (extracted from bulbdef.py)

    with open('granules.txt', 'r') as fi:
        line = fi.readline()
        while line:
            token = line.split()
            # joao modification - the code was returning an error because it was missing a condition to exit the loop when it finishes reading
            if len(line)<2:
                break
            # joao modification-end
            gid = int(token[0])
            pos = (float(token[1]), float(token[2]), float(token[3]))
            ggid2pos.update({ gid:pos })
            pos2ggid.update({ pos:gid })
            line = fi.readline()
        
        netpyne_cell={}
        netpyne_cell =  {'geom':{'diam':1,'nseg':1}}
        for k in ggid2pos.keys():
            # # creating sections with equivalent names in the NetPyNE format to store its properties
            print('generating granule cell: ',k, '\tw/ coord-> x:',ggid2pos[k][0],'y:',ggid2pos[k][1],'z:',ggid2pos[k][2],)
            granule_dict =  {   k:  
                                    {   'pos':  {   'x':ggid2pos[k][0],
                                                    'y':ggid2pos[k][1],
                                                    'z':ggid2pos[k][2],
                                                },
                                    }
                            }
            netpyne_cell.update(granule_dict)

        import json
        with open('./netpyne_cells/netpyne_granule_cells.json', 'w', encoding='utf-8') as f:
            json.dump(netpyne_cell, f, ensure_ascii=False)
            # json.dump(netpyne_cell, f, ensure_ascii=False, indent=4)

def generateMitralGranuleSynapses(dictionary_fileName):
    from bulbdict import BulbDict
    dic = BulbDict(dictionary_fileName)   
    printDicFile=False
    if printDicFile:[print(syn,dic.gid_dict[syn],'\n\n') for syn in dic.gid_dict.keys()]
    '''
    - GID reference for the synapses

    The model creates a GID for each synapse, to work around the fact that 2 compartments share more that 2 synapses

    OBS: THIS FILE ONLY GENERATE SYNAPSES IF THE GLOMERULI HAS BEEN SIMULATED IN THE MODEL, ONCE IT USES THE OUTPUT DATA IN '.dic' FILE

    gid reference:

    {
        synapse_gid:    (m/mt gid,  section,    weight(?)m/mt->g,       granule gid,  0,    weight(?)g->m/mt  )
    }

    e.g.:
    {   
        ...
        198688510:      (959,       35,         0.9791666865348816,     195876,         0,  0.8675885796546936),
        194165278:      (959,       36,         0.2500004470348358,     178068,         0,  0.8759732842445374),
        192622990:      (959,       36,         0.7500013709068298,     171996,         0,  0.9216063022613525),
        195622222:      (959,       37,         0.24999994039535522,    183804,         0,  0.7420779466629028), 
        192684458:      (959,       37,         0.7499998211860657,     172238,         0,  0.7902942895889282),
        ...
    }

    '''
    netpyne_cell={}
    netpyne_cell =  {}
    total_syns=len(dic.gid_dict.keys())
    for ind,syn_gid in enumerate(dic.gid_dict.keys()):
        # # creating sections with equivalent names in the NetPyNE format to store its properties
        print('Converting synapses | syn #: ',ind,'\t synapses remaining: ',total_syns-ind ,'\t|\tsynapse gid: ',syn_gid,)# '\tw/ coord-> x:',dic.gid_dict[k][0],'y:',dic.gid_dict[k][1],'z:',dic.gid_dict[k][2],)
        granule_dict =  {syn_gid:dic.gid_dict[syn_gid]}
        netpyne_cell.update(granule_dict)
    import json
    with open('./netpyne_synapses.json', 'w', encoding='utf-8') as f:
        json.dump(netpyne_cell, f, ensure_ascii=False)
        # json.dump(netpyne_cell, f, ensure_ascii=False, indent=4)


def generateMitralGranuleSynapses2(inputfolder, outputfolder, weightfn, dicfn, connfn, simcellfn, simglomfn):
    shutil.copyfile(os.path.join(inputfolder, "granules.py"), os.path.join(outputfolder, "granules.py"))
    shutil.copyfile(os.path.join(inputfolder, "params.py"), os.path.join(outputfolder, "params.py"))
    shutil.copyfile(os.path.join(inputfolder, "args.py"), os.path.join(outputfolder, "args.py"))
    shutil.copyfile(os.path.join(inputfolder, "realgloms.txt"), os.path.join(outputfolder, "realgloms.txt"))
    shutil.copyfile(os.path.join(inputfolder, "misc.py"), os.path.join(outputfolder, "misc.py"))

    from bulbdict_llb import BulbDict

    dic = BulbDict(dicfn)   
    
    printDicFile = False

    if printDicFile:
        [print(syn, dic.gid_dict[syn], '\n\n') for syn in dic.gid_dict.keys()]
    
    '''
    - GID reference for the synapses

    The model creates a GID for each synapse, to work around the fact that 2 compartments share more that 2 synapses

    OBS: THIS FILE ONLY GENERATE SYNAPSES IF THE GLOMERULI HAS BEEN SIMULATED IN THE MODEL, ONCE IT USES THE OUTPUT DATA IN '.dic' FILE

    gid reference:

    {
        synapse_gid:    (m/mt gid,  section,    weight(?)m/mt->g,       granule gid,  0,    weight(?)g->m/mt  )
    }

    e.g.:
    {   
        ...
        198688510:      (959,       35,         0.9791666865348816,     195876,         0,  0.8675885796546936),
        194165278:      (959,       36,         0.2500004470348358,     178068,         0,  0.8759732842445374),
        192622990:      (959,       36,         0.7500013709068298,     171996,         0,  0.9216063022613525),
        195622222:      (959,       37,         0.24999994039535522,    183804,         0,  0.7420779466629028), 
        192684458:      (959,       37,         0.7499998211860657,     172238,         0,  0.7902942895889282),
        ...
    }

    '''

    with open(weightfn, "r") as f:
        weights = json.load(f)

    simulated_cells = []
    netpyne_cell = {}

    all_connections = {}

    total_syns = len(dic.gid_dict.keys())

    for ind, syn_gid in enumerate(dic.gid_dict.keys()):
        # creating sections with equivalent names in the NetPyNE format to store its properties
        # print('Converting synapses | syn #: ', ind, '\t synapses remaining: ', total_syns-ind, '\t|\tsynapse gid: ', syn_gid,)# '\tw/ coord-> x:',dic.gid_dict[k][0],'y:',dic.gid_dict[k][1],'z:',dic.gid_dict[k][2],)
        granule_dict = {syn_gid: dic.gid_dict[syn_gid]}
        mt_cell = granule_dict[syn_gid][0]
        sec = granule_dict[syn_gid][1]
        gr_cell = granule_dict[syn_gid][3]
        if mt_cell not in all_connections:
            all_connections[mt_cell] = {}
        if sec not in all_connections[mt_cell]:
            all_connections[mt_cell][sec] = {}
        if gr_cell not in all_connections[mt_cell][sec]:
            all_connections[mt_cell][sec][gr_cell] = [0, 0]
        syn_gid_str = str(syn_gid)
        syn_gid_str_plus = str(syn_gid - 1)
        if syn_gid_str in weights.keys():
            all_connections[mt_cell][sec][gr_cell][0] = float(weights[syn_gid_str])
        if syn_gid_str_plus in weights.keys():
            all_connections[mt_cell][sec][gr_cell][1] = float(weights[syn_gid_str_plus])
        if mt_cell not in simulated_cells:
            simulated_cells.append(mt_cell)
    simulated_cells_sorted = sorted(simulated_cells)

    sim_gloms = []

    for i in simulated_cells_sorted:
        for k in glom_cells_dict:
            if i in glom_cells_dict[k] and int(k) not in sim_gloms:
                sim_gloms.append(int(k))
    sim_gloms_sorted = sorted(sim_gloms)

    with open(simglomfn, "w") as sgf:
        json.dump({"sim_gloms": sim_gloms_sorted}, sgf)

    with open(simcellfn, 'w') as sf:
        json.dump({"sim_cells": simulated_cells_sorted}, sf)

    with open(connfn, 'w', encoding='utf-8') as f:
        #json.dump(all_connections, f, ensure_ascii=False, indent=4)
        json.dump(all_connections, f, ensure_ascii=False)


def cat_weights(inputfolder, weightfilename):
    '''
    Concatenate the files containing the weight generated during the simulation
    into a single file named weightfilename
    '''

    # initialized final weight dictionary
    w_dict = {}

    # read individual weights files and insert them into a dictionary
    listdir = os.listdir(inputfolder)
    print(inputfolder)

    for i in listdir:
        if "olfactory_bulb.weight.dat." in i:
            with open(os.path.join(inputfolder, i), "r") as f:
                all_lines = f.readlines()
                for l in all_lines:
                    [gid, weight, cond] = l.split(" ")
                    if gid not in w_dict:
                        w_dict[gid] = float(weight)
    w_dict_sorted = OrderedDict(sorted(w_dict.items()))

    with open(weightfilename, "w") as dd:
        json.dump(w_dict_sorted, dd)

    return 

def cat_dic(inputfolder, dicfilename):

    # read individual .dic files, concatenate them and save the result file to storage
    listdir = os.listdir(inputfolder)
    fext = open(dicfilename, "wb")

    for f in listdir:
        if "olfactory_bulb.dic." in f:
            fo = open(os.path.join(inputfolder, f), "rb")
            shutil.copyfileobj(fo, fext)
            fo.close()
    fext.close()


def main():
    parser = argparse.ArgumentParser(description = '''The sim_dict_to_json.py script  ''' \
        '''convert the dictionary printed as output of the olfactory bulb simulator ''' \
        '''into a .json file containing information of the simulated mitral cells,  ''' \
        '''the relative granule cells and the synapses between these two cell types.''')
    parser.add_argument("--inputfolder", type=str, required=False, default="./", 
        help="path to the folder containing the weight.dat and dictionary files -> e.g.: ./mysim/ - default: './   '")
    parser.add_argument("--outputfolder", type=str, required=False, default="./output", 
        help="path to the output folder (if not present, it will be created) -> e.g.: ./my_results - default: './output'" )
    args = parser.parse_args()

    inputfolder = args.inputfolder
    outputfolder = args.outputfolder

    if not os.path.exists(outputfolder):
        os.makedirs(outputfolder)

    connectionsfilename = os.path.join(outputfolder, "connections.json")
    weightfilename = os.path.join(outputfolder, "weights.json")
    dicfilename = os.path.join(outputfolder, "bulb.dic")
    simcellfilename = os.path.join(outputfolder, "simcells.json")
    simglomfilename = os.path.join(outputfolder, "simgloms.json")

    map_glom_to_cells()
    cat_weights(inputfolder, weightfilename)
    cat_dic(inputfolder, dicfilename)
    generateMitralGranuleSynapses2(inputfolder, outputfolder, weightfilename, 
                                   dicfilename, connectionsfilename, simcellfilename,
                                   simglomfilename)
        
    #gloms=list(range(127))   # generate all mitral cells for a subset of glomeruli
    # gloms=[5,7,30]   # generate all mitral cells for a subset of glomeruli
    #mitrals=[]          # list to generate individual mitral cells (only in case it is not already created in the list of glomeruli)
    # Generates a JSON file with the positions of the glomeruli
    #generateGlomPositions(gloms)
    #generateMitraCells(gloms,mitrals)
    #generateGranulePositions(gloms)
    #dictionary_fileName='olfactory_bulb.dic.00'
    #generateMitralGranuleSynapses(dictionary_fileName)
    

if __name__ == "__main__":
    main()