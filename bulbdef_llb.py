# Olfactory Bulb space and boundary surfaces
import params_llb
#all_dict = {
#    "bulb_center" : (50.0, 1275.0, 0.0),
#    "bulb_axis"   : (2100.0, 2800.0, 2100.0),
#    "glom_axis"   : (2400.0, 3100.0, 2400.0),
#    "gran_bnd_up" : (1400.0, 2100.0, 1400.0),
#    "gran_bnd_dw" : (600.0, 1300.0, 600.0),
#    "mid_tufted_dw" : (1875.0, 2500.0, 1875.0),
#    "glom_radius" : 50.0,
#    "gran_voxel" : 18,
#    "gran_radius" : 5.0,
#}


# granules parameters
# This value changes the number of granules that will be generated
# the GCL volume is approximately 1910088333.383 (um^3)
# so the number of granules filling the GCL is volume/(grid_dim^3)
#gran_voxel  = 18
#gran_radius = 5.0

#gran_connect_radius = 60.0
#gran_type_slope = 10.0

# load glomeruli positions
#glom_coord = []
#with open('./realgloms.txt', 'r') as f:
#  l = f.readline()
#  while l:
#    tk = l.split()
#    p = ()
#    for _tk in tk:
#      p += (float(_tk), )
#    glom_coord.append(p)
#    l = f.readline()


#glom_coord = glom_coord[:127]
# cell numbers    
#Ngloms = len(glom_coord) # glomeruli

# mitral
#gid_mitral_begin = 0
#Nmitral_per_glom = 5 # mitral per glomerolus
#Nmitral = Ngloms * Nmitral_per_glom

# middle tufted
#gid_mtufted_begin = gid_mitral_begin+Nmitral
#Nmtufted_per_glom = 2*Nmitral_per_glom # twice than mitral!
#Nmtufted = Ngloms * Nmtufted_per_glom

#aid_granule_begin = gid_mtufted_begin + Nmtufted
import granules
granules.init(params.bulbCenter, params.granAxisUp, params.granAxisDw, params.gran_voxel, params.gid_granule_begin)
#granules.init()
#Ngranule = len(granules.ggid2pos)

#gid_blanes_begin = gid_granule_begin+Ngranule
#NBlanes_per_glom = 3
#NBlanes = Ngloms*NBlanes_per_glom
#gid_mbs_begin = gid_blanes_begin+NBlanes
#gid_bc2gc_begin = gid_mbs_begin+Nmtufted*NBlanes

# reciprocal synapse
#gid_mgrs_begin = gid_bc2gc_begin+NBlanes*Ngranule
#if gid_mgrs_begin % 2 != 0:
#  gid_mgrs_begin += 1

def gid_is_mitral(gid):
  return gid >= params.gid_mitral_begin and gid < params.gid_mitral_begin+params.Nmitral

def gid_is_mtufted(gid):
  return gid >= params.gid_mtufted_begin and gid < params.gid_mtufted_begin+params.Nmtufted

def gid_is_granule(gid):
  return gid >= params.gid_granule_begin and gid < params.gid_granule_begin+params.Ngranule

def gid_is_blanes(gid):
  return gid >= params.gid_blanes_begin and gid < params.gid_blanes_begin+params.NBlanes

''' useful only for mitral and tufted cells'''
def cellid2glomid(cgid):
  if cgid < params.Nmitral:
    return cgid/params.Nmitral_per_glom
  elif cgid < params.gid_mtufted_begin+params.Nmtufted:
    return (cgid-params.gid_mtufted_begin)/params.Nmtufted_per_glom
  return None

''' reciprocal synapses id to cells '''
def rs2cell_gid(gid):
  gid -= params.gid_mgrs_begin
  if gid % 2:
    gid -= 1
    
    example.py
  gidtarget = int(gid / (2*params.Ngloms))
  gidsrc = gid % (2*params.Ngloms)
  return gidsrc, gidtarget


''' cell id to rs '''
def cell2rs_gid(gidsrc, gidtarget):
  if gidsrc > gidtarget: # ggid and m/t gid inverted
    aux = gidsrc
    gidsrc = gidtarget
    gidtarget = aux

    print('auu')

    slot = -1
  else:
    slot = 0
    
  return gid_mgrs_begin + slot + 2 * (gidtarget * Ngloms + cellid2glomid(gidsrc))
