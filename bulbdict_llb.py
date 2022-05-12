from struct import unpack
import bulbdef_llb
#from .bulbdef import gid_is_mitral as ismitral, gid_is_mtufted as ismtufted

class BulbDict:
  def __init__(self, filename):
    self.gid_dict = {}
    self.mgid_dict = {}
    self.mtgid_dict = {}
    self.ggid_dict = {}
    
    with open(filename, 'rb') as fi:
      rec = fi.read(22)
      while rec:
        
        rsgid, cgid, isec, xc, ggid, xg = unpack('>LLHfLf', rec)
        
        self.gid_dict.update({ rsgid:(cgid, isec, xc, ggid, 0, xg) })

        # by GC 
        if ggid not in self.ggid_dict:
          g2rs_set = set()
          self.ggid_dict[ggid] = g2rs_set
        else:
          g2rs_set = self.ggid_dict[ggid]
        g2rs_set.add(rsgid-1)
       

        # by M/T cells
        if bulbdef_llb.gid_is_mitral(cgid):
          cell_map = self.mgid_dict
        elif bulbdef_llb.gid_is_mtufted(cgid):
          cell_map = self.mtgid_dict


        if cgid not in cell_map:
          rs_set = set()
          g_set = set()
          cell_map[cgid] = (rs_set, g_set)
        else:
          rs_set, g_set = cell_map[cgid]

        rs_set.add(rsgid)
        g_set.add(ggid)
        
        
        rec = fi.read(22)

    print((list(self.mgid_dict.keys())))
    print((list(self.mtgid_dict.keys())))
    
  def __is_rs_of(self, gid, gidfunc):
    try:
      ci = self.gid_dict[gid]
    except KeyError:
      return False

    return gidfunc(ci[0])
  
  def is_mc_rs(self, gid):
    return self.__is_rs_of(gid, bulbdef.gid_is_mitral)

  def is_mtc_rs(self, gid):
    return self.__is_rs_of(gid, gid_is_mtufted)

  def is_gc_rs(self, gid):
    return gid % 2 == 0 and (gid+1) in self.gid_dict

  def mgid2ggid(self, gid):
    return self.mgid_dict[gid][1]
  
  def mtgid2ggid(self, gid):
    return self.mtgid_dict[gid][1]

  def mgid2mgrsid(self, mgid):
    return self.mgid_dict[mgid][0]

  def mtgid2mgrsid(self, mgid):
    return self.mtgid_dict[mgid][0]

  def ggid2mgrsid(self, ggid):
    return self.ggid_dict[ggid]
  
  def query(self, gid):
    if gid % 2 != 0:
      gid += 1
    return self.gid_dict[gid]


if __name__ == '__main__':
  dic = BulbDict('column-control.dic')
  #for gid in dic.gid_dict.keys():
  #  if gid % 2 == 0:
  #    print gid
  # return an information tuple
  # it contains (mitral id, secden index, secden section, granule id, 0, priden section)
  #print dic.query(dic.gid_dict.keys()[0])

  # gc soma gid connected to a mitral cell
  #print dic.gid2ggid(185)

  # all gid connected to a mitral cell
  #print dic.gid2mgrsid(185)
  
