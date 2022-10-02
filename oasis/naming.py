#///////////////////////////////////////////////////////////////
#// naming.py
#// Copyright (C) 2018 Mehdi Najafi (mnuoft@gmail.com)
#//
#// Distribution of this library is not allowed in any form.
#///////////////////////////////////////////////////////////////

__author__ = "Mehdi Najafi <mnuoft@gmail.com>"
__date__ = "2018-04-21"
__copyright__ = "Copyright (C) 2018 " + __author__
__license__  = "GNU GPL version 3 or any later version" 

#///////////////////////////////////////////////////////////////

#from __future__ import print_function
import inspect, sys, os

_srcpath_ = inspect.getsourcefile(lambda:0)
_srcpath_ = _srcpath_[:_srcpath_.rfind('/')]

# Retrive the boundary information from the info file
def read_mesh_info(mesh_info_path, key):
    info = open(mesh_info_path, 'r').read()
    # Looking for the given key
    p1 = info.find(key)
    if p1<0:
        return [], [], []
    p1 += len(key)
    p2 = info.find('<', p1)
    if p2<0:
        buf = info[p1:]
    else:
        buf = info[p1:p2-1]
    lines = buf.split('\n')
    ids = []
    idfr = []
    ida = []
    idr = []
    fcs = []
    # Reaing at the key values
    for line in lines:
        ls = line.split()
        if len(ls) > 1:
            # id
            ids.append(int(ls[0]))# eval(ls[ 0]))
            # radius
            idr.append(eval(ls[-3]))
            # area
            ida.append(eval(ls[-2]))
            # flowrate or arearatio
            s = ls[-1].replace('A[','a[').replace('R[','r[')
            s = s.replace('A',ls[-2]).replace('R',ls[-3])
            idfr.append(s)
            # wave form
            fcs.append(ls[1])
    # evaluate all the expressions in the flowrates and outflow ratios
    for i,expr in enumerate(idfr):
        for j,k in enumerate(ids):
             expr = expr.replace( 'r[%d]'%k, str(idr[j])).replace( 'a[%d]'%k, str(ida[j]))
        idfr[i] = eval(expr)

    # sum of area ratio correction:
    if key == '<OUTLETS>':
        idfr[-1] = 1.0 - sum(idfr[:-1])
        
    return ids, idfr, ida, fcs

# check if the period is mentioned in the fc waveform file
def get_period_from_fcs(fcs):
    periods = [951.0 for f in fcs]
    for i,f in enumerate(fcs):
        fcs_i_filename = f.split(':')[-1]
        if os.path.exists( os.path.join('./data', fcs_i_filename) ):
            fcs_ifname = os.path.join('./data', fcs_i_filename)
        else:
            # print (os.path.dirname(os.path.abspath(__file__))); sys.stdout.flush()
            # srcpath = os.path.dirname(os.path.abspath(__file__))
            fcs_ifname = os.path.join(_srcpath_, 'data', fcs_i_filename)
            if not os.path.exists( fcs_ifname ):
                print ('<!> Cannot find the waveform:', fcs_i_filename)
        for line in open(fcs_ifname,'r').readlines():
            if len(line.strip()) == 0: continue
            if line.strip()[0] in ['#','!','/']:
                p = line.find('period_ms')
                if p < 0: p = line.find('period')
                if p > 0:
                    periods[i] = float(''.join((ch if ch in '0123456789.-e' else ' ') for ch in line[p+9:]).strip().split(' ')[0])
    return periods

def get_period_from_info(mesh_name):
    id_in, Q_means, _, fcs = read_mesh_info('./data/'+mesh_name+'.info', '<INLETS>')
    periods = get_period_from_fcs(fcs)
    if len(periods) > 1:
        period = min(periods)
    else:
        period = periods[0]
    return period

def get_case_fullname(mesh_name, uo, cy, ts):
    id_in, Q_means, _, fcs = read_mesh_info('./data/'+mesh_name+'.info', '<INLETS>')
    periods = get_period_from_fcs(fcs)
    if len(periods) > 1:
        period = min(periods)
    else:
        period = periods[0]
    txt = ''
    for i, id in enumerate(id_in): txt += '_I%d_%s_Q%d'%(id,fcs[i].replace(':','_'),int(Q_means[i]*100))
    txt += '_Per%d'%int(period) + "_Newt370" + "_ts" + str(ts) + "_cy" + str(cy) + "_uO" + str(uo)
    return "art_" + mesh_name + txt

if __name__ == '__main__':
    print ( get_case_fullname(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4]) )
