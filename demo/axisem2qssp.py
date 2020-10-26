#!/usr/bin/env python
# -*- coding: utf-8 -*-


def qka2qp(qka, qmu, cp, cs):
    coefficient = 4./3.
    if cs == 0.0:
        return qka
    else:
        inverse_of_Qp = (1. - coefficient*(cs**2)/(cp**2))/qka + coefficient*(cs**2)/(cp**2)/qmu
        return 1./inverse_of_Qp


# stdata = []

fname = '1dmodel_axisem.bm'
fname2 = '1dmodel_axisem_qssp_noQ.txt'

f2 = open(fname2,"w+")

ii = 0
with open(fname,'r') as f:
    d = f.readlines()

    for i,line in enumerate(d):
        if i <= 5:
            header = line.split()
            print(header)
            continue
        else:
            tmp = line.split()
            if tmp[0] == '#':
                continue
            ii += 1
            sdata = list(map(float, tmp))
            # a_float_m = list(a_float_m)
            # if i % 1000 == 0:
            # AxiSEM: radius [m], rho [kg/m^3], vpv [m/s] vps [m/s] qka , qmu
            # QSSP: NO. depth[km] vp[km/s]  vs[km/s]   ro[g/cm^3] qp         qs
            # qp = qka2qp(sdata[4], sdata[5], sdata[2], sdata[3])
            # qs = qmu
            qp = 1.0e8
            if sdata[5] == 0.0:
                qs = 0.0
            else:
                qs = 1.0e8
            tmp = '{:>4} {:>9.4f} {:>10.5f} {:>10.5f} {:>10.5f} {:>10.5e} {:>10.5e}'.format(ii, 6371.0 - sdata[0]/1000.0, sdata[2]/1000.0, sdata[3]/1000.0, sdata[1]/1000.0, qp, qs)
            f2.write(tmp+'\n')
