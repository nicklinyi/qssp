from obspy import Stream, Trace, UTCDateTime, read

import numpy as np
import os, glob

header = None


fnames = glob.glob('*.dat')

for fname in fnames:
    stdata = []
    with open(fname,'r') as f:
        d = f.readlines()

        for i,line in enumerate(d):
            if i == 0:
                header = line.split()
                print(header)
                continue
            else:
                tmp = line.split()
                sdata = list(map(float, tmp))
                # a_float_m = list(a_float_m)
                stdata.append(sdata)
                # if i % 1000 == 0:
                #     print(sdata)
    print(len(stdata))
    print(len(stdata[0])) 

    array = np.array(stdata)
    print(array.shape)

    dirname = 'sac'

    if not os.path.exists(dirname):
        os.mkdir(dirname)

    suffix = 'BH'+fname.split('_')[2][0].upper()

    for i,sta_name in enumerate(header[1:len(header)]):
        
        trace = Trace()
        # trace.id = sta_name.replace('_','.')+'.'+suffix
        print('Processing '+trace.id)
        trace.data = array[:,i+1]
        trace.stats.sampling_rate = 100
        trace.stats.delta = 1.0 / trace.stats.sampling_rate
        trace.stats.network = sta_name.split('_')[0]
        trace.stats.station = sta_name.split('_')[1]
        trace.stats.location = '00' # TODO
        trace.stats.channel = suffix
        trace.id = trace.stats.network+'.'+trace.stats.station+'.'+trace.stats.location+'.'+trace.stats.channel
        trace.stats.starttime = UTCDateTime(array[0,0])
        trace.stats._format = 'SAC'
        trace.write(dirname+'/'+trace.id+'.SAC', format='SAC') 

        st = read(dirname+'/'+trace.id+'.SAC')
        st[0].stats.sac.b = array[0,0]
        st[0].stats.sac.e = array[len(stdata)-1,0]

        st.write(dirname+'/'+trace.id+'.SAC', format='SAC') 






