import Cmash
import pickle
import numpy as np
from astropy.table import Table


with open('/mnt/rdliang/AGN/BBH_AGN/data/O4_BBH_meta.pkl','rb') as handle:
    data = pickle.load(handle)
    handle.close()

#Calculate spatial and temporal coverage and upper limits
data_output = []
for event_name in list(data.keys()):
    if event_name[1:3] == '23':
        continue
    else:
        gw_info = Cmash.get_event_info(data[event_name])
        search = Cmash.Search()
        maxoverlap, coverage, pnt = search.from_gw(gw_info,event_name)
        if maxoverlap > 0.:
            uplim = float(np.mean(pnt['uplim']))
            data_output.append([event_name,maxoverlap,coverage,uplim])

print(data_output)
data_table = Table(np.array(data_output),names=['Name','t_overlap','coverage','uplim'],dtype=['U20']+['f8']*3)
data_table.write('/mnt/rdliang/AGN/BBH_AGN/data/O4_EP.csv',format='csv',overwrite=True)