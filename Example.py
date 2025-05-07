import Cmash
import pickle


with open('/mnt/rdliang/AGN/BBH_AGN/data/O4_BBH_meta.pkl','rb') as handle:
    data = pickle.load(handle)
    handle.close()

for event_name in data.keys():
    gw_info = Cmash.get_event_info(data[event_name])
    search = Cmash.Search()
    maxoverlap, coverage, pnt = search.from_gw(gw_info,event_name,l23='bb')