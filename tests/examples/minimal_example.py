# -*- coding: utf-8 -*-
"""
for use with test_cli.py, requires files from online source

"""

import hatyan
hatyan.close('all')
import urllib

dir_testdata = "https://raw.githubusercontent.com/Deltares/hatyan/main/tests/data_unitsystemtests/"

current_station = 'VLISSGN'

# retrieve files from online
file_list = [f'{current_station}_ana.txt',f'{current_station}_ext.txt']
for fname in file_list:
    url = dir_testdata+fname
    with urllib.request.urlopen(url) as response:
        data = response.read()
    with open(fname, "w") as f:
        f.write(data.decode('utf-8'))

file_data_comp1 = f'{current_station}_ana.txt'

times_pred = slice("2019-01-01","2019-02-01", "10min")

comp_fromfile = hatyan.read_components(filename=file_data_comp1)

#prediction and validation
ts_prediction = hatyan.prediction(comp=comp_fromfile, times=times_pred)
hatyan.write_dia(ts=ts_prediction, filename='prediction_%s_%s.dia'%(times_pred.step,current_station))
ts_ext_prediction = hatyan.calc_HWLW(ts=ts_prediction)
hatyan.write_dia(ts=ts_ext_prediction, filename='prediction_%s_%s_HWLW.dia'%(times_pred.step,current_station))
fig, (ax1,ax2) = hatyan.plot_timeseries(ts=ts_prediction, ts_ext=ts_ext_prediction)
fig.savefig('prediction_%s_%s'%(times_pred.step, current_station))
