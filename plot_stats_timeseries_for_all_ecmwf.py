import sys,os
import netCDF4
import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime, timedelta
import xarray as xr
from xskillscore import rmse
#from xskillscore import me
import xarray.ufuncs as xu
import cftime
import pandas as pd
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from functools import reduce

import subprocess, atexit
from subprocess import Popen, PIPE
import io   

#from mpl_toolkits.basemap import Basemap

def main(argc, argv):

   if not (argc == 4):
        print ('Usage: %s <fcst hour> <var to plot (rmse)> <level>' % argv[0])
        sys.exit(1)

   fcsh=int(sys.argv[1])
   field=sys.argv[2]
   L=float(sys.argv[3])

   if (field == 'rh'):
       variable_label = "Relative Humidity (%)"
   if (field == 't'):
       variable_label = "Temperature (C)"

   # Fixed parameters
   
   exp = ["m-oper", "m-control","m-foper"]

   DATADIR=[]
   for s in [0,1,2]:
      print('/glade/scratch/cherubin/oper/post/validation/%s/ec_statdata/ts/' % exp[s])
      DATADIR.append('/glade/scratch/cherubin/oper/post/validation/%s/ec_statdata/ts/' % exp[s])

   print (DATADIR)

   BASELINEDIR = '/glade/scratch/cherubin/oper/post/validation/m-conv/ec_statdata/ts/'

   if (fcsh == 3 or fcsh == 9):
      cycles = [3,9,15,21]
   else:
      cycles = [0,6,12,18]

   # open figure before cycle

   plt.figure()
   fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(20, 5))
   axs.set_title('Corrected RMS - WRF vs GFS analysis')
   rootfigname = ('ts_wrfVSec_rms_%02dh_%s_%d' % (fcsh, field, int(L) ))
   print(rootfigname)

   # initialize empty dataframe
#   cols = ["time","cycle","f_hr","ts_bias_m-oper","ts_bias_m-control","ts_bias_m-conv",
#                                 "ts_rmse_m-oper","ts_rmse_m-control","ts_rmse_m-conv",
#                                 "ts_rmse_corr_m-oper","ts_rmse_corr_m-control","ts_rmse_corr_m-conv"]
   df = pd.DataFrame()

   # 0, m-oper, blu; 1, m-control, cyan
   for s in [0,1,2]:
      EXPDIR = DATADIR[s]

      for cyc in cycles:
         analysis_file = ('%s/ecmwf_ts.t%02dz+%02dh.nc' % (BASELINEDIR,cyc,fcsh))
         ref_file = ('%s/wrfprs_ts.t%02dz+%02dh.nc' % (BASELINEDIR,cyc,fcsh))
         exp_file = ('%s/wrfprs_ts.t%02dz+%02dh.nc' % (EXPDIR,cyc,fcsh))
      
         print(analysis_file)
         print(ref_file)
         print(exp_file)
      
         # Load GFS/ECMWF analysis file
         ds_a = xr.open_dataset(analysis_file)
      
         # lon - 360.0 to match coordinates in WRF files. No need in ECMWF?
         # ds_a = ds_a.assign_coords(lon=(ds_a.lon - 360.) )
         
         # Load WRF reference file
         ds_ref = xr.open_dataset(ref_file)
         ds_ref['lon'] = ds_a['lon']   # They are the same but different at the fourth decimal, so bypassing by equalling
         # Load WRF experiment file
         ds_exp = xr.open_dataset(exp_file)
         ds_exp['lon'] = ds_a['lon']   # They are the same but different at the fourth decimal, so bypassing by equalling
         print(ds_exp.coords)
      
         sel_a  =ds_a[field].sel(lev=L).drop('lev')
         sel_ref=ds_ref[field].sel(lev=L).drop('lev')
         sel_exp=ds_exp[field].sel(lev=L).drop('lev')
      
         # Timeserie of RMSE for ocean data only.
         #Skip na so to skip lat longitude raw of Nan
        
         sea_ds_a = ds_a.where(ds_a == 0)
         sea_sel_a  =ds_a[field].sel(lev=L).drop('lev')
         sea_ds_ref = ds_ref.where(ds_a == 0)
         sea_sel_ref  =ds_ref[field].sel(lev=L).drop('lev')
         sea_ds_exp = ds_exp.where(ds_a == 0)
         sea_sel_exp  =ds_exp[field].sel(lev=L).drop('lev')

         # --- Added to mask istances of RH << 0 or >> 100%% in ECMWF analyses!!! Why is this happening?
         if (field == 'rh'):
            sea_sel_a = sea_sel_a.where(sea_sel_a > 0)
            sea_sel_ref = sea_sel_ref.where(sea_sel_a > 0)
            sea_sel_exp = sea_sel_exp.where(sea_sel_a > 0)

            sea_sel_a = sea_sel_a.where(sea_sel_a <= 100)
            sea_sel_ref = sea_sel_ref.where(sea_sel_a <= 100)
            sea_sel_exp = sea_sel_exp.where(sea_sel_a <= 100)

         ts_bias_exp = (sea_sel_exp - sea_sel_a).mean(dim=['lat','lon'], skipna=True)
         ts_bias_ref = (sea_sel_ref - sea_sel_a).mean(dim=['lat','lon'], skipna=True)
         #df_ts_bias_exp = ts_bias_exp.to_pandas()

         label=('ts_bias_%s' % exp[s])
         df_ts_bias_exp = ts_bias_exp.to_dataframe(name=label)
         label_ref=('ts_bias_m-conv')
         df_ts_bias_ref = ts_bias_ref.to_dataframe(name=label_ref)

         ts_rmse_exp = rmse(sea_sel_exp , sea_sel_a, dim=['lat','lon'],skipna=True) 
         ts_rmse_ref = rmse(sea_sel_ref , sea_sel_a, dim=['lat','lon'],skipna=True) 
         label=('ts_rmse_%s' % exp[s])
         df_ts_rmse_exp = ts_rmse_exp.to_dataframe(name=label)
         label_ref=('ts_rsme_m-conv')
         df_ts_rmse_ref = ts_rmse_ref.to_dataframe(name=label_ref)
      
         ts_rmse_corr_exp = xu.sqrt(ts_rmse_exp*ts_rmse_exp - ts_bias_exp*ts_bias_exp)
         ts_rmse_corr_ref = xu.sqrt(ts_rmse_ref*ts_rmse_ref - ts_bias_ref*ts_bias_ref)
         label=('ts_rmse_corr_%s' % exp[s])
         df_ts_rmse_corr_exp = ts_rmse_corr_exp.to_dataframe(name=label)
         label_ref=('ts_rsme_corr_m-conv')
         df_ts_rmse_corr_ref = ts_rmse_corr_ref.to_dataframe(name=label_ref)

         df_local = [df_ts_bias_exp, df_ts_bias_ref, df_ts_rmse_exp, df_ts_rmse_ref, df_ts_rmse_corr_exp, df_ts_rmse_corr_ref]
         df_merged = reduce(lambda  left,right: pd.merge(left,right,on=['time'],how='outer'), df_local)
         #df_merged = reduce(lambda  left,right: pd.merge(left,right,on=['time'],how='outer'), df_local).fillna('void')
         df_merged['cycle']=cyc
         df_merged['f_hr']=fcsh
         print("df_merged")
         print(df_merged)

         df = df.combine_first(df_merged)
         df = df.reindex(df.columns, axis=1)
         print(df)

         # -- Time formatter ---
         #mdates = ts_rmse_corr_exp.coords['time'][:].values
   
         locator = mdates.AutoDateLocator()
         formatter = mdates.ConciseDateFormatter(locator)
         formatter.formats = ['%y',  # ticks are mostly years
                              '%b',       # ticks are mostly months
                              '%d',       # ticks are mostly days
                              '%H:%M',    # hrs
                              '%H:%M',    # min
                              '%S.%f', ]  # secs
         # these are mostly just the level above...
         formatter.zero_formats = [''] + formatter.formats[:-1]
         # ...except for ticks that are mostly hours, then it is nice to have
         # month-day:
         formatter.zero_formats[3] = '%d-%b'
   
         formatter.offset_formats = ['',
                                     '%Y',
                                     '%b %Y',
                                     '%d %b %Y',
                                     '%d %b %Y',
                                     '%d %b %Y %H:%M', ]
   
         axs.xaxis.set_major_locator(locator)
         axs.xaxis.set_major_formatter(formatter)
      
         print(ts_rmse_exp)

         if (s == 0): 
             col = 'blue'
         if (s == 1):
             col = 'cyan'
         if (s == 2):
             col = 'g'
   
         if(field == 'rh'):
             ts_rmse_corr_exp.plot.line(color=col,marker='o', ylim=[12,25])
             ts_rmse_corr_ref.plot.line(color='red',marker='o',ylim=[12,25])
             #ts_rmse_corr_exp.plot.line(color='blue',marker='o', ylim=[5,18], xlim=['20201120 00:00:00','20201128 00:00:00'])
             #ts_rmse_corr_ref.plot.line(color='red',marker='o',ylim=[5,18], xlim=['20201120 00:00:00','20201128 00:00:00'])
         elif (field == 't'):
             ts_rmse_corr_exp.plot.line(color=col,marker='o')
             ts_rmse_corr_ref.plot.line(color='red',marker='o')
             #ts_rmse_corr_exp.plot.line(color='blue',marker='o', xlim=['20201120 00:00:00','20201128 00:00:00'])
             #ts_rmse_corr_ref.plot.line(color='red',marker='o', xlim=['20201120 00:00:00','20201128 00:00:00'])
      
   df.to_csv('ec_stats.csv')

   plt.ylabel(variable_label)
   for label in axs.get_xticklabels():
       label.set_rotation(0)
       label.set_horizontalalignment('center')
   
   figname = ('%s_2.png' % rootfigname)
   plt.savefig(figname)
   figname = ('%s_2.eps' % rootfigname)
   plt.savefig(figname, format='eps')

#   dfn = pd.DataFrame(columns=['n_m-oper','n_m-control','n_m-foper'])
   dfn = pd.DataFrame()
   dfn_list=[]
#   for s in [0,1]:
#      cmd=('/usr/bin/bash ./read_n_obs.sh %s ' % exp[s])
#      call_bash = Popen(cmd, shell=True, stdout=subprocess.PIPE)
#      output = call_bash.communicate()[0].decode('utf-8').strip()
#      df_nobs = pd.DataFrame()
#      df_nobs = pd.read_csv(io.StringIO(output), sep=",").set_index("time")
#      df_nobs.index = pd.to_datetime(df_nobs.index)
#      print(df_nobs)
#      print(df_nobs.size)
#      df_n = df_nobs.groupby('time').sum().reset_index().set_index("time")
#      print(df_n)
#      print(df_n.size)
#
#      # create new dataframe for padding with zeros
#      df_o = pd.DataFrame(index=pd.date_range(start='11/20/2020',
#                                             end=df.index.max(), freq='3H'))
#      df_n = df_o.join(df_n).fillna(0)
#      print(df_n)
#      
#      name=('obs.%s.csv' % exp[s])
#      df_n.to_csv(name)
#
      #sys.exit(1)

if __name__ == "__main__":
    main(len(sys.argv), sys.argv)
