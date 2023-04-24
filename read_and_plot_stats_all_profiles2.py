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

#from mpl_toolkits.basemap import Basemap

def main(argc, argv):

   if not (argc == 6):
        print ('Usage: %s <exp_name=ALL> <cycling hour> <fcst hour> <var to plot (rmse)> <level>' % argv[0])
        sys.exit(1)

   exp=sys.argv[1]
   cyc=int(sys.argv[2])
   fcsh=int(sys.argv[3])

   field=sys.argv[4]
   L=float(sys.argv[5])

   if (field == 'rh'):
       variable_label = "Relative Humidity (%)"
   if (field == 't'):
       variable_label = "Temperature (C)"

   # Fixed parameters


  # Experiments labels
   exp = ["m-oper", "m-foper"]

   EXPDIR=[]
   for s in [0,1]:
      print('/glade/scratch/cherubin/oper/post/validation/%s/statdata/ts/' % exp[s])
      EXPDIR.append('/glade/scratch/cherubin/oper/post/validation/%s/statdata/ts/' % exp[s])
 
   print (EXPDIR)
   BASELINEDIR = '/glade/scratch/cherubin/oper/post/validation/m-conv/statdata/ts/'
 
   #if (fcsh == 3 or fcsh == 9):
   #   cycles = [3,9]
   #else:
   #   cycles = [0,6,12]
 
   # open figures before cycle
 
   fig1, axs1 = plt.subplots(nrows=1, ncols=1, figsize=(5, 10))
   rootfigname_1 = ('profiles_wrfVSgfs_rms_t%02d+%02dh_%s_%d' % (cyc, fcsh, field, int(L) ))
   print(rootfigname_1)

   fig2, axs2 = plt.subplots(nrows=7, ncols=5, figsize=(30, 40))
   rootfigname_2 = ('ts_profiles_wrfVSgfs_rms_t%02d+%02dh_%s_%d' % (cyc, fcsh, field, int(L) ))
   print(rootfigname_2)

#   fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(20, 5))
#   axs.set_title('Corrected RMS - WRF vs GFS analysis')
 
   # 0, m-oper, blu; 1, m-foper, green

   for s in [0,1]:
      DATADIR = EXPDIR[s]
      analysis_file = ('%s/gfs_ts.t%02dz+%02dh.nc' % (BASELINEDIR,cyc,fcsh))
      ref_file = ('%s/wrfprs_ts.t%02dz+%02dh.nc' % (BASELINEDIR,cyc,fcsh))
      exp_file = ('%s/wrfprs_ts.t%02dz+%02dh.nc' % (DATADIR,cyc,fcsh))
   
      print(analysis_file)
      print(ref_file)
      print(exp_file)
   
      # Load GFS/ECMWF analysis file
      ds_a = xr.open_dataset(analysis_file)
   
      # lon - 360.0 to match coordinates in WRF files.
      ds_a = ds_a.assign_coords(lon=(ds_a.lon - 360.) )
      
      # Load WRF reference file
      ds_ref = xr.open_dataset(ref_file)
      # Load WRF experiment file
      ds_exp = xr.open_dataset(exp_file)
      print(ds_exp.coords)
   
      sel_a  =ds_a[field].sel(lev=L).drop('lev')
      sel_ref=ds_ref[field].sel(lev=L).drop('lev')
      sel_exp=ds_exp[field].sel(lev=L).drop('lev')
   
      bias_exp = (sel_exp - sel_a).mean(dim='time')
      bias_ref = (sel_ref - sel_a).mean(dim='time')
   
      rmse_exp = rmse(sel_exp , sel_a, dim='time') 
      rmse_ref = rmse(sel_ref , sel_a, dim='time') 
   
      rmse_exp_corr = xu.sqrt(rmse_exp*rmse_exp - bias_exp*bias_exp) 
      rmse_ref_corr = xu.sqrt(rmse_ref*rmse_ref - bias_ref*bias_ref) 
   
      # Timeserie of RMSE for ocean data only.
      #Skip na so to skip lat longitude raw of Nan
     
      sea_ds_a = ds_a.where(ds_a["xland"] == 0)
      sea_sel_a  =sea_ds_a[field].sel(lev=L).drop('lev')
      sea_ds_ref = ds_ref.where(ds_ref["xland"] == 0)
      sea_sel_ref  =sea_ds_ref[field].sel(lev=L).drop('lev')
      sea_ds_exp = ds_exp.where(ds_exp["xland"] == 0)
      sea_sel_exp  =sea_ds_exp[field].sel(lev=L).drop('lev')
      
      ts_bias_exp = (sea_sel_exp - sea_sel_a).mean(dim=['lat','lon'], skipna=True)
      ts_bias_ref = (sea_sel_ref - sea_sel_a).mean(dim=['lat','lon'], skipna=True)
   
      ts_rmse_exp = rmse(sea_sel_exp , sea_sel_a, dim=['lat','lon'],skipna=True) 
      ts_rmse_ref = rmse(sea_sel_ref , sea_sel_a, dim=['lat','lon'],skipna=True) 
   
      ts_rmse_corr_exp = xu.sqrt(ts_rmse_exp*ts_rmse_exp - ts_bias_exp*ts_bias_exp)
      ts_rmse_corr_ref = xu.sqrt(ts_rmse_ref*ts_rmse_ref - ts_bias_ref*ts_bias_ref)
   
      # Vertical profile of RMSE, averaged over lat, lon (no land), and over time)
      # re-use variable
      sea_sel_a = sea_ds_a[field]
      sea_sel_ref =sea_ds_ref[field]
      sea_sel_exp =sea_ds_exp[field]
   
      vert_bias_exp = (sea_sel_exp - sea_sel_a).mean(dim=['lat','lon','time'], skipna=True)
      vert_bias_ref = (sea_sel_ref - sea_sel_a).mean(dim=['lat','lon','time'], skipna=True)
   
      vert_rmse_exp = rmse(sea_sel_exp , sea_sel_a, dim=['lat','lon','time'],skipna=True)
      vert_rmse_ref = rmse(sea_sel_ref , sea_sel_a, dim=['lat','lon','time'],skipna=True)
      
      vert_rmse_corr_exp = xu.sqrt(vert_rmse_exp*vert_rmse_exp - vert_bias_exp*vert_bias_exp)
      vert_rmse_corr_ref = xu.sqrt(vert_rmse_ref*vert_rmse_ref - vert_bias_ref*vert_bias_ref)
   
      # Vertical profile of RMSE, averaged over lat, lon (no land) - timeseries
      # re-use variable
      sea_sel_a = sea_ds_a[field]
      sea_sel_ref = sea_ds_ref[field]
      sea_sel_exp = sea_ds_exp[field]
   
      ts_vert_bias_exp = (sea_sel_exp - sea_sel_a).mean(dim=['lat','lon'], skipna=True)
      ts_vert_bias_ref = (sea_sel_ref - sea_sel_a).mean(dim=['lat','lon'], skipna=True)
   
      ts_vert_rmse_exp = rmse(sea_sel_exp , sea_sel_a, dim=['lat','lon'],skipna=True)
      ts_vert_rmse_ref = rmse(sea_sel_ref , sea_sel_a, dim=['lat','lon'],skipna=True)
   
      ts_vert_rmse_corr_exp = xu.sqrt(ts_vert_rmse_exp*ts_vert_rmse_exp - ts_vert_bias_exp*ts_vert_bias_exp)
      ts_vert_rmse_corr_ref = xu.sqrt(ts_vert_rmse_ref*ts_vert_rmse_ref - ts_vert_bias_ref*ts_vert_bias_ref)
   
   #   sel=ds_a.sel(lev=500, time=ds_a.coords['time'][0].values)
   #   df_a = ds_a.to_dataframe()
      #print(ds_a.rh.xs(L,level='lev'))
   
      if(field == 'rh'):
         bar_levels = np.arange(0,50,5)
      elif (field == 't'):
         bar_levels = np.arange(0,2,0.2)

      if (s == 0):
         col = 'blue'
      if (s == 1):
         col = 'g'
   
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
   
      plt.figure(1) # Open figure for Average Vertical Profile of rMS
   
      vert_rmse_corr_exp.plot.line(y='lev', ax=axs1, yincrease=False,  color=col)
      vert_rmse_corr_ref.plot.line(y='lev', ax=axs1, yincrease=False, color='red')
   #   vert_rmse_exp.plot.line(y='lev', yincrease=False,  color='blue', linestyle='dashed')
   #   vert_rmse_ref.plot.line(y='lev', yincrease=False, color='red', linestyle='dashed')
      plt.ylabel('Pressure (mb)',fontsize=14)
      plt.xlabel(variable_label,fontsize=14)
      axs1.set_title('Corrected RMS - WRF vs GFS analysis - %02d cycle' % cyc, fontsize=14)
      axs1.xaxis.label.set_size(14)
   
      times = ts_vert_rmse_corr_exp.coords['time'][:].values

      # -- Second Figure
      plt.figure(2) # Open figure for Time Serie of Vertical Profile of rMS
   
      # should plot vertical corrected rms for each day of timeserie
      
      for axis,t in zip(axs2.flat,range(len(times))):  
        if(field == 'rh'):
          ts_vert_rmse_corr_exp[t,:].plot.line(y='lev', ax=axis, xlim=[0,28], yincrease=False,color=col)
          ts_vert_rmse_corr_ref[t,:].plot.line(y='lev', ax=axis, xlim=[0,28], yincrease=False,color='red')
        elif (field == 't'): 
          ts_vert_rmse_corr_exp[t,:].plot.line(y='lev', ax=axis, xlim=[0,2], yincrease=False,color=col)
          ts_vert_rmse_corr_ref[t,:].plot.line(y='lev', ax=axis, xlim=[0,2], yincrease=False,color='red')
#        axis.xaxis.label.set_size(20)
        axis.tick_params(axis='both', labelsize=20)


  # Save Figures
   
   plt.figure(1)
   figname = ('%s_1.png' % rootfigname_1)
   plt.savefig(figname)
   figname = ('%s_1.eps' % rootfigname_1)
   plt.savefig(figname, format='eps')

   plt.figure(2)
   figname = ('%s_2.png' % rootfigname_2)
   plt.savefig(figname)
   figname = ('%s_2.eps' % rootfigname_2)
   plt.savefig(figname, format='eps')

   sys.exit(0)


if __name__ == "__main__":
    main(len(sys.argv), sys.argv)
