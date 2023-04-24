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

   if not (argc == 5):
        print ('Usage: %s <exp_name> <fcst hour> <var to plot (rmse)> <level>' % argv[0])
        sys.exit(1)

# include loop over cycles
   exp=sys.argv[1]
#   cyc=int(sys.argv[2])
   fcsh=int(sys.argv[2])
   field=sys.argv[3]
   L=float(sys.argv[4])

   if (field == 'rh'):
       variable_label = "Relative Humidity (%)"
   if (field == 't'):
       variable_label = "Temperature (C)"

   # Fixed parameters

   DATADIR = ('/glade/scratch/cherubin/oper/post/validation/%s/statdata/ts/' % exp)
   BASELINEDIR = '/glade/scratch/cherubin/oper/post/validation/m-conv/statdata/ts/'

   if (fcsh == 3 or fcsh == 9):
      cycles = [3,9]
   else:
      cycles = [0,6,12]

   # open figure before cycle

   plt.figure()
   fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(20, 5))
   axs.set_title('Corrected RMS - WRF vs GFS analysis')

   for cyc in cycles:
      analysis_file = ('%s/gfs_ts.t%02dz+%02dh.nc' % (BASELINEDIR,cyc,fcsh))
      ref_file = ('%s/wrfprs_ts.t%02dz+%02dh.nc' % (BASELINEDIR,cyc,fcsh))
      exp_file = ('%s/wrfprs_ts.t%02dz+%02dh.nc' % (DATADIR,cyc,fcsh))
   
      print(analysis_file)
      print(ref_file)
      print(exp_file)
   
      rootfigname = ('ts_wrfVSgfs_rms_t%02d+%02dh_%s_%d' % (cyc, fcsh, field, int(L) ))
      print(rootfigname)
   
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
   
      # Timeserie of RMSE for ocean data only.
      #Skip na so to skip lat longitude raw of Nan
     
      sea_ds_a = ds_a.where(ds_a == 0)
      sea_sel_a  =ds_a[field].sel(lev=L).drop('lev')
      sea_ds_ref = ds_ref.where(ds_a == 0)
      sea_sel_ref  =ds_ref[field].sel(lev=L).drop('lev')
      sea_ds_exp = ds_exp.where(ds_a == 0)
      sea_sel_exp  =ds_exp[field].sel(lev=L).drop('lev')
      
      ts_bias_exp = (sea_sel_exp - sea_sel_a).mean(dim=['lat','lon'], skipna=True)
      ts_bias_ref = (sea_sel_ref - sea_sel_a).mean(dim=['lat','lon'], skipna=True)
   
      ts_rmse_exp = rmse(sea_sel_exp , sea_sel_a, dim=['lat','lon'],skipna=True) 
      ts_rmse_ref = rmse(sea_sel_ref , sea_sel_a, dim=['lat','lon'],skipna=True) 
   
      ts_rmse_corr_exp = xu.sqrt(ts_rmse_exp*ts_rmse_exp - ts_bias_exp*ts_bias_exp)
      ts_rmse_corr_ref = xu.sqrt(ts_rmse_ref*ts_rmse_ref - ts_bias_ref*ts_bias_ref)

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

      if(field == 'rh'):
          ts_rmse_corr_exp.plot.line(color='blue',marker='o', ylim=[5,18])
          ts_rmse_corr_ref.plot.line(color='red',marker='o',ylim=[5,18])
          #ts_rmse_corr_exp.plot.line(color='blue',marker='o', ylim=[5,18], xlim=['20201120 00:00:00','20201128 00:00:00'])
          #ts_rmse_corr_ref.plot.line(color='red',marker='o',ylim=[5,18], xlim=['20201120 00:00:00','20201128 00:00:00'])
      elif (field == 't'):
          ts_rmse_corr_exp.plot.line(color='blue',marker='o')
          ts_rmse_corr_ref.plot.line(color='red',marker='o')
          #ts_rmse_corr_exp.plot.line(color='blue',marker='o', xlim=['20201120 00:00:00','20201128 00:00:00'])
          #ts_rmse_corr_ref.plot.line(color='red',marker='o', xlim=['20201120 00:00:00','20201128 00:00:00'])
   
   #   sys.exit(0)
   plt.ylabel(variable_label)
   for label in axs.get_xticklabels():
       label.set_rotation(0)
       label.set_horizontalalignment('center')
   
   
   figname = ('%s_2.png' % rootfigname)
   plt.savefig(figname)
   figname = ('%s_2.eps' % rootfigname)
   plt.savefig(figname, format='eps')
   
if __name__ == "__main__":
    main(len(sys.argv), sys.argv)
