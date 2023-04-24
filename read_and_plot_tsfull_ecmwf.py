import sys,os
import netCDF4
import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
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

   exp=sys.argv[1]
   fcsh=int(sys.argv[2])

   field=sys.argv[3]
   L=float(sys.argv[4])

   # Fixed parameters

   DATADIR = ('/glade/scratch/cherubin/oper/post/validation/%s/ec_statdata/tsfull/' % exp)
   BASELINEDIR = '/glade/scratch/cherubin/oper/post/validation/m-conv/ec_statdata/tsfull/'

   analysis_file = ('%s/ecmwf_tsfullz+%02dh.nc' % (BASELINEDIR,fcsh))
   ref_file = ('%s/wrfprs_tsfullz+%02dh.nc' % (BASELINEDIR,fcsh))
   exp_file = ('%s/wrfprs_tsfullz+%02dh.nc' % (DATADIR,fcsh))

   print(analysis_file)
   print(ref_file)
   print(exp_file)

   rootfigname = ('ec_fullts_rms_%02dh_%s_%d' % (fcsh, field, int(L) ))
   print(rootfigname)

   # Load GFS/ECMWF analysis file
   ds_a = xr.open_dataset(analysis_file)
   print(ds_a)
   sys.exit(1)

   # lon - 360.0 to match coordinates in WRF files.
   #ds_a = ds_a.assign_coords(lon=(ds_a.lon - 360.) )
   
   # Load WRF reference file
   ds_ref = xr.open_dataset(ref_file)
   ds_ref['lon'] = ds_a['lon']   # They are the same but different at the fourth decimal, so bypassing by equalling
   # Load WRF experiment file
   ds_exp = xr.open_dataset(exp_file)
   ds_exp['lon'] = ds_a['lon']   # They are the same but different at the fourth decimal, so bypassing by equalling

   sel_a  =ds_a[field].sel(lev=L).drop('lev')
   sel_ref=ds_ref[field].sel(lev=L).drop('lev')
   sel_exp=ds_exp[field].sel(lev=L).drop('lev')

   print(sel_a["time"])
   print(sel_exp["time"])

   bias_exp = (sel_exp - sel_a).mean(dim='time')
   bias_ref = (sel_ref - sel_a).mean(dim='time')


   rmse_exp = rmse(sel_exp , sel_a, dim='time') 
   rmse_ref = rmse(sel_ref , sel_a, dim='time') 

   rmse_exp_corr = xu.sqrt(rmse_exp*rmse_exp - bias_exp*bias_exp) 
   rmse_ref_corr = xu.sqrt(rmse_ref*rmse_ref - bias_ref*bias_ref) 

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

   # Vertical profile of RMSE, averaged over lat, lon (no land), and over time)
   # re-use variable
   sea_sel_a = ds_a[field]
   sea_sel_ref = ds_ref[field]
   sea_sel_exp = ds_exp[field]

   vert_bias_exp = (sea_sel_exp - sea_sel_a).mean(dim=['lat','lon','time'], skipna=True)
   vert_bias_ref = (sea_sel_ref - sea_sel_a).mean(dim=['lat','lon','time'], skipna=True)

   vert_rmse_exp = rmse(sea_sel_exp , sea_sel_a, dim=['lat','lon','time'],skipna=True)
   vert_rmse_ref = rmse(sea_sel_ref , sea_sel_a, dim=['lat','lon','time'],skipna=True)
   
   vert_rmse_corr_exp = xu.sqrt(vert_rmse_exp*vert_rmse_exp - vert_bias_exp*vert_bias_exp)
   vert_rmse_corr_ref = xu.sqrt(vert_rmse_ref*vert_rmse_ref - vert_bias_ref*vert_bias_ref)

   # Vertical profile of RMSE, averaged over lat, lon (no land) - timeseries
   # re-use variable
   sea_sel_a = ds_a[field]
   sea_sel_ref = ds_ref[field]
   sea_sel_exp = ds_exp[field]

   ts_vert_bias_exp = (sea_sel_exp - sea_sel_a).mean(dim=['lat','lon'], skipna=True)
   ts_vert_bias_ref = (sea_sel_ref - sea_sel_a).mean(dim=['lat','lon'], skipna=True)

   ts_vert_rmse_exp = rmse(sea_sel_exp , sea_sel_a, dim=['lat','lon'],skipna=True)
   ts_vert_rmse_ref = rmse(sea_sel_ref , sea_sel_a, dim=['lat','lon'],skipna=True)

   ts_vert_rmse_corr_exp = xu.sqrt(ts_vert_rmse_exp*ts_vert_rmse_exp - ts_vert_bias_exp*ts_vert_bias_exp)
   ts_vert_rmse_corr_ref = xu.sqrt(ts_vert_rmse_ref*ts_vert_rmse_ref - ts_vert_bias_ref*ts_vert_bias_ref)

#   sel=ds_a.sel(lev=500, time=ds_a.coords['time'][0].values)
#   df_a = ds_a.to_dataframe()
   #print(ds_a.rh.xs(L,level='lev'))


   lon_formatter = LongitudeFormatter(number_format='.1f',
                                      degree_symbol='',
                                      dateline_direction_label=True)
   lat_formatter = LatitudeFormatter(number_format='.1f',
                                     degree_symbol='')

   if(field == 'rh'):
      bar_levels = np.arange(0,50,5)
   elif (field == 't'):
      bar_levels = np.arange(0,2,0.2)

   plt.figure()

   xr.set_options(cmap_sequential='jet')

   fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 8), subplot_kw={'projection': ccrs.PlateCarree()})
   fig.suptitle('BIAS CORRECTED RMSE comparison: Conv vs MW ')

   # First subplot 
   rmse_ref_corr.plot(ax=ax1, levels=bar_levels, cbar_kwargs={'ticks': bar_levels})
   #rmse_ref.plot(ax=ax1, levels=bar_levels, add_colorbar=False)
   ax1.coastlines(color='white')
   ax1.set_xticks(np.arange(-165,-145,5), crs=ccrs.PlateCarree())
   ax1.set_yticks(np.arange(10,30,5), crs=ccrs.PlateCarree())
   ax1.xaxis.set_major_formatter(lon_formatter)
   ax1.yaxis.set_major_formatter(lat_formatter)

   # Second subplot 
   #ax2 = plt.axes(projection=ccrs.PlateCarree())
   rmse_exp_corr.plot(ax=ax2, levels=bar_levels, cbar_kwargs={'ticks': bar_levels})
   ax2.coastlines(color='white')
   ax2.set_xticks(np.arange(-165,-145,5), crs=ccrs.PlateCarree())
   ax2.set_yticks(np.arange(10,30,5), crs=ccrs.PlateCarree())
   ax2.xaxis.set_major_formatter(lon_formatter)
   ax2.yaxis.set_major_formatter(lat_formatter)

   #plt.tight_layout()
   figname = ('%s_0.png' % rootfigname)
   plt.savefig(figname)

   plt.figure()

   xr.set_options(cmap_sequential='jet')

   fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 8), subplot_kw={'projection': ccrs.PlateCarree()})
   fig.suptitle('RMSE comparison: Conv vs MW ')

   # First subplot 
   rmse_ref.plot(ax=ax1, levels=bar_levels, cbar_kwargs={'ticks': bar_levels})
   #rmse_ref.plot(ax=ax1, levels=bar_levels, add_colorbar=False)
   ax1.coastlines(color='white')
   ax1.set_xticks(np.arange(-165,-145,5), crs=ccrs.PlateCarree())
   ax1.set_yticks(np.arange(10,30,5), crs=ccrs.PlateCarree())
   ax1.xaxis.set_major_formatter(lon_formatter)
   ax1.yaxis.set_major_formatter(lat_formatter)

   # Second subplot 
   #ax2 = plt.axes(projection=ccrs.PlateCarree())
   rmse_exp.plot(ax=ax2, levels=bar_levels, cbar_kwargs={'ticks': bar_levels})
   ax2.coastlines(color='white')
   ax2.set_xticks(np.arange(-165,-145,5), crs=ccrs.PlateCarree())
   ax2.set_yticks(np.arange(10,30,5), crs=ccrs.PlateCarree())
   ax2.xaxis.set_major_formatter(lon_formatter)
   ax2.yaxis.set_major_formatter(lat_formatter)

   #plt.tight_layout()
   figname = ('%s_1.png' % rootfigname)
   plt.savefig(figname)

   plt.figure()
   fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(20, 5))

#   ts_rmse_exp.plot.line(color='blue',linestyle='dashed')
#   ts_rmse_ref.plot.line(color='red',linestyle='dashed')
   ts_rmse_corr_exp.plot.line(color='blue')
   ts_rmse_corr_ref.plot.line(color='red')
   figname = ('%s_2.png' % rootfigname)
   plt.savefig(figname)

   plt.figure()
   fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(5, 10))

   vert_rmse_corr_exp.plot.line(y='lev', yincrease=False,  color='blue')
   vert_rmse_corr_ref.plot.line(y='lev', yincrease=False, color='red')
#   vert_rmse_exp.plot.line(y='lev', yincrease=False,  color='blue', linestyle='dashed')
#   vert_rmse_ref.plot.line(y='lev', yincrease=False, color='red', linestyle='dashed')
   figname = ('%s_3.png' % rootfigname)
   plt.savefig(figname)

   times = ts_vert_rmse_corr_exp.coords['time'][:].values

   plt.figure()
   fig, axs = plt.subplots(nrows=13, ncols=10, figsize=(30, 60))

   # should plot vertical corrected rms for each day of timeserie
   
   for axis,t in zip(axs.flat,range(len(times))):  
       ts_vert_rmse_corr_exp[t,:].plot.line(y='lev', ax=axis, yincrease=False,color='blue')
       ts_vert_rmse_corr_ref[t,:].plot.line(y='lev', ax=axis, yincrease=False,color='red')

   figname = ('%s_4.png' % rootfigname)
   plt.savefig(figname)



   sys.exit(0)


if __name__ == "__main__":
    main(len(sys.argv), sys.argv)
