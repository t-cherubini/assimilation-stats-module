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
import seaborn as sns

#from mpl_toolkits.basemap import Basemap

def main(argc, argv):

   if not (argc == 6):
        print ('Usage: %s <exp_name> <cycling hour> <fcst hour> <var to plot (rmse)> <level>' % argv[0])
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

   DATADIR = ('/glade/scratch/cherubin/oper/post/validation/%s/statdata/ts/' % exp)
   BASELINEDIR = '/glade/scratch/cherubin/oper/post/validation/m-conv/statdata/ts/'

   analysis_file = ('%s/gfs_ts.t%02dz+%02dh.nc' % (BASELINEDIR,cyc,fcsh))
   ref_file = ('%s/wrfprs_ts.t%02dz+%02dh.nc' % (BASELINEDIR,cyc,fcsh))
   exp_file = ('%s/wrfprs_ts.t%02dz+%02dh.nc' % (DATADIR,cyc,fcsh))

   print(analysis_file)
   print(ref_file)
   print(exp_file)

   rootfigname = ('wrfVSgfs_rms_t%02d+%02dh_%s_%d' % (cyc, fcsh, field, int(L) ))
   print(rootfigname)

   # Load GFS/ECMWF analysis file
   ds_a = xr.open_dataset(analysis_file)

   # lon - 360.0 to match coordinates in WRF files.
   ds_a = ds_a.assign_coords(lon=(ds_a.lon - 360.) )

   print(ds_a)

   
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
   sea_sel_ref = sea_ds_ref[field]
   sea_sel_exp = sea_ds_exp[field]

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

# convert to dataframe for later plotting
   print(rmse_exp)
   df_rmse_exp=rmse_exp.to_dataframe()
   print(df_rmse_exp)
   df_rmse_exp['label']=('%s' % exp)
   df_rmse_ref=rmse_ref.to_dataframe()
   df_rmse_ref['label']=('%s' % 'm-conv')
   df = pd.concat([df_rmse_exp, df_rmse_ref], ignore_index=True,axis=0)
  
   median_exp=('Median RMSE %s = %.5f' % (exp, round(df_rmse_exp.rh.median(),3)))
   median_ref=('Median RMSE m-conv = %.5f' % round(df_rmse_ref.rh.median(),3))

   plt.figure(figsize=(10,6))
   f, axes=plt.subplots(1,1)
   #sns.histplot(df,x='rh',discrete=True,stat='probability',common_norm=False, shrink=0.8, ax=axes, hue='label',multiple='dodge')
   sns.histplot(df,x='rh',discrete=True,stat='probability',common_norm=False, shrink=0.8, ax=axes, hue='label')
   axes.axvline(df_rmse_exp.rh.median(), color='b', ls='-', lw=1, label=('Median %s'% exp ))
   axes.axvline(df_rmse_ref.rh.median(), color='orange', ls='--', lw=1, label=('Median %s' % 'm-conv'))
   axes.text(20,0.07,median_exp, fontsize='small')
   axes.text(20,0.06,median_ref,fontsize='small')
   axes.set_ylim([0,0.177])
   axes.set_xlim([0,30])
   axes.set_xlabel('RH RMSE')
   pngfile=('rmse_hist.png')
   plt.savefig(pngfile, dpi=700)
   epsfile=('rmse_hist.eps')
   plt.savefig(epsfile, forma='eps')
   plt.close()

   #sys.exit(1)


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

   fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 8), subplot_kw={'projection': ccrs.PlateCarree()})
   fig.suptitle('BIAS CORRECTED RMSE comparison: Conv vs %s ' % exp)

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

   # Third subplot 
   #ax2 = plt.axes(projection=ccrs.PlateCarree())
   (rmse_exp_corr-rmse_ref_corr).plot(ax=ax3, levels=np.arange(-5,5,0.5))
   ax3.coastlines(color='white')
   ax3.set_xticks(np.arange(-165,-145,5), crs=ccrs.PlateCarree())
   ax3.set_yticks(np.arange(10,30,5), crs=ccrs.PlateCarree())
   ax3.xaxis.set_major_formatter(lon_formatter)
   ax3.yaxis.set_major_formatter(lat_formatter)

   #plt.tight_layout()
   figname = ('%s_0.png' % rootfigname)
   plt.savefig(figname)
   figname = ('%s_0.eps' % rootfigname)
   plt.savefig(figname, format='eps')

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
   figname = ('%s_1.eps' % rootfigname)
   plt.savefig(figname, format='eps')

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

   plt.figure()
   fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(20, 5))

   axs.xaxis.set_major_locator(locator)
   axs.xaxis.set_major_formatter(formatter)

#   ts_rmse_exp.plot.line(color='blue',linestyle='dashed')
#   ts_rmse_ref.plot.line(color='red',linestyle='dashed')

   if(field == 'rh'):
       ts_rmse_corr_exp.plot.line(color='blue',marker='o', ylim=[5,18], xlim=['20201120 00:00:00','20201128 00:00:00'])
       ts_rmse_corr_ref.plot.line(color='red',marker='o',ylim=[5,18], xlim=['20201120 00:00:00','20201128 00:00:00'])
   elif (field == 't'):
       ts_rmse_corr_exp.plot.line(color='blue',marker='o', xlim=['20201120 00:00:00','20201128 00:00:00'])
       ts_rmse_corr_ref.plot.line(color='red',marker='o', xlim=['20201120 00:00:00','20201128 00:00:00'])

   plt.ylabel(variable_label)
   for label in axs.get_xticklabels():
       label.set_rotation(0)
       label.set_horizontalalignment('center')

   axs.set_title('Corrected RMS - WRF vs GFS analysis - %02d cycle' % cyc)

   figname = ('%s_2.png' % rootfigname)
   plt.savefig(figname)
   figname = ('%s_2.eps' % rootfigname)
   plt.savefig(figname, format='eps')

   plt.figure()
   fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(5, 10))

   vert_rmse_corr_exp.plot.line(y='lev', yincrease=False,  color='blue')
   vert_rmse_corr_ref.plot.line(y='lev', yincrease=False, color='red')
#   vert_rmse_exp.plot.line(y='lev', yincrease=False,  color='blue', linestyle='dashed')
#   vert_rmse_ref.plot.line(y='lev', yincrease=False, color='red', linestyle='dashed')
   plt.ylabel('Pressure (mb)')
   plt.xlabel(variable_label)
   axs.set_title('Corrected RMS - WRF vs GFS analysis - %02d cycle' % cyc)
   figname = ('%s_3.png' % rootfigname)
   plt.savefig(figname)
   figname = ('%s_3.eps' % rootfigname)
   plt.savefig(figname, format='eps')

   times = ts_vert_rmse_corr_exp.coords['time'][:].values

   plt.figure()
   fig, axs = plt.subplots(nrows=7, ncols=5, figsize=(30, 40))

   # should plot vertical corrected rms for each day of timeserie
   
   for axis,t in zip(axs.flat,range(len(times))):  
     if(field == 'rh'):
       ts_vert_rmse_corr_exp[t,:].plot.line(y='lev', ax=axis, xlim=[0,25], yincrease=False,color='blue')
       ts_vert_rmse_corr_ref[t,:].plot.line(y='lev', ax=axis, xlim=[0,25], yincrease=False,color='red')
     elif (field == 't'): 
       ts_vert_rmse_corr_exp[t,:].plot.line(y='lev', ax=axis, xlim=[0,2], yincrease=False,color='blue')
       ts_vert_rmse_corr_ref[t,:].plot.line(y='lev', ax=axis, xlim=[0,2], yincrease=False,color='red')

   figname = ('%s_4.png' % rootfigname)
   plt.savefig(figname)
   figname = ('%s_4.eps' % rootfigname)
   plt.savefig(figname, format='eps')



   sys.exit(0)


if __name__ == "__main__":
    main(len(sys.argv), sys.argv)
