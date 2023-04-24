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

   if not (argc == 7):
        print ('Usage: %s <exp_name> <cycle> <fcst hour> <n of cycle in the serie> <var to plot (rmse)> <level>' % argv[0])
        sys.exit(1)

   exp=sys.argv[1]
   cyc=int(sys.argv[2])
   fcsh=int(sys.argv[3])
   cyc_step=int(sys.argv[4])

   field=sys.argv[5]
   L=float(sys.argv[6])

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

   rootfigname = ('wrfVSgfs_bias_%02d+%02dh_%s_%d' % (cyc, fcsh, field, int(L) ))
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

   print(sel_exp[cyc_step,:,:])

   bias_exp = (sel_exp[cyc_step,:,:] - sel_a[cyc_step,:,:])
   bias_ref = (sel_ref[cyc_step,:,:] - sel_a[cyc_step,:,:])

   lon_formatter = LongitudeFormatter(number_format='.1f',
                                      degree_symbol='',
                                      dateline_direction_label=True)
   lat_formatter = LatitudeFormatter(number_format='.1f',
                                     degree_symbol='')

   if(field == 'rh'):
      bar_levels = np.arange(-5,5,0.5)
   elif (field == 't'):
      bar_levels = np.arange(2,2,0.1)

   plt.figure()

#   xr.set_options(cmap_sequential='jet')

   #fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 8), subplot_kw={'projection': ccrs.PlateCarree()})
   fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8), subplot_kw={'projection': ccrs.PlateCarree()})
   fig.suptitle('BIAS comparison Conv vs %s'% exp)

   # First subplot 
   bias_ref.plot(ax=ax1, levels=bar_levels, cbar_kwargs={'ticks': bar_levels})
   #rmse_ref.plot(ax=ax1, levels=bar_levels, add_colorbar=False)
   ax1.coastlines(color='white')
   ax1.set_xticks(np.arange(-165,-145,5), crs=ccrs.PlateCarree())
   ax1.set_yticks(np.arange(10,30,5), crs=ccrs.PlateCarree())
   ax1.xaxis.set_major_formatter(lon_formatter)
   ax1.yaxis.set_major_formatter(lat_formatter)
   #ax1.title.set_text('CO')

   # Second subplot 
   #ax2 = plt.axes(projection=ccrs.PlateCarree())
   bias_exp.plot(ax=ax2, levels=bar_levels, cbar_kwargs={'ticks': bar_levels})
   ax2.coastlines(color='white')
   ax2.set_xticks(np.arange(-165,-145,5), crs=ccrs.PlateCarree())
   ax2.set_yticks(np.arange(10,30,5), crs=ccrs.PlateCarree())
   ax2.xaxis.set_major_formatter(lon_formatter)
   ax2.yaxis.set_major_formatter(lat_formatter)
   #ax2.title.set_text('%s ' % exp)

   # Third subplot 
   #ax2 = plt.axes(projection=ccrs.PlateCarree())
#   (rmse_exp_corr-rmse_ref_corr).plot(ax=ax3, levels=np.arange(-5,5,0.5))
#   ax3.coastlines(color='white')
#   ax3.set_xticks(np.arange(-165,-145,5), crs=ccrs.PlateCarree())
#   ax3.set_yticks(np.arange(10,30,5), crs=ccrs.PlateCarree())
#   ax3.xaxis.set_major_formatter(lon_formatter)
#   ax3.yaxis.set_major_formatter(lat_formatter)

   #plt.tight_layout()
   figname = ('%s_0.png' % rootfigname)
   plt.savefig(figname)
   figname = ('%s_0.eps' % rootfigname)
   plt.savefig(figname, format='eps')

   plt.figure()

   sys.exit(0)


if __name__ == "__main__":
    main(len(sys.argv), sys.argv)
