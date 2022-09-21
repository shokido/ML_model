# Python script for constructing input ocean grid file  
import netCDF4 as ncdf
import numpy as np
from write_gdesfile import out_gdesfile_lonlat
# Grid setting========================================
dz=20.0                   # vertical resolution [m]
top_level=-0.0           # Top level           [m]
bottom_level=-1000.0      # Bottom level        [m]
stretch=-5.0
# Calculation region==================================  
lon_w=140.5 
lon_e=180.5
lat_s=30.5
lat_n=50.5
dlon=0.5
dlat=0.5
omega=7.292e-5
grid_name="nwpac"
fname_out="../FILES/grid_"+grid_name+".nc"
fname_gdes="desg_"+grid_name+".txt"
# Coordinate setting
# Lon
if (lon_w != lon_e):
    nlon=int((lon_e-lon_w)/dlon) + 1
    lon=np.linspace(lon_w,lon_e,nlon)
else:
    nlon=1
    lon=np.zeros(1)
    lon[0]=lon_w
# Lat
if (lat_n != lat_s):
    nlat=int((lat_n-lat_s)/dlat) + 1
    lat=np.linspace(lat_s,lat_n,nlat)
else:
    nlat=1
    lat=np.zeros(1)
    lat[0]=lat_s
# Lev
nz = int((top_level - bottom_level) / dz) + 1
xi=np.linspace(0.0,1.0,nz)
lev_q = abs(bottom_level) * np.log(1.0-xi*(1.0-np.exp(stretch)))/ stretch
lev_q=lev_q[::-1]
z_q=lev_q*(-1.0)
z_rho = z_q
z_rho[0]=1.5*z_q[0]-0.5*z_q[1]
for iz in range(1,nz):
    z_rho[iz]=0.5*(z_q[iz-1]+z_q[iz])
#for iz in range(0,nz-2):
#     z_rho[iz]=0.5*(z_q[iz]+z_q[iz+1])
#z_rho[nz-1]=1.5*z_q[nz-1]-0.5*z_q[nz-2]
lev_rho=z_rho*(-1.0)
coriolis=2.0*omega*np.sin(lat * np.pi/180.0)

nlev_out=len(z_rho)
nlat_out=len(lat)
nlon_out=len(lon)

# Output file
nc_out=ncdf.Dataset(fname_out,'w')
nc_out.createDimension('lev',nlev_out)
nc_out.createDimension('lat',nlat_out)
nc_out.createDimension('lon',nlon_out)

lon_out=nc_out.createVariable('lon',lon.dtype, ('lon'))
lon_out.long_name = 'longitude'
lon_out.units = 'degrees_east'

lat_out=nc_out.createVariable('lat',lat.dtype, ('lat'))
lat_out.long_name = 'latitude'
lat_out.units = 'degrees_north'

lev_rho_out=nc_out.createVariable('lev',lev_rho.dtype, ('lev'))
lev_rho_out.long_name = 'level at rho-grid'
lev_rho_out.units = 'm'

z_rho_out=nc_out.createVariable('z_rho',z_rho.dtype, ('lev'))
z_rho_out.long_name = 'z at rho-grid'
z_rho_out.units = 'm'

# lev_q_out=nc_out.createVariable('lev_q',lev_q.dtype, ('lev'))
# lev_q_out.long_name = 'level at q-grid'
# lev_q_out.units = 'm'

z_q_out=nc_out.createVariable('z_q',z_q.dtype, ('lev'))
z_q_out.long_name = 'z at q-grid'
z_q_out.units = 'm'

coriolis_out=nc_out.createVariable('f',coriolis.dtype, ('lat'))
coriolis_out.long_name = 'Coriolis parameter'
coriolis_out.units = '1/s'

lon_out[:]=lon[:]
lat_out[:]=lat[:]
lev_rho_out[:]=lev_rho[:]
z_rho_out[:]=z_rho[:]
z_q_out[:]=z_q[:]
coriolis_out[:]=coriolis[:]
nc_out.close()
out_gdesfile_lonlat(lon,lat,fname_gdes)
