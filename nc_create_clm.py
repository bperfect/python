# -*- coding: utf-8 -*-
"""
Created on Wed Aug 03 15:26:54 2016

@author: bperfect
"""

import numpy as np
from datetime import datetime
try:
  import netCDF4 as netCDF
except:
  import netCDF3 as netCDF

def calcStrat(d):
    #Linear stratification, corresponds to N=1e-3 Hz
    #Use dTdz = (N/.0408)^2 (uses ROMS equation of state)
    dTdz = .0006
    return 17+dTdz*d
    #print(np.exp(d/1000))
    #return 12+7.5*np.exp(d/1000)
    
def nc_create_clm(filename, grd, u_ini, v_ini, ubar_ini, vbar_ini, temp_ini, stratified=True):
    print('Ini file being created')

    # create file
    nc = netCDF.Dataset(filename, 'w', format='NETCDF3_64BIT')
    #nc = netCDF.Dataset(filename, 'w', format='NETCDF4')
    nc.Description = 'ROMS file'
    nc.Author = 'pyroms_toolbox.nc_create_roms_file'
    nc.Created = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    nc.title = 'ROMS file'
    nc.set_auto_mask(True)
    

    nc.createDimension('xi_rho', np.size(grd.hgrid.mask_rho,1))
    nc.createDimension('xi_u', np.size(grd.hgrid.mask_u,1))
    nc.createDimension('xi_v', np.size(grd.hgrid.mask_v,1))
    nc.createDimension('xi_psi', np.size(grd.hgrid.mask_psi,1))
    nc.createDimension('eta_rho', np.size(grd.hgrid.mask_rho,0))
    nc.createDimension('eta_u', np.size(grd.hgrid.mask_u,0))
    nc.createDimension('eta_v', np.size(grd.hgrid.mask_v,0))
    nc.createDimension('eta_psi', np.size(grd.hgrid.mask_psi,0))
    nc.createDimension('s_rho', grd.vgrid.N)
    nc.createDimension('s_w', grd.vgrid.Np)
    nc.createDimension('tracer',1)
    nc.createDimension('ocean_time', None)
                    
# Set u and ubar
    nc.createVariable('u','f8',('s_rho','eta_u','xi_u'))
    nc.createVariable('ubar','f8',('eta_u','xi_u'))
    u = nc.variables['u'][:]
    ubar = nc.variables['ubar'][:]
    for i in range(u.shape[1]):
        for j in range(u.shape[2]):
            ubar[i,j] = ubar_ini
            for k in range(u.shape[0]):
                u[k,i,j] = u_ini
    nc.variables['u'][:] = u
    nc.variables['ubar'][:] = ubar
      
# Set v and vbar                
    nc.createVariable('v','f8',('s_rho','eta_v','xi_v'))
    nc.createVariable('vbar','f8',('eta_v','xi_v'))
    v = nc.variables['v'][:]
    vbar = nc.variables['vbar'][:]
    for i in range(v.shape[1]):
        for j in range(v.shape[2]):
            vbar[i,j] = vbar_ini
            for k in range(v.shape[0]):
                v[k,i,j] = v_ini
    nc.variables['v'][:] = v
    nc.variables['vbar'][:] = vbar  

#Set temperature
    if stratified:
        h = grd.vgrid.h
        s_rho = grd.vgrid.s_rho
        nc.createVariable('temp','f8',('s_rho','eta_rho','xi_rho'))
        temp = nc.variables['temp'][:]
        for i in range(temp.shape[1]):
            for j in range(temp.shape[2]):
                elevation = h[i,j]*s_rho[:]
                #print('elevation', elevation)
                for k in range(temp.shape[0]):
                    temp[k,i,j] = calcStrat(elevation[k])
        nc.variables['temp'][:] = temp
    else:
        nc.createVariable('temp','f8',('s_rho','eta_rho','xi_rho'),fill_value = temp_ini)       
    
    nc.close()

