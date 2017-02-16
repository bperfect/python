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
    
def nc_create_nudge(filename, grd, u_ini, v_ini, ubar_ini, vbar_ini, temp_ini, stratified=True):
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

    
                    
# Create M2 and M3 Nudge variables
    nc.createVariable('M3_NudgeCoef','f8',('s_rho','eta_u','xi_rho'))
    nc.createVariable('M2_NudgeCoef','f8',('eta_rho','xi_rho'))    
    nc.variables['M2_NudgeCoef'].long_name = "2D momentum inverse nudging coefficients" 
    nc.variables['M2_NudgeCoef'].units = "day-1" 
    nc.variables['M2_NudgeCoef'].coordinates = "xi_rho eta_rho " 
    nc.variables['M3_NudgeCoef'].long_name = "3D momentum inverse nudging coefficients" 
    nc.variables['M3_NudgeCoef'].units = "day-1" 
    nc.variables['M3_NudgeCoef'].coordinates = "xi_rho eta_rho s_rho"     
    
 # Fill M2 and M3 nudge matrices. Note: these could be separated, but are together for tercity   
    m3 = nc.variables['M3_NudgeCoef'][:]
    m2 = nc.variables['M2_NudgeCoef'][:]
    for i in range(m3.shape[1]):
        for j in range(m3.shape[2]):
            m2[i,j] = ubar_ini
            for k in range(m3.shape[0]):
                m3[k,i,j] = u_ini
    nc.variables['M3_NudgeCoef'][:] = m3
    nc.variables['M2_NudgeCoef'][:] = m2
    
# Create tracer nudging variables
    tracer_NudgeCoef = nc.createVariable('tracer_NudgeCoef','f8',('s_rho', 'eta_rho', 'xi_rho'))
    temp_NudgeCoef = nc.createVariable('temp_NudgeCoef','f8',('s_rho', 'eta_rho', 'xi_rho'))
    salt_NudgeCoef = nc.createVariable('salt_NudgeCoef','f8',('s_rho', 'eta_rho', 'xi_rho'))
    nc.variables['tracer_NudgeCoef'].long_name = "generic tracer inverse nudging coefficients" 
    nc.variables['tracer_NudgeCoef'].units = "day-1" 
    nc.variables['tracer_NudgeCoef'].coordinates = "xi_rho eta_rho s_rho " 
    nc.variables['temp_NudgeCoef'].long_name = "temp inverse nudging coefficients" 
    nc.variables['temp_NudgeCoef'].units = "day-1" 
    nc.variables['temp_NudgeCoef'].coordinates = "xi_rho eta_rho s_rho " 
    nc.variables['salt_NudgeCoef'].long_name = "salt inverse nudging coefficients" 
    nc.variables['salt_NudgeCoef'].units = "day-1" 
    nc.variables['salt_NudgeCoef'].coordinates = "xi_rho eta_rho s_rho " 

# Fill tracer matrices. Should keep these all together, since I'll keep these the same
    temp = nc.variables['temp_NudgeCoef'][:]
    for i in range(temp.shape[1]):
        for j in range(temp.shape[2]):
            for k in range(temp.shape[0]):
                temp[k,i,j] = u_ini
    nc.variables['temp_NudgeCoef'][:] = temp
    nc.variables['salt_NudgeCoef'][:] = temp
    nc.variables['tracer_NudgeCoef'][:] = temp
      
    nc.close()
