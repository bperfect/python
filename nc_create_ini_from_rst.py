# -*- coding: utf-8 -*-
"""
Created on Wed Aug 03 22:13:03 2016

@author: bperfect
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Aug 03 15:26:54 2016

@author: bperfect
"""

import numpy as np
from matplotlib import pyplot as plt
from datetime import datetime
try:
  import netCDF4 as netCDF
except:
  import netCDF3 as netCDF
from scipy.interpolate import griddata

def calcStrat(d, N):
    #Linear stratification, corresponds to N=1e-3 Hz
    #Use dTdz = (N/.0408)^2 (uses ROMS equation of state)
    dTdz = (N/.0408)**2
    return 17+dTdz*d
    #print(np.exp(d/1000))
    #return 12+7.5*np.exp(d/1000)
    
def nc_create_ini_from_rst(filename, grd, N, geostrophic, stratified, rst_file):
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

    # write time and grid information
    nc.createVariable('theta_s', 'f8', ())
    nc.variables['theta_s'].long_name = 'S-coordinate surface control parameter'
    nc.variables['theta_s'][:] = grd.vgrid.theta_s
    
    nc.createVariable('spherical','f8', (),fill_value=0)
    nc.variables['spherical'].long_name = 'Boolean for spherical or cartesian coords'
    
    #ROMS ignores hraw, so I don't need it
    #nc.createVariable('hraw','f8', ('eta_rho','xi_rho'),fill_value=10)
    
    nc.createVariable('theta_b', 'f8', ())
    nc.variables['theta_b'].long_name = 'S-coordinate bottom control parameter'
    nc.variables['theta_b'][:] = grd.vgrid.theta_b
    
    nc.createVariable('Vtransform','f8',(),fill_value=2)
    
    nc.createVariable('Vstretching','f8',(),fill_value=4)

    nc.createVariable('Tcline', 'f8', ())
    nc.variables['Tcline'].long_name = 'S-cordinate surface/bottom layer width'
    nc.variables['Tcline'].units = 'meter'
    nc.variables['Tcline'][:] = grd.vgrid.Tcline

    nc.createVariable('hc', 'f8', ())
    nc.variables['hc'].long_name = 'S-coordinate parameter, critical depth'
    nc.variables['hc'].units = 'meter'
    nc.variables['hc'][:] = grd.vgrid.hc

    nc.createVariable('s_rho', 'f8', ('s_rho'))
    nc.variables['s_rho'].long_name = 'S-coordinate at RHO-points'
    nc.variables['s_rho'].valid_min = '-1'
    nc.variables['s_rho'].valid_max = '0'
    nc.variables['s_rho'].field = 's_rho,scalar'
    nc.variables['s_rho'][:] = grd.vgrid.s_rho

    nc.createVariable('s_w', 'f8', ('s_w'))
    nc.variables['s_w'].long_name = 'S-coordinate at W-points'
    nc.variables['s_w'].valid_min = '-1'
    nc.variables['s_w'].valid_max = '0'
    nc.variables['s_w'].field = 's_w,scalar'
    nc.variables['s_w'][:] = grd.vgrid.s_w

    nc.createVariable('Cs_r', 'f8', ('s_rho'))
    nc.variables['Cs_r'].long_name = 'S-coordinate stretching curves at RHO-points'
    nc.variables['Cs_r'].valid_min = '-1'
    nc.variables['Cs_r'].valid_max = '0'
    nc.variables['Cs_r'].field = 'Cs_r,scalar'
    nc.variables['Cs_r'][:] = grd.vgrid.Cs_r

    nc.createVariable('Cs_w', 'f8', ('s_w'))
    nc.variables['Cs_w'].long_name = 'S-coordinate stretching curves at W-points'
    nc.variables['Cs_w'].valid_min = '-1'
    nc.variables['Cs_w'].valid_max = '0'
    nc.variables['Cs_w'].field = 'Cs_w,scalar'
    nc.variables['Cs_w'][:] = grd.vgrid.Cs_w
    
    #########################################################################
    nc_coarse=netCDF.Dataset(rst_file,'r')
    # u field
    u_coarse = nc_coarse.variables['u'][:]
    u_shape=u_coarse.shape
    u_coarse=u_coarse[u_shape[0]-1,:,:,:]
    # v field
    v_coarse = nc_coarse.variables['v'][:]
    v_shape=v_coarse.shape
    v_coarse=v_coarse[v_shape[0]-1,:,:,:]
    # temp
    temp_coarse = nc_coarse.variables['temp'][:]
    temp_shape=temp_coarse.shape
    temp_coarse=temp_coarse[temp_shape[0]-1,:,:,:]
    # rho
    rho_coarse = nc_coarse.variables['rho'][:]
    rho_shape=rho_coarse.shape
    rho_coarse=rho_coarse[rho_shape[0]-1,:,:,:]
    # zeta
    zeta_coarse = nc_coarse.variables['zeta'][:]
    zeta_shape=zeta_coarse.shape
    zeta_coarse=zeta_coarse[zeta_shape[0]-1,:,:]
    # ubar
    ubar_coarse = nc_coarse.variables['ubar'][:]
    ubar_shape=ubar_coarse.shape
    ubar_coarse=ubar_coarse[ubar_shape[0]-1,:,:]
    # vbar
    vbar_coarse = nc_coarse.variables['vbar'][:]
    vbar_shape=vbar_coarse.shape
    vbar_coarse=vbar_coarse[vbar_shape[0]-1,:,:]

############# Set u ##################################################
    x_coarse = nc_coarse.variables['x_u'][:]
    y_coarse = nc_coarse.variables['y_u'][:]
    Cs_r_coarse = nc_coarse.variables['Cs_r'][:]
    h_coarse = nc_coarse.variables['h'][:]
    
    x_fine = grd.hgrid.x_u
    y_fine = grd.hgrid.y_u
    Cs_r_fine = grd.vgrid.Cs_r
    h_fine = grd.vgrid.h
    # Coarse coordinates
    x_coarse2 = np.repeat(x_coarse[np.newaxis,:,:],Cs_r_coarse.size,axis=0)
    y_coarse2 = np.repeat(y_coarse[np.newaxis,:,:],Cs_r_coarse.size,axis=0)
    z_coarse = np.ones(x_coarse2.shape) #remember: arrays are objects and don't deep copy automatically
    for i in range(z_coarse.shape[1]):
        for j in range(z_coarse.shape[2]):  
            z_coarse[:,i,j] = h_coarse[i,j]*Cs_r_coarse
    # Fine Coordinates    
    x_fine2 = np.repeat(x_fine[np.newaxis,:,:],Cs_r_fine.size,axis=0)
    y_fine2 = np.repeat(y_fine[np.newaxis,:,:],Cs_r_fine.size,axis=0)
    z_fine = np.ones(x_fine2.shape)
    for i in range(z_fine.shape[1]):
        for j in range(z_fine.shape[2]):  
            z_fine[:,i,j,] = h_fine[i,j]*Cs_r_fine
    # Interpolate and assign variables
    nc.createVariable('u','f8',('s_rho','eta_u','xi_u'))
    nc.createVariable('ubar','f8',('eta_u','xi_u'))

    u_fine = nc.variables['u'][:]
    u_fine=griddata(np.stack((np.ndarray.flatten(x_coarse2),np.ndarray.flatten(y_coarse2),np.ndarray.flatten(z_coarse))).T,np.ndarray.flatten(u_coarse),(x_fine2,y_fine2,z_fine),method='nearest')
    #u_fine=griddata(np.stack((np.ndarray.flatten(x_coarse2),np.ndarray.flatten(y_coarse2),np.ndarray.flatten(z_coarse))).T,np.ndarray.flatten(u_coarse),(x_fine,y_fine,z_fine[:,:,14]),method='nearest')
    ubar_fine = nc.variables['ubar'][:]
    ubar_fine = griddata(np.stack((np.ndarray.flatten(x_coarse),np.ndarray.flatten(y_coarse))).T,np.ndarray.flatten(ubar_coarse),(x_fine,y_fine),method='nearest')   
    nc.variables['u'][:] = u_fine
    nc.variables['ubar'][:] = ubar_fine

    
    #plt.matshow(u_fine[14,:,:]) #appears to be a load bearing print call
    #plt.matshow(u_coarse[14,:,:])
  
############# Set v ##################################################
    x_coarse = nc_coarse.variables['x_v'][:]
    y_coarse = nc_coarse.variables['y_v'][:]  
    x_fine = grd.hgrid.x_v
    y_fine = grd.hgrid.y_v

    # Coarse coordinates
    x_coarse2 = np.repeat(x_coarse[np.newaxis,:,:],Cs_r_coarse.size,axis=0)
    y_coarse2 = np.repeat(y_coarse[np.newaxis,:,:],Cs_r_coarse.size,axis=0)
    z_coarse = np.ones(x_coarse2.shape) #remember: arrays are objects and don't deep copy automatically
    for i in range(z_coarse.shape[1]):
        for j in range(z_coarse.shape[2]):  
            z_coarse[:,i,j] = h_coarse[i,j]*Cs_r_coarse
    # Fine Coordinates    
    x_fine2 = np.repeat(x_fine[np.newaxis,:,:],Cs_r_fine.size,axis=0)
    y_fine2 = np.repeat(y_fine[np.newaxis,:,:],Cs_r_fine.size,axis=0)
    z_fine = np.ones(x_fine2.shape)
    for i in range(z_fine.shape[1]):
        for j in range(z_fine.shape[2]):  
            z_fine[:,i,j,] = h_fine[i,j]*Cs_r_fine
    # Interpolate and assign variables
    nc.createVariable('v','f8',('s_rho','eta_v','xi_v'))
    nc.createVariable('vbar','f8',('eta_v','xi_v'))
    v_fine = nc.variables['v'][:]
    vbar_fine = nc.variables['vbar'][:]
    v_fine=griddata(np.stack((np.ndarray.flatten(x_coarse2),np.ndarray.flatten(y_coarse2),np.ndarray.flatten(z_coarse))).T,np.ndarray.flatten(v_coarse),(x_fine2,y_fine2,z_fine),method='nearest')
    vbar_fine = griddata(np.stack((np.ndarray.flatten(x_coarse),np.ndarray.flatten(y_coarse))).T,np.ndarray.flatten(vbar_coarse),(x_fine,y_fine),method='nearest')
    nc.variables['v'][:] = v_fine
    nc.variables['vbar'][:] = vbar_fine
    #plt.matshow(vbar_fine)
    #plt.matshow(vbar_coarse)
    #plt.matshow(v_fine[:,:,14])    
    #plt.matshow(v_coarse[:,:,14])  
    
############# Set rho,temp,zeta ##################################################
    x_coarse = nc_coarse.variables['x_rho'][:]
    y_coarse = nc_coarse.variables['y_rho'][:]
    x_fine = grd.hgrid.x_rho
    y_fine = grd.hgrid.y_rho
    # Coarse coordinates
    x_coarse2 = np.repeat(x_coarse[np.newaxis,:,:],Cs_r_coarse.size,axis=0)
    y_coarse2 = np.repeat(y_coarse[np.newaxis,:,:],Cs_r_coarse.size,axis=0)
    z_coarse = np.ones(x_coarse2.shape) #remember: arrays are objects and don't deep copy automatically
    for i in range(z_coarse.shape[1]):
        for j in range(z_coarse.shape[2]):  
            z_coarse[:,i,j] = h_coarse[i,j]*Cs_r_coarse
    # Fine Coordinates    
    x_fine2 = np.repeat(x_fine[np.newaxis,:,:],Cs_r_fine.size,axis=0)
    y_fine2 = np.repeat(y_fine[np.newaxis,:,:],Cs_r_fine.size,axis=0)
    z_fine = np.ones(x_fine2.shape)
    for i in range(z_fine.shape[1]):
        for j in range(z_fine.shape[2]):  
            z_fine[:,i,j,] = h_fine[i,j]*Cs_r_fine
    # Interpolate and assign variables
    nc.createVariable('temp','f8',('s_rho','eta_rho','xi_rho'))
    temp_fine = nc.variables['temp'][:]
    nc.createVariable('zeta','f8',('eta_rho','xi_rho'))
    zeta_fine = nc.variables['zeta'][:]  
    nc.createVariable('rho','f8',('s_rho','eta_rho','xi_rho'))
    zeta_fine = nc.variables['zeta'][:] 
    temp_fine=griddata(np.stack((np.ndarray.flatten(x_coarse2),np.ndarray.flatten(y_coarse2),np.ndarray.flatten(z_coarse))).T,np.ndarray.flatten(temp_coarse),(x_fine2,y_fine2,z_fine),method='nearest')
    zeta_fine = griddata(np.stack((np.ndarray.flatten(x_coarse),np.ndarray.flatten(y_coarse))).T,np.ndarray.flatten(zeta_coarse),(x_fine,y_fine),method='nearest')
    rho_fine=griddata(np.stack((np.ndarray.flatten(x_coarse2),np.ndarray.flatten(y_coarse2),np.ndarray.flatten(z_coarse))).T,np.ndarray.flatten(rho_coarse),(x_fine2,y_fine2,z_fine),method='nearest')
    
    nc.variables['temp'][:] = temp_fine
    nc.variables['zeta'][:] = zeta_fine
    nc.variables['rho'][:] = rho_fine

#############################################################
    

    nc.createVariable('h', 'f8', ('eta_rho', 'xi_rho'))
    nc.variables['h'].long_name = 'bathymetry at RHO-points'
    nc.variables['h'].units ='meter'
    nc.variables['h'].coordinates = 'lon_rho y_rho'
    nc.variables['h'].field = 'bath, scalar'
    nc.variables['h'][:] = grd.vgrid.h

    nc.createVariable('x_rho', 'f8', ('eta_rho', 'xi_rho'))
    nc.variables['x_rho'].long_name = 'longitude of RHO-points'
    nc.variables['x_rho'].units = 'degree_east'
    nc.variables['x_rho'].field = 'x_rho, scalar'
    nc.variables['x_rho'][:] = grd.hgrid.x_rho

    nc.createVariable('y_rho', 'f8', ('eta_rho', 'xi_rho'))
    nc.variables['y_rho'].long_name = 'yitude of RHO-points'
    nc.variables['y_rho'].units = 'degree_north'
    nc.variables['y_rho'].field = 'y_rho, scalar'
    nc.variables['y_rho'][:] = grd.hgrid.y_rho

    nc.createVariable('x_u', 'f8', ('eta_u', 'xi_u'))
    nc.variables['x_u'].long_name = 'longitude of U-points'
    nc.variables['x_u'].units = 'degree_east'
    nc.variables['x_u'].field = 'x_u, scalar'
    nc.variables['x_u'][:] = grd.hgrid.x_u

    nc.createVariable('y_u', 'f8', ('eta_u', 'xi_u'))
    nc.variables['y_u'].long_name = 'latitude of U-points'
    nc.variables['y_u'].units = 'degree_north'
    nc.variables['y_u'].field = 'y_u, scalar'
    nc.variables['y_u'][:] = grd.hgrid.y_u

    nc.createVariable('x_v', 'f8', ('eta_v', 'xi_v'))
    nc.variables['x_v'].long_name = 'longitude of V-points'
    nc.variables['x_v'].units = 'degree_east'
    nc.variables['x_v'].field = 'x_v, scalar'
    nc.variables['x_v'][:] = grd.hgrid.x_v

    nc.createVariable('y_v', 'f8', ('eta_v', 'xi_v'))
    nc.variables['y_v'].long_name = 'latitude of V-points'
    nc.variables['y_v'].units = 'degree_north'
    nc.variables['y_v'].field = 'y_v, scalar'
    nc.variables['y_v'][:] = grd.hgrid.y_v

    nc.createVariable('ocean_time', 'f8', ('ocean_time'))
 #   nc.variables['ocean_time'].long_name = ocean_time.long_name
    nc.variables['ocean_time'].units = 'days'
 #   try:
 #       nc.variables['ocean_time'].field = ocean_time.field
 #   except:
 #       nc.variables['ocean_time'].field = 'ocean_time, unlimited'
    nc.variables['ocean_time'][:] = nc_coarse.variables['ocean_time'][:][-1]/24.0/3600.0
    nc.close()
