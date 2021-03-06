# -*- coding: utf-8 -*-
"""
Created on Tue Dec 13 16:26:54 2016

@author: bperfect
"""

import hgrid
import vgrid
from grid import *
import numpy as np
from nc_create_grid import *
from nc_create_ini import *
from nc_create_ini_from_rst import *

use_rst = True
if use_rst==True:
    grid_x = 240 #210
    grid_y = 160 #140
    grid_z = 70 #60
    rst_file = '../../startup/ocean_rst_coarse.nc'
    ini_file = '../../startup/ocean_ini_startup2.nc'
    grd_file = '../../startup/ocean_grd_startup2.nc'
else:
    ini_file = '../../startup/ocean_ini_startup.nc'
    grd_file = '../../startup/ocean_grd_startup.nc'
    grid_x = 81 #210
    grid_y = 54 #140
    grid_z = 35 #60

u_ini = 0.0
v_ini = 0.0
ubar_ini = 0.0
vbar_ini = 0.0
u_clm = 0.1
v_clm = 0.0
ubar_clm = 0.1
vbar_clm = 0.0
temp_clm=13


f = 0.000001
temp_ini=13
N=0.002

tidal = False
geostrophic = True
stratified = True
refined = True

#Seamount geometry
el = 100000       #y-domain width
sigma = 7000    #seamount width
height = 1000    #seamount height
depth = 3000    #depth of ocean
xl = 150000     #x-domain width

Fr = u_clm/N/height
print('Froude number is ' + str(Fr))
Ro = u_clm/f/sigma
print('Rossby number is ' + str(Ro))

theta_b = 3
theta_s = 0.65
clm_file = 'ocean_clm.nc'
nud_file = 'ocean_nud.nc'

#Follow the same code, for the refined grid
if refined:
    grid2_x = 20
    grid2_y = 20
    grid2_z = 30
    refinement_factor = 3.
    

#secify bathymetry
def calcHeight(horGrid1):
    return depth-height*np.exp(-(horGrid1.x_rho-xl/4)**2/sigma**2-(horGrid1.y_rho-el/2)**2/sigma**2)

#Create the original array 
grid1 = np.mgrid[-1:grid_y+1:(grid_y+3)*1j,-1:grid_x+1:(grid_x+3)*1j]
xgrid=grid1[0]
ygrid=grid1[1]

#generate horizontal grid object
horGrid1 = CGrid(np.around(grid1[1])*xl/grid_x,np.around(grid1[0])*el/grid_y)
h=calcHeight(horGrid1) 
#generate vertical grid object  
s4 = s_coordinate_4(h, theta_b,theta_s, 100, grid_z) #theta_b, theta_s, Tcline, N,
#combine the horizontal and vertical grid to make a ROMS_Grid object
grd = ROMS_Grid("Seamount grid", horGrid1, s4)
#write grid object to a netCDF file
nc_create_grid(grd_file, grd, f, True)
#write initial condition netCDF file

#rst_file = 'ocean_rst.nc'
if use_rst == True:
    u_fine = nc_create_ini_from_rst(ini_file,grd,N,geostrophic,stratified,rst_file)
else:
    nc_create_ini(ini_file, grd, u_ini, v_ini, ubar_ini, vbar_ini, temp_ini, N, geostrophic, stratified)

#nc_create_clm(clm_file, grd, u_clm, v_clm, ubar_clm, vbar_clm, temp_clm, stratified)
#nc_create_nudge(nud_file, grd, u_ini, v_ini, ubar_ini, vbar_ini, temp_ini, stratified)
#clm and nudge untested with nesting -10 Oct 16

#Warning: floating point arithmetic errors will arise without use of np.around (rounding)
#Other refinement factors may require the mpgrid command to change as well (untested)
if refined:
    grid2 = np.mgrid[grid_y/2-grid2_y/2-1./refinement_factor:grid_y/2+grid2_y/2+1./refinement_factor:(3+refinement_factor*grid2_y)*1j,
                 grid_x/2-grid2_x/2-1./refinement_factor:grid_x/2+grid2_x/2+1./refinement_factor:(3+refinement_factor*grid2_x)*1j]
    horGrid2 = CGrid(np.around(refinement_factor*grid2[1])*xl/grid_x/refinement_factor, 
                 np.around(grid2[0]*refinement_factor)*el/grid_y/refinement_factor)
    h_nest=calcHeight(horGrid2)
    s4_nest = s_coordinate_4(h_nest, theta_b, theta_s, 100, grid2_z)
    grd_nest = ROMS_Grid("nested grid", horGrid2, s4_nest)
    nc_create_grid('ocean_grd_nest.nc', grd_nest, f, True)
    #add in magic values necessary to create the contact file with contact.m
    nc = netCDF.Dataset('ocean_grd_nest.nc','a')
    nc.type = "GRID file"
    nc.title = "Refined grid area"
    nc.parent_grid = "ocean_grd.nc"
    #Aligned for the center of the grid
    nc.parent_Imin = grid_x/2-grid2_x/2
    nc.parent_Imax = grid_x/2+grid2_x/2
    nc.parent_Jmin = grid_y/2-grid2_y/2
    nc.parent_Jmax = grid_y/2+grid2_y/2
    nc.refine_factor = 3
    nc.history = "Magic value file trying to do a grid refinement run"
    nc.close()
    nc_create_ini('ocean_ini_nest.nc', grd_nest, u_ini, v_ini, ubar_ini, vbar_ini, temp_ini, geostrophic, stratified)
    nc = netCDF.Dataset('ocean_ini_nest.nc','a')
    

