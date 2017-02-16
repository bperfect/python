# -*- coding: utf-8 -*-
"""
Created on Tue Aug 02 22:59:51 2016

@author: bperfect
"""

import hgrid
import vgrid
from grid import *
import numpy as np
from nc_create_grid import *
from nc_create_ini import *

grid_x = 120
grid_y = 120
grid_z = 30
u_ini = 0.5
v_ini = 0
ubar_ini = 0.5
vbar_ini = 0
f = 0.0001
temp_ini=13
geostrophic = True

a=np.mgrid[0:32000:123j,0:320000:123j] #increment this value by 3 from what you want in the .in file
#generate horizontal grid
tempGrid = CGrid(a[0],a[1])

#psecify bathymetry
def calcHeight(tempGrid):
    return 5000-800*np.exp(-(tempGrid.x_rho-320000/2)**2/20000**2-(tempGrid.y_rho-320000/2)**2/20000**2)
 
h=calcHeight(tempGrid) 
#generate vertical grid  
s4 = s_coordinate_4(h, 6.5, 2, 100, grid_z) #theta_b, theta_s, Tcline, N,
#generate grid object
grd = ROMS_Grid("bradGrid", tempGrid, s4)
#write grid netCDF file
nc_create_grid('ocean_grd.nc', grd, f, True)
#write initial condition netCDF file


nc_create_ini('ocean_ini.nc', grd, u_ini, v_ini, ubar_ini, vbar_ini, temp_ini, geostrophic)
# need to solve an 'inconsistency between theta_s and theta_b






