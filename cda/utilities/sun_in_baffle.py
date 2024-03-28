#!/usr/bin/env python
#coding : utf-8

###############################################################################
#
# *** sun_in_baffle.py ***
#
# DESCRIPTION
# The following code defines the function 'sun_in_baffle', which is used to 
# check if, at a given time and position of the scanner of PHEBUS, the 
# direct sunlight might enter the baffle (it must not happen). 
# Therefore, the function checks if the Sun is at more than 92° away from 
# the line of sight of PHEBUS (sufficient condition) and if not, it checks
# if the Sun is occulted by Mercury or the radiator of MPO.
#
# Version 1.0
# Marco Merusi (marco.merusi@latmos.ipsl.fr), LATMOS (Guyancourt, France)
#
###############################################################################
#### IMPORT LIBRARIES #########################################################

# import spiceypy
if "spy" not in dir() : 
    import spiceypy as spy

# import numpy and linalg
if "np" not in dir() :
    import numpy as np
if "lal" not in dir() :
    from numpy import linalg as lal

# Manage directories
import sys
import os
# Get the parent directory
parent_dir = os.path.dirname(os.path.realpath(__file__))

# Add the parent directory to sys.path
sys.path.append(parent_dir)

# import useful functions
from phestar_functions_spice import spos2vzecl
from phestar_functions import angular_sep

#### START FUNCTIONS ##########################################################

### sun_in_baffle
def sun_in_baffle(time, scanpos, mpo_sun_ecl, mpo_mer_ecl) :
    '''
    This function checks if the direct sunlight can hazardously enter the baffle
    of PHEBUS. First it verifies if the Sun is at least 92° away from the line 
    of sight of PHEBUS or if it's behind Mercury with respect to MPO. If not, 
    it checks if the Sun is behind the radiator of MPO according to some angles 
    and scanner positions.

    Parameters
    ----------
    time : (float) ephemeris time of the observation, in seconds past J2000.
    scanpos : (integer) scanner position of PHEBUS, from 100 to 3990.
    mpo_sun_ecl : (float array) 3-element array representing the vector 
        MPO -> Sun in ecliptic reference frame.
    mpo_mer_ecl : (float array) 3-element array representing the vector 
        MPO -> Mercury in ecliptic reference frame.

    Returns
    -------
    If the Sun is not a hazard, the observation is feasible and the function 
    returns 0, otherwise it returns 1.

    '''    
    # If PHEBUS is in parking position, there is no problem.
    if (0 <= scanpos < 100) or (3990 < scanpos <= 4096) :
        return(0)
    
    # Check the two sufficient conditions: the angle with Z_FOV is greater than
    # 90° (92°, for safety), or the Sun is behind Mercury.
    
    # Compute the line of sight of PHEBUS
    phebfov = spos2vzecl(scanpos, time)
    vz_ecl = phebfov[0]
    
    # Angle > 92° (1.6 rad)?
    if angular_sep(vz_ecl,mpo_sun_ecl) >= 1.6 :
        return(0)
    
    # Sun completely behind Mercury?
    
    r_mpo_sun = lal.norm(mpo_sun_ecl)
    r_ang_sun = np.arctan2(700000,r_mpo_sun) 
    r_mpo_mer = lal.norm(mpo_mer_ecl)
    r_ang_mer = np.arctan2(2440,r_mpo_mer)
    ang_sun_mpo_mer = angular_sep(mpo_sun_ecl, mpo_mer_ecl)
    if (r_mpo_sun > r_mpo_mer) and (np.degrees(ang_sun_mpo_mer) + np.degrees(r_ang_sun) < np.degrees(r_ang_mer)) :
        return(0)
    
    # If the function didn't return with the previous conditions, it means that
    # the Sun illuminates directly MPO and the angle with the LoS of PHEBUS is
    # less than 90°. Therefore, the observation is possible only if the Sun is 
    # occulted by the radiator of MPO.
        
    # Compute the vector MPO -> Sun in MPO spacecraft RF by changing base
    # Generate the matrix of change of basis
    matrix_change = spy.pxform("ECLIPJ2000", "MPO_SPACECRAFT", time)
    # Change basis and compute the new vector
    mpo_sun_mpo = np.dot(matrix_change, mpo_sun_ecl)
    
    # Convert MPO -> Sun from MPO spacecraft to URF
    mpo_sun_urf = [-mpo_sun_mpo[2], mpo_sun_mpo[1], mpo_sun_mpo[0]]
    
    # If the Sun is on PHEBUS' side of the radiator, PHEBUS can't observe.
    # In that case, the Yurf component of the Sun is negative.
    if mpo_sun_urf[1] <= 0 :
        return(1)
    
    # If the Sun is at positive values of Yurf, more checks are needed.
    
    scanang = scanpos * 2*np.pi / 4096 # absolute scanner angle in degrees 
    # Compute the angle of the scanner with respect to the axis Zurf (plane Zurf,Xurf).
    scanang_z = scanang + np.pi/4
    # Compute the position of the baffle aperture in the plane (Zurf, Xurf)
    baffle_urf_zx_mm = [246+246*np.cos(scanang_z), 246+246*np.sin(scanang_z)]
    # convert it from mm to km
    baffle_urf_zx_km = [(246+246*np.cos(scanang_z)) / (1e-6), (246+246*np.sin(scanang_z)) / (1e-6)]
    
    # Position of the sun in the plane (Zurf, Xurf)
    sun_urf_zx = [mpo_sun_urf[2], mpo_sun_urf[0]]
    
    # Compute the position of the Sun in the plane (Zurf, Xurf) centered on the baffle aperture
    baf_sun_urf_zx = np.zeros(2)
    baf_sun_urf_zx[0] = sun_urf_zx[0] - baffle_urf_zx_km[0]
    baf_sun_urf_zx[1] = sun_urf_zx[1] - baffle_urf_zx_km[1]
    
    # Distance baffle aperture - Sun in the plane (Zurf, Xurf)
    r_baf_sun_zx = np.sqrt(baf_sun_urf_zx[0]**2 + baf_sun_urf_zx[1]**2)
    # Unit components of the vector baf_sun_urf_zx
    baf_sun_urf_z = baf_sun_urf_zx[0] / r_baf_sun_zx
    baf_sun_urf_x = baf_sun_urf_zx[1] / r_baf_sun_zx
    
    # Now divide the radiator in 4 rectangular quadrants with the baffle aperture
    # at the origin. According to the angle between the Sun and the 
    # Zurf or Xurf axis of the baffle:
    #    1. Compute the two sides of the corresponding quadrant and the diagonal;
    #    2. Compute the angle between one side and the diagonal;
    #    3. Compute the length of the segment going from the baffle to the edge 
    #       of the radiator in the direction of the Sun, through simple trigonometry.
    
    # The sides of the quadrants are computed knowing the coordinates of the baffle
    # and the distances from the real origin of the URF system to the sides of the radiator.
    # Distance URF real origin -> side at -Zurf = 1220 mm
    # Distance URF real origin -> side at +Zurf = 2412 mm
    # Distance URF real origin -> side at +Xurf = 1703 mm
    # Distance URF real origin -> side at -Xurf = 0 mm
    
    # In the following, 'pz' and 'nz' correspond to '+Z' and '-Z'. Same goes for 'px' and 'nx'
    
    # First or third quadrant
    if baf_sun_urf_z * baf_sun_urf_x > 0 :
        
        # First quadrant (0° < angle < 90°)
        if baf_sun_urf_x > 0 :
            
            dist_baf_pz_side = 2412 - baffle_urf_zx_mm[0]
            dist_baf_px_side = 1703 - baffle_urf_zx_mm[1]
            dist_corn_pzpx = np.sqrt(dist_baf_pz_side**2 + dist_baf_px_side**2)
            
            ang_pz_baf_cornpzpx = np.arccos(dist_baf_pz_side/dist_corn_pzpx)
            ang_pz_baf_sun = np.arccos(baf_sun_urf_z)
            
            if ang_pz_baf_sun <= ang_pz_baf_cornpzpx :
                
                angmod = ang_pz_baf_sun
                dist_baf_sidesun = dist_baf_pz_side / np.cos(angmod)
            
            if ang_pz_baf_sun > ang_pz_baf_cornpzpx :
                
                angmod = np.pi/2 - ang_pz_baf_sun
                dist_baf_sidesun = dist_baf_px_side / np.cos(angmod)
            
        # Third quadrant (180° < angle < 270°)
        if baf_sun_urf_x < 0 :
            
            dist_baf_nz_side = 1220 + baffle_urf_zx_mm[0]
            dist_baf_nx_side = baffle_urf_zx_mm[1]
            dist_corn_nznx = np.sqrt(dist_baf_nz_side**2 + dist_baf_nx_side**2)
            
            ang_nz_baf_cornnznx = np.arccos(dist_baf_nz_side / dist_corn_nznx)
            ang_nz_baf_sun = np.arccos(abs(baf_sun_urf_z))
            
            if ang_nz_baf_sun <= ang_nz_baf_cornnznx :
                
                angmod = ang_nz_baf_sun
                dist_baf_sidesun = dist_baf_nz_side / np.cos(angmod)
            
            if ang_nz_baf_sun > ang_nz_baf_cornnznx :
                
                angmod = np.pi/2 - ang_nz_baf_sun
                dist_baf_sidesun = dist_baf_nx_side / np.cos(angmod)
    
    # Second or fourth quadrant
    if baf_sun_urf_z * baf_sun_urf_x < 0 :
        
        # Second quadrant
        if baf_sun_urf_x > 0 :
            
            dist_baf_nz_side = 1220 + baffle_urf_zx_mm[0]
            dist_baf_px_side = 1703 - baffle_urf_zx_mm[1]
            dist_corn_nzpx = np.sqrt(dist_baf_nz_side**2 + dist_baf_px_side**2)
            
            ang_nz_baf_cornnzpx = np.arccos(dist_baf_nz_side / dist_corn_nzpx)
            ang_nz_baf_sun = np.arccos(abs(baf_sun_urf_z))
            
            if ang_nz_baf_sun <= ang_nz_baf_cornnzpx :
                
                angmod = ang_nz_baf_sun
                dist_baf_sidesun = dist_baf_nz_side / np.cos(angmod)
            
            if ang_nz_baf_sun > ang_nz_baf_cornnzpx :
                
                angmod = np.pi/2 - ang_nz_baf_sun
                dist_baf_sidesun = dist_baf_px_side / np.cos(angmod)
        
        # Fourth quadrant
        if baf_sun_urf_x < 0 :
            
            dist_baf_pz_side = 2412 - baffle_urf_zx_mm[0]
            dist_baf_nx_side = baffle_urf_zx_mm[1]
            dist_corn_pznx = np.sqrt(dist_baf_pz_side**2 + dist_baf_nx_side**2)
            
            ang_pz_baf_cornpznx = np.arccos(dist_baf_pz_side / dist_corn_pznx)
            ang_pz_baf_sun = np.arccos(baf_sun_urf_z)
            
            if ang_pz_baf_sun <= ang_pz_baf_cornpznx :
                
                angmod = ang_pz_baf_sun
                dist_baf_sidesun = dist_baf_pz_side / np.cos(angmod)
            
            if ang_pz_baf_sun > ang_pz_baf_cornpznx :
                
                angmod = np.pi/2 - ang_pz_baf_sun
                dist_baf_sidesun = dist_baf_nx_side / np.cos(angmod)
    
    # Sun aligned with one of the axes
    if baf_sun_urf_z == 1 and baf_sun_urf_x == 0 :
        dist_baf_sidesun = 2412 - baffle_urf_zx_mm[0]
    if baf_sun_urf_z == 0 and baf_sun_urf_x == 1 :
        dist_baf_sidesun = 1703 - baffle_urf_zx_mm[1]
    if baf_sun_urf_z == -1 and baf_sun_urf_x == 0 :
        dist_baf_sidesun = 1220 + baffle_urf_zx_mm[0]
    if baf_sun_urf_z == 0 and baf_sun_urf_x == -1 :
        dist_baf_sidesun = 0
        
    
    height_baffle = 110 # height of the baffle from the radiator in mm
    
    # Compute the angle between the direction MPO->Sun in urf system and the +Yurf axis
    ang_sun_mpo_yurf = angular_sep(mpo_sun_urf, [0,1,0])
    
    # Compute the angle at which the baffle sees the radiator edge in direction of the Sun
    ang_baffle_to_side = np.arctan2(dist_baf_sidesun, height_baffle)
    
    margin_angle = np.radians(5) # margin of degrees
    
    # If the angle of the Sun with +Yurf (plus a small margin) is smaller than
    # tha angle at which the baffle sees the edge of the radiator, the Sun is not a hazard.
    if ang_baffle_to_side > ang_sun_mpo_yurf + margin_angle :
        return(0)
    else :
        return(1)  
    
    
#### END FUNCTIONS ############################################################
