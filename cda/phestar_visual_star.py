#!/usr/bin/env python
#coding : utf-8

###############################################################################
#
# *** phestar_visual_star.py ***
#
# DESCRIPTION
# The following code defines the function 'phestar_visual_star', which is used
# in the PHESTAR software. The aim of this function is to generate a plot in 
# 2 dimensions (longitude on X from 0° to 360°, latitude on Y from -90° to 90°) 
# representing the 'celestial sphere' seen from MPO in ecliptic reference frame.
# According to the parameters given in input (inserted by the user in the main 
# program), the plot shows the stars from a catalog, the Sun and some planets 
# (Mercury and its disk, Venus, the Earth) with their names, the projection of 
# the line of sight of PHEBUS and the four corners of its slit. In addition, if 
# a star is arbitrarily close to the line of sight, its name is shown in the plot.
# For any details, visit my Github (https://github.com/marcomeru/phestar_phebus)!
#
# Version 1.0
# Marco Merusi (marco.merusi@latmos.ipsl.fr), LATMOS (Guyancourt, France)
#
###############################################################################

#### IMPORT LIBRARIES #########################################################

# import SPICE
if "spy" not in dir() :
    import spiceypy as spy
 
# numpy
if "np" not in dir() :
    import numpy as np
if "lal" not in dir() :
    from numpy import linalg as lal

# import plot
if "plt" not in dir() :
    import matplotlib.pyplot as plt

# import pandas
if "pd" not in dir() :
    import pandas as pd

# pysimpleGUI
if "sg" not in dir() :
    import PySimpleGUI as sg

# Manage directories
import sys
import os
# Get the parent directory
parent_dir = os.path.dirname(os.path.realpath(__file__))

# Add the parent directory to sys.path
sys.path.append(parent_dir)

# Useful functions
from utilities.phestar_functions import angular_sep, circumf_3d_vec, convert3to2deg
from utilities.phestar_functions_spice import pos_body_rf, spos2vzecl
from utilities.sun_in_baffle import sun_in_baffle 


#### START FUNCTION ###########################################################

def phestar_visual_star(kernels_path, star_data_file, etdate, scannerposs, slitx, slity, max_ang, showsun, showmer, showven, showear) :
    '''
    This function generates a plot in 2 dimensions (longitude, latitude) 
    representing the 'celestial sphere' seen from MPO in ecliptic reference frame.
    According to the parameters given in input, the plot shows the stars from 
    a catalog, the Sun and some planets the projection of the line of sight of 
    PHEBUS and the four corners of its slit. If a star is arbitrarily close to 
    the line of sight, its name is shown in the plot.
    
    Parameters
    ----------
    kernels_path : (string) complete path to the kernel list file (bc_plan.tm).
    star_data_file : (integer array) indices of the stars from the catalog to be
        plotted in the sky map.
    etdate : (string) complete date and time in the form 'YYYY-MM-DDThh:mm:ss'.
    scannerposs : (integer) position of the scanner of PHEBUS.
    slitx : (float) full angular width of the slit in degrees.
    slity : (float) full angular height of the slit in degrees.
    max_ang : (float) angular radius from the line of sight of PHEBUS within 
        which any star's name is shown.
    showsun : (bool) should the Sun be plotted in the map?
    showmer : (bool) should Mercury be plotted in the map?
    showven : (bool) should Venus be plotted in the map?
    showear : (bool) should the Earth be plotted in the map?
    
    Returns
    -------
    It returns 0 at the end.
    
    ''' 
    ###########################################################################
    #### Open and load SPICE kernels ##########################################
    
    # Load the spice kernels if necessary
    try :
        _ = spy.str2et("2025-01-01T20:00:00")
    except :
        spy.furnsh(kernels_path)
    
    ###########################################################################
    
    #### Process the star catalog to find the ideal stars to observe ##########
    # Open the catalog of stars and take only those selected by the user
    starcat = pd.read_csv(parent_dir+r"\utilities\starcat_phebus.csv")
    nbst = len(star_data_file) 
    starnm = []
    starhr = np.zeros(nbst)
    vec_ecli = np.zeros((3,nbst))
    for uu in range(nbst) :
        ind = star_data_file[uu]
        starnm.append(starcat['Star name'][ind])
        starhr[uu] = starcat['HR'][ind]
        vec_ecli[0,uu] = starcat['Star X position'][ind]
        vec_ecli[1,uu] = starcat['Star Y position'][ind]
        vec_ecli[2,uu] = starcat['Star Z position'][ind]
      
    # Convert coordinates of each star to (longitude, latitude)
    lonlatstars = np.zeros((2,nbst))
    for k in range(nbst) :
        lonlatstars[0:2,k] = convert3to2deg(vec_ecli[0:3,k])
    
    #### Convert the date/time to ephemeris time
    ettime_ok = spy.str2et(etdate)
        
    #### Projection of the slit on the sky ####################################
    
    #### Compute the FOV components of PHEBUS in ecliptic RF
    phebfov = spos2vzecl(scannerposs, ettime_ok)
    vz_ecl = phebfov[0]
    vy_ecl = phebfov[1]
    vx_ecl = phebfov[2]
    
    # Compute the longitude and latitude of the LoS (Z_FOV) from MPO
    lonlatfov = convert3to2deg(vz_ecl)
    
    #### Set the horizontal and vertical "radii" of the slit
    xslit = np.radians(slitx / 2) # Half width of the slit 
    yslit = np.radians(slity / 2) # Half height of the slit
    
    #### Identify the corners of the slit and obtain their latitude and longitude
    # Top right corner
    vis1 = (np.cos(yslit) * vz_ecl + np.sin(yslit) * vy_ecl) * np.cos(xslit) + np.sin(xslit) * vx_ecl
    lonlat1 = convert3to2deg(vis1)
    # Bottom right corner
    vis2 = (np.cos(yslit) * vz_ecl - np.sin(yslit) * vy_ecl) * np.cos(xslit) + np.sin(xslit) * vx_ecl
    lonlat2 = convert3to2deg(vis2)
    # Bottom left corner
    vis3 = (np.cos(yslit) * vz_ecl - np.sin(yslit) * vy_ecl) * np.cos(xslit) - np.sin(xslit) * vx_ecl
    lonlat3 = convert3to2deg(vis3)
    # Top left corner
    vis4 = (np.cos(yslit) * vz_ecl + np.sin(yslit) * vy_ecl) * np.cos(xslit) - np.sin(xslit) * vx_ecl
    lonlat4 = convert3to2deg(vis4)

    
    #### Compute positions of Sun and planets #################################
    # Compute vector MPO -> Mercury in ecliptic RF
    mpo_mer_ecl = pos_body_rf("-121", "1", "ECLIPJ2000", ettime_ok)
    # Compute distance MPO-Mercury
    r_mpo_mer = lal.norm(mpo_mer_ecl)
    # Compute angular radius of Mercury seen from MPO
    r_ang_mer = np.arctan2(2440,r_mpo_mer)
    # Compute vectors MPO -> Sun, Venus, Earth in ecliptic RF
    mpo_sun_ecl = pos_body_rf("-121", "10", "ECLIPJ2000", ettime_ok)
    mpo_ven_ecl = pos_body_rf("-121", "2", "ECLIPJ2000", ettime_ok)
    mpo_ear_ecl = pos_body_rf("-121", "3", "ECLIPJ2000", ettime_ok)
    # Compute the distance MPO-Sun
    r_mpo_sun = lal.norm(mpo_sun_ecl)
    # Compute angular radius of the Sun seen from MPO
    r_ang_sun = np.arctan2(700000,r_mpo_sun)
    
    ###########################################################################
       
    #### START PLOT ###########################################################
    
    ax = plt.axes()
    ax.grid()
    # Plot PHEBUS line of sight (Z_FOV) and the corners of the slit
    ax.scatter(lonlatfov[0], lonlatfov[1], c = 'black', label='PHEBUS line of sight', s=5)
    ax.scatter(lonlat1[0], lonlat1[1], label='Slit (Top right)', color='red')
    ax.scatter(lonlat2[0], lonlat2[1], label='Slit (Bottom right)', color='green')
    ax.scatter(lonlat3[0], lonlat3[1], label='Slit (Bottom left)', color='blue')
    ax.scatter(lonlat4[0], lonlat4[1], label='Slit (Top left)', color='goldenrod')
    
    # Plot stars
    vis = 0
    occ = 0
    for j in range(nbst) :
        # Use different colors for the stars if they are visible or occulted
        if angular_sep([vec_ecli[0,j],vec_ecli[1,j],vec_ecli[2,j]], mpo_mer_ecl) > r_ang_mer and vis == 0 :
            ax.scatter(lonlatstars[0,j], lonlatstars[1,j], color='saddlebrown', marker="*", label='Visible star')
            vis = 1
        elif angular_sep([vec_ecli[0,j],vec_ecli[1,j],vec_ecli[2,j]], mpo_mer_ecl) > r_ang_mer and vis != 0 :
            ax.scatter(lonlatstars[0,j], lonlatstars[1,j], color='saddlebrown', marker="*")
        if angular_sep([vec_ecli[0,j],vec_ecli[1,j],vec_ecli[2,j]], mpo_mer_ecl) <= r_ang_mer and occ == 0 :
            ax.scatter(lonlatstars[0,j], lonlatstars[1,j], color='grey', marker="*", label='Occulted star')
            occ = 1
        elif angular_sep([vec_ecli[0,j],vec_ecli[1,j],vec_ecli[2,j]], mpo_mer_ecl) <= r_ang_mer and occ != 0 :
            ax.scatter(lonlatstars[0,j], lonlatstars[1,j], color='grey', marker="*")
        # If some star is close to the LoS, show its name
        if angular_sep([vz_ecl[0],vz_ecl[1],vz_ecl[2]], [vec_ecli[0,j],vec_ecli[1,j],vec_ecli[2,j]]) * 180/np.pi <= max_ang :
            ax.text(*lonlatstars[0:2,j], starnm[j])
        
    # If Sun/planets have to be plotted, check if they are occulted by Mercury
    # and if so, write it with their names
    if showsun == 1 :
        mpo_sun_ecl_lonlat = convert3to2deg(mpo_sun_ecl)
        ax.scatter(*mpo_sun_ecl_lonlat, color='yellow', edgecolors='black')
        ax.text(mpo_sun_ecl_lonlat[0]+0.5, mpo_sun_ecl_lonlat[1]+0.5, "Sun")
        if angular_sep(mpo_sun_ecl, mpo_mer_ecl)+r_ang_sun < r_ang_mer :
            ax.text(mpo_sun_ecl_lonlat[0]+0.5, mpo_sun_ecl_lonlat[1]+0.5, "---")    
           
    if showmer == 1 :
        mpo_mer_ecl_lonlat = convert3to2deg(mpo_mer_ecl)
        ax.scatter(*mpo_mer_ecl_lonlat, color='brown')
        ax.text(mpo_mer_ecl_lonlat[0]+0.5, mpo_mer_ecl_lonlat[1]+0.5, "Mercury")
        
        # Draw the disk of Mercury
        npoints = 200
        circum_mer_lonlat = np.zeros((2,npoints))
        circ_3d = circumf_3d_vec(npoints, mpo_mer_ecl, 2440)
        for qq in range(npoints) :
            circum_mer_lonlat[0:2,qq] = convert3to2deg(circ_3d[0:3,qq])
            if circum_mer_lonlat[0,qq] > 360 : circum_mer_lonlat[0,qq] = circum_mer_lonlat[0,qq] - 360
            if circum_mer_lonlat[0,qq] < 0 : circum_mer_lonlat[0,qq] = circum_mer_lonlat[0,qq] + 360
            if circum_mer_lonlat[1,qq] > 90 : circum_mer_lonlat[1,qq] = circum_mer_lonlat[1,qq] - 180
            if circum_mer_lonlat[1,qq] < -90 : circum_mer_lonlat[1,qq] = circum_mer_lonlat[1,qq] + 180
        ax.scatter(circum_mer_lonlat[0,0], circum_mer_lonlat[1,0], s=2, color='purple', label='Mercury disk') 
        ax.scatter(circum_mer_lonlat[0,1:], circum_mer_lonlat[1,1:], color='purple', s=2)  
        
    if showven == 1 :
        mpo_ven_ecl_lonlat = convert3to2deg(mpo_ven_ecl)
        ax.scatter(*mpo_ven_ecl_lonlat, color='orange')
        ax.text(mpo_ven_ecl_lonlat[0]+0.5, mpo_ven_ecl_lonlat[1]+0.5, "Venus")
        if angular_sep(mpo_ven_ecl, mpo_mer_ecl) < r_ang_mer :
            ax.text(mpo_ven_ecl_lonlat[0]+0.5, mpo_ven_ecl_lonlat[1]+0.5, "-----")
        
    if showear == 1 :
        mpo_ear_ecl_lonlat = convert3to2deg(mpo_ear_ecl)
        ax.scatter(*mpo_ear_ecl_lonlat, color='cyan', edgecolor='black')
        ax.text(mpo_ear_ecl_lonlat[0]+0.5, mpo_ear_ecl_lonlat[1]+0.5, "Earth")
        if angular_sep(mpo_ear_ecl, mpo_mer_ecl) < r_ang_mer :
            ax.text(mpo_ear_ecl_lonlat[0]+0.5, mpo_ear_ecl_lonlat[1]+0.5, "-----")
        
        
    # Set axes limits
    plt.xlim([0, 360])
    plt.ylim([-90, 90])
    
    # Set axes labels
    ax.set_title('Date (Y-M-D): '+etdate[0:4]+'-'+etdate[5:7]+'-'+etdate[8:10]+', time (h:m:s): '+etdate[11:13]+':'+etdate[14:16]+':'+etdate[17:19])
    ax.set_xlabel('Longitude [°]', labelpad=20, fontsize=16)
    ax.set_ylabel('Latitude [°]', labelpad=20, fontsize=16)
    ax.legend() 
    plt.show() # Make plot!
    
    # If line of sight is at less than 95° from the Sun, inform the user
    suns = sun_in_baffle(ettime_ok, scannerposs, mpo_sun_ecl, mpo_mer_ecl)
    
    if suns == 0 :
        textsun = "The Sun does not constitute a hazard for PHEBUS."
    else :
        textsun = "The Sun constitutes a hazard for PHEBUS! The observation is not possible."

    if ettime_ok < spy.str2et('2025-12-01T00:00:00') and (2100 < scannerposs < 2850) :
        textmtm = "The observation is not feasible as the MTM module obstructs the view of PHEBUS!"
    else :
        textmtm = " "
    
    layout_visual = [
        [sg.Text(textsun)],
        [sg.Text(textmtm)],
        [sg.Button("Close")]
        ]
    version = "1.0"
    window_visual = sg.Window(title="PHESTAR "+version+" - Visualize observation", layout=layout_visual)
    while True :
        event6, values9 = window_visual.read()
        if event6 in ("Close", sg.WIN_CLOSED) : # if click "close", close the window
            break
    window_visual.close()
    
    return(0)
    
        
#### END FUNCTION #############################################################
