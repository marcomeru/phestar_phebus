#!/usr/bin/env python
#coding : utf-8

###############################################################################
#
# *** phestar_search_by_star.py ***
#
# DESCRIPTION
# This function is employed in the PHESTAR software and searches when some star
# is visible from PHEBUS and potentially observable because there are no 
# constraints (Sun, obstructions, etc.). 
# The search is made for some stars from a catalog between two dates and times 
# (and time step) chosen by the user, and a star is observable if its separation
# from the line of sight of PHEBUS is smaller than a threshold given in input.
# Knowing the position of a star, the function performs an inverse search to find
# the expected scanner position to see it, and then it tests this scanner position
# to check if the star is in the field of view. 
# Eventually, all the results are listed in a dataframe.
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

# import pandas
if "pd" not in dir() :
    import pandas as pd
 
# import simpleGUI
if "sg" not in dir() :
    import PySimpleGUI as sg

# Manage directories
import sys
import os
# Get the parent directory
parent_dir = os.path.dirname(os.path.realpath(__file__))

# Add the parent directory to sys.path
sys.path.append(parent_dir)

# Useful computing functions
from utilities.phestar_functions import angular_sep
from utilities.phestar_functions_spice import pos_body_rf, spos2vzecl
from utilities.sun_in_baffle import sun_in_baffle


#### START FUNCTION ###########################################################

def phestar_search_by_star(kernels_path, starset, time_start, time_end, dt, max_sep, probar) :
    '''
    This function searches when a star from a catalog can be observed by PHEBUS
    within a time interval and time step given in input. The observations are
    found when the star is perfectly visible to PHEBUS, the Sun is not a hazard
    for PHEBUS, and the angular separation between the star and the line of sight
    of PHEBUS is smaller than a threshold given in input.
    
    Parameters
    ----------
    kernels_path : (string) complete path to the kernel list file (bc_plan.tm).
    starset : (int array) array of the indices of the stars from the catalog 
        for which the observations are searched.
    time_start : (string) start time and date of the interval, in the form
        'YYYY-MM-DDThr:mn:sec'.
    time_end : (string) end time and date of the interval, in the form
        'YYYY-MM-DDThr:mn:sec'.
    dt : (int) time step in seconds on which to run the search.
    max_sep : (float) maximum separation angle between the star and the line of
        sight of PHEBUS, within which the star is considered visible.
    probar : (bool) whether to show the progress bar or not (it can slow the search).
    
    Returns
    -------
    If the search gives no results, it returns the tuple (0, 0). If at least one
    observation is found, it returns a tuple containing the number of observations
    found and the dataframe.
    
    ''' 
    ###########################################################################
    
    #### Open and load SPICE kernels ##########################################
    
    # Load the spice kernels if necessary
    try :
        tim = spy.str2et("2025-01-01T20:00:00")
    except :
        spy.furnsh(kernels_path)
    
    ###########################################################################
    
    #### Process the star catalog to find the stars to observe ################
    
    # Open the catalog of stars and take only those selected by the user
    starcat = pd.read_csv(parent_dir+r"\utilities\starcat_phebus.csv")
    nbst = len(starset) 
    starnm = []
    starhr = np.zeros(nbst)
    vec_ecli = np.zeros((3,nbst))
    for uu in range(nbst) :
        ind = starset[uu]
        starnm.append(starcat['Star name'][ind])
        starhr[uu] = starcat['HR'][ind]
        vec_ecli[0,uu] = starcat['Star X position'][ind]
        vec_ecli[1,uu] = starcat['Star Y position'][ind]
        vec_ecli[2,uu] = starcat['Star Z position'][ind]
            
    ###########################################################################
    
    #### Check the input times ################################################
         
    # Transform the start and end time to seconds from J2000 (ephemeris time)
    # If the times are equal, use it
    if time_start == time_end :
        ettime = [spy.str2et(time_start)]
        dimet = 1
    else :
        
        # Otherwise, compute the times in the time interval in which to check the star observability
        etstart = spy.str2et(time_start)
        etend = spy.str2et(time_end)

        etdurat = etend - etstart + dt
        dimet = round(etdurat/dt)
        ettime = np.zeros(dimet)
        for mm in range(dimet): ettime[mm] = etstart + mm * dt
        
    #### Prepare some arrays and counters for the final storing of results ####
    time_obs = [] # ephemeris time
    time_date_obs = [] # time in date format  
    star_id = [] # Name of star
    star_hr = [] # HR number of star
    ang_los_star = [] # Angle between star and LoS
    scannerpos = [] # Scanner position
    
    D = np.radians(100) # Convert the mirror deviation to radians
    n=0 # Counter of the observable stars in the time interval
    count = 0
    mtm_separation = spy.str2et("2025-12-01T00:00:00")
    niters = dimet * nbst # Expected total number of iterations
    
    ###########################################################################
    
    #### Start iterations on each time value ##################################
    for ii in range(dimet) : # for each time step
        
        # Show progress bar if requested
        if probar == 1 :
            count = count + 1
            if not sg.one_line_progress_meter('PHESTAR (v1.0) - Search by star', count, niters,'Potential observations found: '+str(n)) :
                break
            
        
        for jj in range(nbst) : # for each selected star
            
            ## Start the inverse search
            
            # Obtain the matrix of change of basis ecl -> MPO 
            mtx_ecl_mpo = spy.pxform("ECLIPJ2000", "MPO_SPACECRAFT", ettime[ii])
            # Change the basis and find the vector MPO -> star in MPO RF
            vzmpo_star = np.dot(mtx_ecl_mpo, vec_ecli[0:3,jj])
            # Compute the vector MPO -> star in URF RF
            vzurf_star = [-vzmpo_star[2], vzmpo_star[1], vzmpo_star[0]]
            
            # We obtain D = arccos(vzurf_star[1]), which should be around 100°
            # with a margin = max_sep
            if D - 1.5*np.radians(max_sep) <= np.arccos(vzurf_star[1]) <= D + 1.5*np.radians(max_sep) :
                       
                # Compute the angle smp
                smp_star = np.arccos(vzurf_star[0] / np.sin(D))
                # Compute the scanner angle S
                S_star = smp_star + np.pi/4
                # Obtain the expected scanner position 
                scanpos_star = np.round(S_star * 2048 / np.pi)
                
                # Round scanner position
                scanpos_min = np.floor(scanpos_star)
                scanpos_max = np.ceil(scanpos_star)
                
                ## Test the scanner positions found
                
                ## Floor position
                # Compute the line of sight of PHEBUS in ecliptic RF
                phebfov1 = spos2vzecl(scanpos_min, ettime[ii])
                vz_ecl1 = phebfov1[0]
                
                # Check if the separation between the star and LoS is less than threshold
                if angular_sep(vec_ecli[0:3,jj], vz_ecl1) * 180/np.pi <= max_sep :
                    
                    # Compute the vector MPO -> Mercury in ecliptic RF
                    mpo_mer_ecl = pos_body_rf("-121", "1", "ECLIPJ2000", ettime[ii])
                    r_mpo_mer = lal.norm(mpo_mer_ecl)
                    r_ang_mer = np.arctan2(2440, r_mpo_mer) # Angular radius of Mercury
                    
                    # Check if the star is not behind Mercury
                    if np.degrees(angular_sep(mpo_mer_ecl, vec_ecli[0:3,jj])) > np.degrees(r_ang_mer) :
                        
                        # Before December 2025, the MTM obstructs the view for scanner positions
                        # between 185° (2100) and 250° (2850)
                        if ettime[ii] >= mtm_separation or (ettime[ii] < mtm_separation and (115 <= scanpos_min <= 2100 or 2850 <= scanpos_min <= 3980)) :
                            
                            mpo_sun_ecl = pos_body_rf("-121", "10", "ECLIPJ2000", ettime[ii])
                            # Check if the Sun is not a hazard for PHEBUS  
                            suns = sun_in_baffle(ettime[ii], scanpos_min, mpo_sun_ecl, mpo_mer_ecl)
                    
                            if suns == 0 :
                                                                       
                                timme = spy.timout(ettime[ii], 'YYYY-MM-DDTHR:MN:SC.### ::RND', len('YYYY-MM-DDTHR:MN:SC.### ::RND'))
                                
                                # Store the data of these star in arrays for output
                                time_obs.append(float(ettime[ii])) # ephemeris time
                                time_date_obs.append(timme) # time in date format  
                                star_id.append(str(starnm[jj])) # Name of star
                                star_hr.append(float(starhr[jj])) # HR number of star
                                ang_los_star.append(float(np.degrees(angular_sep(vec_ecli[0:3,jj], vz_ecl1)))) # Angle between star and LoS
                                scannerpos.append(float(scanpos_min)) # Scanner position
                                        
                                n = n + 1 # Update number of observable stars
                    
                    
                    if scanpos_min != scanpos_max :
                    
                        # Ceil position
                        # Compute the line of sight of PHEBUS in ecliptic RF
                        phebfov2 = spos2vzecl(scanpos_max, ettime[ii])
                        vz_ecl2 = phebfov2[0]
                        
                        # Check if the separation between the star and LoS is less than threshold
                        if angular_sep(vec_ecli[0:3,jj], vz_ecl2) * 180/np.pi <= max_sep :
                            
                            # Compute the vector MPO -> Mercury in ecliptic RF
                            mpo_mer_ecl = pos_body_rf("-121", "1", "ECLIPJ2000", ettime[ii])
                            r_mpo_mer = lal.norm(mpo_mer_ecl)
                            r_ang_mer = np.arctan2(2440, r_mpo_mer) # Angular radius of Mercury
                            
                            # Check if the star is not behind Mercury
                            if np.degrees(angular_sep(mpo_mer_ecl, vec_ecli[0:3,jj])) > np.degrees(r_ang_mer) :
                                
                                # Before December 2025, the MTM obstructs the view for scanner positions
                                # between 185° (2100) and 250° (2850)
                                if ettime[ii] >= mtm_separation or (ettime[ii] < mtm_separation and (115 <= scanpos_max <= 2100 or 2850 <= scanpos_max <= 3980)) :
                                    
                                    mpo_sun_ecl = pos_body_rf("-121", "10", "ECLIPJ2000", ettime[ii])
                                    # Check if the Sun is not a hazard for PHEBUS  
                                    suns = sun_in_baffle(ettime[ii], scanpos_max, mpo_sun_ecl, mpo_mer_ecl)
                            
                                    if suns == 0 :
                                                                               
                                        timme = spy.timout(ettime[ii], 'YYYY-MM-DDTHR:MN:SC.### ::RND', len('YYYY-MM-DDTHR:MN:SC.### ::RND'))
                                        
                                        # Store the data of these star in arrays for output
                                        time_obs.append(float(ettime[ii])) # ephemeris time
                                        time_date_obs.append(timme) # time in date format  
                                        star_id.append(str(starnm[jj])) # Name of star
                                        star_hr.append(float(starhr[jj])) # HR number of star
                                        ang_los_star.append(float(np.degrees(angular_sep(vec_ecli[0:3,jj], vz_ecl2)))) # Angle between star and LoS
                                        scannerpos.append(float(scanpos_max)) # Scanner position
                                                
                                        n = n + 1 # Update number of observable stars
    
    #### End of iterations ####################################################

    
    if probar == 1 :
        sg.one_line_progress_meter_cancel() # Close progress bar window
    
    # If no stars found, return 0
    if n == 0 :
        return(0, 0)
    
    # If the program found at least one potential observation, store the data in dataframe
    if n > 0:
        # Format arrays
        # time_obs = np.array(time_obs, dtype=np.float32)
        # star_hr = np.array(star_hr, dtype=np.float32)
        # ang_los_star = np.array(ang_los_star, dtype=np.float32)
        # scannerpos = np.array(scannerpos, dtype=np.float32)
    
        # Store the arrays in a pandas dataframe
        df = { 'time_obs' : time_obs ,  'time_date_obs': time_date_obs , 
              'star_name': star_id, 'star_HR': star_hr, 
              'ang_los_star': ang_los_star, 'scanner_pos': scannerpos}
        obs_dataframe = pd.DataFrame(df)
        
        # Return the number of stars and the dataframe
        return(n, obs_dataframe)
    
#### END FUNCTION #############################################################
