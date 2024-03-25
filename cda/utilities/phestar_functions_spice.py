#!/usr/bin/env python
#coding : utf-8

###############################################################################
#
# *** phestar_functions_spice.py ***
#
# DESCRIPTION
# The following code contains two functions called by the routines of
# the PHESTAR software for their correct functioning and a smoother readability.
# These functions make use of the SPICE package, directly or indirectly.
# Each function is provided with detailed description and comments.
# Here is a tweet-like description of each function:
#     - spos2vzecl : it reads a PHEBUS scanner position and a time, and computes
#         the corresponding vectors of the Field Of View reference frame.
#     - pos_body_rf : it computes the position vector of a body with respect to
#         another one in some reference frame and time.
# For any details, visit my Github (https://github.com/marcomeru/phestar_phebus)!
#
# Version 1.0
# Marco Merusi (marco.merusi@latmos.ipsl.fr), LATMOS (Guyancourt, France)
#
###############################################################################

#### IMPORT LIBRARIES #########################################################

# import numpy
if "np" not in dir() :
    import numpy as np

# import Spice
if "spy" not in dir() : 
    import spiceypy as spy

###############################################################################

#### START FUNCTIONS ##########################################################

### pos_body_rf
def pos_body_rf(observer, target, rframe, time) :
    '''
    This function computes the position vector of a celestial body or spacecraft
    (target) seen from another similar body (observer) in some reference frame
    at a given ephemeris time. 
    The function uses the 'numpy' and 'spiceypy' libraries.
    
    Parameters
    ----------
    observer : (string) the observing body, according to the Spice
        naming convention.
    target : (string) the target body, according to the Spice
        naming convention, of which we want the position with respect to the 
        observer.
    rframe : (string) the reference frame, according to the Spice
        naming convention, on which the observer is centered.
    time : (float) the ephemeris time of the observation
        (in seconds after J2000).

    Returns
    -------
    pos_body : (3-element float array) 3D position of the target body. 

    '''
    pos_body = np.zeros(3)
    # Call spkpos to compute 3D position and light-time distance 
    pos_light_body = spy.spkpos(target, time, rframe, "NONE", observer)
    pos_body[0:3] = pos_light_body[0] # Only interested in the 3D position!
    
    return(pos_body)


###############################################################################


### spos2vzecl
def spos2vzecl(scanpos, time):
    '''
    This function computes the vectors of the FoV of PHEBUS (X_FOV, Y_FOV and Z_FOV),
    in ecliptic reference frame, given a scanner position of the instrument and 
    an epoch (expressed as ephemeris time in seconds after J2000).    
    The function uses the 'numpy' and 'spiceypy' libraries.
    
    Parameters
    ----------
    scanpos : (integer) a position of the scanner of PHEBUS, from 0 
        to 4096. If a number greater than 4096 is given in input, the function 
        computes the modulus of that number by 4096.
    time : (float) an ephemeris time (seconds from J2000).

    Returns
    -------
    vz_ecl, vy_ecl, vx_ecl : (tuple of three arrays) Each one is a 3-element float 
        array describing the position of the axes of the reference frame of the
        field of view of PHEBUS' slit (X_FOV, Y_FOV, Z_FOV) in ecliptic reference frame.

    '''    
    D = np.radians(100) # Fixed deviation angle of the scanner mirror of PHEBUS

    S = scanpos * 2*np.pi / 4096 # Convert scanner position to radians
    smp = S - np.pi/4 # Fixed reference with respect to the parking position

    # Compute the Z_FOV component (line of sight) in URF reference frame
    vz_urf = [np.sin(D) * np.cos(smp), np.cos(D), -np.sin(D) * np.sin(smp)]

    # Transform Z_FOV from URF to MPO reference frame
    vz_mpo = [vz_urf[2], vz_urf[1], -vz_urf[0]]

    # Transform Z_FOV from MPO to ecliptic frame
    mtx_mpo_ecl = spy.pxform("MPO_SPACECRAFT", "ECLIPJ2000", time)
    vz_ecl = np.dot(mtx_mpo_ecl, vz_mpo)
    
    # Compute the Y_FOV component in URF reference frame
    vy_urf = np.zeros(3)
    vy_urf[0] =  np.sin(D/2) * np.sin(D/2) * np.cos(2*S)
    vy_urf[1] = -np.sin(D) * np.sin(smp)
    vy_urf[2] = -np.sin(D/2) * np.sin(D/2) * np.sin(2*S) - np.cos(D/2) * np.cos(D/2)
        
    # Compute the Y component of the FoV system in the MPO reference frame
    vy_mpo = [vy_urf[2], vy_urf[1], -vy_urf[0]]

    # Use matrix to convert Y_FOV from MPO to ecliptic reference frame
    vy_ecl = np.dot(mtx_mpo_ecl, vy_mpo)

    # Compute the X component of the FoV system in the ecliptic reference frame
    vx_ecl = np.cross(vy_ecl, vz_ecl)
    
    return(vz_ecl, vy_ecl, vx_ecl)


#### END FUNCTIONS ############################################################
