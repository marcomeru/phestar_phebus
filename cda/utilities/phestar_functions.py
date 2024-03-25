#!/usr/bin/env python
#coding : utf-8

###############################################################################
#
# *** phestar_functions.py ***
#
# DESCRIPTION
# The following code contains a series of functions called by the routines of
# the PHESTAR software for their correct functioning and a smoother readability.
# These functions don't make use of the SPICE package.
# Each function is provided with detailed description and comments.
# Here is a tweet-like description of each function:
#     - check_date : it verifies if all the numbers of a date and time are
#         correct and the date exists.
#     - circumf_3d_vec : it computes the coordinates of a circumference in 3D
#         space with the center perpendicular to a given vector.
#     - convert3to2deg : it reads a position vector in 3D space and converts it
#         to a direction on the 2D (longitude,latitude) plan.
#     - angular_sep : it reads two position vectors and computes the
#         angle between them with respect to the origin.
# For any details, visit my Github (https://github.com/marcomeru/phestar_phebus)!
#
# Version 1.0
# Marco Merusi (marco.merusi@latmos.ipsl.fr), LATMOS (Guyancourt, France)
#
###############################################################################

#### IMPORT LIBRARIES #########################################################

# import calendar to check leap years
from calendar import isleap

# import numpy and linalg
if "np" not in dir() :
    import numpy as np
if "lal" not in dir() :
    from numpy import linalg as lal

###############################################################################

#### START FUNCTIONS ##########################################################

### check_date
def check_date(date_str) :
    '''
    This function receives a string of date and time-of-day and checks that all
    its values are correct (and the date and time exist).
    The function uses the 'isleap' module from the 'calendar' library.
    
    Parameters
    ----------
    date_str : (string) a complete date and time (in UTC) in the form
        "YYYY-MM-DDThh:mm:ss".

    Returns
    -------
    It returns 2 if the date does not exist, 3 if the time does not exist, 
    4 if the values in the date or time are in invalid format (e.g., letters 
    instead of numbers), or 0 if the date/time is correct.

    '''
    # Remove any opening or closing quotation marks to the string
    if date_str[0] == "'" or date_str[0] == '"' : date_str = date_str[1:-1]
    if date_str[-1] == "'" or date_str[-1] == '"' : date_str = date_str[0:-2]
    
    # Read the substrings corresponding to the date and time and try to convert
    # them to integers. If not possible, return error
    try :
        year = int(date_str[0:4])
    except :
        return(4)
    try :
        month = int(date_str[5:7])
    except :
        return(4)
    try :
        day = int(date_str[8:10])
    except :
        return(4)
    try :
        hour = int(date_str[11:13])
    except :
        return(4)
    try :
        minute = int(date_str[14:16])
    except :
        return(4)
    try :
        second = int(date_str[17:19])
    except :
        return(4)
    
    # Check that day and month exist
    if month < 1 or month > 12 or day < 1 :
        return(2)
    
    # Check that the day exists according to the month and the leap years
    if month == 2 :
        if (isleap(year) == True and day > 29) or (isleap(year) == False and day > 28) :
            return(2)
    if (month == 4 or month == 6 or month == 9 or month == 11) and day > 30 :
        return(2)
    if (month == 1 or month == 3 or month == 5 or month == 7 or month == 8 or month == 10 or month == 12) and day > 31 :
        return(2)
    
    # Check that the time exists
    if hour < 0 or hour > 23 or minute < 0 or minute > 59 or second < 0 or second > 59 :
        return(3)
    
    # Otherwise, if everything is valid and correct, return 0
    return(0)


###############################################################################


### circumf_3d_vec
def circumf_3d_vec(npoints, body_pos, r_body) :
    '''
    This function computes the coordinates of the points of a circumference in
    3D space lying on a plane that is perpendicular to the vector going from 
    the origin to the center of the circle.
    It is used to compute the 3D coordinates of the disk of a planet (or the Sun)
    seen from another body. For example the disk of Mercury seen from MPO.
    The function uses 'numpy' library. 
    
    Parameters
    ----------
    npoints : (float) value stating the number of points to represent on the 
        circumference.
    body_pos : (3-element float array) the position of the center of the
        circumference.
    r_body : (float) the radius of the circumference in km.

    Returns
    -------
    circ : (3Xnpoints float array) points of the circumference in 3D space.

    '''
    # Define the points to generate the circumference from 0 to 2pi
    theta = np.linspace(0, 2*np.pi, npoints)
    circ = np.zeros((3,len(theta)))
    
    # Compute the distance to the body's center
    r_bodypos = lal.norm(body_pos)
    
    # Compute the unit vector to the body's center
    body_pos_dir = [body_pos[0]/r_bodypos, body_pos[1]/r_bodypos, body_pos[2]/r_bodypos]
    
    # Compute the absolute values of the elements of the unit vector
    bodypos_xxyyzz = [abs(body_pos_dir[0]), abs(body_pos_dir[1]), abs(body_pos_dir[2]) ]
    indmax = np.argmax(bodypos_xxyyzz) # Find the index of the largest element
    indmin = np.argmin(bodypos_xxyyzz) # Find the index of the minimum element
    bodypos_xxyyzz[indmax] = -1 # "Remove" the largest element
    secmax = np.argmax(bodypos_xxyyzz) # Find the index of the second largest element
    
    # Form the vector perpendicular to the unit direction to the body
    p = np.zeros(3)
    p[indmax] = body_pos_dir[secmax]
    p[secmax] = -body_pos_dir[indmax]
    p[indmin] = 0
    # Compute the norm of this vector
    p_norm = lal.norm(p)
    # Compute its unit vector
    p_dir = [p[0]/p_norm, p[1]/p_norm, p[2]/p_norm]
    # Compute the unit vector perpendicular to body_pos and p
    b = np.cross(body_pos_dir, p_dir)
    
    # Compute the points of the circumference for each value of the angle theta
    for qq in range(len(theta)) :
        circ[0,qq] = body_pos[0] + r_body * p_dir[0] * np.cos(theta[qq]) + r_body * b[0] * np.sin(theta[qq])
        circ[1,qq] = body_pos[1] + r_body * p_dir[1] * np.cos(theta[qq]) + r_body * b[1] * np.sin(theta[qq])
        circ[2,qq] = body_pos[2] + r_body * p_dir[2] * np.cos(theta[qq]) + r_body * b[2] * np.sin(theta[qq])

    return(circ)


###############################################################################


### convert3to2deg ####
def convert3to2deg(coordinate):
    """
    This function transforms a set of cartesian coordinates XYZ
    to a couple of geographic coordinates (latitude, longitude).
    The function uses the 'numpy' library.

    Parameters
    ----------
    coordinate : (3-element float array) position vector of the body.

    Returns
    -------
    lonlat : (2-element float array) (longitude, latitude) in degrees.

    """
    lonlat = np.zeros(2)
    # Compute the norm of the array defined by the 3 coordinates
    norm = lal.norm(coordinate)
    # Normalize the components of the vector
    xecli = coordinate[0] / norm 
    yecli = coordinate[1] / norm
    zecli = coordinate[2] / norm
    # Compute the longitude 
    lon = np.degrees(np.arctan2(float(yecli), float(xecli))) 
    if lon < 0:
        lon += 360
    # Compute the latitude
    lat = np.degrees(np.arcsin(zecli))
    lonlat[0] = lon
    lonlat[1] = lat
    return(lonlat)


###############################################################################


### angular_sep ####
def angular_sep(fir_pos, sec_pos):
    '''
    This function computes the smaller angular separation between two position
    vectors with respect to their common origin. The two vectors are in the same 
    reference frame and unit.
    The function uses the 'numpy' library and its module 'linalg'.
    
    Parameters
    ----------
    fir_pos : (3-element float array) position vector of the first point.
    tar_pos : (3-element float array) position vector of the second point.

    Returns
    -------
    angle : (float) smaller anglular separation between the two points, in radians.

    '''
    # Normalize the two vectors if their norms are not =1
    if lal.norm(fir_pos) != 1:
        r_fir_pos = lal.norm(fir_pos)
        fir_pos = [ fir_pos[0] / r_fir_pos , fir_pos[1] / r_fir_pos , fir_pos[2] / r_fir_pos ]
    if lal.norm(sec_pos) != 1:
        r_sec_pos = lal.norm(sec_pos)
        sec_pos = [ sec_pos[0] / r_sec_pos , sec_pos[1] / r_sec_pos , sec_pos[2] / r_sec_pos ]
    
    # Compute the angle by scalar product 
    cosbeta = np.dot(fir_pos, sec_pos)
    angle = np.arccos(cosbeta)
    
    return(angle)
        

#### END FUNCTIONS ############################################################
