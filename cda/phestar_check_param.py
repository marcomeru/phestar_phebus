#!/usr/bin/env python
#coding : utf-8

###############################################################################
#
# *** check_param_phebus_stars.py ***
#
# DESCRIPTION
# The following code defines the function 'check_param_phebus_stars', which is 
# used in the PHESTAR software. The aim of this function is to check if all the
# parameters given by the user in input in the search window are valid and correct.
# If they are not, the function returns to the main program showing a popup 
# that informs the user on the error. If everything is correct, it gives green 
# light to the main program to continue. 
# According to the first input parameter, it can check the values for the star
# search or the visualizer. 
# For any details, visit my Github (https://github.com/marcomeru/phestar_phebus)!
#
# Version 1.0
# Marco Merusi (marco.merusi@latmos.ipsl.fr), LATMOS (Guyancourt, France)
#
###############################################################################

#### IMPORT LIBRARIES #########################################################

# Import SPICE
if "spy" not in dir() :
    import spiceypy as spy

# Import pysimpleGUI
if "sg" not in dir() :
    import PySimpleGUI as sg

# Manage directories
import sys
import os
# Get the parent directory
parent_dir = os.path.dirname(os.path.realpath(__file__))

# Add the parent directory to sys.path
sys.path.append(parent_dir)

# Import some useful functions
from utilities.phestar_functions import check_date
from utilities.phestar_functions_spice import pos_body_rf

#### START FUNCTION ###########################################################
def phestar_check_param(whatto, kernel, timest, timend, timestep, scanp, maxang) :
    '''
    This function checks if all the parameters given in input are valid and correct
    (strings are valid, numbers are numbers, etc.).    
    According to the first input parameter, it can check the values for the star
    search or the observation visualizer. In the second case some parameters
    will be ignored.
    In both cases, it uses a nested IF-ELSE check from the first to the last 
    parameter, in order to save time.

    Parameters
    ----------
    whatto : (bool) If it's 0 it's the case of star search, if 1 it's the
        observation visualizer.
    kernel : (string) Complete path to the kernel list file ('bc_plan.tm').
    timest : (string) For the star search, this is the starting date and time. 
        For the visualizer, this is the chosen date and time. It's in the form
        'YYYY-MM-DDThr:mm:ss'.
    timend : (string) For the star search, this is the end date and time. 
        For the visualizer, this is ignored.
    timestep : (float) Time step for the star search, in seconds. For the 
        visualizer, this is ignored.
    scanp : (integer) For the star search, this is ignored. For the visualizer, 
        this is the chosen scanner position. It can range from 100 to 3990.
    maxang : (float) For the star search, this is the maximum angle within which 
        stars are considered observable. For the visualizer, it's the maximum angle
        within which star names are shown in the map. It's expressed in degrees.

    Returns
    -------
    It returns 1 if there is any invalid value that must be corrected. If 
    everything is fine, it returns 0.

    '''
    ###########################################################################
    
    # Case of the search of stars
    if whatto == 0 :
        
        # Check if the kernel list file is valid
        try :
            _ = spy.str2et("2025-01-01T20:00:00")
        except :
            try :
                spy.furnsh(kernel)
            except :
                sg.PopupError("The kernel list file is invalid. Please insert a valid file (e.g., \"bc_plan.tm\")!")
                return(1)
        else :          
            
            # Check the start date
            ret = check_date(timest)
            if ret == 2 : 
                sg.PopupError("The start date does not exist!")
                return(1)
            elif ret == 3 : 
                sg.PopupError("The start time does not exist!")
                return(1)
            elif ret == 4 : 
                sg.PopupError("The start time is in invalid format. Please insert a date/time in the form 'YYYY-MM-DDThh:mm:ss'.")
                return(1)
            elif ret == 0 :                         
                
                # Check if the kernel file has data for MPO at the start date
                try :
                    #spy.furnsh(kernel)
                    tim2 = spy.str2et(timest)
                    _ = pos_body_rf("-121", "1", "ECLIPJ2000", tim2)
                except :
                    sg.PopupError("There are no data on BepiColombo on the inserted start date. It may be too far in the past or in the future. Please insert a closer date!")
                    return(1)
                else :
                
                    # Check the end date
                    ret2 = check_date(timend)
                    if ret2 == 2 : 
                        sg.PopupError("The end date does not exist!")
                        return(1)
                    elif ret2 == 3 : 
                        sg.PopupError("The end time does not exist!")
                        return(1)
                    elif ret2 == 4 : 
                        sg.PopupError("The end time is in invalid format. Please insert a date/time in the form 'YYYY-MM-DDThh:mm:ss'.")
                        return(1)
                    elif ret2 == 0 : 
                        
                        # Check if the kernel file has data for MPO at the end date
                        try :
                            #spy.furnsh(kernel)
                            tim3 = spy.str2et(timend)
                            _ = pos_body_rf("-121", "1", "ECLIPJ2000", tim3)
                        except :
                            sg.PopupError("There are no data on BepiColombo on the inserted end date. It may be too far in the past or in the future. Please insert a closer date!")
                            return(1)
                        else :
                                                            
                            # Check if the end date is after the start date
                            # Check if the end date is after the start date
                            t1 = spy.str2et(timest)
                            t2 = spy.str2et(timend)
                            if t1 > t2 : 
                                sg.PopupError("The start time can't be after the end time!")
                                return(1)
                            else :                            
                               
                                # Check the time step
                                try :
                                    tmp = int(timestep)
                                except :
                                    sg.PopupError("The time step must be a positive integer!")
                                    return(1)
                                else :
                                    
                                    if tmp <= 0 :
                                        sg.PopupError("The time step must be a positive integer!")
                                        return(1)
                                    else :                                    
                                        
                                        # Check the maximum angle
                                        try :
                                            tmp = float(maxang)
                                        except :
                                            sg.PopupError("The maximum angular distance must be a positive number!")
                                            return(1)
                                        else :
                                            if tmp <= 0 :
                                                sg.PopupError("The maximum angular distance must be a positive number!")
                                                return(1)
                                            else :
                                                
                                                # If everything is good, return 0
                                                return(0)
                                                                    
    ###########################################################################
    
    # Case of the observation visualizer
    if whatto == 1 :
        
        # Check if the kernel list file is valid
        try :
            _ = spy.str2et("2025-01-01T20:00:00")
        except :
            try :
                spy.furnsh(kernel)
            except :
                sg.PopupError("The kernel list file is invalid. Please insert a valid file (e.g., \"bc_plan.tm\")!")
                return(1)
        else :              
                                
            rets = check_date(timest)
            if rets == 2 : 
                sg.PopupError("The date does not exist!")
                return(1)
            elif rets == 3 : 
                sg.PopupError("The time does not exist!")
                return(1)
            elif rets == 4 : 
                sg.PopupError("The start time is in invalid format. Please insert a date/time in the form 'YYYY-MM-DDThh:mm:ss'.")
                return(1)
            elif rets == 0 : 
                
                # Check if the kernel file has data for MPO at the input date
                try :
                    #spy.furnsh(kernel)
                    tim5 = spy.str2et(timest)
                    _ = pos_body_rf("-121", "1", "ECLIPJ2000", tim5)
                except :
                    sg.PopupError("There are no data on BepiColombo on the inserted date. It may be too far in the past or in the future. Please insert a closer date!")
                    return(1)
                else :
                
                    # Check the scanner position
                    try :
                        tmp = int(scanp)
                    except :
                        sg.PopupError("The scanner position must be a positive integer between 100 and 3990!")
                        return(1)
                    else :
                        if tmp < 100 or tmp > 3990 :
                            sg.PopupError("The scanner position must be a positive integer between 100 and 3990!")
                            return(1)
                        else :
                            
                            # Check the maximum angle
                            try :
                                tmp = float(maxang)
                            except :
                                sg.PopupError("The maximum angular distance must be a positive number!")
                                return(1)
                            else :
                                if tmp <= 0 :
                                    sg.PopupError("The maximum angular distance must be a positive number!")
                                    return(1)
                                else :
                                    
                                    # If everything is good, return 0
                                    return(0)
    
#### END FUNCTION #############################################################
