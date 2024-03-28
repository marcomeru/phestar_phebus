#!/usr/bin/env python
#coding : utf-8

###############################################################################
#
# *** phestar_main.py ***
#
# DESCRIPTION
# This code is the main program of the small software PHEBUS Stellar Targeting,
# Alignment and Reconnaissance (PHESTAR). It allows the user to search for stars 
# that can potentially be observable over some time interval by PHEBUS, 
# and to visualize those observations in a simple celestial map. 
# Originally, the program was made of a pair of python routines which were 
# suitably modified and connected through this main routine, which provides a 
# nice and simple interface.
# All errors should have been taken care of. If you notice a mistake, or if the
# program crashes due to some unexpected error, please let me know!
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

# import pandas
if "pd" not in dir() :
    import pandas as pd

# import pysimplegui
if "sg" not in dir() :
    import PySimpleGUI as sg

# import the function to search stars
from cda.phestar_search_by_star import phestar_search_by_star

# import the function to visualize
from cda.phestar_visual_star import phestar_visual_star

# import some functions
from cda.phestar_check_param import phestar_check_param

# import timeit
import time

###############################################################################
######## WINDOW: MAIN MENU

version = "(v1.0)"

# Buttons and layout of main menu
search_stars = "Search stars"
visual_obs = "Visualize observation"
layout_menu = [ [sg.Image(r'cda\utilities\logo_phebus.png')],
    [sg.Column([[sg.Button(search_stars), sg.VSeparator(), sg.Button(visual_obs)]], justification="center")],
    [sg.HSeparator()],
    [sg.Button("Close")]
    ]

# Open window
window_menu = sg.Window(title="PHESTAR "+version, layout=layout_menu)
which_com = 0
while True :
    event0, values0 = window_menu.read() # read action on mouse click
    if event0 == search_stars : # search stars
        which_com = 1
        break
    if event0 == visual_obs : # make plot
        which_com = 2
        break
    if event0 in ("Close", sg.WIN_CLOSED) : # Click "close" or X
        break
    
window_menu.close() # Finally, close window

###############################################################################
######## WINDOW: SEARCH STARS

search = 0
if which_com == 1 :
    
    # Open the catalog of stars
    starcat = pd.read_csv(r"cda\utilities\starcat_phebus.csv")
    star_heads = ['', 'Star name', 'HR', 'Star vector X','Star vector Y', 'Star vector Z', 'Sum EUV+FUV']
    star_vals = starcat.values.tolist()
    
    # In the left column of the window there is the star list (selectable)
    layout_lcol = [
        [sg.Table(values = star_vals, headings = star_heads, vertical_scroll_only=False, visible_column_map=[False, True, True, False, False, False, True], key='-STARTABLE-', select_mode=sg.TABLE_SELECT_MODE_EXTENDED)],
        [sg.Button("Select all", key='-SELALL-'), sg.Button("Deselect all", key='-DESELALL-')]
        ]
    
    # In the right column of the window there are the other fields to fill in
    layout_rcol = [
        [sg.Text("Kernel list file ('bc_plan.tm) : "), sg.Input(key='-KERNEL-'), sg.FileBrowse()],
        [sg.Text("Time start in the form 'YYYY-MM-DDThh:mm:ss' : "), sg.Input(key='-TIMEST-')],
        [sg.Text("Time end in the form 'YYYY-MM-DDThh:mm:ss' : "), sg.Input(key='-TIMEND-')],
        [sg.Text("Time step in seconds : "), sg.Input(size=(8,1), key='-TIMESTEP-')],
        [sg.Text("Maximum angular distance from Line of Sight in degrees : "), sg.Input(size=(5,1), key='-MAXANG-')],
        [sg.Checkbox('Show progress bar', default=True, key='-PROGBAR-')]
        ]
    
    # Main window of the star search form
    layout_search_stars = [
        [sg.Text("Welcome to the PHEBUS Stellar Targeting, Alignment and Reconnaissance (PHESTAR) search engine!", font='bold')],
        [sg.Text("This app will search for all the potential observations of the stars from a catalog according to some input parameters. \n"+
                 "They include the file with the list of SPICE kernels, the star catalog, the dates/times of the preferred interval and \n"+
                 "the time step, and the maximum angular distance of the star from the line of sight of PHEBUS.")],
        [sg.Text("Please be sure to select at least one star from the list on the left, and to fill in all the fields on the right correctly before launching the search! \n"+
                 "N.B.! The time step must be a positive integer. However, if the start and end dates/times are equal, the program will (logically) ignore the time step.")],
        [sg.Text("Until the separation of MTM from MPO (December 2025) the field of view of PHEBUS is obstructed for scanner positions from \n"+
                 "2100 to 2850. In addition, PHEBUS can observe when it's out of the parking position, namely for scanner positions from 100 to 3990.\n"+
                 "The baffle of PHEBUS can only be illuminated sideways by the Sun. Direct sunlight must not enter the aperture of the baffle, \n"+
                 "which means that the observation can be done if the Sun is occulted or it's at least 90° away from the line of sight of PHEBUS.")],
        [sg.Column(layout_lcol), sg.VSeparator(), sg.Column(layout_rcol)],
        [sg.HSeparator()],
        [sg.Button("Launch"), sg.Button("Close")]
        ]
    
    
    # Open window
    window_searchstars = sg.Window(title="PHESTAR "+version+" - Search stars", layout=layout_search_stars)
    while True :
        event1, values1 = window_searchstars.read() # Read actions and values
        
        if event1 in ("Close", sg.WIN_CLOSED) : # Click "close" or X
            break
        
        # If the user clicks 'select all', all the stars are selected
        if event1 == '-SELALL-' :
            window_searchstars['-STARTABLE-'].update(select_rows=list(range(len(starcat['Star name']))))    
        # If the user clicks 'deselect all', all the stars are deselected
        if event1 == '-DESELALL-' :
            window_searchstars['-STARTABLE-'].update(select_rows=[])
                
        if event1 == "Launch" : # Click "launch"
            
            # Check that at least one star is selected
            if values1['-STARTABLE-'] == []:
                sg.PopupError("At least one star must be selected for the search!")
            else :
        
                # Check that all required fields are filled in
                if values1['-KERNEL-'] == '' or values1['-TIMEST-'] == '' or values1['-TIMEND-'] == '' or values1['-TIMESTEP-'] == '' or values1['-MAXANG-'] == '' :
                    sg.PopupError("Fill in all the required fields before launching the search!")
                else :
                                 
                    # Call the function to check all the parameters  
                    res = phestar_check_param(0, values1['-KERNEL-'], values1['-TIMEST-'], values1['-TIMEND-'], values1['-TIMESTEP-'], 0, values1['-MAXANG-'])        
                    
                    if res == 0 :
                        # If all values are good and no errors are raised, the search can start
                        kernel_path = values1['-KERNEL-']
                        pheb_cata = values1['-STARTABLE-']
                        tstart = values1['-TIMEST-']
                        tend = values1['-TIMEND-']
                        if tstart == tend :
                            tstep = 0
                        else :
                            tstep = int(values1['-TIMESTEP-'])
                        angdist = float(values1['-MAXANG-'])
                        progbar = 0 if values1['-PROGBAR-'] == False else 1
                        
                        search = 1 # Give it a GO for the search!
                        
                        break # Close the loop
                                                                        
    window_searchstars.close() # close the search window

# If the search was launched successfully, let's search stars
if search == 1 :
    
    start_t = time.time()
    # Call the function that searches stars based on the input parameters
    ressear = phestar_search_by_star(kernel_path, pheb_cata, tstart, tend, tstep, angdist, progbar)
    # If the search gave 0 results, inform the user
    end_t = time.time()
    
    if ressear[0] == 0 :
        layout_reszero = [ 
            [sg.Text("The search gave zero results!")],
            [sg.Text("The search was performed in "+str(np.round(end_t-start_t, 2))+" seconds.")],
            [sg.Button("Close")]
            ]
        
        window_ressearzero = sg.Window(title="PHESTAR "+version+" - Search results", layout=layout_reszero)
        while True :
            event2, values2 = window_ressearzero.read()
            
            if event2 in ("Close", sg.WIN_CLOSED) :
                break
        
        window_ressearzero.close()
    
    # If the search gave results
    if ressear[0] != 0 :
        
        datafr = ressear[1]
        heads = ['time_obs',  'time_date_obs' ,                  
                 'star_name', 'star_HR', 
                 'ang_los_star', 'scanner_pos']
        vals = datafr.values.tolist()
        
        layout_resmore = [ 
            [sg.Text("The search gave "+str(ressear[0])+" result(s)!")],
            [sg.Text("The search was performed in "+str(np.round(end_t-start_t, 2))+" seconds.")],
            [sg.Table(values = vals, headings=heads, auto_size_columns=True, vertical_scroll_only=False, select_mode=sg.TABLE_SELECT_MODE_BROWSE, key='-TABLE-', right_click_selects=True, right_click_menu=['',['Visualize observation']]) ], 
            [sg.Input(visible=False, enable_events=True, key='-IN-'), sg.FileSaveAs("Export table", file_types=(("CSV Files","*.csv"),)), sg.Button("Close") ] ]
        
        # Show the window with the results from the dataframe
        window_ressearmany = sg.Window(title="PHESTAR "+version+" - Search results", layout=layout_resmore, resizable=True)
    
        while True :
            event3, values3 = window_ressearmany.read()
            
            if event3 in ("Close", sg.WIN_CLOSED) :
                break
            
            # The user can right-click the results to plot them
            if event3 == "Visualize observation" :
                roww = values3['-TABLE-']
                if len(roww) == 0 :
                    sg.Popup("No rows selected!")
                    break
                else :
                    nrow = roww[0]
                    myrow = vals[nrow]
                    tt = myrow[1]
                    posscan = int(myrow[5])
                    
                    layout_show_planets = [
                        [sg.Text("Slit size : "), sg.Radio('Across (0.1°x2°)', "RADIO1", key='-ACROSS-', default=True), sg.Radio('Removed (1.4°x3.1°)', "RADIO1", key='-REMOVED-')],
                        [sg.Text("Show Sun/planets : "), sg.Checkbox('Sun', key='-SUN-'), sg.Checkbox('Mercury', key='-MER-'), sg.Checkbox('Venus', key='-VEN-'), sg.Checkbox('Earth', key='-EAR-')],
                        [sg.HSeparator()],
                        [sg.Button("Launch"), sg.Button("Close")]
                        ]
                    window_show_planets = sg.Window(title="PHESTAR "+version+" - Search results", layout=layout_show_planets)
                    while True :
                        event4, values4 = window_show_planets.read()
                        
                        if event4 in ("Close", sg.WIN_CLOSED) :
                            break
                        
                        # If the user wants to see the plot, ask for more information
                        # on the slit size and if they want to plot Sun/planets
                        if event4 == 'Launch' :
                            if values4['-ACROSS-'] == True :
                                slitx = 0.1 # deg
                                slity = 2.0 # deg
                            else :
                                slitx = 1.4 # deg
                                slity = 3.1 # deg
                            sunz = 0 if values4['-SUN-'] == False else 1
                            merz = 0 if values4['-MER-'] == False else 1
                            venz = 0 if values4['-VEN-'] == False else 1
                            earz = 0 if values4['-EAR-'] == False else 1
                            
                            rr = phestar_visual_star(kernel_path, pheb_cata, tt, posscan, slitx, slity, angdist, sunz, merz, venz, earz)
                            
                    window_show_planets.close()
                    
            # The user can export the table as CSV file
            if event3 == "-IN-" :
                filename = values3['-IN-']
                datafr.to_csv(filename)
                continue
    
        window_ressearmany.close()



###############################################################################



######### Visualize observation

# If the user chooses to visualize an observation
if which_com == 2 :
    
    # Open the catalog of stars
    starcat = pd.read_csv(r"cda\utilities\starcat_phebus.csv")
    star_heads = ['', 'Star name', 'HR', 'Star vector X','Star vector Y', 'Star vector Z', 'Sum EUV+FUV']
    star_vals = starcat.values.tolist()
    
    # In the left column of the window there is the star list (selectable)
    layout_visual_lcol = [
        [sg.Table(values = star_vals, headings = star_heads, vertical_scroll_only=False, visible_column_map=[False, True, True, False, False, False, True], key='-STARTABLE2-', select_mode=sg.TABLE_SELECT_MODE_EXTENDED)],
        [sg.Button("Select all", key='-SELALL2-'), sg.Button("Deselect all", key='-DESELALL2-')]
        ]
    
    layout_visual_rcol = [
        [sg.Text("Kernel list file ('bc_plan.tm') : "), sg.Input(key='-KERNEL2-'), sg.FileBrowse()],
        [sg.Text("Time in the form 'YYYY-MM-DDThh:mm:ss' : "), sg.Input(key='-TIME2-')],
        [sg.Text("Scanner position : "), sg.Input(size=(5,1), key='-SCANPOS2-')],
        [sg.Text("Slit size : "), sg.Radio('Across (0.1°x2°)', "RADIO1", default=True, key='-ACROSS2-'), sg.Radio('Removed (1.4°x3.1°)', "RADIO1", key='-REMOVED2-')],
        [sg.Text("Show Sun/planets : "), sg.Checkbox('Sun', key='-SUN2-'), sg.Checkbox('Mercury', key='-MER2-'), sg.Checkbox('Venus', key='-VEN2-'), sg.Checkbox('Earth', key='-EAR2-')],
        [sg.Text("Show names of stars at less than : "), sg.Input(size=(5,1), key='-MAXANG2-'), sg.Text(" degrees from the line of sight of PHEBUS")]
        ]
    
    layout_visual_star = [ 
        [sg.Text("Welcome to the PHEBUS Stellar Targeting, Alignment and Reconnaissance (PHESTAR) visualizer!", font='bold')],
        [sg.Text("This app will allow you to visualize a potential observation of a star from a catalog according to some input parameters. \n"+
                 "The latter include the stars from a catalog, the SPICE kernel list file, the ephemeris time, the slit size, and the maximum \n"+
                 "angular radius from the line of sight of PHEBUS within which to show the name of stars. In addition, you can show the Sun \n "+
                 "and some planets in the plot.")],
        [sg.Text("Please be sure to choose at least one star from the list on the left and to fill in all the fields on the right correctly before launching the search! \n"+
                 "N.B.! The plot is a bidimensional map (longitude, latitude) obtained from the tridimensional ecliptic reference frame \n"+
                 "centered on MPO. For the Sun and planets (except for Mercury), the sizes of their dots in the plot do not account for their apparent \n"+
                 "angular size.")],
        [sg.Text("Until the separation of MTM from MPO (December 2025) the field of view of PHEBUS is obstructed for scanner positions from \n"+
                 "2100 to 2850. In addition, PHEBUS can observe when it's out of the parking position, namely for scanner positions from 100 to 3990.\n"+
                 "The baffle of PHEBUS can only be illuminated sideways by the Sun. Direct sunlight must not enter the aperture of the baffle, \n"+
                 "which means that the observation can be done if the Sun is occulted or it's at least 90° away from the line of sight of PHEBUS.")],
        [sg.Column(layout_visual_lcol), sg.VSeparator(), sg.Column(layout_visual_rcol)],
        [sg.HSeparator()],
        [sg.Button("Launch"), sg.Button("Close")],
        ]
    
    
    # Show window where the user can choose the parameters
    window_visstar = sg.Window(title="PHESTAR "+version+" - Visualize observation", layout=layout_visual_star)
    
    while True :
        event5, values5 = window_visstar.read()
        
        if event5 in ("Close", sg.WIN_CLOSED) : # if click "close", close the window
            break
        
        # If the user clicks 'select all', all the stars are selected
        if event5 == '-SELALL2-' :
            window_visstar['-STARTABLE2-'].update(select_rows=list(range(len(starcat['Star name']))))    
        # If the user clicks 'deselect all', all the stars are deselected
        if event5 == '-DESELALL2-' :
            window_visstar['-STARTABLE2-'].update(select_rows=[])
            
        if event5 == "Launch" : # Launch the plot
            
            # Check that at least one star is selected
            if values5['-STARTABLE2-'] == []:
                sg.PopupError("At least one star must be selected for the search!")
            else :
                
                if values5['-KERNEL2-'] == '' or values5['-TIME2-'] == '' or values5['-SCANPOS2-'] == '' or values5['-MAXANG2-'] == '' :
                    sg.PopupError("Fill in all the fields before launching the search!")  
                else :
                    
                    # Call the function to check all the parameters
                    res1 = phestar_check_param(1, values5['-KERNEL2-'], values5['-TIME2-'], 0, 0, values5['-SCANPOS2-'], values5['-MAXANG2-'])        

                    if res1 == 0 :   
                        # If everything seems all right, we are good to go
                        kernel_path = values5['-KERNEL2-']
                        starcat_sel = values5['-STARTABLE2-']
                        tt = values5['-TIME2-']
                        posscan = int(values5['-SCANPOS2-'])
                        if values5['-ACROSS2-'] == True :
                            slitx = 0.1 # deg
                            slity = 2.0 # deg
                        else :
                            slitx = 1.4 # deg
                            slity = 3.1 # deg
                        sunz = 0 if values5['-SUN2-'] == False else 1
                        merz = 0 if values5['-MER2-'] == False else 1
                        venz = 0 if values5['-VEN2-'] == False else 1
                        earz = 0 if values5['-EAR2-'] == False else 1
                        max_ang_dist = float(values5['-MAXANG2-'])
                        
                        # Call the function that makes the plot
                        rr = phestar_visual_star(kernel_path, starcat_sel, tt, posscan, slitx, slity, max_ang_dist, sunz, merz, venz, earz)
                        
                        
                    
    window_visstar.close()           
            
#### END PROGRAM ##############################################################
