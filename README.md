# PHESTAR (version 1.0)
## PHEBUS Stellar Targeting, Alignment and Reconnaissance

Marco Merusi (marco.merusi@latmos.ipsl.fr)<br>
LATMOS (Guyancourt, France)

### Introduction
In the cruise phase and the Mercury orbit phase of ESA BepiColombo spacecraft, PHEBUS instrument is calibrated by pointing stars with a known spectrum. The choice of the stars to observe depends not only on scientific factors, such as their physical properties, but also orbital constraints, that is when some star crosses the slit of the instrument at some scanner position. While the orbit and attitude of Bepi was already determined in advance (with some margin for slight changes) for the next few years, the rotation of PHEBUS can be planned with relatively shorter notice based on the scanner position that allows the observation of a star. In this scope, the PHESTAR software assists the user with the planning of potential observations of stars from a catalog within an arbitrary time interval. In addition, the software can realize plots of the celestial sphere relative to BepiColombo showing the stars from the catalog and the field of view of PHEBUS.

### Installation
The software is based on some files to download, a couple of python packages to install, and a large dataset to have available locally, on cloud or on the network. 

The files in this repository have to be downloaded into a folder on the local machine, and they include:
- `phestar_main.py`: it's the main program managing the graphic interface and showing the results;
- in the `/cda/` directory:
  - `phestar_check_param.py`: a function checking all the parameters given in input by the user in the interface;
  - `phestar_search_by_star.py`: a function searching for all potential observations of some stars from a catalog within a time interval;
  - `phestar_visual_star.py`: a function making the graphic visualization of an observation in a 2D plot.
- in the `/cda/utilities/` directory:
  - `phestar_functions.py` and `phestar_functions_spice.py`: two modules containing several useful functions (the latter makes use of the Spiceypy package);
  - `sun_in_baffle.py`: a function verifying whether the Sun is a hazard for PHEBUS;
  - `starcat_phebus.csv`: the star catalog with positions, names and brightness;
  - `logo_phebus.png`: the PHEBUS logo.

These codes are based on two python libraries that need to be installed on the machine, [PySimpleGUI](https://www.pysimplegui.org/en/latest/) (for the graphic interface) and [Spiceypy](https://spiceypy.readthedocs.io/en/v2.3.1/index.html) (for the computation of orbits and ephemeris). They can be installed through the commands:
```
pip install pysimplegui
```
and:
```
pip install spiceypy
```
Finally, the user must have direct access to the [SPICE](https://naif.jpl.nasa.gov/naif/index.html) kernel dataset for BepiColombo, and in particular the file `bc_plan.tm`. The SPICE kernels can be downloaded from [this link](https://s2e2.cosmos.esa.int/bitbucket/projects/SPICE_KERNELS/repos/bepicolombo/browse). The whole kernel bundle is more than 3 gigabytes and is a huge dataset of binary files that are accessed through the Spiceypy package.

The PHESTAR software can be opened by running the file `phestar_main.py`. 

### Functionalities
In the planning of the observation of stars with PHEBUS, PHESTAR carries out two tasks: **Search stars** and **Visualize observation**.

#### 1) Search stars
Through a form, the user can insert some parameters, which include:
- the interesting stars from a catalog. The user is free to select all the stars or only some of them; 
- the files with the list of kernels that must be loaded for the computation of the orbits. The file is typically called `bc_plan.tm` and is provided with the freeware [SPICE package](https://s2e2.cosmos.esa.int/bitbucket/projects/SPICE_KERNELS/repos/bepicolombo/browse) in the directory `/bepicolombo/kernels/mk/`. The package for BepiColombo takes up a few gigabytes of memory, so the user is expected to have it stored locally or on cloud;   
- the date and time corresponding to the start and end of the time interval in which to search for observations, and a time step (in seconds) on which to search. The date/time must be written in the form `YYYY-MM-DDThh:mm:ss`, though the characters `-`, `:` and the `T` can be replaced by any other character. The dates can be any, as long as the SPICE package finds data on the orbit of BepiColombo, and they are transformed to ephemeris times (seconds from J2000);
- the maximum angle from the line of sight of PHEBUS within which to consider the observations as "potential". The best configuration is when the star crosses the shape of the slit perpendicularly at half height, and therefore a logical value to give for the maximum angle would be half width of the slit;
- whether to show an informative progress bar during the search. Although it is useful to keep track of the elapsed time, several tests showed that it slows down the program by a factor 9. Therefore, in case of a large expected number of iterations, it would be better to untick this box.

The program iterates over each ephemeris time in the interval, defined by the time step, and for each selected star from the catalog it makes an inverse computation of the corresponding scanner position that would allow PHEBUS to see the star. If such a scanner position exists (due to some angular constraints related to the orientation of the spacecraft), the program tests this position through a direct computation to find whether the line of sight is within some input maximum separation from the star. At the end of the search, it shows a table with some essential information on the observation(s) (date and time, ephemeris time, name and HR number of the star, scanner position, angular separation from the line of sight). If the user is interested, they can save the table as a `.csv` file or, by right-clicking on one row of the table, plot the celestial map of that observation in terms of longitude and latitude in ecliptic reference frame centered on MPO. In this case, the user is required to choose whether the slit is present ("across", with a field of view of 0.1°x2°) or not ("removed", with field of view of 1.4°x3.1°), and whether to show the positions of the Sun and the planets (Mercury and its disk, Venus and the Earth). If Mercury is occulting the Sun or the planets, their name labels are crossed out, while stars will have a different color.

The search is launched successfully if all the fields are filled in correctly.
For any combination of date/time and scanner position, the observation is potentially feasible if:
- i<sub>a</sub>) the direction MPO-Sun is more than 90° (92° for safety) away from the line of sight of PHEBUS, or 
- i<sub>b</sub>) the Sun is occulted to PHEBUS by the disk of Mercury or by the radiator of the MPO, 
- ii) the star is perfectly visible (not occulted by planets)
- iii) the star is within the angular distance from the line of sight specified by the user (typically, half width of the slit). 

The conditions i<sub>a</sub>/i<sub>b</sub> are due to the fact that the sunlight can illuminate the baffle of PHEBUS sideways, but the sunlight must never enter the baffle aperture directly.

In addition, until its separation from MPO (September 2026), the MTM module completely obstructs the field of view of PHEBUS for scanner positions between 2100 and 2850. After the separation, that interval of scanner positions will be perfectly viable. Also the positions 0-109 and 3991-4096 are not considered because they correspond to PHEBUS being in parking position, where it can't observe.



#### 2) Visualize observation 
The user is required to choose the following parameters:
- the interesting stars from a catalog. The user is free to select all the stars or only some of them; 
- the files with the list of kernels that must be loaded for the computation of the orbits. The file is called `bc_plan.tm` and is provided with the freeware [SPICE package](https://s2e2.cosmos.esa.int/bitbucket/projects/SPICE_KERNELS/repos/bepicolombo/browse) in the directory `/bepicolombo/kernels/mk/`. The package for BepiColombo takes up a few gigabytes of memory, so the user is expected to have it stored locally;   
- the date and time corresponding to the celestial map that the user wants to show. The date/time must be written in the form `YYYY-MM-DDThh:mm:ss`, though the characters `-`, `:` and the `T` can be replaced by any other character. Any date is valid, as long as the SPICE package finds data on the orbit of BepiColombo (since its launch and for the next few years), and it is transformed to ephemeris times (seconds from J2000);
- the value of the scanner position under which the line of sight of PHEBUS and the projection of the slit on the sky should be plotted. The position can span integers between 100 and 3990;
- whether the slit is present ("across", with a field of view of 0.1°x2°) or not ("removed", with field of view of 1.4°x3.1°);
- whether to show the Sun and the planets (Mercury and its disk, Venus, the Earth) in the plot;   	  
- the maximum angle, measured radially from the line of sight of PHEBUS, within which to show the name of a star. 
    
This part of the software simply makes a 2D plot (longitude and latitude in ecliptic reference frame centered on MPO) of the sky showing the stars, the Sun and planets if required, and the projection of the slit. The python function that makes this plot is the same that makes the plot mentioned at the end of part 1). However, whereas in part 1) the user could make a plot-map only for the potential observations that resulted from the star search, in part 2) the user is free to choose any single date/time and scanner position they want. A popup window also informs the user whether the Sun is a hazard for PHEBUS or the observation is not possible due to obstructions by the MTM.


### Future work
- Develop an executable version of this software that doesn't require the download of all the folder and the installation of python packages;
- Develop a similar software for the observation of Mercury and its exosphere. It will allow to search for observations of the planet according to a series of input parameters and it will plot 2D maps of the observation on the planet and 3D models of the relative geometry.


### Acronyms
- ESA : European Space Agency
- LATMOS : Laboratoire Atmosphères, Observations Spatiales
- MPO : Mercury Planetary Orbiter
- MTM : Mercury Transfer Module
- PHEBUS : Probing of Hermean Exosphere By Ultraviolet Spectroscopy
- PHESTAR : PHEBUS Stellar Targeting, Alignment and Reconnaissance
- SPICE : Spacecraft, Planet, Instrument, C-matrix, Events
