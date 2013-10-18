EJ-lab-code
===========

Collection of code I have written for EJ's lab.

The main file to run is /Matlab/motion_script.m

The options at the beginning of the file should be able to recreate any figure I ever created.

To start I would recommend setting run_opt.auto_set = true and all t/f run_opts to true.  To get things running on the new filesystem, the new paths to the data files I used will have to be hardcoded in.  The script and functions it calls are well documented and should be fairly self explanatory.  The plots generated should all be labeled well and should be able to be interpreted without description.  Also, note the sliders at the bottom of the plots which allow the user to vary different parameters depending on the plot.

For the purpose of this distribution I am leaving several files such as export_fig/* and loop_slider_n.m in the repository which I have refactored into my repository MarvinT/Useful-Matlab-Functions, however, if I do any work on this repository later I will remove the duplicate files.

Another file of note is /scripting/make_script.py.  This is simply a python file that generates a bash script to download the data files to a local directory, perform the proper whitening and sorting, and upload the analyzed files to the network drives.  I used this to really quickly generate scripts which would process whole datasets that I could run over night.  The script generation is fairly robust and checks for existence of files and folders to try and resume where it previously left off.  This file will of course be broken in the new data structure, but could prove extremely useful if fixed and expanded to include more functionality than just the analysis I was performing.
