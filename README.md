# Balance_Function_Jet
This code is for analysing data(p-p) inside the jet using the Balance Function Method. The data are from CERN CMS' run of 2016,2017, 2018.
The set up is usually implemented in CERN's computing platform(LXPLUS). Root Framework is used to analyse the data.
Below is the brief description of the codes.
1) Tree_Analyzer_NewNew_Condor_CMSpythia.C-> It helps to process  and  filter  the data and and stored in a 3d Vector(event->jet->tracks).
2) read_tree_New_CMSpythia.h-> This header file is used to retreive the data from a tree(a data structure in Root framework).
3) histogram_definition_CMSpythia.h->It declares all the histogram required for the analysis.
4) function_definition_CMSpythia.h->This is one of the main code that have different functions for correlation measurement.
5) input_variables_CMSpythia.h->It declares the important variables for the analysis such as the jet's radius, jet's momentum etc.
6) coordinateTools.h-> This code is used to transform the tracks(particle) coordinates from detector's coordinates to jet's coordinates. This is not used in this datasets as they already have the transformed coordinates,nevertheless it will be useful in understanding how they transform the coordinates.
7) Condor_Local_Run.sh->This shell file runs the program locally, it helps to check the results of the code before submitting to the Cern's computing grid.
8) Condor_Submit.sh->After getting the correct result from local run, the task is submitted to the Cern's Computing grid.
9) File_Write.C->This code keeps track of the submitted tasks.
10) .txt file-> This are the root file(source file), from which data are retrieved.Each of them are labelled according to the CMS's run year.


