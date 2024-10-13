#!/bin/bash

rm -rf *.C~
rm -rf *.h~
rm -rf *.sh~
rm -rf *.txt~

# 1st Arguments : if you want to test run, true. It will take only 4 root files. For full jobs, it will be false
# 2nd Arguments : Your input txt file, where all the root files are there
# 3rd Arguments : Which step we are doing or some name

root -l -b -q "Create_subtxtFile.C(true, \"CMSDATA_pp_13TeV_Run2016_MiniAODV2_quickcheck.txt\", \"Jets_TuneCP5_CMSdata\")"

echo "Again enter the number of jobs you have given:"
read njobs

# 1st Arguments : Divided txt files are stored in the path "Divided_txt_files/Step1/". "Outfile_Step1" is the name of recreated txt files.
# Note that the 3rd argument of the macro "Create_subtxtFile.C", which is "Step1" in this case, will be the same in the path after "Divided_txt_files/", e.g., Divided_txt_files/Step1.
# 2nd Arguments : Number of txt files created. You will enter this when running this script. It will be taken automatically
# 4th Arguments : Output root file name. Better to give the same name as in the 3rd argument name of the macro "Create_subtxtFile.C", which is "Step1"

root -l -b -q "Tree_Analyzer_NewNew_Condor_CMSpythia.C(\"Divided_txt_files/Jets_TuneCP5_CMSdata/Outfile_Jets_TuneCP5_CMSdata_0.txt\", $njobs, \"/eos/user/s/slaishra/CMSpythia_injet_OutPutFiles\")" # new mac
# root -l -b -q "Tree_Analyzer_NewNew_Local.C(\"Divided_txt_files/Jets_TuneCP5/Outfile_Jets_TuneCP5\", $njobs, \"/eos/user/r/rpradhan/Research/With_NirbhayBhai/Out_RootFile\")" # new mac
# root -l -b -q "Tree_Analyzer_NewNew.C(\"Divided_txt_files/Jets_TuneCP5/Outfile_Jets_TuneCP5\", $njobs, \"/Users/raghumatapita/cernbox/Research/With_NirbhayBhai/Out_RootFile\")" # new mac
# root -l -b -q "Tree_Analyzer_NewNew.C(\"Divided_txt_files/Jets_TuneCP5/Outfile_Jets_TuneCP5\", $njobs, \"/Users/PIITAMATA/cernbox/Research/With_NirbhayBhai/Out_RootFile\")" # new mac

