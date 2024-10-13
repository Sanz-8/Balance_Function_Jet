#!/bin/bash

rm -rf *.C~
rm -rf *.h~
rm -rf *.h~
rm -rf *.sh~

# Run this as --> sh Cond_Submit_CMSpythia.sh step3b 4 0

# 1st Arguments : which step (step1(TuneCP5), step2(TuneCP5_RH_CR2), step3a, step3b, step3c) 
# Output folder name automatic created in the destination path
# 2nd arguments : which file (1->TuneCP5.txt, 2->Step_2.txt, 3->Step_3a.txt, 4->Step_3b.txt, 5->Step_3c.txt)
# 3rd arguments : is test run (1 or 0)

which_step=$2
file="Step_1.txt"

if [ $which_step == 1 ]; then
    file="Pythia_pThat15.txt"
    echo ${file}
elif [ $which_step == 2 ]; then
    file="CMSDATA_pp_13TeV_Run2016_MiniAODV2.txt"
    echo ${file}
elif [ $which_step == 3 ]; then
    file="CMSDATA_pp_13TeV_Run2017_MiniAODV2.txt"
    echo ${file}
elif [ $which_step == 4 ]; then
    file="CMSDATA_pp_13TeV_Run2018_MiniAODV2.txt"
    echo ${file}
elif [ $which_step == 5 ]; then
    file="Default_WOCR.txt"
    echo ${file}
else
    echo "There are not step included"
fi

istestrun=$3
testRun_condition=true

if [ $istestrun == 1 ]; then
    testRun_condition=true
else
    testRun_condition=false
fi

root -l -b -q "File_Write.C(${testRun_condition}, \"${file}\", \"$1\", \"/eos/user/s/slaishra/CMSpythia_injet_OutPutFiles\", 4, \"workday/\")"

condor_submit Submission_Files/$1/Submission_$1.sub

