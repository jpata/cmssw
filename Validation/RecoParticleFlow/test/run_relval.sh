#!/bin/bash
#Script to run RECO and DQM sequences on existing files using cmsDriver.py
#More background information: 
#https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideCmsDriver
#https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookDataFormats

hostname
env
voms-proxy-info

#abort on error
set +e


#number of events to process per job
#passed through condor, but here is a default value
if [ -z "$PERJOB" ]; then
    PERJOB=200
fi

#set default conditions
CONDITIONS=auto:phase1_2017_realistic
ERA=Run2_2017,run2_nanoAOD_94XMiniAODv1
NTHREADS=1

#Argument parsing
if [ "$#" -ne 3 ]; then
    echo "Must pass exactly 3 arguments: run_relval.sh [QCD|QCDPU] [reco|dqm] [njob]"
    exit 0
fi

#index of the job is used to keep track of which events / files to process in the reco step
NJOB=$(($3 + 1))

#set CMSSW environment and go to condor work dir
LAUNCHDIR=`pwd`
source /cvmfs/cms.cern.ch/cmsset_default.sh

#this environment variable comes from the condor submit script
cd $CMSSW_BASE
eval `scram runtime -sh`

#if the _CONDOR_SCRATCH_DIR is not defined, we are not inside a condor batch job
if [ -z "$_CONDOR_SCRATCH_DIR" ]; then
    cd $LAUNCHDIR
else
    cd $_CONDOR_SCRATCH_DIR
fi

##RelVal samples
if [ "$1" == "QCD" ]; then
    INPUT_FILELIST=${CMSSW_BASE}/src/Validation/RecoParticleFlow/tmp/das_cache/QCD_FlatPt_noPU/RelValQCD_FlatPt_15_3000HS_13__CMSSW_10_4_0_pre4-103X_mc2017_realistic_v2-v1__GEN-SIM-DIGI-RAW.txt
    NAME=QCD
elif [ "$1" == "QCDPU" ]; then
    INPUT_FILELIST=${CMSSW_BASE}/src/Validation/RecoParticleFlow/tmp/das_cache/QCD_FlatPt_PU25ns/RelValQCD_FlatPt_15_3000HS_13__CMSSW_10_4_0_pre4-PU25ns_103X_mc2017_realistic_v2-v1__GEN-SIM-DIGI-RAW.txt
    NAME=QCDPU
elif [ "$1" == "ZMM" ]; then
    INPUT_FILELIST=${CMSSW_BASE}/src/Validation/RecoParticleFlow/tmp/das_cache/ZMM/RelValZMM_13__CMSSW_10_4_0_pre4-103X_mc2017_realistic_v2-v1__GEN-SIM-DIGI-RAW.txt
    NAME=ZMM
elif [ "$1" == "MinBias" ]; then
    INPUT_FILELIST=${CMSSW_BASE}/src/Validation/RecoParticleFlow/tmp/das_cache/MinBias/RelValMinBias_13__CMSSW_10_4_0_pre4-103X_mc2017_realistic_v2-v1__GEN-SIM-DIGI-RAW.txt
    NAME=MinBias
else
    echo "Argument 1 must be [QCD|QCDPU|ZMM|MinBias] but was $1"
    exit 1
fi

##Which step to do
if [ "$2" == "reco" ]; then
    STEP="RECO"
elif [ "$2" == "dqm" ]; then
    STEP="DQM"
else
    echo "Argument 2 must be [reco|dqm] but was $2"
    exit 1
fi

#skip njob*perjob events
SKIPEVENTS=$(($NJOB * $PERJOB))

#Just print out environment last time for debugging
echo $INPUT_FILELIST $NAME $STEP $SKIPEVENTS
#env

if [ $STEP == "RECO" ]; then
    #Start of workflow
    echo "Making subdirectory $NAME"

    if [ -e $NAME ]; then
        echo "directory $NAME exists, aborting"
        exit 1
    fi

    mkdir $NAME
    cd $NAME

    FILENAME=`sed -n "${NJOB}p" $INPUT_FILELIST`
    echo "FILENAME="$FILENAME
    #Run the actual CMS reco with particle flow.
    echo "Running step RECO" 
    cmsDriver.py step3 --runUnscheduled  --conditions $CONDITIONS -s RAW2DIGI,L1Reco,RECO,RECOSIM,EI,PAT --datatier RECOSIM,AODSIM,MINIAODSIM --nThreads $NTHREADS -n -1 --era $ERA --eventcontent RECOSIM,AODSIM,MINIAODSIM --filein file:$FILENAME --fileout file:step3.root | tee step3.log  2>&1
   
    #NanoAOD
    #On lxplus, this step takes about 1 minute / 1000 events
    #Can be skipped if doing DQM directly from RECO
    #cmsDriver.py step4 --conditions $CONDITIONS -s NANO --datatier NANOAODSIM --nThreads $NTHREADS -n $N --era $ERA --eventcontent NANOAODSIM --filein file:step3_inMINIAODSIM.root --fileout file:step4.root > step4.log 2>&1
elif [ $STEP == "DQM" ]; then
    echo "Running step DQM" 

    cd $NAME
    
    #get all the filenames and make them into a python-compatible list of strings
    #STEP3FNS=`ls -1 step3*MINIAODSIM*.root | sed 's/^/"file:/;s/$/",/' | tr '\n' ' '`
    du step3*MINIAODSIM*.root | grep -v "^0" | awk '{print $2}' | sed 's/^/file:/' > step3_filelist.txt
    cat step3_filelist.txt 

    #Run the DQM sequences (PF DQM only)
    #override the filenames here as cmsDriver does not allow multiple input files and there is no easy way to merge EDM files
    cmsDriver.py step5 --conditions $CONDITIONS -s DQM:@pfDQM --datatier DQMIO --nThreads $NTHREADS --era $ERA --eventcontent DQM --filein filelist:step3_filelist.txt --fileout file:step5.root -n -1 | tee step5.log 2>&1

    #Harvesting converts the histograms stored in TTrees to be stored in folders by run etc
    cmsDriver.py step6 --conditions $CONDITIONS -s HARVESTING:@pfDQM --era $ERA --filetype DQM --filein file:step5.root --fileout file:step6.root | tee step6.log 2>&1
fi

#echo "Exit code was $?"
#tail -n3 *.log

cd ..

find . -name "*"
