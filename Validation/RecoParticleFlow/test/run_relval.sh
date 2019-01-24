#!/bin/bash
#Script to run RECO and DQM sequences on existing files using cmsDriver.py
#More background information: 
#https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideCmsDriver
#https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookDataFormats

#env

#abort on error
set -e


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
    echo "Must pass exactly 3 arguments: run_relval.sh [TTbar|QCD|ZMM] [reco|dqm] [njob]"
    exit 0
fi

#index of the job is used to keep track of which events / files to process in the reco step
NJOB=$3

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
if [ "$1" == "TTbar" ]; then
    INPUT_FILE=root://cms-xrd-global.cern.ch///store/relval/CMSSW_9_4_11_cand2/RelValTTbar_13/GEN-SIM-DIGI-RAW/94X_mc2017_realistic_v15-v1/10000/FA0F368D-D4D2-E811-A159-0025905B8600.root
    NAME=TTbar
    echo "TTbar sample needs to be implemented"
    exit 1
elif [ "$1" == "QCD" ]; then
    INPUT_FILELIST=${CMSSW_BASE}/src/Validation/RecoParticleFlow/test/filelist_QCD.txt
    NAME=QCD
elif [ "$1" == "ZMM" ]; then
    INPUT_FILE=root://cms-xrd-global.cern.ch///store/relval/CMSSW_9_4_11_cand2/RelValZMM_13/GEN-SIM-DIGI-RAW/PU25ns_94X_mc2017_realistic_v15-v1/10000/E81E507D-5DD4-E811-9512-0025905A497A.root 
    NAME=ZMM
    echo "ZMM sample needs to be implemented"
    exit 1
else
    echo "Argument 1 must be [TTbar|QCD|ZMM] but was $1"
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


    #Run the actual CMS reco with particle flow.
    #On lxplus, this step takes about 15 minutes / 1000 events
    echo "Running step RECO" 
    cmsDriver.py step3  --customise_commands="process.source.skipEvents=cms.untracked.uint32($SKIPEVENTS)" --runUnscheduled  --conditions $CONDITIONS -s RAW2DIGI,L1Reco,RECO,RECOSIM,EI,PAT --datatier RECOSIM --nThreads $NTHREADS -n $PERJOB --era $ERA --eventcontent RECOSIM --filein filelist:$INPUT_FILELIST --fileout file:step3.root | tee step3.log  2>&1
   
    #NanoAOD
    #On lxplus, this step takes about 1 minute / 1000 events
    #Can be skipped if doing DQM directly from RECO
    #cmsDriver.py step4 --conditions $CONDITIONS -s NANO --datatier NANOAODSIM --nThreads $NTHREADS -n $N --era $ERA --eventcontent NANOAODSIM --filein file:step3_inMINIAODSIM.root --fileout file:step4.root > step4.log 2>&1
elif [ $STEP == "DQM" ]; then
    echo "Running step DQM" 

    cd $NAME
    
    #get all the filenames and make them into a python-compatible list of strings
    STEP3FNS=`ls -1 step3*.root | sed 's/^/"file:/;s/$/",/' | tr '\n' ' '`
    echo "step3 filenames for DQM: "$STEP3FNS

    #Run the DQM sequences (PF DQM only)
    #override the filenames here as cmsDriver does not allow multiple input files and there is no easy way to merge EDM files
    cmsDriver.py step5 --customise_commands="process.source.fileNames=cms.untracked.vstring($STEP3FNS)" --conditions $CONDITIONS -s DQM:@pfDQM --datatier DQMIO --nThreads $NTHREADS --era $ERA --eventcontent DQM --filein step3*.root --fileout file:step5.root | tee step5.log 2>&1

    #Harvesting converts the histograms stored in TTrees to be stored in folders by run etc
    cmsDriver.py step6 --conditions $CONDITIONS -s HARVESTING:@pfDQM --era $ERA --filetype DQM --filein file:step5.root --fileout file:step6.root | tee step6.log 2>&1
fi

#echo "Exit code was $?"
#tail -n3 *.log

cd ..

find . -name "*"
