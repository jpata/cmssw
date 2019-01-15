#!/bin/bash
#https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideCmsDriver
#https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookDataFormats
#abort on error
set -e

N=5000
CONDITIONS=auto:phase1_2017_realistic
ERA=Run2_2017,run2_nanoAOD_94XMiniAODv1
NTHREADS=10

#Argument parsing
if [ "$#" -ne 2 ]; then
    echo "Must pass exactly 2 arguments: run_relval.sh [TTbar|QCD|ZMM] [reco|dqm]"
    exit 0
fi

##RelVal samples
if [ "$1" == "TTbar" ]; then
    INPUT_FILE=root://cms-xrd-global.cern.ch///store/relval/CMSSW_9_4_11_cand2/RelValTTbar_13/GEN-SIM-DIGI-RAW/94X_mc2017_realistic_v15-v1/10000/FA0F368D-D4D2-E811-A159-0025905B8600.root
    NAME=TTbar
elif [ "$1" == "QCD" ]; then
    INPUT_FILE=root://cms-xrd-global.cern.ch///store/relval/CMSSW_9_4_11_cand2/RelValQCD_FlatPt_15_3000HS_13/GEN-SIM-DIGI-RAW/PU25ns_94X_mc2017_realistic_v15-v1/10000/F2FEBFDF-69D4-E811-ADF1-0CC47A4C8F0C.root 
    NAME=QCD
elif [ "$1" == "ZMM" ]; then
    INPUT_FILE=root://cms-xrd-global.cern.ch///store/relval/CMSSW_9_4_11_cand2/RelValZMM_13/GEN-SIM-DIGI-RAW/PU25ns_94X_mc2017_realistic_v15-v1/10000/E81E507D-5DD4-E811-9512-0025905A497A.root 
    NAME=ZMM
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
    cmsDriver.py step3  --runUnscheduled  --conditions $CONDITIONS -s RAW2DIGI,L1Reco,RECO,RECOSIM,EI,PAT --datatier RECOSIM,AODSIM,MINIAODSIM --nThreads $NTHREADS -n $N --era $ERA --eventcontent RECOSIM,AODSIM,MINIAODSIM --filein file:$INPUT_FILE --fileout file:step3.root > step3.log  2>&1
   
    #NanoAOD
    #On lxplus, this step takes about 1 minute / 1000 events
    #Can be skipped if doing DQM directly from RECO
    #cmsDriver.py step4 --conditions $CONDITIONS -s NANO --datatier NANOAODSIM --nThreads $NTHREADS -n $N --era $ERA --eventcontent NANOAODSIM --filein file:step3_inMINIAODSIM.root --fileout file:step4.root > step4.log 2>&1
elif [ $STEP == "DQM" ]; then
    echo "Running step DQM" 

    cd $NAME

    #Run the DQM sequences (PF only)
    cmsDriver.py step5 --conditions $CONDITIONS -s DQM:@pfDQM --datatier DQMIO --nThreads $NTHREADS -n $N --era $ERA --eventcontent DQM --filein file:step3.root --fileout file:step5.root > step5.log 2>&1

    #Harvesting converts the histograms stored in TTrees to be stored in folders by run etc
    cmsDriver.py step6 --conditions $CONDITIONS -s HARVESTING:@pfDQM --era $ERA --filetype DQM --filein file:step5.root --fileout file:step6.root > step6.log 2>&1
fi

echo "Exit code was $?"
tail -n3 *.log

cd ..
