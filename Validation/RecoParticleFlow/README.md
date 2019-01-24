
# Quickstart

1. Set up using `source test/setup_cmssw.sh`
2. Merge the pfvalidation branch using `git cms-merge-topic jpata:pfvalidation`
3. Make a temporary directory in e.g. `$CMSSW_BASE/src/Validation/RecoParticleFlow/tmp`
4. Run the RECO sequence using `$CMSSW_BASE/src/Validation/RecoParticleFlow/test/run_relval.sh QCD reco $NJOB` where $NJOB is an integer in [0,249] that is used to select the events for this job, assuming 200 events per job and a RelVal dataset of 50k events
5. Run the DQM sequence using `test/run_relval.sh QCD dqm 0`, in this case `$NJOB` has no relevance 
6. Produce the plots using `test/compare.py`

# Running on condor
The reco sequence takes about 1-2 hours / 100 events on batch. We have prepared condor scripts to facilitate this on lxbatch. In order to submit the condor jobs, go to a temporary working directory e.g. `$CMSSW_BASE/src/Validation/RecoParticleFlow/tmp/QCD`, create a log directory with `mkdir log` and submit the condor jobs using `condor_submit $CMSSW_BASE/src/Validation/RecoParticleFlow/test/condor_sub.jdl`. The input files for the jobs are configured in the `run_relval.sh` script. Once the jobs are done, the DQM sequence is able to use all step3 EDM files in a single job.
