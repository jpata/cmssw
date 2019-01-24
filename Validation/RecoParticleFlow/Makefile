TMPDIR=tmp
RELVALCMD=${CMSSW_BASE}/src/Validation/RecoParticleFlow/test/run_relval.sh

all: TTbar QCD ZMM
dqm: TTbar_dqm QCD_dqm ZMM_dqm

TTbar: TTbar_reco TTbar_dqm
QCD: QCD_reco QCD_dqm
ZMM: ZMM_reco ZMM_dqm

TTbar_reco:
	cd ${TMPDIR} && ${RELVALCMD} TTbar reco 0
TTbar_dqm:
	cd ${TMPDIR} && ${RELVALCMD} TTbar dqm 0

QCD_reco:
	cd ${TMPDIR} && ${RELVALCMD} QCD reco 0
QCD_dqm:
	cd ${TMPDIR} && ${RELVALCMD} QCD dqm 0

ZMM_reco:
	cd ${TMPDIR} && ${RELVALCMD} ZMM reco 0
ZMM_dqm:
	cd ${TMPDIR} && ${RELVALCMD} ZMM dqm 0

.PHONY: TTbar_reco QCD_reco ZMM_reco TTbar_dqm QCD_dqm ZMM_dqm
