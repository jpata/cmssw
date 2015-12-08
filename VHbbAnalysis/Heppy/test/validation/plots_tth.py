import ROOT
ROOT.gROOT.SetBatch(True)
import sys, os
sys.path.append(os.environ.get("CMSSW_BASE") + "/src/VHbbAnalysis/Heppy/test")
from vhbb_combined import *
from PhysicsTools.HeppyCore.framework.looper import Looper

dirs = ["Loop_validation_tth_sl_dl_tth_hbb"]

vars_to_plot = [
    "Jet_pt",
    "Jet_eta",
    "Jet_btagCSV",
    "Jet_btagBDT",
    "Jet_mcFlavour",
    "Jet_bTagWeight",
    "Jet_bTagWeightJESUp",
    "Jet_bTagWeightJESDown",
    "Jet_bTagWeightHFUp",
    "Jet_bTagWeightHFDown",
    "Jet_bTagWeightLFUp",
    "Jet_bTagWeightLFDown",
    "Jet_bTagWeightStats1Up",
    "Jet_bTagWeightStats1Down",
    "Jet_bTagWeightStats2Up",
    "Jet_bTagWeightStats2Down",
    "nGenJet",
    "GenJet_pt",
    "GenJet_numBHadrons",
    "GenJet_numCHadrons",
]


components = [
    cfg.MCComponent(
        files = [
            "root://xrootd-cms.infn.it///store/mc/RunIISpring15DR74/ttHTobb_M125_13TeV_powheg_pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/00000/141B9915-1F08-E511-B9FF-001E675A6AB3.root"
        ],
        name = "tth_hbb",
        isMC = True
    ),
    cfg.MCComponent(
        files = [
            "root://xrootd-cms.infn.it///store/mc/RunIISpring15DR74/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v2/00000/0AB045B5-BB0C-E511-81FD-0025905A60B8.root"
        ],
        name = "ttjets",
        isMC = True
    )
]

def process_dir(d):
    print "Processing",d
    tf = ROOT.TFile(d + "/tree.root")
    tt = tf.Get("tree")
    if tt.GetEntries() <= 100:
        print "WARN: low efficiency", d
   
    npos = tf.Get("CountPosWeight").GetBinContent(1)
    nneg = tf.Get("CountNegWeight").GetBinContent(1)
    ntot = npos + nneg
    print "Ngen", ntot, npos, nneg

    for v in vars_to_plot:
        tt.Draw(v + " >> h")
        h = tf.Get("h")
        print v, round(h.Integral(), 2), round(h.GetMean(), 2), round(h.GetRMS(), 2) 

if __name__ == '__main__':
    for comp in components:
        print "processing",comp
        config.components = [comp] 
        looper = Looper( 'Loop_validation_tth_sl_dl_' + comp.name, config, nPrint = 0, nEvents = 10)
        looper.loop()
        looper.write()

    for d in dirs:
        process_dir(d)
