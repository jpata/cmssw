#! /usr/bin/env python
import sys

from PhysicsTools.Heppy.utils.cmsswPreprocessor import CmsswPreprocessor

from vhbb import *
from VHbbAnalysis.Heppy.AdditionalBTag import AdditionalBTag
from VHbbAnalysis.Heppy.AdditionalBoost import AdditionalBoost
from VHbbAnalysis.Heppy.GenHFHadronMatcher import GenHFHadronMatcher


# Add Boosted Information

boostana=cfg.Analyzer(
    verbose=False,
    class_object=AdditionalBoost,
)
sequence.insert(sequence.index(VHbb),boostana)

genhfana=cfg.Analyzer(
    verbose=False,
    class_object=GenHFHadronMatcher,
)
sequence.insert(sequence.index(VHbb),genhfana)

treeProducer.collections["ungroomedFatjets"] = NTupleCollection("ungroomedFatjets",
                                                                fatjetType,
                                                                10,
                                                                help = "CA, R=1.5, pT > 200 GeV, no grooming")

treeProducer.collections["trimmedFatjets"] = NTupleCollection("trimmedFatjets",
                                                              fourVectorType,
                                                              10,
                                                              help="CA, R=1.5, pT > 200 GeV, trimmed with R=0.2 and f=0.06")

treeProducer.collections["httCandidates"] = NTupleCollection("httCandidates",
                                                             httType,
                                                             10,
                                                             help="MultiR HEPTopTagger Candidates")

# Add b-Tagging Information

btagana=cfg.Analyzer(
    verbose=False,
    class_object=AdditionalBTag,
)
sequence.insert(sequence.index(VHbb),btagana)

# Add Information on generator level hadronic tau decays

from VHbbAnalysis.Heppy.TauGenJetAnalyzer import TauGenJetAnalyzer
TauGenJet = cfg.Analyzer(
    verbose = False,
    class_object = TauGenJetAnalyzer,
)
sequence.insert(sequence.index(VHbb),TauGenJet)

treeProducer.collections["tauGenJets"] = NTupleCollection("GenHadTaus", genTauJetType, 15, help="Generator level hadronic tau decays")

# Run Everything

import sys
if len(sys.argv) == 3:
    out = sys.argv[2]
else:
    out = "Loop"

preprocessor = CmsswPreprocessor("combined_cmssw.py")
config.preprocessor=preprocessor
if __name__ == '__main__':
    from PhysicsTools.HeppyCore.framework.looper import Looper
    looper = Looper( out, config, nPrint = 1, nEvents = 1000000)
    import time
    import cProfile
    p = cProfile.Profile(time.clock)
    p.runcall(looper.loop)
    p.print_stats()
    looper.write()
