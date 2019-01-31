import FWCore.ParameterSet.Config as cms
from Validation.RecoParticleFlow.defaults_cfi import ptbins, etabins, response_distribution_name

def ResponsePlot(name, title, responseNbins, responseLow, responseHigh, ptBinLow, ptBinHigh, etaBinLow, etaBinHigh):
    return cms.PSet(
        name = cms.string(name),
        title = cms.string(title),
        responseNbins = cms.uint32(responseNbins),
        responseLow = cms.double(responseLow),
        responseHigh = cms.double(responseHigh),
        ptBinLow = cms.double(ptBinLow),
        ptBinHigh = cms.double(ptBinHigh),
        etaBinLow = cms.double(etaBinLow),
        etaBinHigh = cms.double(etaBinHigh),
    )

#Jet response is plotted in histograms which can be subdivided by pt and |eta| of the genjet.
#To minimize the amount of logic on the C++ side, we define all response plots here.
#Each plot has low and high pt and |eta| edges, the plot is filled only if the genjet
#is in the bin defined by the edges.
#It is your job here to make sure you define the bins in a non-overlapping way if
#you want to emulate a 2D map over (pT, |eta|) of 1D histograms.
def createResponsePlots():
    response_plots = []
    #we always use a range [ibin, ibin+1) 
    for ietabin in range(len(etabins)-1):
        for iptbin in range(len(ptbins)-1):

            response_plots += [ResponsePlot(
                response_distribution_name(iptbin, ietabin),
                "Jet response (pT/pTgen) in {0} <= pt < {1}, {2} <= |eta| < {3}".format(ptbins[iptbin], ptbins[iptbin+1], etabins[ietabin], etabins[ietabin+1]),
                100, 0.0, 3.0, ptbins[iptbin], ptbins[iptbin+1], etabins[ietabin], etabins[ietabin+1]
            )]
    return response_plots

#matchRecoJetToGenJet = cms.EDProducer('MatchRecToGen',
#        srcGen = cms.InputTag('ak4PFJets'),
#        srcRec = cms.InputTag('ak4GenJets')
#    )

pfDQMAnalyzer = cms.EDProducer("ParticleFlowDQM",

    #match these reco-jets to the gen-jets and compute jet response
    recoJetCollection = cms.InputTag('slimmedJets'),
    genJetCollection = cms.InputTag('slimmedGenJets'),
    jetDeltaR = cms.double(0.2),

    responsePlots = cms.VPSet(createResponsePlots())

)

pfDQM = cms.Sequence(
#    matchRecoJetToGenJet *
    pfDQMAnalyzer
)
