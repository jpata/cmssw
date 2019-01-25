import FWCore.ParameterSet.Config as cms

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
    ptbins = [
        10,24,32,43,56,74,97,133,174,245,300,362,430,507,592,686,846,1032,1248,1588,
        2000,2500,3000,4000,6000
    ]
    etabins = [0.0, 0.5, 1.3, 2.1, 2.5, 3.0]

    response_plots = []
    #we always use a range [ibin, ibin+1) 
    for ietabin in range(len(etabins)-1):
        for iptbin in range(len(ptbins)-1):
            #convert 0.5 -> "05"
            eta_string = "{0:.1f}".format(etabins[ietabin+1]).replace(".", "")

            response_plots += [ResponsePlot(
                "reso_dist_{0:.0f}_{1:.0f}_eta{2}".format(ptbins[iptbin], ptbins[iptbin+1], eta_string),
                "Jet response (pT/pTgen) in {0} <= pt < {1}, {2} <= |eta| < {3}".format(ptbins[iptbin], ptbins[iptbin+1], etabins[ietabin], etabins[ietabin+1]),
                100, 0.0, 3.0, ptbins[iptbin], ptbins[iptbin+1], etabins[ietabin], etabins[ietabin+1]
            )]
    return response_plots

pfDQM = cms.EDAnalyzer("ParticleFlowDQM",

    #match these reco-jets to the gen-jets and compute jet response
    recoJetCollection = cms.InputTag('ak4PFJets'),
    genJetCollection = cms.InputTag('ak4GenJets'),
    jetDeltaR = cms.double(0.25),

    responsePlots = cms.VPSet(createResponsePlots())

)
