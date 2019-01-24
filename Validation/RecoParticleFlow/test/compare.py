#!/usr/bin/env python
import sys
import os
import ROOT

# This is an example of plotting the standard tracking validation
# plots from an explicit set of DQM root files.

from Validation.RecoTrack.plotting.validation import SimpleValidation, SimpleSample

from Validation.RecoTrack.plotting.plotting import Subtract, FakeDuplicate, CutEfficiency, Transform, AggregateBins, ROC, Plot, PlotEmpty, PlotGroup, PlotOnSideGroup, PlotFolder, Plotter
from Validation.RecoTrack.plotting.html import PlotPurpose

outputDir = "plots" # Plot output directory
description = "Simple ParticleFlow comparison"

plotterDrawArgs = dict(
    separate=False, # Set to true if you want each plot in it's own canvas
#    ratio=False,   # Uncomment to disable ratio pad
)

def parse_sample_string(ss):
    spl = ss.split(":")
    if not (len(spl) >= 3):
        raise Exception("Sample must be in the format name:DQMfile1.root:DQMfile2.root:...")
    
    name = spl[0]
    files = spl[1:]
    for fi in files:
        if not os.path.isfile(fi):
            raise FileError("Could not read DQM file {0}".format(fi))
    return name, files
  
def parse_args():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--sample", type=str, action='append')
    parser.add_argument("-p", "--plots", type=str, action='append', default=None)
    args = parser.parse_args()

    #collect all the SimpleSample objects    
    samples = []
    plots = []
 
    sample_strings = args.sample
    for ss in sample_strings:
        name, files = parse_sample_string(ss)
        samp = SimpleSample(name, name, [(fn, "Option {0}".format(i)) for fn, i in zip(files, range(len(files)))])
        samples += [samp]
    
    if not (args.plots is None):
        pass
 
    return samples, plots

samples, plots = parse_args()

def getall(d, basepath="/"):
    "Generator function to recurse into a ROOT file/dir and yield (path, obj) pairs"
    for key in d.GetListOfKeys():
        kname = key.GetName()
        if key.IsFolder():
            for i in getall(d.Get(kname), basepath+kname+"/"):
                yield i
        else:
            yield basepath+kname, d.Get(kname).ClassName()

if len(plots) == 0:
    for samp in samples:
        for fn, label in samp._fileLegends:
            print fn
            tf = ROOT.TFile(fn)
            tf.ReadAll()
            for name, objtype in getall(tf):
                if objtype.startswith("TH1"):
                    print name

 
def addPlots(plotter, folder, name, section, bin_range):
    folders = [folder]
    plots = [PlotGroup(name, [Plot("Bin{0}".format(ibin)) for ibin in bin_range])]
    plotter.append("ParticleFlow", folders, PlotFolder(*plots, loopSubFolders=False, page="pf", section=section))

plotter = Plotter()

addPlots(plotter, "DQMData/Run 1/Physics/Run summary/JetResponse/ByGenJetPt", "ByGenJetPt1", "JetResponse", range(0,6))
addPlots(plotter, "DQMData/Run 1/Physics/Run summary/JetResponse/ByGenJetPt", "ByGenJetPt2", "JetResponse", range(6,12))

addPlots(plotter, "DQMData/Run 1/Physics/Run summary/JetResponse/ByGenJetEta", "ByGenJetEta1", "JetResponse", range(0,6))
addPlots(plotter, "DQMData/Run 1/Physics/Run summary/JetResponse/ByGenJetEta", "ByGenJetEta2", "JetResponse", range(6,12))
addPlots(plotter, "DQMData/Run 1/Physics/Run summary/JetResponse/ByGenJetEta", "ByGenJetEta3", "JetResponse", range(12,13))

val = SimpleValidation(samples, outputDir)
report = val.createHtmlReport(validationName=description)
val.doPlots([plotter],
    plotterDrawArgs=plotterDrawArgs,
)
report.write()
