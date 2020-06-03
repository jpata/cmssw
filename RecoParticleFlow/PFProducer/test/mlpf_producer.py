import FWCore.ParameterSet.Config as cms

process = cms.Process('MLPF')

process.task = cms.Task()

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load(
    'Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

process.maxEvents = cms.untracked.PSet(
    input=cms.untracked.int32(-1)
)

# Input source
process.source = cms.Source("PoolSource",
    secondaryFileNames=cms.untracked.vstring(),
    fileNames=cms.untracked.vstring("http://jpata.web.cern.ch/jpata/edm_pf_2ev.root"),
    skipEvents=cms.untracked.uint32(0)
)

process.mlpfproducer = cms.EDProducer("MLPFProducer",
    src=cms.InputTag("particleFlowBlock"),
    min_batch_size=cms.uint32(1000),
    model_path=cms.string("RecoParticleFlow/PFProducer/data/mlpf/mlpf_2020_05_19.pb")
) 

process.sequence = cms.Sequence(process.mlpfproducer)
process.p = cms.Path(process.sequence)
process.output = cms.OutputModule("PoolOutputModule",
    outputCommands=cms.untracked.vstring('keep *'),
    fileName=cms.untracked.string('test.root')
)

process.outpath  = cms.EndPath(process.output)
