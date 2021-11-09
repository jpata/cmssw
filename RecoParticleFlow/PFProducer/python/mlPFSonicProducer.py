import FWCore.ParameterSet.Config as cms

mlPFSonicProducer = cms.EDProducer("MLPFSonicProducer",
    Client = cms.PSet(
        timeout = cms.untracked.uint32(300),
        mode = cms.string("Async"),
        modelName = cms.string("mlpf_acat2021"),
        modelConfigPath = cms.FileInPath("HeterogeneousCore/SonicTriton/data/models/mlpf_acat2021/config.pbtxt"),
        verbose = cms.untracked.bool(False),
        allowedTries = cms.untracked.uint32(0),
        useSharedMemory = cms.untracked.bool(True),
        compression = cms.untracked.string(""),
    )
)
