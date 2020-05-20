import FWCore.ParameterSet.Config as cms

def customize_step3(process):
    process.mlpfproducer = cms.EDProducer("MLPFProducer",
        src=cms.InputTag("particleFlowBlock"),
        min_batch_size=cms.uint32(1000),
        model_path=cms.string("RecoParticleFlow/PFProducer/data/mlpf/mlpf_2020_05_19.pb")
    ) 
    process.ak4MLPFJets = process.ak4PFJets.clone()
    process.ak4MLPFJets.src = cms.InputTag("mlpfproducer")

    process.MINIAODSIMoutput.outputCommands.append('keep recoPFCandidates_*_*_*')
    process.MINIAODSIMoutput.outputCommands.append('keep *_ak4MLPFJets_*_*')
    process.MINIAODSIMoutput.outputCommands.append('keep *_ak4PFJets_*_*')
    process.mlpf_path = cms.Path(process.mlpfproducer*process.ak4MLPFJets)

    process.schedule.insert(-1, process.mlpf_path)  
    return process
