autoDQM = { 'common' : ['DQMOfflineCommon',
                        'PostDQMOffline',
                        'DQMHarvestCommon+DQMCertCommon'],
            'commonSiStripZeroBias' : ['DQMOfflineCommonSiStripZeroBias',
                                       'PostDQMOffline',
                                       'DQMHarvestCommonSiStripZeroBias+DQMCertCommon'],
            'trackingOnlyDQM' : ["DQMOfflineTracking",
                                 "PostDQMOffline",
                                 "DQMHarvestTracking"],
            'outerTracker': ['DQMOuterTracker',
                             'PostDQMOffline',
                             'DQMHarvestOuterTracker'],
            'lumi': ['DQMOfflineLumi',
                     'PostDQMOffline',
                     'DQMHarvestLumi'],
            'muon': ['DQMOfflineMuon',
                     'PostDQMOffline',
                     'DQMHarvestMuon+DQMCertMuon'],
            'hcal':     ['DQMOfflineHcal',
                         'PostDQMOffline',
                         'DQMHarvestHcal'],
            'hcal2': ['HcalDQMOfflineSequence',
                      'PostDQMOffline',
                      'HcalDQMOfflinePostProcessor'],
            'jetmet':  ['DQMOfflineJetMET',
                        'PostDQMOffline',
                        'DQMHarvestJetMET+DQMCertJetMET'],
            'ecal':       ['DQMOfflineEcal',
                           'PostDQMOffline',
                           'DQMHarvestEcal+DQMCertEcal'],
            'egamma':       ['DQMOfflineEGamma',
                             'PostDQMOffline',
                           'DQMHarvestEGamma'],
            'btag':       ['DQMOfflineBTag',
                           'PostDQMOffline',
                           'DQMHarvestBTag'],
            'HLTMon':     ['HLTMonitoring',
                           'PostDQMOffline',
                           'HLTMonitoringClient'],
            'HLTMonPA' :  ['HLTMonitoringPA', 'PostDQMOffline', 'HLTMonitoringClientPA'],
            'express':       ['@commonSiStripZeroBias+@muon+@hcal+@jetmet+@ecal',
                              'PostDQMOffline',
                              '@commonSiStripZeroBias+@muon+@hcal+@jetmet+@ecal'],
            'allForPrompt':  ['@common+@muon+@hcal+@jetmet+@ecal+@egamma',
                              'PostDQMOffline',
                              '@common+@muon+@hcal+@jetmet+@ecal+@egamma'],
            'miniAODDQM': ['DQMOfflineMiniAOD',
                           'PostDQMOfflineMiniAOD',
                           'DQMHarvestMiniAOD'],
            'nanoAODDQM': ['DQMOfflineNanoAOD',
                           'PostDQMOffline',
                           'DQMHarvestNanoAOD'],
            'pfDQM': ['DQMOfflinePF',
                           'PostDQMOffline',
                           'DQMHarvestCommon'],
            'standardDQM': ['DQMOffline',
                            'PostDQMOffline',
                            'dqmHarvesting'],
            'ExtraHLT': ['DQMOfflineExtraHLT',
                         'PostDQMOffline',
                         'dqmHarvestingExtraHLT'],
            'standardDQMFakeHLT': ['DQMOfflineFakeHLT',
                                   'PostDQMOffline',
                                   'dqmHarvestingFakeHLT'],
            'liteDQMHI': ['liteDQMOfflineHeavyIons',
                          'PostDQMOffline',
                          'dqmHarvesting']
            }

_phase2_allowed = ['trackingOnlyDQM','outerTracker','muon','hcal','hcal2','egamma']
autoDQM['phase2'] = ['','','']
for i in [0,2]:
    autoDQM['phase2'][i] = '+'.join([autoDQM[m][i] for m in _phase2_allowed])
autoDQM['phase2'][1] = 'PostDQMOffline'
