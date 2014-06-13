import FWCore.ParameterSet.Config as cms
hgcalAnalyzer = cms.EDAnalyzer('HGCSimHitsAnalyzer',
                                       DetectorName = cms.string("HGCalEESensitive"),
                                       CaloHitSource = cms.string("HGCHitsEE"),
                                       Verbosity     = cms.untracked.int32(1)
)
