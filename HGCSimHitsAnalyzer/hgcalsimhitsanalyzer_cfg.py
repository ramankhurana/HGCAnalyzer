import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("IOMC.EventVertexGenerators.VtxSmearedGauss_cfi")
process.load("Geometry.HGCalCommonData.testHGCXML_cfi")
process.load("Geometry.TrackerNumberingBuilder.trackerNumberingGeometry_cfi")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.EventContent.EventContent_cff")
process.load("SimG4Core.Application.g4SimHits_cfi")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( 
    input = cms.untracked.int32(10)
    )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
    'file:gensimeElevent100GeV.root'
    #'file:EventGenerator/MIPStudy/GenSimEvents_Eta_175_PT_50.root'
    #'file:electron_sample_50_3.root'
    #'file:/afs/cern.ch/work/k/kchatter/public/EventGen/electron_sample_70_3.root'
    #'file:/afs/cern.ch/work/k/kchatter/public/EventGen/gensimPionPlusevent100GeVEta20_NoEE.root'
    #'file:/afs/cern.ch/work/k/kchatter/public/EventGen/gensimPionPlusevent100GeVEta20_NoEE_NoHEF.root'
    #'file:gensimPionPlusevent100GeVEta20_NoEE.root'
    #'file:gensimPionPlusevent100GeVEta20_NoEE_NoHEF.root'
    )
                            )

process.demo = cms.EDAnalyzer('HGCSimHitsAnalyzer',
                              CaloHitSource = cms.string("HGCHitsEE"),
                              DetectorName = cms.string("HGCalEESensitive"),
                              
                              #CaloHitSource = cms.string("HGCHitsHEfront"), 
                              #DetectorName = cms.string("HGCalHESiliconSensitive"),
                              
                              #CaloHitSource = cms.string("HGCHitsHEback"), 
                              #DetectorName = cms.string("HGCalHEScintillatorSensitive"),
                              
                              OutFile      = cms.string("OUTPUTFILE.root"),
                              )


process.EE = process.demo.clone()
process.EE.CaloHitSource = cms.string("HGCHitsEE")
process.EE.DetectorName = cms.string("HGCalEESensitive")
process.EE.OutFile      = cms.string("OUTPUTFILEEE.root")

process.HEfront = process.demo.clone()
process.HEfront.CaloHitSource = cms.string("HGCHitsHEfront")
process.HEfront.DetectorName = cms.string("HGCalHESiliconSensitive")
process.HEfront.OutFile = cms.string("OUTPUTFILEEEHEf.root")

process.HEback = process.demo.clone()
process.HEback.CaloHitSource = cms.string("HGCHitsHEback")
process.HEback.DetectorName = cms.string("HGCalHEScintillatorSensitive")
process.HEback.OutFile = cms.string("OUTPUTFILEEEHEb.root")

process.p = cms.Path(process.EE
                     #process.HEfront+
                     #process.HEback
                     )
