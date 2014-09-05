import FWCore.ParameterSet.Config as cms

process = cms.Process("Digianalyzer")

process.load("FWCore.MessageService.MessageLogger_cfi")
                                                                                               
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2023HGCalMuonReco_cff')
process.load ('Configuration.Geometry.GeometryExtended2023HGCalMuon_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedGauss_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.DigiToRaw_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                            # replace 'myfile.root' with the source file you want to use
                            fileNames = cms.untracked.vstring(
        #'file:/afs/cern.ch/work/k/khurana/HGCAL/Digi/CMSSW_6_2_0_SLHC16/src/Validation/HGCalValidation/Muon_500Events.root'
        'file:../../Samples/Muon_500Events.root'
        )
                            )

process.digianalyzerEE = cms.EDAnalyzer('DigiAnalyzer',
                                        DetectorName = cms.string("HGCalEESensitive"),
                                        DigiSource = cms.string("HGCDigisEE"),
                                        Verbosity     = cms.untracked.int32(1)                                      
                                        )


process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("DigiTree.root"),
                                   closeFileFast = cms.untracked.bool(True)
                                   )



process.p = cms.Path(process.digianalyzerEE )
