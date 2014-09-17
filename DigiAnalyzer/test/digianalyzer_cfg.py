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
        #'file:/afs/cern.ch/work/k/khurana/HGCAL/Digi/CMSSW_6_2_0_SLHC16/src/EventGeneration/HGC/FourMuEvents_140PU_200GenSimDigi.root'
        #'file:/afs/cern.ch/work/k/khurana/HGCAL/Digi/CMSSW_6_2_0_SLHC16/src/Validation/HGCalValidation/Muon_500Events.root'
        #'file:../../Samples/Muon_500Events.root'
        #'file:/afs/cern.ch/work/k/khurana/HGCAL/Digi/CMSSW_6_2_0_SLHC16/src/Samples/12200_FourMuPt1_200+FourMuPt_1_200_Extended2023HGCalMuon_GenSimFull+DigiFull_Extended2023HGCalMuon+RecoFull_Extended2023HGCalMuon+HARVESTFull_Extended2023HGCalMuon/step2_1event.root'
        #'file:/afs/cern.ch/work/k/khurana/HGCAL/Digi/CMSSW_6_2_0_SLHC16/src/Samples/step2.root'
        #'file:/afs/cern.ch/work/k/khurana/HGCAL/Digi/CMSSW_6_2_0_SLHC16/src/EventGeneration/HGC/step2_10FourMuEvents_140PU.root'
        
        #'file:/afs/cern.ch/work/k/khurana/HGCAL/Digi/CMSSW_6_2_0_SLHC16/src/EventGeneration/HGC/step2_SingleMuonEvents_140PU_NoShaping.root'
        #'file:/afs/cern.ch/work/k/khurana/HGCAL/Digi/CMSSW_6_2_0_SLHC16/src/EventGeneration/HGC/gen-sim-digi_Pt-10_Eta-1.75_ID-13.root'
        #'file:/tmp/khurana/gen-sim-digi-Pt-10.0_Eta-2.0__TauValue-20.0.root'
        #'file:/tmp/khurana/gen-sim-digi-Pt-10.0_Eta-2.0_140PU__TauValue-5.0.root'
        'root://eoscms.cern.ch//store/user/khurana/HGCAL/Digi/SingleParticle/INPUTFILENAME'
        )
                            )

process.digianalyzerEE = cms.EDAnalyzer('DigiAnalyzer',
                                        DetectorName = cms.string("HGCalEESensitive"),
                                        DigiSource = cms.string("HGCDigisEE"),
                                        Verbosity     = cms.untracked.int32(1)                                      
                                        )


process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("DigiHisto-INPUTFILENAME.root"),
                                   closeFileFast = cms.untracked.bool(True)
                                   )



process.p = cms.Path(process.digianalyzerEE )
