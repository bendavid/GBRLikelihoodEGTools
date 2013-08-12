import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration/StandardSequences/MagneticField_38T_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('Configuration/EventContent/EventContent_cff')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

#process.source = cms.Source("PoolSource",
    ## replace 'myfile.root' with the source file you want to use
    #fileNames = cms.untracked.vstring(
        #'file:myfile.root'
    #)
#)




process.source = cms.Source("PoolSource",
   fileNames = cms.untracked.vstring('/store/relval/CMSSW_5_3_6-START53_V14/RelValH130GGgluonfusion/GEN-SIM-RECO/v2/00000/202DD4DB-F929-E211-8F53-001A92810AF2.root'),
)

        

process.egenergyanalyzer = cms.EDAnalyzer('EGEnergyAnalyzerSemiParm'
)


process.GlobalTag.globaltag = 'START53_V14::All'    

                    
process.p = cms.Path(process.egenergyanalyzer)