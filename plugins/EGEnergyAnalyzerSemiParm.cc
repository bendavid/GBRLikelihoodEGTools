// -*- C++ -*-
//
// Package:    EGEnergyAnalyzerSemiParm
// Class:      EGEnergyAnalyzerSemiParm
// 
/**\class EGEnergyAnalyzerSemiParm EGEnergyAnalyzerSemiParm.cc GBRWrap/EGEnergyAnalyzerSemiParm/src/EGEnergyAnalyzerSemiParm.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Josh Bendavid
//         Created:  Tue Nov  8 22:26:45 CET 2011
// $Id: EGEnergyAnalyzerSemiParm.cc,v 1.3 2011/12/14 21:08:11 bendavid Exp $
//
//
 

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "TFile.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
//#include "CondCore/DBOutputService/interface/PoolDBOutputService.h"
//#include "CondCore/DBCommon/interface/CoralServiceManager.h"

#include "CondCore/DBOutputService/interface/PoolDBOutputService.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "HiggsAnalysis/GBRLikelihoodEGTools/interface/EGEnergyCorrectorSemiParm.h"
#include "HiggsAnalysis/GBRLikelihoodEGTools/interface/EGEnergyCorrectorTraditional.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"


//
// class declaration
//

class EGEnergyAnalyzerSemiParm : public edm::EDAnalyzer {
   public:
      explicit EGEnergyAnalyzerSemiParm(const edm::ParameterSet&);
      ~EGEnergyAnalyzerSemiParm();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      EGEnergyCorrectorSemiParm corV4;
      EGEnergyCorrectorSemiParm corV5;
      EGEnergyCorrectorSemiParm corV6;
      EGEnergyCorrectorSemiParm corV7;      
      EGEnergyCorrectorSemiParm corV8;      
      
      EGEnergyCorrectorTraditional corV4T;

      // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
EGEnergyAnalyzerSemiParm::EGEnergyAnalyzerSemiParm(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed

}


EGEnergyAnalyzerSemiParm::~EGEnergyAnalyzerSemiParm()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
EGEnergyAnalyzerSemiParm::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

  if (!corV4.IsInitialized()) {
    corV4.Initialize("/afs/cern.ch/user/b/bendavid/CMSSWshapessl6/CMSSW_6_2_0_pre7/src/HiggsAnalysis/GBRLikelihoodEGTools/data/regweights_v4_forest_ph.root",4);
  }
  
  if (!corV5.IsInitialized()) {
    corV5.Initialize("/afs/cern.ch/user/b/bendavid/CMSSWshapessl6/CMSSW_6_2_0_pre7/src/HiggsAnalysis/GBRLikelihoodEGTools/data/regweights_v5_forest_ph.root",5);
  }
  
  if (!corV6.IsInitialized()) {
    corV6.Initialize("/afs/cern.ch/user/b/bendavid/CMSSWshapessl6/CMSSW_6_2_0_pre7/src/HiggsAnalysis/GBRLikelihoodEGTools/data/regweights_v6_8TeV_forest_ph.root",6);
  }  
  
  if (!corV7.IsInitialized()) {
    corV7.Initialize("/afs/cern.ch/user/b/bendavid/CMSSWshapessl6/CMSSW_6_2_0_pre7/src/HiggsAnalysis/GBRLikelihoodEGTools/data/regweights_v7_8TeV_forest_ph.root",7);
  } 
  
  if (!corV8.IsInitialized()) {
    corV8.Initialize("/afs/cern.ch/user/b/bendavid/CMSSWshapessl6/CMSSW_6_2_0_pre7/src/HiggsAnalysis/GBRLikelihoodEGTools/data/regweights_v8_8TeV_forest_ph.root",8);
  }    
  
  if (!corV4T.IsInitialized()) {
    corV4T.Initialize("/afs/cern.ch/user/b/bendavid/CMSSWshapessl6/CMSSW_6_2_0_pre7/src/HiggsAnalysis/GBRLikelihoodEGTools/data/regweights_v4_traditional_ph.root");
  }  

  // get photon collection
  Handle<reco::PhotonCollection> hPhotonProduct;
  iEvent.getByLabel("photons",hPhotonProduct);
  
  EcalClusterLazyTools lazyTools(iEvent, iSetup, edm::InputTag("reducedEcalRecHitsEB"), 
                                 edm::InputTag("reducedEcalRecHitsEE"));  
  
  Handle<reco::VertexCollection> hVertexProduct;
  iEvent.getByLabel("offlinePrimaryVerticesWithBS", hVertexProduct);      
  
  Handle<double> hRho;
  iEvent.getByLabel(edm::InputTag("kt6PFJets","rho"), hRho);  
  
  for (reco::PhotonCollection::const_iterator it = hPhotonProduct->begin(); it!=hPhotonProduct->end(); ++it) {
    double ecor, sigeovere, mean, sigma, alpha1, n1, alpha2, n2, pdfval;

    
    corV4.CorrectedEnergyWithErrorV4(*it, *hVertexProduct, *hRho, lazyTools, iSetup,ecor, sigma, alpha1, n1, alpha2, n2, pdfval);
    printf("V4:  sceta = %5f, default = %5f, corrected = %5f, sigmaE/E = %5f, alpha1 = %5f, n1 = %5f, alpha2 = %5f, n2 = %5f, pdfval = %5f\n", it->superCluster()->eta(), it->energy(),ecor,sigma,alpha1,n1,alpha2,n2,pdfval);
//     
    corV5.CorrectedEnergyWithErrorV5(*it, *hVertexProduct, *hRho, lazyTools, iSetup,ecor, sigma, alpha1, n1, alpha2, n2, pdfval);
    printf("V5:  sceta = %5f, default = %5f, corrected = %5f, sigmaE/E = %5f, alpha1 = %5f, n1 = %5f, alpha2 = %5f, n2 = %5f, pdfval = %5f\n", it->superCluster()->eta(), it->energy(),ecor,sigma,alpha1,n1,alpha2,n2,pdfval);
    
    corV6.CorrectedEnergyWithErrorV6(*it, *hVertexProduct, *hRho, lazyTools, iSetup,ecor, sigeovere, mean, sigma, alpha1, n1, alpha2, n2, pdfval);
    printf("V6:  sceta = %5f, default = %5f, corrected = %5f, sigmaE/E = %5f, alpha1 = %5f, n1 = %5f, alpha2 = %5f, n2 = %5f, pdfval = %5f, meancb = %5f, sigmacb = %5f\n", it->superCluster()->eta(), it->energy(),ecor,sigeovere,alpha1,n1,alpha2,n2,pdfval,mean,sigma);    
    
    corV7.CorrectedEnergyWithErrorV7(*it, *hVertexProduct, *hRho, lazyTools, iSetup,ecor, sigeovere, mean, sigma, alpha1, n1, alpha2, n2, pdfval);
    printf("V7:  sceta = %5f, default = %5f, corrected = %5f, sigmaE/E = %5f, alpha1 = %5f, n1 = %5f, alpha2 = %5f, n2 = %5f, pdfval = %5f, meancb = %5f, sigmacb = %5f\n", it->superCluster()->eta(), it->energy(),ecor,sigeovere,alpha1,n1,alpha2,n2,pdfval,mean,sigma);        
    
    corV8.CorrectedEnergyWithErrorV8(*it, *hVertexProduct, *hRho, lazyTools, iSetup,ecor, sigeovere, mean, sigma, alpha1, n1, alpha2, n2, pdfval);
    printf("V8:  sceta = %5f, default = %5f, corrected = %5f, sigmaE/E = %5f, alpha1 = %5f, n1 = %5f, alpha2 = %5f, n2 = %5f, pdfval = %5f, meancb = %5f, sigmacb = %5f\n", it->superCluster()->eta(), it->energy(),ecor,sigeovere,alpha1,n1,alpha2,n2,pdfval,mean,sigma);            
    
    corV4T.CorrectedEnergyWithErrorV4Traditional(*it, *hVertexProduct, *hRho, lazyTools, iSetup,ecor, sigma);
    printf("V4T: sceta = %5f, default = %5f, corrected = %5f\n", it->superCluster()->eta(), it->energy(),ecor);
   

  }  




}


// ------------ method called once each job just before starting event loop  ------------
void 
EGEnergyAnalyzerSemiParm::beginJob()
{



}

// ------------ method called once each job just after ending the event loop  ------------
void 
EGEnergyAnalyzerSemiParm::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
EGEnergyAnalyzerSemiParm::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
EGEnergyAnalyzerSemiParm::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
EGEnergyAnalyzerSemiParm::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
EGEnergyAnalyzerSemiParm::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
EGEnergyAnalyzerSemiParm::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(EGEnergyAnalyzerSemiParm);
