//--------------------------------------------------------------------------------------------------
// $Id $
//
// EGEnergyCorrectorSemiParm
//
// Helper Class for applying regression-based energy corrections with optimized BDT implementation
//
// Authors: J.Bendavid
//--------------------------------------------------------------------------------------------------

#ifndef EGAMMATOOLS_EGEnergyCorrectorSemiParm_H
#define EGAMMATOOLS_EGEnergyCorrectorSemiParm_H
    
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "RecoEgamma/EgammaTools/interface/EcalClusterLocal.h"
#include "RooArgList.h"

class RooWorkspace;
class RooRealVar;
class RooAbsPdf;
class RooAbsReal;
class HybridGBRForest;
class HybridGBRForestD;
class EcalClusterLazyTools;

class EGEnergyCorrectorSemiParm {
  public:
    EGEnergyCorrectorSemiParm();
    ~EGEnergyCorrectorSemiParm(); 

    void Initialize(std::string regweights, int version);
    Bool_t IsInitialized() const { return _isInitialized; }
        
    void CorrectedEnergyWithErrorV4(const reco::Photon &p, const reco::VertexCollection& vtxcol, double rho, EcalClusterLazyTools &clustertools, const edm::EventSetup &es, double &ecor, double &cbsigma, double &cbalpha1, double &cbn1, double &cbalpha2, double &cbn2, double &pdfpeakval);
    void CorrectedEnergyWithErrorV4(const reco::GsfElectron &e, const reco::VertexCollection& vtxcol, double rho, EcalClusterLazyTools &clustertools, const edm::EventSetup &es, double &ecor, double &cbsigma, double &cbalpha1, double &cbn1, double &cbalpha2, double &cbn2, double &pdfpeakval);    

    //V5 version folding detector along phi symmetries to avoid bias from MC intercalibration
    void CorrectedEnergyWithErrorV5(const reco::Photon &p, const reco::VertexCollection& vtxcol, double rho, EcalClusterLazyTools &clustertools, const edm::EventSetup &es, double &ecor, double &cbsigma, double &cbalpha1, double &cbn1, double &cbalpha2, double &cbn2, double &pdfpeakval);
    void CorrectedEnergyWithErrorV5(const reco::GsfElectron &e, const reco::VertexCollection& vtxcol, double rho, EcalClusterLazyTools &clustertools, const edm::EventSetup &es, double &ecor, double &cbsigma, double &cbalpha1, double &cbn1, double &cbalpha2, double &cbn2, double &pdfpeakval); 
    
    //V6: Updated training
    void CorrectedEnergyWithErrorV6(const reco::Photon &p, const reco::VertexCollection& vtxcol, double rho, EcalClusterLazyTools &clustertools, const edm::EventSetup &es, double &ecor, double &sigEoverE, double &cbmean, double &cbsigma, double &cbalpha1, double &cbn1, double &cbalpha2, double &cbn2, double &pdfpeakval);    
    void CorrectedEnergyWithErrorV6(const reco::GsfElectron &e, const reco::VertexCollection& vtxcol, double rho, EcalClusterLazyTools &clustertools, const edm::EventSetup &es, double &ecor, double &sigEoverE, double &cbmean, double &cbsigma, double &cbalpha1, double &cbn1, double &cbalpha2, double &cbn2, double &pdfpeakval);    
    
    //V7: Updated training, phi-blind (folding detector along phi symmetries to avoid bias from MC intercalibration)
    void CorrectedEnergyWithErrorV7(const reco::Photon &p, const reco::VertexCollection& vtxcol, double rho, EcalClusterLazyTools &clustertools, const edm::EventSetup &es, double &ecor, double &sigEoverE, double &cbmean, double &cbsigma, double &cbalpha1, double &cbn1, double &cbalpha2, double &cbn2, double &pdfpeakval);        
    void CorrectedEnergyWithErrorV7(const reco::GsfElectron &e, const reco::VertexCollection& vtxcol, double rho, EcalClusterLazyTools &clustertools, const edm::EventSetup &es, double &ecor, double &sigEoverE, double &cbmean, double &cbsigma, double &cbalpha1, double &cbn1, double &cbalpha2, double &cbn2, double &pdfpeakval);      
    
    //V8: Wrapper function to call V6 for barrel and V7 for endcap
    void CorrectedEnergyWithErrorV8(const reco::Photon &p, const reco::VertexCollection& vtxcol, double rho, EcalClusterLazyTools &clustertools, const edm::EventSetup &es, double &ecor, double &sigEoverE, double &cbmean, double &cbsigma, double &cbalpha1, double &cbn1, double &cbalpha2, double &cbn2, double &pdfpeakval);     
    void CorrectedEnergyWithErrorV8(const reco::GsfElectron &e, const reco::VertexCollection& vtxcol, double rho, EcalClusterLazyTools &clustertools, const edm::EventSetup &es, double &ecor, double &sigEoverE, double &cbmean, double &cbsigma, double &cbalpha1, double &cbn1, double &cbalpha2, double &cbn2, double &pdfpeakval);      
    
    
  protected:
    std::vector<float> _vals;
    
    HybridGBRForest *_foresteb;
    HybridGBRForest *_forestee;

    HybridGBRForestD *_forestDeb;
    HybridGBRForestD *_forestDee;    
    
    RooRealVar *_mean;
    RooRealVar *_tgt;
    RooRealVar *_sigma;
    RooRealVar *_n1;
    RooRealVar *_n2;  
    
    RooAbsReal *_meanlim;
    RooAbsReal *_sigmalim;
    RooAbsReal *_n1lim;
    RooAbsReal *_n2lim;        
    
    RooAbsPdf *_pdf;
    
    RooArgList _args;
    
    Bool_t _isInitialized;
    
    EcalClusterLocal _ecalLocal;
    
    };


#endif
