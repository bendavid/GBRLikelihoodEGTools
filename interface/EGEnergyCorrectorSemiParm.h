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

class RooWorkspace;
class RooArgList;
class RooRealVar;
class RooAbsPdf;
class RooAbsReal;
class EcalClusterLazyTools;

class EGEnergyCorrectorSemiParm {
  public:
    EGEnergyCorrectorSemiParm();
    ~EGEnergyCorrectorSemiParm(); 

    void Initialize(std::string regweights);
    Bool_t IsInitialized() const { return _isInitialized; }
        
    void CorrectedEnergyWithErrorV4(const reco::Photon &p, const reco::VertexCollection& vtxcol, double rho, EcalClusterLazyTools &clustertools, const edm::EventSetup &es, double &ecor, double &cbsigma, double &cbalpha1, double &cbn1, double &cbalpha2, double &cbn2, double &pdfpeakval);
    void CorrectedEnergyWithErrorV4(const reco::GsfElectron &e, const reco::VertexCollection& vtxcol, double rho, EcalClusterLazyTools &clustertools, const edm::EventSetup &es, double &ecor, double &cbsigma, double &cbalpha1, double &cbn1, double &cbalpha2, double &cbn2, double &pdfpeakval);    
    
  protected:
    RooWorkspace *_ws; 
    const RooArgList *_varseb;   
    const RooArgList *_varsee;    
    RooRealVar *_tgteb;
    RooRealVar *_tgtee;    
    RooAbsPdf *_pdfeb;
    RooAbsPdf *_pdfee;
    
    RooAbsReal *_meaneb;
    RooAbsReal *_sigmaeb;
    RooAbsReal *_n1eb;
    RooAbsReal *_n2eb;
    
    RooAbsReal *_meanee;
    RooAbsReal *_sigmaee;
    RooAbsReal *_n1ee;
    RooAbsReal *_n2ee;    
    
    Bool_t _isInitialized;
    
    EcalClusterLocal _ecalLocal;
    
    };


#endif