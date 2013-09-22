//--------------------------------------------------------------------------------------------------
// $Id $
//
// EGEnergyCorrectorTraditional
//
// Helper Class for applying regression-based energy corrections with optimized BDT implementation
//
// Authors: J.Bendavid
//--------------------------------------------------------------------------------------------------

#ifndef EGAMMATOOLS_EGEnergyCorrectorTraditional_H
#define EGAMMATOOLS_EGEnergyCorrectorTraditional_H
    
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "RecoEgamma/EgammaTools/interface/EcalClusterLocal.h"
#include "RooArgList.h"

class GBRForest;
class EcalClusterLazyTools;

class EGEnergyCorrectorTraditional {
  public:
    EGEnergyCorrectorTraditional();
    ~EGEnergyCorrectorTraditional(); 

    void Initialize(std::string regweights);
    Bool_t IsInitialized() const { return _isInitialized; }
        
    void CorrectedEnergyWithErrorV4Traditional(const reco::Photon &p, const reco::VertexCollection& vtxcol, double rho, EcalClusterLazyTools &clustertools, const edm::EventSetup &es, double &ecor, double &cbsigma);
    void CorrectedEnergyWithErrorV4Traditional(const reco::GsfElectron &e, const reco::VertexCollection& vtxcol, double rho, EcalClusterLazyTools &clustertools, const edm::EventSetup &es, double &ecor, double &cbsigma);    
    
    
  protected:
    std::vector<float> _vals;
    
    GBRForest *_foresteb;
    //GBRForest *_forestee;
    
    Bool_t _isInitialized;
    
    EcalClusterLocal _ecalLocal;
    
    };


#endif
