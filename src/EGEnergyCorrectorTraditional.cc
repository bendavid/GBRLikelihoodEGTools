#include <TFile.h>
#include "../interface/EGEnergyCorrectorTraditional.h"
#include "CondFormats/EgammaObjects/interface/GBRForest.h"
#include "CondFormats/DataRecord/interface/GBRWrapperRcd.h"
#include "FWCore/Framework/interface/ESHandle.h" 
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

using namespace reco;

//--------------------------------------------------------------------------------------------------
EGEnergyCorrectorTraditional::EGEnergyCorrectorTraditional() :
_foresteb(0),
//_forestee(0),
_isInitialized(kFALSE)
{
  // Constructor.
}


//--------------------------------------------------------------------------------------------------
EGEnergyCorrectorTraditional::~EGEnergyCorrectorTraditional()
{
  
  if (_isInitialized) {
    delete _foresteb;
    //delete _forestee;
  }
 
}

//--------------------------------------------------------------------------------------------------
void EGEnergyCorrectorTraditional::Initialize(std::string regweights) {
    _isInitialized = kTRUE;

    //initialize eval vector
    _vals.resize(37);
    
    //load forests from file
    TFile *fgbr = TFile::Open(regweights.c_str(),"READ");    
    fgbr->GetObject("EBCorrection", _foresteb);
    //fgbr->GetObject("EGRegressionForest_EE", _forestee);
    fgbr->Close();

}


//--------------------------------------------------------------------------------------------------
void EGEnergyCorrectorTraditional::CorrectedEnergyWithErrorV4Traditional(const reco::Photon &p, const reco::VertexCollection& vtxcol, double rho, EcalClusterLazyTools &clustertools, const edm::EventSetup &es, double &ecor, double &cbsigma) {
  
  //not implemented at the moment
  cbsigma = 0.;
  
  const SuperClusterRef s = p.superCluster();
  const CaloClusterPtr b = s->seed(); //seed  basic cluster

  Bool_t isbarrel =  b->hitsAndFractions().at(0).first.subdetId()==EcalBarrel;
  
  if (!isbarrel) {
    ecor = s->rawEnergy() + s->preshowerEnergy();
    return;
  }
    
//   //basic supercluster variables
  _vals[0]  = s->rawEnergy();
  _vals[1]  = s->eta();
  _vals[2]  = s->phi();
  _vals[3]  = p.r9();
  _vals[4] = s->etaWidth();
  _vals[5] = s->phiWidth();
  _vals[6] = double(s->clustersSize());
  _vals[7] = p.hadTowOverEm();
  _vals[8] = rho;
  _vals[9] = double(vtxcol.size());

  //seed basic cluster variables
  double bemax = clustertools.eMax(*b);
  double be2nd = clustertools.e2nd(*b);
  double betop = clustertools.eTop(*b);
  double bebottom = clustertools.eBottom(*b);
  double beleft = clustertools.eLeft(*b);
  double beright = clustertools.eRight(*b);

  double be2x5max = clustertools.e2x5Max(*b);
  double be2x5top = clustertools.e2x5Top(*b);
  double be2x5bottom = clustertools.e2x5Bottom(*b);
  double be2x5left = clustertools.e2x5Left(*b);
  double be2x5right = clustertools.e2x5Right(*b);
  
  double be5x5 = clustertools.e5x5(*b);

  _vals[10] = b->eta()-s->eta();
  _vals[11] = reco::deltaPhi(b->phi(),s->phi());
  _vals[12] = b->energy()/s->rawEnergy();
  _vals[13] = clustertools.e3x3(*b)/be5x5;
  _vals[14] = sqrt(clustertools.localCovariances(*b)[0]); //sigietaieta
  _vals[15] = sqrt(clustertools.localCovariances(*b)[2]); //sigiphiiphi
  _vals[16] = clustertools.localCovariances(*b)[1];       //sigietaiphi
  _vals[17] = bemax/be5x5;                       //crystal energy ratio gap variables   
  _vals[18] = be2nd/be5x5;
  _vals[19] = betop/be5x5;
  _vals[20] = bebottom/be5x5;
  _vals[21] = beleft/be5x5;
  _vals[22] = beright/be5x5;
  _vals[23] = be2x5max/be5x5;                       //crystal energy ratio gap variables   
  _vals[24] = be2x5top/be5x5;
  _vals[25] = be2x5bottom/be5x5;
  _vals[26] = be2x5left/be5x5;
  _vals[27] = be2x5right/be5x5;

  if (isbarrel) {
    //additional energy ratio (always ~1 for endcap, therefore only included for barrel)
    _vals[28] = be5x5/b->energy();
    
    //local coordinates and crystal indices (barrel only)    
    
    //seed cluster
    float betacry, bphicry, bthetatilt, bphitilt;
    int bieta, biphi;
    _ecalLocal.localCoordsEB(*b,es,betacry,bphicry,bieta,biphi,bthetatilt,bphitilt);
    
    _vals[29] = bieta; //crystal ieta
    _vals[30] = biphi; //crystal iphi
    _vals[31] = bieta%5; //submodule boundary eta symmetry
    _vals[32] = biphi%2; //submodule boundary phi symmetry
    _vals[33] = (TMath::Abs(bieta)<=25)*(bieta%25) + (TMath::Abs(bieta)>25)*((bieta-25*TMath::Abs(bieta)/bieta)%20);  //module boundary eta approximate symmetry
    _vals[34] = biphi%20; //module boundary phi symmetry
    _vals[35] = betacry; //local coordinates with respect to closest crystal center at nominal shower depth
    _vals[36] = bphicry;

  }
  else {
    //preshower energy ratio (endcap only)
    _vals[28]  = s->preshowerEnergy()/s->rawEnergy();
  }
  
  double den = s->rawEnergy();
  GBRForest *forest = _foresteb;  
//   if (isbarrel) {
//     den = s->rawEnergy();
//     forest = _foresteb;
//   }
//   else {
//     den = s->rawEnergy() + s->preshowerEnergy();
//     forest = _forestee;
//   }
    
  ecor = den*forest->GetResponse(&_vals[0]);

  
  return;
  
}

//--------------------------------------------------------------------------------------------------
void EGEnergyCorrectorTraditional::CorrectedEnergyWithErrorV4Traditional(const reco::GsfElectron &e, const reco::VertexCollection& vtxcol, double rho, EcalClusterLazyTools &clustertools, const edm::EventSetup &es, double &ecor, double &cbsigma) {
  
  //not implemented at the moment
  cbsigma = 0.;
  
  const SuperClusterRef s = e.superCluster();
  const CaloClusterPtr b = s->seed(); //seed  basic cluster

  Bool_t isbarrel =  b->hitsAndFractions().at(0).first.subdetId()==EcalBarrel;
  
  if (!isbarrel) {
    ecor = s->rawEnergy() + s->preshowerEnergy();
    return;
  }
    
//   //basic supercluster variables
  _vals[0]  = s->rawEnergy();
  _vals[1]  = s->eta();
  _vals[2]  = s->phi();
  _vals[3]  = e.r9();
  _vals[4] = s->etaWidth();
  _vals[5] = s->phiWidth();
  _vals[6] = double(s->clustersSize());
  _vals[7] = e.hcalOverEcalBc();
  _vals[8] = rho;
  _vals[9] = double(vtxcol.size());

  //seed basic cluster variables
  double bemax = clustertools.eMax(*b);
  double be2nd = clustertools.e2nd(*b);
  double betop = clustertools.eTop(*b);
  double bebottom = clustertools.eBottom(*b);
  double beleft = clustertools.eLeft(*b);
  double beright = clustertools.eRight(*b);

  double be2x5max = clustertools.e2x5Max(*b);
  double be2x5top = clustertools.e2x5Top(*b);
  double be2x5bottom = clustertools.e2x5Bottom(*b);
  double be2x5left = clustertools.e2x5Left(*b);
  double be2x5right = clustertools.e2x5Right(*b);
  
  double be5x5 = clustertools.e5x5(*b);

  _vals[10] = b->eta()-s->eta();
  _vals[11] = reco::deltaPhi(b->phi(),s->phi());
  _vals[12] = b->energy()/s->rawEnergy();
  _vals[13] = clustertools.e3x3(*b)/be5x5;
  _vals[14] = sqrt(clustertools.localCovariances(*b)[0]); //sigietaieta
  _vals[15] = sqrt(clustertools.localCovariances(*b)[2]); //sigiphiiphi
  _vals[16] = clustertools.localCovariances(*b)[1];       //sigietaiphi
  _vals[17] = bemax/be5x5;                       //crystal energy ratio gap variables   
  _vals[18] = be2nd/be5x5;
  _vals[19] = betop/be5x5;
  _vals[20] = bebottom/be5x5;
  _vals[21] = beleft/be5x5;
  _vals[22] = beright/be5x5;
  _vals[23] = be2x5max/be5x5;                       //crystal energy ratio gap variables   
  _vals[24] = be2x5top/be5x5;
  _vals[25] = be2x5bottom/be5x5;
  _vals[26] = be2x5left/be5x5;
  _vals[27] = be2x5right/be5x5;

  if (isbarrel) {
    //additional energy ratio (always ~1 for endcap, therefore only included for barrel)
    _vals[28] = be5x5/b->energy();
    
    //local coordinates and crystal indices (barrel only)    
    
    //seed cluster
    float betacry, bphicry, bthetatilt, bphitilt;
    int bieta, biphi;
    _ecalLocal.localCoordsEB(*b,es,betacry,bphicry,bieta,biphi,bthetatilt,bphitilt);
    
    _vals[29] = bieta; //crystal ieta
    _vals[30] = biphi; //crystal iphi
    _vals[31] = bieta%5; //submodule boundary eta symmetry
    _vals[32] = biphi%2; //submodule boundary phi symmetry
    _vals[33] = (TMath::Abs(bieta)<=25)*(bieta%25) + (TMath::Abs(bieta)>25)*((bieta-25*TMath::Abs(bieta)/bieta)%20);  //module boundary eta approximate symmetry
    _vals[34] = biphi%20; //module boundary phi symmetry
    _vals[35] = betacry; //local coordinates with respect to closest crystal center at nominal shower depth
    _vals[36] = bphicry;

  }
  else {
    //preshower energy ratio (endcap only)
    _vals[28]  = s->preshowerEnergy()/s->rawEnergy();
  }
  
  double den = s->rawEnergy();
  GBRForest *forest = _foresteb;  
//   if (isbarrel) {
//     den = s->rawEnergy();
//     forest = _foresteb;
//   }
//   else {
//     den = s->rawEnergy() + s->preshowerEnergy();
//     forest = _forestee;
//   }
    
  ecor = den*forest->GetResponse(&_vals[0]);

  
  return;
  
}


