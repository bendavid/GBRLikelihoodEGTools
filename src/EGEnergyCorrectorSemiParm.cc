#include <TFile.h>
#include "../interface/EGEnergyCorrectorSemiParm.h"
#include "CondFormats/EgammaObjects/interface/GBRForest.h"
#include "CondFormats/DataRecord/interface/GBRWrapperRcd.h"
#include "FWCore/Framework/interface/ESHandle.h" 
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "RooWorkspace.h"
#include "RooArgList.h"
#include "RooRealVar.h"
#include "RooAbsReal.h"
#include "RooAbsPdf.h"
#include "RooConstVar.h"
#include "TStreamerInfo.h"
#include "HiggsAnalysis/GBRLikelihood/interface/RooHybridBDTAutoPdf.h"
#include "HiggsAnalysis/GBRLikelihood/interface/RooDoubleCBFast.h"
#include "HiggsAnalysis/GBRLikelihood/interface/HybridGBRForest.h"


using namespace reco;

//--------------------------------------------------------------------------------------------------
EGEnergyCorrectorSemiParm::EGEnergyCorrectorSemiParm() :
_foresteb(0),
_forestee(0),
_mean(0),
_tgt(0),
_sigma(0),
_n1(0),
_n2(0),
_meanlim(0),
_sigmalim(0),
_n1lim(0),
_n2lim(0),
_pdf(0),
_isInitialized(kFALSE)
{
  // Constructor.
}


//--------------------------------------------------------------------------------------------------
EGEnergyCorrectorSemiParm::~EGEnergyCorrectorSemiParm()
{
  
  if (_isInitialized) {
    delete _foresteb;
    delete _forestee;
  }
 
}

//--------------------------------------------------------------------------------------------------
void EGEnergyCorrectorSemiParm::Initialize(std::string regweights) {
    _isInitialized = kTRUE;

    //initialize eval vector
    _vals.resize(37);
    
    //load forests from file
    TFile *fgbr = TFile::Open(regweights.c_str(),"READ");    
    fgbr->GetObject("EGRegressionForest_EB", _foresteb);
    fgbr->GetObject("EGRegressionForest_EE", _forestee);
    fgbr->Close();

    //recreate pdf with constraint transformations (can't load directly from file due to weird RooWorkspace IO features)
    
    _tgt = new RooRealVar("tgt","",1.);
    _mean = new RooRealVar("mean","",1.);
    _sigma = new RooRealVar("sigma","",1.);
    _n1 = new RooRealVar("n1","",2.);
    _n2 = new RooRealVar("n2","",2.);
    
    _sigmalim = new RooRealConstraint("sigmalim","",*_sigma,0.0002,0.5);
    _meanlim = new RooRealConstraint("meanlim","",*_mean,0.2,2.0);
    _n1lim = new RooRealConstraint("n1lim","",*_n1,1.01,110.);
    _n2lim = new RooRealConstraint("n2lim","",*_n2,1.01,110.);     
    
    _pdf = new RooDoubleCBFast("sigpdf","",*_tgt,RooFit::RooConst(1.),*_sigmalim,RooFit::RooConst(2.0),*_n1lim,RooFit::RooConst(1.0),*_n2lim);
    
    //add to RooArgList for proper garbage collection
    _args.addOwned(*_tgt);
    _args.addOwned(*_mean);
    _args.addOwned(*_sigma);
    _args.addOwned(*_n1);
    _args.addOwned(*_n2);
    _args.addOwned(*_sigmalim);
    _args.addOwned(*_meanlim);
    _args.addOwned(*_n1lim);
    _args.addOwned(*_n2lim);
    _args.addOwned(*_pdf);    

}


//--------------------------------------------------------------------------------------------------
void EGEnergyCorrectorSemiParm::CorrectedEnergyWithErrorV4(const reco::Photon &p, const reco::VertexCollection& vtxcol, double rho, EcalClusterLazyTools &clustertools, const edm::EventSetup &es, double &ecor, double &cbsigma, double &cbalpha1, double &cbn1, double &cbalpha2, double &cbn2, double &pdfpeakval) {
  
  const SuperClusterRef s = p.superCluster();
  const CaloClusterPtr b = s->seed(); //seed  basic cluster

  Bool_t isbarrel =  b->hitsAndFractions().at(0).first.subdetId()==EcalBarrel;
    
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
    _vals[29]  = s->preshowerEnergy()/s->rawEnergy();
  }
  
  double den;
  HybridGBRForest *forest;  
  if (isbarrel) {
    den = s->rawEnergy();
    forest = _foresteb;
  }
  else {
    den = s->rawEnergy() + s->preshowerEnergy();
    forest = _forestee;
  }
  
  _tgt->setVal(1.0); //evaluate pdf at peak position
  
  //set raw response variables from GBRForest
  _sigma->setVal(forest->GetResponse(&_vals[0],0));
  _mean->setVal(forest->GetResponse(&_vals[0],1));
  _n1->setVal(forest->GetResponse(&_vals[0],2));
  _n2->setVal(forest->GetResponse(&_vals[0],3));
  
  //retrieve final pdf parameter values from transformed forest outputs
  ecor = den/_meanlim->getVal();
  cbsigma = _sigmalim->getVal();
  cbalpha1 = 2.0;  //alpha hardcoded in this version of the regression
  cbn1 = _n1lim->getVal();
  cbalpha2 = 1.0;  //alpha hardcoded in this version of the regression
  cbn2 = _n2lim->getVal();
  
  pdfpeakval = _pdf->getVal(*_tgt);
  
  return;
  
}

//--------------------------------------------------------------------------------------------------
void EGEnergyCorrectorSemiParm::CorrectedEnergyWithErrorV4(const reco::GsfElectron &e, const reco::VertexCollection& vtxcol, double rho, EcalClusterLazyTools &clustertools, const edm::EventSetup &es, double &ecor, double &cbsigma, double &cbalpha1, double &cbn1, double &cbalpha2, double &cbn2, double &pdfpeakval) {
  
  const SuperClusterRef s = e.superCluster();
  const CaloClusterPtr b = s->seed(); //seed  basic cluster

  Bool_t isbarrel =  b->hitsAndFractions().at(0).first.subdetId()==EcalBarrel;
    
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
    _vals[29]  = s->preshowerEnergy()/s->rawEnergy();
  }
  
  double den;
  HybridGBRForest *forest;  
  if (isbarrel) {
    den = s->rawEnergy();
    forest = _foresteb;
  }
  else {
    den = s->rawEnergy() + s->preshowerEnergy();
    forest = _forestee;
  }
  
  _tgt->setVal(1.0); //evaluate pdf at peak position
  
  //set raw response variables from GBRForest
  _sigma->setVal(forest->GetResponse(&_vals[0],0));
  _mean->setVal(forest->GetResponse(&_vals[0],1));
  _n1->setVal(forest->GetResponse(&_vals[0],2));
  _n2->setVal(forest->GetResponse(&_vals[0],3));
  
  //retrieve final pdf parameter values from transformed forest outputs
  ecor = den/_meanlim->getVal();
  cbsigma = _sigmalim->getVal();
  cbalpha1 = 2.0;  //alpha hardcoded in this version of the regression
  cbn1 = _n1lim->getVal();
  cbalpha2 = 1.0;  //alpha hardcoded in this version of the regression
  cbn2 = _n2lim->getVal();
  
  pdfpeakval = _pdf->getVal(*_tgt);
  
  return;
  
}
