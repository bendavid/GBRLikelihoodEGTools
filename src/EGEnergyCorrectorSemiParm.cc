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
#include "TStreamerInfo.h"
#include "HiggsAnalysis/GBRLikelihood/interface/RooHybridBDTAutoPdf.h"

using namespace reco;

//--------------------------------------------------------------------------------------------------
EGEnergyCorrectorSemiParm::EGEnergyCorrectorSemiParm() :
_ws(0),
_isInitialized(kFALSE)
{
  // Constructor.
}


//--------------------------------------------------------------------------------------------------
EGEnergyCorrectorSemiParm::~EGEnergyCorrectorSemiParm()
{
  
  if (_ws) delete _ws;
 
}

//--------------------------------------------------------------------------------------------------
void EGEnergyCorrectorSemiParm::Initialize(std::string regweights) {
    _isInitialized = kTRUE;

    //retrieve RooWorkspace from file
    //TFile *fgbr = TFile::Open(regweights.c_str(),"READ");
    //TFile::SetReadStreamerInfo(true);
    
    //SetErrorHandler(DefaultErrorHandler);    
    
    TFile *fgbr = new TFile(regweights.c_str(),"READ");    
    _ws = static_cast<RooWorkspace*>(fgbr->Get("EGRegressionWorkspace"));
    fgbr->Close();

    //retrieve and cache required pointers from workspace
    
    RooGBRFunction *func_EB = static_cast<RooGBRFunction*>(_ws->arg("func_EB"));
    RooGBRFunction *func_EE = static_cast<RooGBRFunction*>(_ws->arg("func_EE"));
    
    _varseb = &func_EB->Vars();
    _varsee = &func_EE->Vars();
    
    _tgteb = _ws->var("tgtvar_EB");
    _tgtee = _ws->var("tgtvar_EE");

    _meaneb = _ws->var("sigmeanlim_EB");
    _sigmaeb = _ws->var("sigwidthlim_EB");
    _n1eb = _ws->var("signlim_EB");
    _n2eb = _ws->var("sign2lim_EB");
    
    _meanee = _ws->var("sigmeanlim_EE");
    _sigmaee = _ws->var("sigwidthlim_EE");
    _n1ee = _ws->var("signlim_EE");
    _n2ee = _ws->var("sign2lim_EE");    
    
    _pdfeb = _ws->pdf("sigpdf_EB");
    _pdfee = _ws->pdf("sigpdf_EE");
        

}


//--------------------------------------------------------------------------------------------------
void EGEnergyCorrectorSemiParm::CorrectedEnergyWithErrorV4(const reco::Photon &p, const reco::VertexCollection& vtxcol, double rho, EcalClusterLazyTools &clustertools, const edm::EventSetup &es, double &ecor, double &cbsigma, double &cbalpha1, double &cbn1, double &cbalpha2, double &cbn2, double &pdfpeakval) {
  
  const SuperClusterRef s = p.superCluster();
  const CaloClusterPtr b = s->seed(); //seed  basic cluster

  Bool_t isbarrel =  b->hitsAndFractions().at(0).first.subdetId()==EcalBarrel;
  
  const RooArgList *vars;
  if (isbarrel) vars = _varseb;
  else vars = _varsee;
    
//   //basic supercluster variables
  static_cast<RooRealVar&>((*vars)[0])  = s->rawEnergy();
  static_cast<RooRealVar&>((*vars)[1])  = s->eta();
  static_cast<RooRealVar&>((*vars)[2])  = s->phi();
  static_cast<RooRealVar&>((*vars)[3])  = p.r9();
  static_cast<RooRealVar&>((*vars)[4]) = s->etaWidth();
  static_cast<RooRealVar&>((*vars)[5]) = s->phiWidth();
  static_cast<RooRealVar&>((*vars)[6]) = double(s->clustersSize());
  static_cast<RooRealVar&>((*vars)[7]) = p.hadTowOverEm();
  static_cast<RooRealVar&>((*vars)[8]) = rho;
  static_cast<RooRealVar&>((*vars)[9]) = double(vtxcol.size());

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

  static_cast<RooRealVar&>((*vars)[10]) = b->eta()-s->eta();
  static_cast<RooRealVar&>((*vars)[11]) = reco::deltaPhi(b->phi(),s->phi());
  static_cast<RooRealVar&>((*vars)[12]) = b->energy()/s->rawEnergy();
  static_cast<RooRealVar&>((*vars)[13]) = clustertools.e3x3(*b)/be5x5;
  static_cast<RooRealVar&>((*vars)[14]) = sqrt(clustertools.localCovariances(*b)[0]); //sigietaieta
  static_cast<RooRealVar&>((*vars)[15]) = sqrt(clustertools.localCovariances(*b)[2]); //sigiphiiphi
  static_cast<RooRealVar&>((*vars)[16]) = clustertools.localCovariances(*b)[1];       //sigietaiphi
  static_cast<RooRealVar&>((*vars)[17]) = bemax/be5x5;                       //crystal energy ratio gap variables   
  static_cast<RooRealVar&>((*vars)[18]) = be2nd/be5x5;
  static_cast<RooRealVar&>((*vars)[19]) = betop/be5x5;
  static_cast<RooRealVar&>((*vars)[20]) = bebottom/be5x5;
  static_cast<RooRealVar&>((*vars)[21]) = beleft/be5x5;
  static_cast<RooRealVar&>((*vars)[22]) = beright/be5x5;
  static_cast<RooRealVar&>((*vars)[23]) = be2x5max/be5x5;                       //crystal energy ratio gap variables   
  static_cast<RooRealVar&>((*vars)[24]) = be2x5top/be5x5;
  static_cast<RooRealVar&>((*vars)[25]) = be2x5bottom/be5x5;
  static_cast<RooRealVar&>((*vars)[26]) = be2x5left/be5x5;
  static_cast<RooRealVar&>((*vars)[27]) = be2x5right/be5x5;

  if (isbarrel) {
    //additional energy ratio (always ~1 for endcap, therefore only included for barrel)
    static_cast<RooRealVar&>((*vars)[28]) = be5x5/b->energy();
    
    //local coordinates and crystal indices (barrel only)    
    
    //seed cluster
    float betacry, bphicry, bthetatilt, bphitilt;
    int bieta, biphi;
    _ecalLocal.localCoordsEB(*b,es,betacry,bphicry,bieta,biphi,bthetatilt,bphitilt);
    
    static_cast<RooRealVar&>((*vars)[29]) = bieta; //crystal ieta
    static_cast<RooRealVar&>((*vars)[30]) = biphi; //crystal iphi
    static_cast<RooRealVar&>((*vars)[31]) = bieta%5; //submodule boundary eta symmetry
    static_cast<RooRealVar&>((*vars)[32]) = biphi%2; //submodule boundary phi symmetry
    static_cast<RooRealVar&>((*vars)[33]) = (TMath::Abs(bieta)<=25)*(bieta%25) + (TMath::Abs(bieta)>25)*((bieta-25*TMath::Abs(bieta)/bieta)%20);  //module boundary eta approximate symmetry
    static_cast<RooRealVar&>((*vars)[34]) = biphi%20; //module boundary phi symmetry
    static_cast<RooRealVar&>((*vars)[35]) = betacry; //local coordinates with respect to closest crystal center at nominal shower depth
    static_cast<RooRealVar&>((*vars)[36]) = bphicry;

  }
  else {
    //preshower energy ratio (endcap only)
    static_cast<RooRealVar&>((*vars)[29])  = s->preshowerEnergy()/s->rawEnergy();
  }
  
  double den;
  if (isbarrel) {
    den = s->rawEnergy();
  }
  else {
    den = s->rawEnergy() + s->preshowerEnergy();
  }
  
  RooRealVar *tgtvar;
  RooAbsReal *mean;
  RooAbsReal *sigma;
  RooAbsReal *n1;
  RooAbsReal *n2;
  RooAbsPdf *pdf;
  
  if (isbarrel) {
    tgtvar = _tgteb;
    mean = _meaneb;
    sigma = _sigmaeb;
    n1 = _n1eb;
    n2 = _n2eb;
    pdf = _pdfeb;
  }
  else {
    tgtvar = _tgtee;
    mean = _meanee;
    sigma = _sigmaee;
    n1 = _n1ee;
    n2 = _n2ee;    
    pdf = _pdfee;
  }
  
  tgtvar->setVal(1.0); //evaluate pdf at peak position

  ecor = den/mean->getVal();
  cbsigma = sigma->getVal();
  cbalpha1 = 2.0;  //alpha hardcoded in this version of the regression
  cbn1 = n1->getVal();
  cbalpha2 = 1.0;  //alpha hardcoded in this version of the regression
  cbn2 = n2->getVal();
  
  pdfpeakval = pdf->getVal(*tgtvar);
  
  return;
  
}

//--------------------------------------------------------------------------------------------------
void EGEnergyCorrectorSemiParm::CorrectedEnergyWithErrorV4(const reco::GsfElectron &e, const reco::VertexCollection& vtxcol, double rho, EcalClusterLazyTools &clustertools, const edm::EventSetup &es, double &ecor, double &cbsigma, double &cbalpha1, double &cbn1, double &cbalpha2, double &cbn2, double &pdfpeakval) {
  
  const SuperClusterRef s = e.superCluster();
  const CaloClusterPtr b = s->seed(); //seed  basic cluster

  Bool_t isbarrel =  b->hitsAndFractions().at(0).first.subdetId()==EcalBarrel;
  
  const RooArgList *vars;
  if (isbarrel) vars = _varseb;
  else vars = _varsee;
    
//   //basic supercluster variables
  static_cast<RooRealVar&>((*vars)[0])  = s->rawEnergy();
  static_cast<RooRealVar&>((*vars)[1])  = s->eta();
  static_cast<RooRealVar&>((*vars)[2])  = s->phi();
  static_cast<RooRealVar&>((*vars)[3])  = e.r9();
  static_cast<RooRealVar&>((*vars)[4]) = s->etaWidth();
  static_cast<RooRealVar&>((*vars)[5]) = s->phiWidth();
  static_cast<RooRealVar&>((*vars)[6]) = double(s->clustersSize());
  static_cast<RooRealVar&>((*vars)[7]) = e.hcalOverEcalBc();
  static_cast<RooRealVar&>((*vars)[8]) = rho;
  static_cast<RooRealVar&>((*vars)[9]) = double(vtxcol.size());

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

  static_cast<RooRealVar&>((*vars)[10]) = b->eta()-s->eta();
  static_cast<RooRealVar&>((*vars)[11]) = reco::deltaPhi(b->phi(),s->phi());
  static_cast<RooRealVar&>((*vars)[12]) = b->energy()/s->rawEnergy();
  static_cast<RooRealVar&>((*vars)[13]) = clustertools.e3x3(*b)/be5x5;
  static_cast<RooRealVar&>((*vars)[14]) = sqrt(clustertools.localCovariances(*b)[0]); //sigietaieta
  static_cast<RooRealVar&>((*vars)[15]) = sqrt(clustertools.localCovariances(*b)[2]); //sigiphiiphi
  static_cast<RooRealVar&>((*vars)[16]) = clustertools.localCovariances(*b)[1];       //sigietaiphi
  static_cast<RooRealVar&>((*vars)[17]) = bemax/be5x5;                       //crystal energy ratio gap variables   
  static_cast<RooRealVar&>((*vars)[18]) = be2nd/be5x5;
  static_cast<RooRealVar&>((*vars)[19]) = betop/be5x5;
  static_cast<RooRealVar&>((*vars)[20]) = bebottom/be5x5;
  static_cast<RooRealVar&>((*vars)[21]) = beleft/be5x5;
  static_cast<RooRealVar&>((*vars)[22]) = beright/be5x5;
  static_cast<RooRealVar&>((*vars)[23]) = be2x5max/be5x5;                       //crystal energy ratio gap variables   
  static_cast<RooRealVar&>((*vars)[24]) = be2x5top/be5x5;
  static_cast<RooRealVar&>((*vars)[25]) = be2x5bottom/be5x5;
  static_cast<RooRealVar&>((*vars)[26]) = be2x5left/be5x5;
  static_cast<RooRealVar&>((*vars)[27]) = be2x5right/be5x5;

  if (isbarrel) {
    //additional energy ratio (always ~1 for endcap, therefore only included for barrel)
    static_cast<RooRealVar&>((*vars)[28]) = be5x5/b->energy();
    
    //local coordinates and crystal indices (barrel only)    
    
    //seed cluster
    float betacry, bphicry, bthetatilt, bphitilt;
    int bieta, biphi;
    _ecalLocal.localCoordsEB(*b,es,betacry,bphicry,bieta,biphi,bthetatilt,bphitilt);
    
    static_cast<RooRealVar&>((*vars)[29]) = bieta; //crystal ieta
    static_cast<RooRealVar&>((*vars)[30]) = biphi; //crystal iphi
    static_cast<RooRealVar&>((*vars)[31]) = bieta%5; //submodule boundary eta symmetry
    static_cast<RooRealVar&>((*vars)[32]) = biphi%2; //submodule boundary phi symmetry
    static_cast<RooRealVar&>((*vars)[33]) = (TMath::Abs(bieta)<=25)*(bieta%25) + (TMath::Abs(bieta)>25)*((bieta-25*TMath::Abs(bieta)/bieta)%20);  //module boundary eta approximate symmetry
    static_cast<RooRealVar&>((*vars)[34]) = biphi%20; //module boundary phi symmetry
    static_cast<RooRealVar&>((*vars)[35]) = betacry; //local coordinates with respect to closest crystal center at nominal shower depth
    static_cast<RooRealVar&>((*vars)[36]) = bphicry;

  }
  else {
    //preshower energy ratio (endcap only)
    static_cast<RooRealVar&>((*vars)[29])  = s->preshowerEnergy()/s->rawEnergy();
  }
  
  double den;
  if (isbarrel) {
    den = s->rawEnergy();
  }
  else {
    den = s->rawEnergy() + s->preshowerEnergy();
  }
  
  RooRealVar *tgtvar;
  RooAbsReal *mean;
  RooAbsReal *sigma;
  RooAbsReal *n1;
  RooAbsReal *n2;
  RooAbsPdf *pdf;
  
  if (isbarrel) {
    tgtvar = _tgteb;
    mean = _meaneb;
    sigma = _sigmaeb;
    n1 = _n1eb;
    n2 = _n2eb;
    pdf = _pdfeb;
  }
  else {
    tgtvar = _tgtee;
    mean = _meanee;
    sigma = _sigmaee;
    n1 = _n1ee;
    n2 = _n2ee;    
    pdf = _pdfee;
  }
  
  tgtvar->setVal(1.0); //evaluate pdf at peak position

  ecor = den/mean->getVal();
  cbsigma = sigma->getVal();
  cbalpha1 = 2.0;  //alpha hardcoded in this version of the regression
  cbn1 = n1->getVal();
  cbalpha2 = 1.0;  //alpha hardcoded in this version of the regression
  cbn2 = n2->getVal();
  
  pdfpeakval = pdf->getVal(*tgtvar);
  
  return;
  
}
