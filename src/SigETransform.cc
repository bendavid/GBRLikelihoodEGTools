#include <TFile.h>
#include "../interface/SigETransform.h"
#include "RooWorkspace.h"
#include "RooArgList.h"
#include "RooRealVar.h"
#include "RooAbsReal.h"
#include "RooAbsPdf.h"
#include "RooConstVar.h"
#include "TStreamerInfo.h"
#include "HiggsAnalysis/GBRLikelihood/interface/RooHybridBDTAutoPdf.h"
#include <cmath>

//--------------------------------------------------------------------------------------------------
SigETransform::SigETransform() :
_eerrvarEB(0),
_eerrvarEE(0),
_eerrvarinvEB(0),
_eerrvarinvEE(0),
_scetavarEB(0),
_scetavarEE(0),
_scetavarinvEB(0),
_scetavarinvEE(0),
_energyvarEB(0),
_energyvarEE(0),
_cdfEB(0),
_cdfEE(0),
_cdfinvEB(0),
_cdfinvEE(0),
_isInitialized(kFALSE)
{
  // Constructor.
}


//--------------------------------------------------------------------------------------------------
SigETransform::~SigETransform()
{
 
}

//--------------------------------------------------------------------------------------------------
void SigETransform::Initialize(std::string regweights, int version) {
    
    if (version==1) {
      TFile *fin = TFile::Open(regweights.c_str(),"READ");
      RooWorkspace *win = static_cast<RooWorkspace*>(fin->Get("wssigetrans"));
      
      RooCondAddPdf *pdfEB = static_cast<RooCondAddPdf*>(win->pdf("eerrpdf_EB"));
      RooCondAddPdf *pdfEE = static_cast<RooCondAddPdf*>(win->pdf("eerrpdf_EE"));
      RooCondAddPdf *pdfinvEB = static_cast<RooCondAddPdf*>(win->pdf("eerrpdf_invEB"));
      RooCondAddPdf *pdfinvEE = static_cast<RooCondAddPdf*>(win->pdf("eerrpdf_invEE"));
      
      _eerrvarEB = win->var("eerr_EB");
      _eerrvarEE = win->var("eerr_EE");
      _eerrvarinvEB = win->var("eerr_invEB");
      _eerrvarinvEE = win->var("eerr_invEE");
      
      _scetavarEB = win->var("sceta_EB");
      _scetavarEE = win->var("sceta_EE");
      _scetavarinvEB = win->var("sceta_invEB");
      _scetavarinvEE = win->var("sceta_invEE");      
      
      _energyvarEB = win->var("energy_EB");
      _energyvarEE = win->var("energy_EE");      
      
      _cdfEB = pdfEB->createCDF(*_eerrvarEB);
      _cdfEE = pdfEE->createCDF(*_eerrvarEE);
      _cdfinvEB = pdfinvEB->createCDF(*_eerrvarinvEB);
      _cdfinvEE = pdfinvEE->createCDF(*_eerrvarinvEE);
      
      _args.addOwned(*_cdfEB);
      _args.addOwned(*_cdfEE);
      _args.addOwned(*_cdfinvEB);
      _args.addOwned(*_cdfinvEE);
      
      _isInitialized = true;
    }    

}


//--------------------------------------------------------------------------------------------------
double SigETransform::sigEoverETranformed(double sigEoverE, double sceta, double energy) {
  
  bool isbarrel = std::abs(sceta)<1.479;
  
  if (isbarrel) {
    _eerrvarEB->setVal(-999.);
    _eerrvarinvEB->setVal(-999.);
    
    _scetavarEB->setVal(-999.);
    _scetavarinvEB->setVal(-999.);
    
    _energyvarEB->setVal(-999.);

    _eerrvarEB->setVal(sigEoverE);
    _eerrvarinvEB->setVal(sigEoverE);
    
    _scetavarEB->setVal(sceta);
    _scetavarinvEB->setVal(sceta);
    
    _energyvarEB->setVal(energy);
    
    double cdfval = _cdfEB->getVal();
    double sigEoverEtransformed = _cdfinvEB->findRoot(*_eerrvarinvEB,0.,1.,cdfval);
    return sigEoverEtransformed;
  }
  else {
    _eerrvarEE->setVal(-999.);
    _eerrvarinvEE->setVal(-999.);
    
    _scetavarEE->setVal(-999.);
    _scetavarinvEE->setVal(-999.);
    
    _energyvarEE->setVal(-999.);
    
    _eerrvarEE->setVal(sigEoverE);
    _eerrvarinvEE->setVal(sigEoverE);
    
    _scetavarEE->setVal(sceta);
    _scetavarinvEE->setVal(sceta);
    
    _energyvarEE->setVal(energy);
    
    double cdfval = _cdfEE->getVal();
    double sigEoverEtransformed = _cdfinvEE->findRoot(*_eerrvarinvEE,0.,1.,cdfval);    
    return sigEoverEtransformed;
  }
  
}
