//--------------------------------------------------------------------------------------------------
// $Id $
//
// SigETransform
//
// Helper Class for applying transformation to SigmaE/E to decorrelated dependence on photon energy
//
// Authors: J.Bendavid
//--------------------------------------------------------------------------------------------------

#ifndef EGAMMATOOLS_SigETransform_H
#define EGAMMATOOLS_SigETransform_H
    
#include "RooArgList.h"

class RooWorkspace;
class RooRealVar;
class RooAbsPdf;
class RooAbsReal;
class HybridGBRForest;
class HybridGBRForestD;
class EcalClusterLazyTools;

class SigETransform {
  public:
    SigETransform();
    ~SigETransform(); 

    void Initialize(std::string regweights, int version);
    Bool_t IsInitialized() const { return _isInitialized; }
        
    double sigEoverETranformed(double sigEoverE, double sceta, double energy);
        
    
  protected:
      
    RooRealVar *_eerrvarEB;
    RooRealVar *_eerrvarEE;
    RooRealVar *_eerrvarinvEB;
    RooRealVar *_eerrvarinvEE;

    RooRealVar *_scetavarEB;
    RooRealVar *_scetavarEE;
    RooRealVar *_scetavarinvEB;
    RooRealVar *_scetavarinvEE;

    RooRealVar *_energyvarEB;
    RooRealVar *_energyvarEE;
    
    RooAbsReal *_cdfEB;
    RooAbsReal *_cdfEE;
    RooAbsReal *_cdfinvEB;
    RooAbsReal *_cdfinvEE;
    
    RooArgList _args;
    
    Bool_t _isInitialized;
    
    
    };


#endif
