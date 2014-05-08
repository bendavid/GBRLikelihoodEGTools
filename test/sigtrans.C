#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooExponential.h"
#include "RooGaussian.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "RooConstVar.h"
#include "RooDataSet.h"
#include "RooHybridBDTAutoPdf.h"
#include "RooFormulaVar.h"
#include "RooProdPdf.h"
#include "RooUniform.h"
#include "TRandom.h"
#include "TGraph.h"
#include "RooAddPdf.h"
#include "RooNDKeysPdf.h"
#include "RooExtendPdf.h"
#include "RooMinimizer.h"
#include "TFile.h"
#include "TNtuple.h"
#include "HybridGBRForest.h"
#include "RooProduct.h"
#include "RooChebychev.h"
#include "RooBernstein.h"
#include "RooPolynomial.h"
#include "RooGenericPdf.h"
//#include "HZZ2L2QRooPdfs.h"
#include "RooArgSet.h"
#include "RooArgList.h"
#include "RooCBShape.h"
//#include "RooCBShapeModified.h"
//#include "RooBernsteinFast.h"
#include "RooWorkspace.h"
#include "TH1D.h"
#include "TChain.h"
#include "TCut.h"
#include "TLine.h"
#include "TLegend.h"
#include "RooRandom.h"
#include "RooAddition.h"
#include "TSystem.h"
#include "RooCBFast.h"
#include "RooDoubleCBFast.h"
#include "RooGaussianFast.h"
//#include "RooBernsteinFast.h"
#include "TProfile.h"
#include "RooRevCBFast.h"
#include "RooBifurGauss.h"
#include "SigETransform.h"

SigETransform *g_sigtrans;
double sigetransformed(double sige,double sceta,double energy) {
 
  return g_sigtrans->sigEoverETranformed(sige,sceta,energy);
  
}

void sigtrans() {
 
TString dirname = "/scratch/bendavid/root/bare/sigtransplotsMay7/";
gSystem->mkdir(dirname,true);
gSystem->cd(dirname);  
  
 g_sigtrans = new SigETransform;
 g_sigtrans->Initialize("/home/bendavid/CMSSW_6_1_2/src/HiggsAnalysis/GBRLikelihoodEGTools/data/sigetrans_v1.root",1);
 
 TChain *tree = new TChain("RunLumiSelectionMod/MCProcessSelectionMod/HLTModP/GoodPVFilterMod/PhotonMvaMod/JetPub/JetCorrectionMod/SeparatePileUpMod/ElectronIDMod/MuonIDMod/PhotonPairSelectorPresel/PhotonTreeWriterPresel/hPhotonTree");
 tree->Add("/home/mingyang/cms/hist/hgg-2013Final8TeV_train/merged/hgg-2013Final8TeV_s12-diphoj-m60-v7n_noskim.root");
 
 //tree->Print();
 
 TCut selcut = "(ph1.pt > (mass/3.0) && ph2.pt > (mass/4.0) && mass>100. && mass<180. && ph1.idmva>-0.2 && ph2.idmva>-0.2)";
 TCut selweight = "mcweight";
 TCut prescale100 = "(evt%10==0)";
 TCut ebcentral = "abs(ph1.sceta)<1. && abs(ph2.sceta)<1.";
 
 tree->SetAlias("ph1.sigeoetrans","sigetransformed(ph1.eerr/ph1.e,ph1.sceta,ph1.e)");
 tree->SetAlias("ph2.sigeoetrans","sigetransformed(ph2.eerr/ph2.e,ph2.sceta,ph2.e)");
 tree->SetAlias("sigmom","0.5*sqrt( (ph1.eerr/ph1.e)^2 + (ph1.esmearing/ph1.e)^2 + (ph2.eerr/ph2.e)^2 + (ph2.esmearing/ph2.e)^2)");
 tree->SetAlias("sigmomtrans","0.5*sqrt( (ph1.sigeoetrans)^2 + (ph1.esmearing/ph1.e)^2 + (ph2.sigeoetrans)^2 + (ph2.esmearing/ph2.e)^2)");
 
 //tree->Draw("sigmom/(masserrsmeared/mass)",selcut*selweight*prescale100);
 
 TCanvas *c0 = new TCanvas;
 tree->Draw("sigmom:sigmomtrans>>htmp0(100,0.,0.03,100,0.,0.03)",selcut*selweight*prescale100,"COLZ");
 c0->SaveAs("sigmvtrans.eps");
 
 TCanvas *c1 = new TCanvas;
 tree->Draw("sigmom:mass>>htmp1(100,100.,180.,100,0.,0.03)",selcut*selweight*prescale100,"COLZ"); 
 c1->SaveAs("sigmvmass.eps");
 
 TCanvas *c2 = new TCanvas;
 tree->Draw("sigmomtrans:mass>>htmp2(100,100.,180.,100,0.,0.03)",selcut*selweight*prescale100,"COLZ"); 
 c2->SaveAs("sigmtransvmass.eps");
 
 TCanvas *c3 = new TCanvas;
 tree->Draw("sigmom:TMath::Max(abs(ph1.sceta),abs(ph2.sceta))>>htmp3(100,0.,2.5,100,0.,0.03)",selcut*selweight*prescale100,"COLZ");   
 c3->SaveAs("sigmvmaxeta.eps");
 
 TCanvas *c4 = new TCanvas;
 tree->Draw("sigmomtrans:TMath::Max(abs(ph1.sceta),abs(ph2.sceta))>>htmp4(100,0.,2.5,100,0.,0.03)",selcut*selweight*prescale100,"COLZ");  
 c4->SaveAs("sigmtransvmaxeta.eps");
 
 TCanvas *c5 = new TCanvas;
 tree->Draw("sigmom:mass>>htmp5(100,100.,180.,100,0.,0.03)",selcut*selweight*prescale100*ebcentral,"COLZ"); 
 c5->SaveAs("sigmvmassEBcentral.eps");
 
 TCanvas *c6 = new TCanvas;
 tree->Draw("sigmomtrans:mass>>htmp6(100,100.,180.,100,0.,0.03)",selcut*selweight*prescale100*ebcentral,"COLZ"); 
 c6->SaveAs("sigmtransvmassEBcentral.eps"); 
 
 //tree->Draw("sigmom>>(100,0.,0.07)",selcut*selweight*prescale100,"HIST");
 
 
 
 //new TCanvas;
 //tree->Draw("sigmomtrans",selcut*selweight*prescale100,"ESAME");
 
 return; 
  
}