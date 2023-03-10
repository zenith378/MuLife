#include "TTree.h"
#include "TString.h"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include <stdio.h>
#include "math.h"
#include "RooRealVar.h"
#include "RooGaussModel.h"
#include "RooTruthModel.h"
#include "RooAddPdf.h"
#include "RooDecay.h"
#include "RooPlot.h"
#include "TH1.h"
#include "RooFitResult.h"
#include "RooDataHist.h"
#include "RooHist.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "RooConstVar.h"
#include "RooPolynomial.h"
#include "TAxis.h"
#include "RooArgSet.h"
#include "TRootCanvas.h"
#include "TF1.h"
#include "TPaveText.h"
#include <string.h>
#include <filesystem>
#include "RooNumConvPdf.h"
#include "TLatex.h"

using namespace RooFit;

void energy()
{
   //----------------------BLOCK 1------------------------//
   //------------------ Data Reading ---------------------//

   //---------- Define string for data handling----------//
   TString path_to_file = "../Dati/";

   TString fname = path_to_file + "spettro_elett_9mar23.dat";
   //TString fname1 = path_to_file + "cal_9mar23.dat";


   TString hname = "h";
   TString info = "Electron energy";
   TString date = "08/03/23";
   TString authors = "G. Cordova, A. Giani";
   TString acqtime = "15 min FAULTY";
   //---------histogram name for fit plot----------//

   // auto df = ROOT::RDF::MakeCsvDataFrame(fname,false,'\t');

   //---------- Define Tree ---------//
   Int_t currentIgnoreLevel = gErrorIgnoreLevel;
   gErrorIgnoreLevel = kError;
   TTree *tree = new TTree("tree", "tree");
   tree->ReadFile(fname, "x/D:y/D:z");
   //tree->ReadFile(fname1);
   gErrorIgnoreLevel = currentIgnoreLevel;

   Int_t N = tree->GetEntries();

   //---------- Tree Branches -------//
   Double_t time;
   Double_t ch0;
   Double_t ch1;

   tree->SetBranchAddress("x", &time);
   tree->SetBranchAddress("y", &ch0);
   tree->SetBranchAddress("z", &ch1);

   auto min = 0.001;
   auto max = 2.5;
   auto bins = 50;
   TString ffit;
   ffit.Form("MuLife, " + info + " %.2f-%.0f", min, max);

   //--------- Define Histrogram -----//
   TH1D *h = new TH1D("ch0",hname,bins,min,max);

   TH1D *h1 = new TH1D("ch1",hname,bins,min,max);

   //-- Fill Histogram with stop-start signal--//

   for (Int_t i = 0; i < N; i++)
   {
      tree->GetEntry(i);
      h->Fill(ch0);
      h1->Fill(ch1);
   }
   auto c = new TCanvas("c", "rawhist", 950, 800);
   //gPad->SetLogy();
   h->SetBinErrorOption(TH1::EBinErrorOpt::kPoisson);
   h->GetYaxis()->SetTitle("Counts");
   h->GetXaxis()->SetTitle("Voltage [V]");
   h->SetTitle("Raw Counts " + info+" CH0 (5 #mus delayed)");
   h->Draw();
   auto tp = new TPaveText(0.15, 0.7, 0.35, 0.85, "NDC");
   tp->AddText("MuLife");
   tp->AddText(authors);
   tp->AddText("09/03/23 15.5h");
   //tp->AddText("Run3 07/03/23 14.5h");
   tp->Draw();

   auto c1 = new TCanvas("c1", "rawhist", 950, 800);
   //gPad->SetLogy();
   h1->SetBinErrorOption(TH1::EBinErrorOpt::kPoisson);
   h1->GetYaxis()->SetTitle("Counts");
   h1->GetXaxis()->SetTitle("Voltage [V]");
   h1->SetTitle("Raw Counts " + info+" CH1");
   h1->Draw();
   auto tp1 = new TPaveText(0.15, 0.7, 0.35, 0.85, "NDC");
   tp1->AddText("MuLife");
   tp1->AddText(authors);
   tp1->AddText("09/03/23 15.5h");
   //tp->AddText("Run3 07/03/23 14.5h");
   tp1->Draw();

   return;
}