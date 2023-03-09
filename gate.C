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

void gate()
{
   //----------------------BLOCK 1------------------------//
   //------------------ Data Reading ---------------------//

   //---------- Define string for data handling----------//
   TString path_to_file = "Dati/Gate/";

   TString fname = path_to_file + "1680_9mar23.dat";


   TString hname = "Run1long_tm_backcost_0830_mediumbin";
   TString info = "calibration at 1.68 #mus";
   //TString date = "/03/23";
   TString authors = "G. Cordova, A. Giani";
   TString acqtime = "14.5h FAULTY";
   //---------histogram name for fit plot----------//

   // auto df = ROOT::RDF::MakeCsvDataFrame(fname,false,'\t');

   //---------- Define Tree ---------//
   Int_t currentIgnoreLevel = gErrorIgnoreLevel;
   gErrorIgnoreLevel = kError;
   TTree *tree = new TTree("tree", "tree");
   tree->ReadFile(fname, "x/D:y");
   //tree->ReadFile(fname1);
   //tree->ReadFile(fname2);
   gErrorIgnoreLevel = currentIgnoreLevel;

   TTree *data_tree = new TTree("data tree", "tree of acquired data");

   Int_t N = tree->GetEntries();

   //---------- Tree Branches -------//
   Double_t time;
   Double_t channel;
   Double_t time_start;
   Double_t time_stop;
   Double_t effective_time;
   tree->SetBranchAddress("x", &channel);
   tree->SetBranchAddress("y", &time);
   data_tree->Branch("start", &time_start);
   data_tree->Branch("stop", &time_stop);
   data_tree->Branch("eff_time", &effective_time);
   auto min = 1.64;
   auto max = 1.74;
   auto bins = 5;
   TString ffit;
   ffit.Form("MuLife, " + info + " (tm+ background) %.2f-%.0f", min, max);

   //--------- Define Histrogram -----//
   TH1D *h = new TH1D("h", hname, bins, min, max);

   TH1D *h1 = new TH1D("h1", "t1 distribution", 100, 0, 700);

   TH1D *h2 = new TH1D("h", "t2 distribution", 100, 0, 700);

   //-- Fill Histogram with stop-start signal--//
   auto tmax = 85.8485;
   auto tmaxx = 687.144;
   for (Int_t i = 0; i < N; i++)
   {
      tree->GetEntry(i);
      if (channel == 1)
      {
         auto ttemp = time;
         h1->Fill(ttemp);
         // data_tree->Fill();
      }
      if (channel == 2)
      {
         auto t2 = time;

         auto channel2 = channel;
         tree->GetEntry(i - 1);
         auto t1 = time;
         auto channel1 = channel;
         /*if(channel2-channel1!=1){
         std::cout<< "Channel1: "<< channel1 << std::endl;
         std::cout<< "Channel2: "<< channel2 << std::endl;
         }*/
         // h1->Fill(t1);
         h2->Fill(t2);
         if (channel2 == 2 && channel1 == 1)
         {
            time_start = t1;
            time_stop = t2;
            auto decaytime = (t2 - t1)+8*pow(10,-9);
            effective_time = decaytime * pow(10, 6);
            if (decaytime < 0)
            {
               decaytime = decaytime + tmax;
               if (t1 > 620 && t2 < 70)
                  decaytime = decaytime - tmax + tmaxx;
            }
            if (decaytime * pow(10, 6) < max && decaytime * pow(10, 6) > min)
            {
               effective_time = decaytime * pow(10, 6);
               h->Fill(effective_time);
               data_tree->Fill();
            }
         }
      }
   }

   auto c = new TCanvas("c", "rawhist", 950, 800);
   gPad->SetLogy();
   h->SetBinErrorOption(TH1::EBinErrorOpt::kPoisson);
   h->GetYaxis()->SetTitle("Counts");
   h->GetXaxis()->SetTitle("Time [#mus]");
   h->SetTitle("Raw Counts " + info);
   h->Draw();

   return;
}