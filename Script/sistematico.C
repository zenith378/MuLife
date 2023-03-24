#include "TTree.h"
#include "TH1.h"
#include <iostream>
#include "TCanvas.h"
#include "TGraph.h"
#include "TPaveText.h"
#include "TArrow.h"
#include "RooGaussian.h"
#include "RooAddPdf.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooHist.h"
#include "RooArgSet.h"
#include "RooDataSet.h"

using namespace RooFit;


void sistematico()
{
   //----------------------BLOCK 1------------------------//
   //------------------ Data Reading ---------------------//

   //---------- Define string for data handling----------//
   TString path_to_file="../Dati/Sistematico/";

   TString fname = path_to_file+"sist_1mar23.dat";

   TString hname = "Risoluzione";
   TString date = "02/03/23";
   TString authors = "G. Cordova, A. Giani";
   TString acqtime = "10 min";
   //---------histogram name for fit plot----------//
   TString ffit = "Sistematico, "+date+" 50 ms";

   // auto df = ROOT::RDF::MakeCsvDataFrame(fname,false,'\t');

   //---------- Define Tree ---------//
   TTree *tree = new TTree("tree", "tree");
   tree->ReadFile(fname, "x/D:y");

   Int_t N = tree->GetEntries();

   //---------- Tree Branches -------//
   Double_t time;
   Double_t channel;
   Double_t effective_time;
   tree->SetBranchAddress("x", &channel);
   tree->SetBranchAddress("y", &time);
   auto min=160;
   auto max=260;
   auto bins=5;
   //--------- Define Histrogram -----//
   TH1D *h = new TH1D("h", hname, bins, min, max);






   //-- Fill Histogram with stop-start signal--//
   auto tmax=85.8485;
   auto tmaxx=687.144;
   for (Int_t i = 0; i < N; i++)
   {
      tree->GetEntry(i);
      auto channel2 = channel;
      if (channel2 == 2)
      {
         auto t2 = time;
         tree->GetEntry(i - 1);
         auto t1 = time;
         auto channel1 = channel;
         if (channel1 == 1)
         {
            auto decaytime = (t2 - t1);

            if (decaytime < -tmaxx) decaytime = decaytime + tmaxx;
            if (-tmaxx<decaytime < -tmax) decaytime = decaytime + tmax;
            effective_time = decaytime*pow(10,9);
            //if (effective_time < max && effective_time > min)
            //{

               h->Fill(effective_time);
            //}
         }
      }
   }
      
   std::cout <<"Overflow: " << h->GetBinContent(N+1) << std::endl;
   std::cout <<"Underflow: " << h->GetBinContent(0) << std::endl;

   auto c = new TCanvas("c", "raw", 950, 800);
   //gPad->SetLogy();
   h->GetYaxis()->SetTitle("Counts");
   h->GetXaxis()->SetTitle("Time [#mus]");
   //h->GetXaxis()->SetNdivisions(-502);

   h->SetTitle("Preliminary histogram for an expected signal of 180 ns");
   h->Draw();
   auto tp = new TPaveText(183,5000,194,6300);
   tp->AddText("MuLife");
   tp->AddText(authors);
   tp->AddText("02/03/23 4 min");
   tp->Draw();

   
   return;
}