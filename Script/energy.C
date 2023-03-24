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
#include "TGraph.h"
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
   TString path_to_file = "../Dati/Energy/";

   TString run0 = path_to_file + "spettro_elett_9mar23.dat";
   TString run1 = path_to_file + "spettro_elet_14mar23.dat";
   TString run11 = path_to_file + "spettro2_14mar23.dat";
   TString run2 = path_to_file + "spettro_elettrone_15mar.dat";
   TString run3 = path_to_file + "spettro_16mar23.dat";

   TString hname = "h";
   TString info = "Electron energy";
   TString date = "08/03/23";
   TString authors = "G. Cordova, A. Giani";
   TString acqtime = "15 min FAULTY";
   //---------histogram name for fit plot----------//

   // auto df = ROOT::RDF::MakeCsvDataFrame(fname,false,'\t');

   //---------- Define Tree ---------//
   // Int_t currentIgnoreLevel = gErrorIgnoreLevel;
   // gErrorIgnoreLevel = kError;
   TTree *tree = new TTree("tree", "tree");
   tree->ReadFile(run2, "x/D:y/D:z");
   // tree->ReadFile(run11);
   //tree->ReadFile(run2);
   // gErrorIgnoreLevel = currentIgnoreLevel;

   Int_t N = tree->GetEntries();
   //---------- Tree Branches -------//
   Double_t time;
   Double_t ch0;
   Double_t ch1;

   tree->SetBranchAddress("x", &time);
   tree->SetBranchAddress("y", &ch1);
   tree->SetBranchAddress("z", &ch0);


   auto min = -2.5;
   auto max = 2.5;
   auto bins = 100;
   TString ffit;
   ffit.Form("MuLife, " + info + " %.2f-%.0f", min, max);

   //--------- Define Histrogram -----//
   TH1D *h = new TH1D("h", "ch0-ch1_energy", 500, -0.2, 0.2);

   TH1D *h1 = new TH1D("h1", "raw ch1", bins, -1, 2.5);

   TH1D *h2 = new TH1D("h2", "deltav/v", 100, -2, 20);

   Double_t temp1 = 0;
   Double_t temp2 = 0;
   Double_t arr_v1[N];
   Double_t arr_v2[N];
   Double_t arr_deltav[N];

   for (Int_t i = 0; i < N; i++)
   {
      tree->GetEntry(i);

      temp2 = ch0;
      temp1 = ch1;
      if (ch0 > 0.1 && ch1>0.0001 && ch1 <0.7)
      {
         temp2=ch0+0.078*ch1;
         arr_deltav[i] = temp2 - temp1;
         arr_v1[i] = temp1;
         arr_v2[i] = temp2;
         Double_t rap = arr_deltav[i] / arr_v1[i];
         h->Fill(temp2 - temp1);
         h1->Fill(temp1);
         h2->Fill(rap);
      }
   }

    std::cout << h2->GetBinContent(101) << std::endl;

   auto c = new TCanvas("c", "Electron Energy Spectrum", 950, 800);
   //gPad->SetLogy();
   h->SetBinErrorOption(TH1::EBinErrorOpt::kPoisson);
   h->GetYaxis()->SetTitle("Counts");
   h->GetXaxis()->SetTitle("Voltage [V]");
   h->SetTitle("Ch0-Ch1 " + info + " CH0>0.1 && 0.0001<CH1<0.7");
   h->Draw();
   auto tp = new TPaveText(0.15, 0.7, 0.35, 0.85, "NDC");
   tp->AddText("MuLife");
   tp->AddText(authors);
   // tp->AddText("Run0 9/03/23 17h");
   //tp->AddText("Run1 14/03/23 14h");
   //tp->AddText("Run2 15/03/23 23h");
   tp->AddText("Run3 16/03/23 86h");
   tp->Draw();

   auto c1 = new TCanvas("c1", "rawhist", 950, 800);

   // c1->cd();

    gPad->SetLogy();
   h1->SetBinErrorOption(TH1::EBinErrorOpt::kPoisson);
   h1->GetYaxis()->SetTitle("Counts");
   h1->GetXaxis()->SetTitle("Voltage [V]");
   h1->SetTitle("Raw " + info + " from CH1 (biased from Muon Energy)");
   h1->Draw();
   auto tp1 = new TPaveText(0.15, 0.7, 0.35, 0.85, "NDC");
   tp1->AddText("MuLife");
   tp1->AddText(authors);
   // tp1->AddText("Run0 9/03/23 17h");
   tp1->AddText("Run1 14/03/23 14h");
   tp1->AddText("Run2 15/03/23 23h");
   // tp->AddText("Run3 07/03/23 14.5h");
   tp1->Draw();

   auto c2 = new TCanvas("c2", "deltav/v", 950, 800);

   auto g = new TGraph(N, arr_v1, arr_deltav);
   TAxis *axis = g->GetXaxis();

   axis->SetLimits(-0.1,2.5);
   g->SetTitle("#Deltav vs v;CH1 [V]; CH0-CH1 [V]");
   g->Draw("AP");

   auto c3 = new TCanvas("c3", "CH1 vs CH0", 950, 800);

   auto g1 = new TGraph(N, arr_v1, arr_v2);
   TAxis *axis1 = g1->GetXaxis();

   axis1->SetLimits(-0.1,2.5);
   g1->SetTitle("Scatter plot of read signal (CH0 is read 5 #mus after CH1);CH1 [V]; CH0 [V]");
   g1->Draw("AP");

   auto c4 = new TCanvas("c4", "deltav/v histogram", 950, 800);
   gPad->SetLogy();

   auto tp2 = new TPaveText(0.15, 0.7, 0.35, 0.85, "NDC");
   tp2->AddText("MuLife");
   tp2->AddText(authors);
   // tp1->AddText("Run0 9/03/23 17h");
   //tp2->AddText("Run1 14/03/23 14h");
   //tp2->AddText("Run2 15/03/23 23h");
   //tp2->AddText("Overflow: 6");

   h2->GetYaxis()->SetTitle("Counts");
   h2->GetXaxis()->SetTitle("#Deltav/v");
   h2->SetTitle("CH0-CH1 / CH1 for Electron Energy");
   h2->Draw();
   tp2->Draw();

   return;
   
}