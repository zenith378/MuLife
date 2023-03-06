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
#include "TLatex.h"
#include <vector>
#include "TGraphErrors.h"
#include "TLine.h"

using namespace RooFit;

void binning()
{

   //-----------------------------BLOCK 2-------------------------------//
   //-------- Constructing the models to be used for fitting -----------//

   //-----Define variable (time) and import histogram ---------//
   RooRealVar t("t", "time [microsec]", 0, 20);

   //-------lifetime variable---------//
   RooRealVar tau("tau", "mean life of muon", 2.2, 0.1, 10.);

   //----------Resolution function for signal--------------//
   RooRealVar fsig("fsig", "signal component", 0.6, 0.01, 0.99);
   //--------------delta----------------//
   // Build a truth resolution model (delta function)

   //--------------gauss1---------------//
   RooRealVar bias1("bias1", "bias1", 0.039);
   RooRealVar sigma1("sigma1", "sigma1", 0.005);
   RooGaussModel gm1("gm1", "gauss model 1", t, bias1, sigma1);

   // Construct decay(t) (x) gauss1(t)
   RooDecay decay_gm1("decay_gm1", "decay", t, tau, gm1, RooDecay::SingleSided);

   //-------Resolution function for background-------------//
   // RooDecay background("background","background exponential",t,tauback,tm1, RooDecay::SingleSided);
   // RooConstVar a0("a0","constant background",5);
   RooPolynomial background("background", "background polynomial", t, RooArgList());

   //-----------final pdf-------------//
   // RooDecay model = decay_tm;
   RooAddPdf model("model", "signal and background", RooArgList(decay_gm1, background), RooArgList(fsig));

   model.fixAddCoefNormalization(RooArgSet(t));
   model.fixCoefNormalization(RooArgSet(t));

   auto max = 20.;
   auto bins = 5;
   auto min = 1.;
   std::vector<double> bins_vec;
   std::vector<double> tau_fit;
   std::vector<double> err;
   std::vector<double> err_x;
   Int_t entries;
   while (bins < 30)
   {
      //----------------------BLOCK 1------------------------//
      //------------------ Data Reading ---------------------//

      //---------- Define string for data handling----------//
      TString path_to_file = "Dati/MuLife/";

      TString fname = path_to_file + "22febbraio2023.dat";
      TString fname1 = path_to_file + "run1_23feb23.dat";

      //---------histogram name for fit plot----------//

      // auto df = ROOT::RDF::MakeCsvDataFrame(fname,false,'\t');

      //---------- Define Tree ---------//

      // Silence warnings during TTree::ReadFile
      Int_t currentIgnoreLevel = gErrorIgnoreLevel;
      gErrorIgnoreLevel = kError;
      TTree *tree = new TTree("tree", "tree");
      tree->ReadFile(fname, "x/D:y");
      tree->ReadFile(fname1);
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
      //--------- Define Histrogram -----//
      TH1D *h = new TH1D("h", "hname", bins, min, max);

      //-- Fill Histogram with stop-start signal--//
      auto tmax = 85.8485;
      auto tmaxx = 687.144;
      for (Int_t i = 0; i < N; i++)
      {
         tree->GetEntry(i);
         if (channel == 1)
         {
            auto ttemp = time;
         }
         if (channel == 2)
         {
            auto t2 = time;

            auto channel2 = channel;
            tree->GetEntry(i - 1);
            auto t1 = time;
            auto channel1 = channel;

            if (channel2 == 2 && channel1 == 1)
            {
               time_start = t1;
               time_stop = t2;
               auto decaytime = (t2 - t1);
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
      entries=h->GetEntries();
      //------------------------BLOCK 3---------------------------//
      //------------------ Fiting and drawing --------------------//
      RooDataHist rh("rh", "rh", t, Import(*h));

      RooFitResult *fitResult = model.fitTo(rh, IntegrateBins(0.000001), // Extended(true),
                                            Verbose(false), Warnings(false), Save(),
                                            PrintEvalErrors(-1), PrintLevel(-1));
      //--------- print result on terminal -------//
      // fitResult->Print("v");

   // list of the parameter values of the fitted function
   RooArgList lf = fitResult->floatParsFinal();
   // covariant matrix
   TMatrixDSym cov = fitResult->covarianceMatrix();
      Double_t tauu = static_cast<RooAbsReal &>(lf[1]).getVal();
      Double_t err_tau = std::sqrt(cov[1][1]);
   //std::cout << nsig1 << std::endl;

      // compute Baker-Cousins chi2:
      // converting model to a TF1 and data to a TH1

      auto fn = model.asTF(t, RooArgList(), t);
      // auto h1 = data->createHistogram("x");

      // pdf is normalized, need to scale it to number of events * bin width
      auto f = new TF1(
          "f1", [&](double *t, double *)
          { return h->Integral() * h->GetBinWidth(1) * fn->EvalPar(t, nullptr); },
          0, 30, 0);

      double chi2_BC = h->Chisquare(f, "L"); // Baker-Cousins chi2

      bins_vec.push_back(bins);
      tau_fit.push_back(tauu);
      err.push_back(err_tau);
      err_x.push_back(0.5);
      bins = bins + 1;
   }
   Int_t nvc = bins_vec.size();
   Double_t *tt = &bins_vec[0];
   Double_t *cc = &tau_fit[0];
   Double_t *tt_err = &err[0];
   Double_t *x_err = &err_x[0];

   auto c = new TCanvas("chi", "chi", 950, 800);
   auto g1 = new TGraphErrors(nvc, tt, cc,x_err,tt_err);
   g1->SetTitle("Stability of fitted parameter for binning change;Number of bins;tau [#mus]");
   //g1->SetMinimum(0.);
   auto tp = new TPaveText(0.6, 0.6, 0.85, 0.85, "NDC");
   tp->AddText("MuLife");
   tp->AddText("G. Cordova, A. Giani");
   tp->AddText("Run0 22/02/23 26h");
   tp->AddText("Run1 23/02/23 116h (FAULTY)");
   TString entr_str;
   entr_str.Form("Entries: %d",entries);
   tp->AddText(entr_str);
   tp->AddText("Cuts: 1 #mus < t < 20 #mus");
   tp->Draw();

   //c->SetLogy();
   g1->Draw("AP");

   tp->Draw();

   return;
}