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
#include "TGraph.h"
#include "TLine.h"

using namespace RooFit;

void cut()
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
   RooTruthModel tm("tm", "truth model", t);

   //--------------gauss1---------------//
   RooRealVar bias1("bias1", "bias1", 0.039);
   RooRealVar sigma1("sigma1", "sigma1", 0.005);
   RooGaussModel gm1("gm1", "gauss model 1", t, bias1, sigma1);

   // Construct decay(t) (x) gauss1(t)
   RooDecay decay_gm1("decay_gm1", "decay", t, tau, tm, RooDecay::SingleSided);

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
   auto bins = 20;
   auto min = 0.1;
   std::vector<double> chi_vec;
   std::vector<double> t_min;
   std::vector<double> tau_fit;
   while (min < 1)
   {
      //----------------------BLOCK 1------------------------//
      //------------------ Data Reading ---------------------//

      //---------- Define string for data handling----------//
      TString path_to_file = "../Dati/MuLife/";

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
               auto decaytime = (t2 - t1)+40*pow(10,-9);
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
      Double_t nsig1 = static_cast<RooAbsReal &>(lf[1]).getVal();
   std::cout << nsig1 << std::endl;

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

      chi_vec.push_back(chi2_BC);
      t_min.push_back(min);
      min = min + 0.01;
   }
   Int_t nvc = t_min.size();
   Double_t *tt = &t_min[0];
   Double_t *cc = &chi_vec[0];
   auto c = new TCanvas("chi", "chi", 950, 800);
   auto g1 = new TGraph(nvc, tt, cc);
   g1->SetTitle("#chi_{BC}^{2} vs lower time cut;time cut [#mus];#chi_{BC}^{2}");
   //g1->SetMinimum(0.);
   auto l1 = new TLine(0.03, 25.989, 1.05, 25.989);
   auto l2 = new TLine(0.03, 28.869, 1.05, 28.869);
   auto l3 = new TLine(0.03, 18, 1.05, 18);
   l1->SetLineColor(3);
   l1->SetLineWidth(1);
   l2->SetLineColor(kOrange);
   l2->SetLineWidth(1);
   l3->SetLineColor(kRed);
   l3->SetLineWidth(1);
   auto cl1 = new TText(0.8, 23, "significativity 0.05");
   cl1->SetTextSize(21);
   cl1->SetTextFont(43);
   cl1->SetTextColor(3);
   auto cl2 = new TText(0.6, 30, "significativity 0.01");
   cl2->SetTextSize(21);
   cl2->SetTextFont(43);
   cl2->SetTextColor(kOrange);
   auto cl3 = new TText(0.2, 19, "ndof=18");
   cl3->SetTextSize(21);
   cl3->SetTextFont(43);
   cl3->SetTextColor(kRed);
   auto tp = new TPaveText(0.7, 90., 1., 160.);
   tp->AddText("MuLife");
   TString date = "06/03/23";
   TString authors = "G. Cordova, A. Giani";
   tp->AddText(authors);
   tp->AddText("Run0 22/02/23 26h");
   tp->AddText("Run1 23/02/23 115h FAULTY");

   c->SetLogy();
   g1->Draw("AC*");
   l1->Draw();
   l2->Draw();
   l3->Draw();
   tp->Draw();
   cl1->Draw();
   cl2->Draw();
   cl3->Draw();
   return;
}