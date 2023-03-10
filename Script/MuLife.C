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

void MuLife()
{
   //----------------------BLOCK 1------------------------//
   //------------------ Data Reading ---------------------//

   //---------- Define string for data handling----------//
   TString path_to_file = "../Dati/MuLife/";

   TString fname = path_to_file + "22febbraio2023.dat";
   TString fname1 = path_to_file + "run1_23feb23.dat";
   TString fname2 = path_to_file + "run4_8mar.dat";

   TString hname = "Run1long_tm_backcost_0830_mediumbin";
   TString info = "FAKE";
   TString date = "23/03/23";
   TString authors = "G. Cordova, A. Giani";
   TString acqtime = "14.5h FAULTY";
   //---------histogram name for fit plot----------//

   // auto df = ROOT::RDF::MakeCsvDataFrame(fname,false,'\t');

   //---------- Define Tree ---------//
   Int_t currentIgnoreLevel = gErrorIgnoreLevel;
   gErrorIgnoreLevel = kError;
   TTree *tree = new TTree("tree", "tree");
   tree->ReadFile(fname, "x/D:y");
   tree->ReadFile(fname1);
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
   auto min = 0.45;
   auto max = 20.;
   auto bins = 20;
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
   auto ch1 = new TCanvas("ch1", "t1 distribution", 950, 800);
   h1->GetYaxis()->SetTitle("Counts");
   h1->GetXaxis()->SetTitle("Time [s]");
   h1->SetTitle("Counts CH1 " + info);
   h1->Draw();

   auto ch2 = new TCanvas("ch2", "t2 distribution", 950, 800);
   h2->GetYaxis()->SetTitle("Counts");
   h2->GetXaxis()->SetTitle("Time [s]");
   h2->SetTitle("Counts CH2 " + info);
   h2->Draw();

   //-----------------------------BLOCK 2-------------------------------//
   //-------- Constructing the models to be used for fitting -----------//

   RooRealVar start("start", "start time", 0, 690);
   RooRealVar stop("stop", "stop time", 0, 690);

   //-----Define variable (time) and import histogram ---------//
   RooRealVar t("t", "time [microsec]", 0, 10000);
   RooDataHist rh("rh", "rh", t, Import(*h));
   RooDataSet raw_data("raw_data", "columnar dataset of t1 and t2", data_tree, RooArgSet(start, stop));
   // RooDataSet data();

   //-------lifetime variable---------//
   RooRealVar tau("tau", "mean life of muon", 2.2, 0.001, 10.);

   //----------Resolution function for signal--------------//
   RooRealVar fsig("fsig", "signal component", 0.6, 0.01, 0.99);
   //--------------delta----------------//
   // Build a truth resolution model (delta function)
   RooTruthModel tm1("tm", "truth model", t);

   // Construct decay(t) (x) delta(t)
   RooDecay decay_tm("decay_tm", "decay", t, tau, tm1, RooDecay::SingleSided);

   //--------------gauss1---------------//
   RooRealVar bias1("bias1", "bias1", 0.039);
   RooRealVar sigma1("sigma1", "sigma1", 0.005);
   RooGaussModel gm1("gm1", "gauss model 1", t, bias1, sigma1);

   // Construct decay(t) (x) gauss1(t)
   RooDecay decay_gm1("decay_gm1", "decay", t, tau, gm1, RooDecay::SingleSided);

   //-------backgorund variable-------//
   RooRealVar tauback("tauback", "background constant", 100, 10, 100000);

   //-------Resolution function for background-------------//
   // RooDecay background("background","background exponential",t,tauback,tm1, RooDecay::SingleSided);
   // RooConstVar a0("a0","constant background",5);
   RooPolynomial background("background", "background polynomial", t, RooArgList());

   RooNumConvPdf back_conv("back_conv","background with gaus",t,background,gm1);
   //-----------final pdf-------------//
   // RooDecay model = decay_tm;
   RooAddPdf model("model", "signal and background", RooArgList(decay_tm, background), RooArgList(fsig));

   model.fixAddCoefNormalization(RooArgSet(t));
   model.fixCoefNormalization(RooArgSet(t));

   //------------------------BLOCK 3---------------------------//
   //------------------ Fiting and drawing --------------------//

   RooFitResult *fitResult = model.fitTo(rh, IntegrateBins(0.000001), // Extended(true),
                                         Verbose(false), Warnings(false), Save(),
                                         PrintEvalErrors(-1), PrintLevel(-1));
   //--------- print result on terminal -------//
   fitResult->Print("v");

   //------------ Plot data and PDF overlaid----------------//
   RooPlot *xframe = t.frame(Title(ffit)); // define frame

   rh.plotOn(xframe, MarkerStyle(6), MarkerSize(1)); // plot data

   model.plotOn(xframe, Components(background), LineColor(41), LineStyle(kDashed));
   model.plotOn(xframe, Components(decay_tm), LineColor(30), LineStyle(9));

   model.plotOn(xframe, LineWidth(2), LineColor(kRed)); // plot fitted pdf

   // compute Baker-Cousins chi2:
   // converting model to a TF1 and data to a TH1

   auto fn = model.asTF(t, RooArgList(), t);
   // auto h1 = data->createHistogram("x");

   // pdf is normalized, need to scale it to number of events * bin width
   auto f = new TF1(
       "f1", [&](double *t, double *)
       { return h->Integral() * h->GetBinWidth(1) * fn->EvalPar(t, nullptr); },
       0, 30, 0);

   new TCanvas();
   // need to use DrawClone since RooFit object are deleted at the end of the macro
   h->DrawClone(); // SISTEMA ERRORI ISTOGRAMMA
   f->DrawClone("Same");

   double chi2_N = h->Chisquare(f);       // Neyman chi2
   double chi2_BC = h->Chisquare(f, "L"); // Baker-Cousins chi2
   std::cout << "Baker-Cousins chi2 = " << chi2_BC << " Neyman chi2 " << chi2_N << std::endl;

   //-------------plot parameters on figure----------------//
   
   RooArgSet display(tau, tauback, fsig); // parameters to display on figure
   model.paramOn(xframe,
                 Parameters(display),
                 Layout(0.45, 0.6, 0.9), // position
                 Format("NE", AutoPrecision()));
   rh.statOn(xframe, Layout(0.8, 0.95, 0.92),
             What("N"));

   //-----------Pull Plot--------------//
   RooHist *hpull = xframe->pullHist();

   hpull->SetMarkerStyle(6);
   hpull->SetLineWidth(0); // no line

   //-----------Final Canvas----------//
   auto c1 = new TCanvas("Fit", ffit, 950, 800);

   TPad *pad1 = new TPad("pad1", "The pad 80 of the height", 0.0, 0.2, 1.0, 1.0);  // divide canvas in 2 (fit)
   TPad *pad2 = new TPad("pad2", "The pad 20 of the height", 0.0, 0.005, 1, 0.25); //(residuals)
   pad1->Draw();
   pad2->Draw();
   pad1->cd();

   //------------Bellurie (Fuso cc)------------//
   xframe->GetYaxis()->SetTitleOffset(1.5);
   xframe->GetXaxis()->SetTitleSize(0);
   pad1->SetLogy();
   // xframe->SetMinimum(0.001);
   xframe->Draw();

   pad2->cd();
   pad2->SetBottomMargin(0.4);
   hpull->SetMinimum(-4);
   hpull->SetMaximum(4);
   hpull->GetYaxis()->SetNdivisions(4);
   hpull->GetXaxis()->SetTitleOffset(1.3);
   hpull->GetYaxis()->SetTitle("Pull");
   hpull->GetXaxis()->SetTitle("Lifetime [#mus]");
   hpull->GetXaxis()->SetLabelFont(43);
   hpull->GetXaxis()->SetLabelSize(21);
   hpull->GetYaxis()->SetLabelFont(43);
   hpull->GetYaxis()->SetLabelSize(21);
   hpull->GetXaxis()->SetTitleSize(21);
   hpull->GetXaxis()->SetTitleFont(43);
   hpull->GetYaxis()->SetTitleSize(21);
   hpull->GetYaxis()->SetTitleFont(43);
   hpull->SetTitle("");
   hpull->Draw();

   pad1->cd();

   xframe->getAttText()->SetTextSize(0.031); // parameters size and font
   xframe->getAttText()->SetTextFont(43);
   xframe->getAttText()->SetTextSize(21);

   auto tp = new TPaveText(0.15, 0.15, 0.4, 0.35, "NDC");
   tp->AddText("MuLife");
   tp->AddText(authors);
   tp->AddText("Run0 22/02/23 26h");
   tp->AddText("Run1 23/02/23 115h FAULTY");
   //tp->AddText("Run4 FAKE 08/03/23 26h");
   tp->Draw();

   TLatex chi;
   TString chi_string;
   chi_string.Form("#chi_{BC}^{2}: %.2f", chi2_BC);
   // chi.SetTextFont(43);
   // chi.SetTextSize(91);
   chi.DrawLatexNDC(0.75, 0.8, chi_string);

   c1->Update();

   //---------- Save Canvas ---------------//
   // TString ffit_plot = "./Plots/MuLife/Run1/" + hname + "_Fit" + ".pdf";

   // c1->SaveAs(ffit_plot);
   return;
}