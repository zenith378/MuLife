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

   TString run0 = path_to_file + "22febbraio2023.dat";
   TString run1 = path_to_file + "run1_23feb23.dat";
   TString run2 = path_to_file + "run2_2023.dat";
   TString run3 = path_to_file + "run3_7mar23.dat";
   TString run4 = path_to_file + "run4_8mar.dat";
   TString run5 = path_to_file + "run5_10mar23.dat";

   TString hname = "Run1long_tm_backcost_0830_mediumbin";
   TString info = "";
   TString date = "23/03/23";
   TString authors = "G. Cordova, A. Giani";
   TString acqtime = "14.5h FAULTY";
   //---------histogram name for fit plot----------//

   // auto df = ROOT::RDF::MakeCsvDataFrame(fname,false,'\t');

   //---------- Define Tree ---------//
   Int_t currentIgnoreLevel = gErrorIgnoreLevel;
   gErrorIgnoreLevel = kError;
   TTree *tree = new TTree("tree", "tree");
   tree->ReadFile(run0, "x/D:y");
   tree->ReadFile(run1);
   tree->ReadFile(run2);
   tree->ReadFile(run5);
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
   auto min = 0.1;
   auto max = 20.;
   auto bins = 20;
   TString ffit;
   ffit.Form("MuLife " + info + "%.2f-%.0f", min, max);

   //--------- Define Histrogram -----//
   TH1D *h = new TH1D("h", hname, bins, min, max);



   //-- Fill Histogram with stop-start signal--//
   auto tmax = 85.8485;
   auto tmaxx = 687.144;
   auto tmax_c=-tmax;
   auto tmaxx_c=-tmaxx;
   for (Int_t i = 0; i < N; i++)
   {
      tree->GetEntry(i);
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
            if (decaytime < tmaxx_c) decaytime = decaytime + tmaxx;
            if (tmaxx_c<decaytime < tmax_c) decaytime = decaytime + tmax;
            
            effective_time = (decaytime+0*pow(10,-9)) * pow(10, 6);
            if (effective_time < max && effective_time > min)
            {
               h->Fill(effective_time);
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

   //-----------------------------BLOCK 2-------------------------------//
   //-------- Constructing the models to be used for fitting -----------//

   //-----Define variable (time) and import histogram ---------//
   RooRealVar t("t", "time [microsec]", 0, 30);
   RooDataHist rh("rh", "rh", t, Import(*h));


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
   RooRealVar bias1("bias1", "bias1", 0.03);
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
   RooAddPdf model("model", "signal and background", RooArgList(decay_gm1, background), RooArgList(fsig));

   model.fixAddCoefNormalization(RooArgSet(t));
   model.fixCoefNormalization(RooArgSet(t));

   //------------------------BLOCK 3---------------------------//
   //------------------ Fiting and drawing --------------------//

   RooFitResult *fitResult = model.fitTo(rh, IntegrateBins(0.0001), // Extended(true),
                                         Verbose(false), Warnings(false), Save(),
                                         PrintEvalErrors(-1), PrintLevel(-1));
   //--------- print result on terminal -------//
   fitResult->Print("v");

   //------------ Plot data and PDF overlaid----------------//
   RooPlot *xframe = t.frame(Title(ffit)); // define frame

   rh.plotOn(xframe, MarkerStyle(6), MarkerSize(1)); // plot data

   model.plotOn(xframe, Components(background), LineColor(41), LineStyle(kDashed));
   model.plotOn(xframe, Components(decay_gm1), LineColor(30), LineStyle(9));

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

   TPad *pad1 = new TPad("pad1", "The pad 80 of the height", 0.0, 0.2, .9999, 1.0);  // divide canvas in 2 (fit)
   TPad *pad2 = new TPad("pad2", "The pad 20 of the height", 0.0, 0.005, 1., 0.25); //(residuals)
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

   auto tp = new TPaveText(0.15, 0.15, 0.4, 0.4, "NDC");
   tp->AddText("MuLife");
   tp->AddText(authors);
   tp->AddText("Run0 22/02/23 26h");
   tp->AddText("Run1 23/02/23 115h FAULTY");
   tp->AddText("Run2 30/02/23 115h FAULTY");
   //tp->AddText("Run4 08/03/23 26h FAKE");
   tp->AddText("Run5 10/03/23 72h");
   tp->Draw();
   //TLatex chi;
   TString chi_string;
   chi_string.Form("#chi_{BC}^{2}: %.2f", chi2_BC);
   TString  ndof;
   ndof.Form("ndof: %d",bins-2);
   TLatex chi;
   TLatex nd;
   chi.SetTextFont(43);
   chi.SetTextSize(25);
   chi.DrawLatexNDC(0.75, 0.7, chi_string);
   nd.SetTextFont(43);
   nd.SetTextSize(25);
   nd.DrawLatexNDC(0.75, 0.63, ndof);
   c1->Update();


   return;
}