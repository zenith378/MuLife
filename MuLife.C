#include "TTree.h"
#include "TString.h"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include <stdio.h>
#include "math.h"
#include "RooRealVar.h"
#include "RooGaussModel.h"
#include "RooTruthModel.h"
#include "RooDecay.h"
#include "RooPlot.h"
#include "TH1.h"
#include "RooFitResult.h"
#include "RooDataHist.h"
#include "RooHist.h"
#include "TCanvas.h"
#include "TStyle.h"

#include "TAxis.h"

#include "TRootCanvas.h"
#include <string.h>
#include <filesystem>

using namespace RooFit;

void MuLife()
{
   //----------------------BLOCK 1------------------------//
   //------------------ Data Reading ---------------------//

   //---------- Define string for data handling----------//
   TString path_to_file="Dati/MuLife/";

   TString fname = path_to_file+"22febbraio2023.dat";

   TString hname = "Run0long_08-30_background_esponenziale";

   //---------histogram name for fit plot----------//
   TString ffit = "MuLife, run0 22/02/23 26h (truth model + background) 0.8-30";

   // auto df = ROOT::RDF::MakeCsvDataFrame(fname,false,'\t');

   //---------- Define Tree ---------//
   TTree *tree = new TTree("tree", "tree");
   tree->ReadFile(fname, "x/D:y");

   Int_t N = tree->GetEntries();

   //---------- Tree Branches -------//
   Double_t time;
   Double_t channel;
   tree->SetBranchAddress("x", &channel);
   tree->SetBranchAddress("y", &time);
   auto min=0.8;
   auto max=30.;
   auto bins=30;
   //--------- Define Histrogram -----//
   TH1D *h = new TH1D("h", hname, bins, min, max);

   TH1D *h1 = new TH1D("h1", "t1 distribution", 100, 0, 750);

   TH1D *h2 = new TH1D("h", "t2 distribution", 100, 0, 750);





   //-- Fill Histogram with stop-start signal--//
   auto tmax=70;
   auto tmaxx=690;
   for (Int_t i = 0; i < N; i++)
   {
      tree->GetEntry(i);
      if(channel ==1){
         auto ttemp=time;
         h1->Fill(ttemp);
      }
      if (channel == 2)
      {
         auto t2 = time;
         auto channel2= channel;
         tree->GetEntry(i - 1);
         auto t1 = time;
         auto channel1=channel;
         /*if(channel2-channel1!=1){
         std::cout<< "Channel1: "<< channel1 << std::endl;
         std::cout<< "Channel2: "<< channel2 << std::endl;
         }*/
         //h1->Fill(t1);
         h2->Fill(t2);
         if(channel2==2&&channel1==1) {
         auto decaytime = (t2 - t1) * pow(10, 6);
         if(decaytime<0) {
            decaytime=decaytime+tmax;
            if(t1>680&&t2<70) decaytime=decaytime-tmax+tmaxx;
         }
         if(decaytime<max&&decaytime>min) h->Fill(decaytime);
         }
      }
   }
   auto c = new TCanvas("c", "rawhist", 950, 800);
   gPad->SetLogy();
   h->Draw();
   auto ch1 = new TCanvas("ch1", "t1 distribution", 950, 800);

   h1->Draw();
   
   auto ch2 = new TCanvas("ch2", "t2 distribution", 950, 800);

   h2->Draw();






   //-----------------------------BLOCK 2-------------------------------//
   //-------- Constructing the models to be used for fitting -----------//

   //-----Define variable (time) and import histogram ---------//
   RooRealVar t("t", "time [microsec]", 0, 30);
   RooDataHist rh("rh", "rh" ,t , Import(*h));

   //-------lifetime variable---------//
   RooRealVar tau("tau", "mean life of muon", 2.2, 1.5, 3.5);

   //----------Resolution function for signal--------------//
   RooRealVar fsig("fsig", "signal component",0.5,0.2,0.99);
   //--------------delta----------------//
   // Build a truth resolution model (delta function)
   RooTruthModel tm1("tm", "truth model", t);

   // Construct decay(t) (x) delta(t)
   RooDecay decay_tm("decay_tm", "decay", t, tau, tm1, RooDecay::SingleSided);

   //--------------gauss1---------------//
   RooRealVar bias1("bias1", "bias1", 0);
   RooRealVar sigma1("sigma1", "sigma1", 1);
   RooGaussModel gm1("gm1", "gauss model 1", t, bias1, sigma1);

   // Construct decay(t) (x) gauss1(t)
   RooDecay decay_gm1("decay_gm1", "decay", t, tau, gm1, RooDecay::SingleSided);

   //-------backgorund variable-------//
   RooRealVar tauback("tauback","background constant",100,10,100000);

   //-------Resolution function for background-------------//
   RooDecay background("background","background exponential",t,tauback,tm1, RooDecay::SingleSided);


   //-----------final pdf-------------//
   //RooDecay model = decay_tm;
   RooAddPdf model("model","signal and background",RooArgList(decay_tm,background),RooArgList(fsig));





   //------------------------BLOCK 3---------------------------//
   //------------------ Fiting and drawing --------------------//

   RooFitResult *fitResult = model.fitTo(rh, RecoverFromUndefinedRegions(1), 
                                             Verbose(false), Warnings(false), Save(), 
                                             PrintEvalErrors(-1), PrintLevel(-1));
   //--------- print result on terminal -------//
   fitResult->Print("v");


   //------------ Plot data and PDF overlaid----------------//
   RooPlot *xframe = t.frame(Title(ffit)); //define frame

   rh.plotOn(xframe, MarkerStyle(6), MarkerSize(1)); //plot data


   model.plotOn(xframe, Components(background), LineColor(41), LineStyle(kDashed)); 
   model.plotOn(xframe, Components(decay_tm), LineColor(30), LineStyle(9));
   
   model.plotOn(xframe, LineWidth(2), LineColor(kRed)); //plot fitted pdf

   //-------------plot parameters on figure----------------//
   RooArgSet display(tau,tauback,fsig); //parameters to display on figure
   model.paramOn(xframe,
                 Parameters(display),
                 Layout(0.45, 0.6, 0.9), //position
                 Format("NE", AutoPrecision()));
   rh.statOn(xframe,Layout(0.8,0.99,0.9));

   xframe->getAttText()->SetTextSize(0.031); //parameters size and font
   xframe->getAttText()->SetTextFont(42);

   //-----------Pull Plot--------------//
   RooHist *hpull = xframe->pullHist();

   hpull->SetMarkerStyle(6);
   hpull->SetLineWidth(0); //no line

   //-----------Final Canvas----------//
   auto c1 = new TCanvas("Fit", ffit, 950, 800);


   TPad *pad1 = new TPad("pad1", "The pad 80 of the height", 0.0, 0.2, 1.0, 1.0); //divide canvas in 2 (fit)
   TPad *pad2 = new TPad("pad2", "The pad 20 of the height", 0.0, 0.005, 1, 0.25); //(residuals)
   pad1->Draw();
   pad2->Draw();
   pad1->cd();

   //------------Bellurie (Fuso cc)------------//
   xframe->GetYaxis()->SetTitleOffset(1.5);
   xframe->GetXaxis()->SetTitleSize(0);
   //pad1->SetLogy();
   //xframe->SetMinimum(0.001);
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
   c1->Update();


   //---------- Save Canvas ---------------//
   TString ffit_plot = "./Plots/MuLife/" + hname + "_Fit" + ".pdf";

   c1->SaveAs(ffit_plot);
   return;
}