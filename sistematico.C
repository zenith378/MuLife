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
   TString path_to_file="Dati/Sistematico/";

   TString fname = path_to_file+"sist_1mar23.dat";

   TString hname = "Risoluzione";
   TString date = "28/02/23";
   TString authors = "G. Cordova, A. Giani";
   TString acqtime = "10 min";
   //---------histogram name for fit plot----------//
   TString ffit = "Sistematico, "+date+" 50 ms";

   // auto df = ROOT::RDF::MakeCsvDataFrame(fname,false,'\t');

   //---------- Define Tree ---------//
   TTree *tree = new TTree("tree", "tree");
   tree->ReadFile(fname, "x/D:y");

   TTree *data_tree = new TTree("data tree", "tree of acquired data");

   Int_t N = tree->GetEntries();

   //---------- Tree Branches -------//
   Double_t time;
   Double_t channel;
   Double_t effective_time;
   tree->SetBranchAddress("x", &channel);
   tree->SetBranchAddress("y", &time);
   data_tree->Branch("eff_time",&effective_time);
   auto min=195;
   auto max=230;
   auto bins=35;
   //--------- Define Histrogram -----//
   TH1D *h = new TH1D("h", hname, bins, min, max);






   //-- Fill Histogram with stop-start signal--//
   auto tmax=85.8485;
   auto tmaxx=687.144;
   for (Int_t i = 0; i < N; i++)
   {
      tree->GetEntry(i);
      if(channel==2){
         auto t2 = time;
         tree->GetEntry(i - 1);
         if(channel==1){
         auto t1 = time; 
         auto decaytime = (t2 - t1);
         effective_time=decaytime*pow(10,9);
         if(decaytime<0) {
            effective_time=(decaytime+tmax)*pow(10,9);
         }
         if(decaytime<0&&t2<90&&t1>620) {
            effective_time=(decaytime-tmax+tmaxx)*pow(10,9);
         }
         data_tree->Fill();
         h->Fill(effective_time);
         }
      }
   }
      
   std::cout <<"Overflow: " << h->GetBinContent(N+1) << std::endl;
   std::cout <<"Underflow: " << h->GetBinContent(0) << std::endl;

   auto c = new TCanvas("c", "raw", 950, 800);
   //gPad->SetLogy();
   h->GetYaxis()->SetTitle("Counts");
   h->GetXaxis()->SetTitle("Time [ns]");
   //h->GetXaxis()->SetNdivisions(-502);

   h->SetTitle("Preliminary histogram for an expected signal of 180 ns");
   h->Draw();
   auto tp = new TPaveText(197,5000,205,6300);
   tp->AddText("MuLife");
   tp->AddText(authors);
   tp->AddText(date+" "+acqtime);
   tp->Draw();
   TArrow *ar = new TArrow(200,400,200,1500,0.02,"<|");

   auto tm = new TText(198,1700,"6% of signal");
   tm->SetTextSize(21);
   tm->SetTextFont(43);
   ar->Draw();
   tm->Draw();

   
   //-----------------------------BLOCK 2-------------------------------//
   //-------- Constructing the models to be used for fitting -----------//

   RooRealVar start("start","start time",0,690);
   RooRealVar stop("stop","stop time",0,690);

   //-----Define variable (time) and import histogram ---------//
   RooRealVar t("t", "time [nanosec]", min, max);
   RooDataHist rh("rh", "rh" ,t , Import(*h));
   RooDataSet data_set("data_set","data",t,Import(*data_tree));


   //-------lifetime variable---------//
   RooRealVar mu1("mu1", "mean of first gaussian", 220., 219., 222.);
   RooRealVar mu2("mu2", "mean of second gaussian", 200., 199., 201.);
   RooRealVar sigma1("sigma1", "sigma of first gaussian", 0.1, 0.000001, 1.);
   RooRealVar sigma2("sigma2", "sigma of second gaussian", 0.1, 0.000001, 1.);

   RooGaussian gaus1("gaus1","first signal",t,mu1,sigma1);
   RooGaussian gaus2("gaus2","second signal",t,mu2,sigma2);
   //----------Resolution function for signal--------------//
   RooRealVar fsig("fsig", "signal component",0.94,0.001,1.);
 
   //-----------final pdf-------------//
   //RooDecay model = decay_tm;
   RooAddPdf model("model","signals",RooArgList(gaus1,gaus2),RooArgList(fsig));





   //------------------------BLOCK 3---------------------------//
   //------------------ Fiting and drawing --------------------//

   RooFitResult *fitResult = model.fitTo(data_set,
                                             Verbose(false), Warnings(false), Save(), 
                                             PrintEvalErrors(-1), PrintLevel(-1));
   //--------- print result on terminal -------//
   fitResult->Print("v");


   //------------ Plot data and PDF overlaid----------------//
   RooPlot *xframe = t.frame(Title("Guassian Fit for an expected signal of 180 ns")); //define frame

   data_set.plotOn(xframe, MarkerStyle(6), MarkerSize(1)); //plot data


   model.plotOn(xframe, Components(gaus1), LineColor(41), LineStyle(kDashed)); 
   model.plotOn(xframe, Components(gaus2), LineColor(30), LineStyle(9));
   
   model.plotOn(xframe, LineWidth(2), LineColor(kRed)); //plot fitted pdf

   //-------------plot parameters on figure----------------//
   RooArgSet display(mu1,mu2,sigma1,sigma2,fsig); //parameters to display on figure
   model.paramOn(xframe,
                 Parameters(display),
                 Layout(0.35, 0.6, 0.9), //position
                 Format("NE", AutoPrecision()));
   //data_set.statOn(xframe,Layout(0.75,0.95,0.95));

   xframe->getAttText()->SetTextSize(0.031); //parameters size and font
   xframe->getAttText()->SetTextFont(42);

   //-----------Pull Plot--------------//
   RooHist *hpull = xframe->pullHist();

   hpull->SetMarkerStyle(6);
   hpull->SetLineWidth(0); //no line

   //-----------Final Canvas----------//
   auto c1 = new TCanvas("Fit", "Expected signal of 180 ns",950,800);


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
   
   auto tp1 = new TPaveText(197,6000,205,7500);
   tp1->AddText("MuLife");
   tp1->AddText(authors);
   tp1->AddText(date+" "+acqtime);
   tp1->Draw();

   pad2->cd();
   pad2->SetBottomMargin(0.4);
   hpull->SetMinimum(-10);
   hpull->SetMaximum(10);
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
   //TString ffit_plot = "./Plots/MuLife/Run1/" + hname + "_Fit" + ".pdf";

   //c1->SaveAs(ffit_plot);
   
   return;
}