#include "TTree.h"
#include "TH1.h"
#include <iostream>
#include "TCanvas.h"
#include "TGraph.h"
#include "TPaveText.h"
#include "TArrow.h"

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
         if(decaytime<0&&t1<620) {
            effective_time=(decaytime+tmax)*pow(10,9);
         }
         if(decaytime<0&&t2<90&&t1>620) {
            effective_time=(decaytime+tmaxx)*pow(10,9);
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

   h->SetTitle("Preliminary histogram for a pulse of 180 ns");
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

   /*
   //-----------------------------BLOCK 2-------------------------------//
   //-------- Constructing the models to be used for fitting -----------//

   RooRealVar start("start","start time",0,690);
   RooRealVar stop("stop","stop time",0,690);

   //-----Define variable (time) and import histogram ---------//
   RooRealVar t("t", "time [microsec]", 0, 30);
   RooDataHist rh("rh", "rh" ,t , Import(*h));
   RooDataSet raw_data("raw_data","columnar dataset of t1 and t2",data_tree,RooArgSet(start,stop));
   //RooDataSet data();

   //-------lifetime variable---------//
   RooRealVar tau("tau", "mean life of muon", 2.2, 1.5, 3.);

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
   //RooDecay background("background","background exponential",t,tauback,tm1, RooDecay::SingleSided);
   //RooConstVar a0("a0","constant background",5);
   RooPolynomial background("background","background polynomial",t,RooArgList());

   //-----------final pdf-------------//
   //RooDecay model = decay_tm;
   RooAddPdf model("model","signal and background",RooArgList(decay_tm,background),RooArgList(fsig));





   //------------------------BLOCK 3---------------------------//
   //------------------ Fiting and drawing --------------------//

   RooFitResult *fitResult = model.fitTo(rh,
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
   pad1->SetLogy();
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
   TString ffit_plot = "./Plots/MuLife/Run1/" + hname + "_Fit" + ".pdf";

   c1->SaveAs(ffit_plot);
   */
   return;
}