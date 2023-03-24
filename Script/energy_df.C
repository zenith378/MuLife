#include "TString.h"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include <stdio.h>
#include "math.h"

#include "TH1.h"
#include "TCanvas.h"
#include "TStyle.h"

#include "TAxis.h"
#include "RooPolynomial.h"
#include "RooHist.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooFitResult.h"
#include "RooAddPdf.h"
#include "TPaveText.h"
#include "TGraph.h"
#include <string.h>
#include "RooPlot.h"
#include <filesystem>
#include "TLatex.h"
#include "TF1.h"
using namespace RooFit;

void energy_df()
{
   //ROOT::EnableImplicitMT();
   //----------------------BLOCK 1------------------------//
   //------------------ Data Reading ---------------------//

   //---------- Define string for data handling----------//
   TString path_to_file = "../Dati/Energy/";

   TString run0 = path_to_file + "spettro_elett_9mar23.dat";
   TString run1 = path_to_file + "spettro_elet_14mar23.dat";
   TString run11 = path_to_file + "spettro2_14mar23.dat";
   TString run2 = path_to_file + "spettro_elettrone_15mar.dat";
   TString run3 = path_to_file + "spettro_16mar23.dat";


   //Define tree
   TTree *tree1 = new TTree("energy_spectrum", "Energy spectrum");
   tree1->ReadFile(run1, "time/D:E_in/D:E_fin");
   tree1->ReadFile(run0);
   tree1->ReadFile(run11);
   tree1->ReadFile(run2);
   tree1->ReadFile(run3);

   //binning for histograms
   auto min = -2.5;
   auto max = 2.5;
   auto bins = 100;



   //Dataframe
   ROOT::RDataFrame df(*tree1);
   auto df_new = df.Define("Energy", "E_fin- E_in")
                   .Define("DeltaVsuV","(E_fin-E_in)/E_in");
                   //.Filter("E_fin>0.", "Final energy above 0.");
                   //.Filter("E_in<0.7", "Initial energy below 0.7")
                   //.Filter("Energy>0.05","Overall Energy above 0.5")
                   //.Filter("E_fin>0.", "Final energy above 0.");

   auto df_cut = df_new.Filter("E_fin>0.", "Final energy above 0.")
                       .Filter("E_in<0.7", "Initial energy below 0.7")
                       .Filter("E_in>0.","Initial energy above 0");
                       //.Filter("(E_fin-E_in)>-0.083*E_in-0.02","Diagonal cut");


   auto df_cor = df_cut.Define("E_fin_cor","E_fin+0.0838*E_in")
                       .Define("E_cor","E_fin_cor-E_in")
                       .Filter("E_cor>0.05","Overall corrected energy above 0");



   auto report = df.Report(); // request cut report

   //auto c = new TCanvas("c", "", 800, 700);




   TString authors = "G. Cordova, A. Giani";

   auto c2 = new TCanvas("c2", "deltav/v", 950, 800);

   auto g = df_cor.Graph("E_in","E_cor");
   TAxis *axis = g->GetXaxis();

   axis->SetLimits(-0.1,2.5);
   g->SetTitle("#Deltav vs v;E_in [V]; E_fin-E_in [V]");
   g->DrawClone("AP");

   auto c3 = new TCanvas("c3", "CH1 vs CH0", 950, 800);

   auto g1 = df_new.Graph("E_in","E_fin");
   TAxis *axis1 = g1->GetXaxis();

   axis1->SetLimits(-0.1,2.5);
   g1->SetTitle("Scatter plot of read signal (CH0 is read 5 #mus after CH1);CH1 [V]; CH0 [V]");
   g1->DrawClone("AP");




   auto c = new TCanvas("c", "Electron Energy Spectrum", 950, 800);
   //gPad->SetLogy();
   auto h = df_cor.Histo1D({"h", "Electron Energy Spectrum", 50, 0.1, 1.5}, "E_cor");

   h->SetBinErrorOption(TH1::EBinErrorOpt::kPoisson);
   h->GetYaxis()->SetTitle("Counts");
   h->GetXaxis()->SetTitle("Voltage [V]");
   h->SetTitle("Energy in target E_{fin}>0. && 0.<E_{in}<0.7 and E_{tot}>0.1");
   auto tp = new TPaveText(0.15, 0.7, 0.35, 0.85, "NDC");
   tp->AddText("MuLife");
   tp->AddText(authors);
   tp->AddText("Run0 9/03/23 17h");
   tp->AddText("Run1 14/03/23 14h");
   tp->AddText("Run2 15/03/23 23h");
   tp->AddText("Run3 16/03/23 86h");


   TF1 *f = new TF1("f", "[0] * x +[1]",0.28,0.7); 
   f->SetParameters(-500,300);
   f->SetParLimits(0,-600,-300);
   f->SetParLimits(1,200,500);
   h->Fit(f,"ILR");
   gPad->Update();
   //TPaveStats *p = g5->GetListOfFunctions()->FindObject("stats");
   //gPad->Modified(); gPad->Update();
   h->DrawClone("E1");
   tp->Draw();

   f->Draw("SAME");

//Fattore di conversione 77.44


   auto c1 = new TCanvas("c1", "rawhist", 950, 800);

   // c1->cd();
   auto h1 = df_new.Histo1D({"h1","Energy from CH1 (JIT acquired)", 50, -0.1, 0.8}, "E_in");
    gPad->SetLogy();
   h1->SetBinErrorOption(TH1::EBinErrorOpt::kPoisson);
   h1->GetYaxis()->SetTitle("Counts");
   h1->GetXaxis()->SetTitle("Voltage [V]");
   h1->SetTitle("Raw target energy from Muon");
   h1->DrawClone();
   auto tp1 = new TPaveText(0.15, 0.7, 0.35, 0.85, "NDC");
   tp1->AddText("MuLife");
   tp1->AddText(authors);
   tp1->AddText("Run0 9/03/23 17h");
   tp1->AddText("Run1 14/03/23 14h");
   tp1->AddText("Run2 15/03/23 23h");
   tp1->AddText("Run3 16/03/23 86h");
   tp1->Draw();




   auto c4 = new TCanvas("c4", "deltav/v histogram", 950, 800);
   gPad->SetLogy();

   auto tp2 = new TPaveText(0.15, 0.7, 0.35, 0.85, "NDC");
   tp2->AddText("MuLife");
   tp2->AddText(authors);
   tp2->AddText("Run0 9/03/23 17h");
   tp2->AddText("Run1 14/03/23 14h");
   tp2->AddText("Run2 15/03/23 23h");
   tp2->AddText("Run3 16/03/23 86h");
   //tp2->AddText("Overflow: 6");






   auto h2 = df_new.Histo1D({"h2","deltav/v",bins,min,max},"DeltaVsuV");
   h2->GetYaxis()->SetTitle("Counts");
   h2->GetXaxis()->SetTitle("#Deltav/v");
   h2->SetTitle("CH0-CH1 / CH1 for Electron Energy");
   h2->DrawClone();
   tp2->Draw();





   report->Print();           // print cut report on the terminal




/*

   //-----Define variable (time) and import histogram ---------//
   RooRealVar E("E", "Energy [V]", 0.25, .7);
   RooDataHist rh("rh", "rh", E, Import(*h));


   //----------Resolution function for signal--------------//
   RooRealVar fsig("fsig", "signal component", 0.01, 0.001, 0.1);

   RooRealVar a0("a0","offset signal",300,100,500);
   RooRealVar a1("a1","slope of signal",-500,-300,-800);
   RooPolynomial signal("signal","signal of decay",E,RooArgList(a1));
   //-------Resolution function for background-------------//
   // RooDecay background("background","background exponential",t,tauback,tm1, RooDecay::SingleSided);
   // RooConstVar a0("a0","constant background",5);
   //RooPolynomial background("background", "background polynomial", E, RooArgList());

   //-----------final pdf-------------//
   // RooDecay model = decay_tm;
   RooAddPdf model("model", "signal and background", RooArgList(signal));

   //E.setRange("sig", 0.3, 1.5);

   //model.fixAddCoefNormalization(E);
   //model.fixCoefNormalization(E);

   //------------------------BLOCK 3---------------------------//
   //------------------ Fiting and drawing --------------------//

   RooFitResult *fitResult = signal.fitTo(rh, //IntegrateBins(0.0001), // Extended(true),
                                         Verbose(false), Warnings(false), Save(),
                                         PrintEvalErrors(-1), PrintLevel(-1));
   //--------- print result on terminal -------//
   fitResult->Print("v");

   //------------ Plot data and PDF overlaid----------------//
   RooPlot *xframe = E.frame(Title("Fit")); // define frame

   rh.plotOn(xframe, MarkerStyle(6), MarkerSize(1)); // plot data

   //model.plotOn(xframe, Components(background), LineColor(kYellow), LineStyle(kDashed));
   signal.plotOn(xframe, Components(signal), LineColor(kGreen), LineStyle(9));

   //model.plotOn(xframe, LineWidth(2), LineColor(kRed)); // plot fitted pdf

   // compute Baker-Cousins chi2:
   // converting model to a TF1 and data to a TH1

   auto fn = signal.asTF(E, RooArgList(), E);
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
   
   RooArgSet display(a0, a1, fsig); // parameters to display on figure
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
   auto c9 = new TCanvas("Fit", "ffit", 950, 800);

   TPad *pad1 = new TPad("pad1", "The pad 80 of the height", 0.0, 0.2, .9999, 1.0);  // divide canvas in 2 (fit)
   TPad *pad2 = new TPad("pad2", "The pad 20 of the height", 0.0, 0.005, 1., 0.25); //(residuals)
   pad1->Draw();
   pad2->Draw();
   pad1->cd();

   //------------Bellurie (Fuso cc)------------//
   xframe->GetYaxis()->SetTitleOffset(1.5);
   xframe->GetXaxis()->SetTitleSize(0);
   //pad1->SetLogy();
   // xframe->SetMinimum(0.001);
   xframe->Draw();

   pad2->cd();
   pad2->SetBottomMargin(0.4);
   hpull->SetMinimum(-4);
   hpull->SetMaximum(4);
   hpull->GetYaxis()->SetNdivisions(4);
   hpull->GetXaxis()->SetTitleOffset(1.3);
   hpull->GetYaxis()->SetTitle("Pull");
   hpull->GetXaxis()->SetTitle("Energy [V]");
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

   auto tp7 = new TPaveText(0.15, 0.15, 0.4, 0.4, "NDC");
   tp7->AddText("MuLife");
   tp7->AddText(authors);
   tp7->AddText("Run0 9/03/23 17h");
   tp7->AddText("Run1 14/03/23 14h");
   tp7->AddText("Run2 15/03/23 23h");
   tp7->AddText("Run3 16/03/23 86h");
   tp7->Draw();
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



*/

}