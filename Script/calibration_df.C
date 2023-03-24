#include "TString.h"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include <stdio.h>
#include "math.h"
#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TF1.h"
#include "TAxis.h"
#include "TLine.h"
#include "TArrow.h"
#include "TPaveText.h"
#include "TGraph.h"
#include <string.h>
#include <filesystem>
#include "TLatex.h"
#include "TLegend.h"

void calibration_df()
{
   //ROOT::EnableImplicitMT();
   //----------------------BLOCK 1------------------------//
   //------------------ Data Reading ---------------------//

   //---------- Define string for data handling----------//
   TString path_to_file = "../Dati/Energy/";


   TString fname4 = path_to_file + "cal_16mar23.dat";
   TString cal_sotto = path_to_file + "cal_sotto_16mar23.dat";
   TString cal_sopra = path_to_file + "cal_sopra_16mar23.dat";

   TString hname = "h";
   TString info = "Calibration at MIP for Muon Energy higher PMTs";
   TString authors = "G. Cordova, A. Giani";
   //Define tree
   TTree *tree1 = new TTree("energy_spectrum", "Energy spectrum");
   tree1->ReadFile(fname4, "time/D:E_in_v/D:E_fin_v");

   //TTree *mc_tree = new TTree("cal", "monte carlo calibration");
   //mc_tree->ReadFile("../cal.root");

   TFile * f0 = new TFile("../cal.root");
   TTree * t_mc = (TTree *)f0->Get("In the target");
   ROOT::RDataFrame df_mc(*t_mc);


   //binning for histograms
   auto min = -2.5;
   auto max = 2.5;
   auto bins = 100;

   //Dataframe
   ROOT::RDataFrame df(*tree1);
   auto df_new = df.Define("E_in","E_in_v*77.44")
                   .Define("E_fin","E_fin_v*77.44")
                   .Define("Energy", "E_fin- E_in")
                   .Define("DeltaV","(E_fin-E_in)/E_in");
                   //.Filter("E_fin>0.", "Final energy above 0.");
                   //.Filter("E_in<0.7", "Initial energy below 0.7")
                   //.Filter("Energy>0.05","Overall Energy above 0.5")
                   //.Filter("E_fin>0.", "Final energy above 0.");

   auto df_cut = df_new.Filter("E_fin>0.", "Final energy above 0.")
                       //.Filter("E_in<0.7*77.4", "Initial energy below 0.7")
                       .Filter("E_in>0.","Initial energy above 0");
                       //.Filter("(E_fin-E_in)>-0.083*E_in-0.02","Diagonal cut");


   auto df_cor = df_cut.Define("E_fin_cor","E_fin+0.1031*E_in")
                        .Define("E_cor","(E_fin_cor-E_in)");



   auto c = new TCanvas("c", "rawhist", 950, 800);

   auto h1 = df_cut.Histo1D({"h1","Data", 200, 0.1*77.44, 1.7*77.44}, "E_in");
   h1->SetBinErrorOption(TH1::EBinErrorOpt::kPoisson);
   h1->SetMarkerStyle(8);
   h1->SetMarkerSize(1);
   h1->SetLineColor(kBlack);
   h1->SetMarkerColor(kBlack);
   h1->SetStats(0);

   auto h0 = df_mc.Histo1D({"h0","Monte Carlo",200,0.1*77.44,1.7*77.44},"EdepTarget");
   h0->SetTitle("Energy for Muons at MIP (calibration)");
   h0->GetYaxis()->SetTitle("Counts");
   h0->GetXaxis()->SetTitle("Energy [MeV]");
   h0->Scale(0.7343209036);
   h0->SetFillColor(38);
   h0->SetStats(0);
   h0->DrawClone("hist");
   h1->DrawClone("E1 SAME");

   auto tp0 = new TPaveText(0.15, 0.7, 0.35, 0.85, "NDC");
   tp0->AddText("MuLife");
   tp0->AddText(authors);
   tp0->AddText("Calibration 16/03/23 5h");
   tp0->AddText("Cut on E>7.5 MeV to avoid pedestal");
   TArrow *ar0 = new TArrow(0.5*77.44,320,0.75*77.44,320,0.02,"|>");
   TText *t0 = new TText(0.12*77.44,320,"Peak at 62 MeV");
   t0->SetTextSize(21);
   t0->SetTextFont(43);
   ar0->Draw();
   t0->Draw();
   TArrow *ar11 = new TArrow(1.*77.44,240,1.2*77.44,240,0.02,"<|");
   TText *t11 = new TText(1.2*77.44,240,"Horn at 74 MeV due to unequalized PMT");
   t11->SetTextSize(21);
   t11->SetTextFont(43);
   ar11->Draw();
   t11->Draw();
   tp0->Draw();
   auto legend = new TLegend(0.1,0.7,0.48,0.9);
   legend->AddEntry("h1","Data");
   legend->AddEntry("h0","Monte Carlo");
   legend->Draw();
/*

   auto c1 = new TCanvas("c1", "rawhist", 950, 800);

   // c1->cd();

   h1->DrawClone("E1");
   auto tp1 = new TPaveText(0.15, 0.7, 0.35, 0.85, "NDC");
   tp1->AddText("MuLife");
   tp1->AddText(authors);
   tp1->AddText("Calibration 16/03/23 5h");
   tp1->AddText("Cut on CH1>0.1 V to avoid pedestrial");
   TArrow *ar1 = new TArrow(0.5,600,0.75,600,0.02,"|>");
   TText *t1 = new TText(0.12,600,"Peak at 0.8 V");
   t1->SetTextSize(21);
   t1->SetTextFont(43);
   ar1->Draw();
   t1->Draw();
   TArrow *ar11 = new TArrow(1.,400,1.2,400,0.02,"<|");
   TText *t11 = new TText(1.2,400,"Horn at 0.95 V due to unequalized PMT");
   //t11->SetTextSize(21);
   //t11->SetTextFont(43);
   //ar11->Draw();
   //t11->Draw();
   tp1->Draw();

   auto c2 = new TCanvas("c2", "deltav/v", 950, 800);

   auto g = df_new.Graph("E_in","Energy");
   TAxis *axis = g->GetXaxis();

   axis->SetLimits(-0.1*77.44,2.5*77.44);
   g->SetTitle("Charge difference vs Initial Charge;E_{in} [MeV]; E_{fin}-E_{in} [MeV]");
   g->DrawClone("AP");
   auto tp2 = new TPaveText(0.15, 0.7, 0.35, 0.85, "NDC");
   tp2->AddText("MuLife");
   tp2->AddText(authors);
   tp2->AddText("Calibration 16/03/23 5h"); 
   tp2->AddText("E_{fin} is read 5 #mus after E_{in}");
   tp2->Draw();
   TArrow *ar21 = new TArrow(0.02*77.44,0.5*77.44,0.2*77.44,0.5*77.44,0.02,"<|");
   TText *t21 = new TText(0.21*77.44,0.5*77.44,"Trigger before stored event");
   TLine *l2 = new TLine(0.7*77.44,-1.2*77.44,0.7*77.44,0.45*77.44);
   TText *t22 = new TText(0.75*77.44,0.2*77.44,"Threshold for trust region");
   TArrow *ar23 = new TArrow(0.5*77.44,-1.2*77.44,0.5*77.44,-0.7*77.44,0.02,"|>");
   TText *t23 = new TText(0.4*77.44,-1.3*77.44,"Very low decay constant (outliers)");
   t21->SetTextSize(21);
   t21->SetTextFont(43);
   t22->SetTextSize(21);
   t22->SetTextFont(43);
   t23->SetTextSize(21);
   t23->SetTextFont(43);
   ar23->Draw();
   t23->Draw();
   l2->Draw();
   t22->Draw();
   ar21->Draw();
   t21->Draw();





   auto c3 = new TCanvas("c3", "CH1 vs CH0", 950, 800);

   auto g1 = df_new.Graph("E_in","E_fin");
   TAxis *axis1 = g1->GetXaxis();

   axis1->SetLimits(0.,2.5);
   g1->SetTitle("Read signal in CH1 and CH0;E_{in} [MeV]; E_{fin} [MeV]");
   g1->DrawClone("AP");
   auto tp3 = new TPaveText(0.15, 0.7, 0.35, 0.85, "NDC");
   tp3->AddText("MuLife");
   tp3->AddText(authors);
   tp3->AddText("Calibration 16/03/23 5h"); 
   tp3->AddText("CH0 is read 5 #mus after CH1");
   tp3->Draw();

   auto c4 = new TCanvas("c4", "deltav/v histogram", 950, 800);
   gPad->SetLogy();

   auto tp4 = new TPaveText(0.15, 0.7, 0.35, 0.85, "NDC");
   tp4->AddText("MuLife");
   tp4->AddText(authors);
   tp4->AddText("Calibration 16/03/23 5h");
   //tp2->AddText("Overflow: 6");

   auto h2 = df_new.Histo1D({"h2","deltav/v",bins,min*77.44,max*77.44},"DeltaV");
   h2->GetYaxis()->SetTitle("Counts");
   h2->GetXaxis()->SetTitle("(E_{fin}-E_{in})/E_{in}");
   h2->SetTitle("CH0-CH1 / CH1 for Muon Energy");
   h2->DrawClone();
   tp4->Draw();




   auto c5 = new TCanvas("c5", "deltav/v", 950, 800);

   auto g5 = df_cut.Graph("E_in","Energy");
   TAxis *axis5 = g5->GetXaxis();

   axis5->SetLimits(0.,0.8*77.44);
   g5->SetTitle("Charge difference vs Initial Charge;CH1 [V]; CH0-CH1 [V]");
   auto tp5 = new TPaveText(0.15, 0.7, 0.35, 0.85, "NDC");
   tp5->AddText("MuLife");
   tp5->AddText(authors);
   tp5->AddText("Calibration 16/03/23 5h"); 
   tp5->AddText("CH0 is read 5 #mus after CH1");
   tp5->AddText("Cut on CH0>0 and 0<CH1<0.7");

   gStyle->SetOptFit(1111);

   TF1 *f = new TF1("f", "[0] * x",0.,0.7*77.44); 
   g5->Fit(f,"R");
   g5->DrawClone("AP");
   tp5->Draw();
   gPad->Update();
   //TPaveStats *p = g5->GetListOfFunctions()->FindObject("stats");
   //gPad->Modified(); gPad->Update();

   f->Draw("SAME");


   auto c6 = new TCanvas("c6","resolution fit", 950, 800);

   auto h6 = df_cor.Histo1D({"h6","Energy corrected fit",500,-10,10},"E_cor");
   h6->SetBinErrorOption(TH1::EBinErrorOpt::kPoisson);
   h6->GetYaxis()->SetTitle("Counts");
   h6->GetXaxis()->SetTitle("E_{fin}-E_{in} corrected [MeV]");
   h6->SetTitle("E_{fin}-E_{in} corrected for slope");
   auto tp6 = new TPaveText(0.15, 0.7, 0.35, 0.85, "NDC");
   tp6->AddText("MuLife");
   tp6->AddText(authors);
   tp6->AddText("Calibration 16/03/23 5h"); 
   tp6->AddText("E_{fin} is read 5 #mus after E_{in}");
   tp6->AddText("Cut on E_{fin}>0 and 0<E_{in}<0.7");
   h6->DrawClone("E1");
   tp6->Draw();


   auto c7 = new TCanvas("c7", "v", 950, 800);

   auto g7 = df_cor.Graph("E_in","E_cor");
   TAxis *axis7 = g7->GetXaxis();

   axis7->SetLimits(0.,0.8*77.44);
   g7->SetTitle("Charge difference vs Initial Charge corrected for slope;CH1 [V]; CH0-CH1 [V]");
   g7->DrawClone("AP");
   auto tp7 = new TPaveText(0.15, 0.7, 0.35, 0.85, "NDC");
   tp7->AddText("MuLife");
   tp7->AddText(authors);
   tp7->AddText("Calibration 16/03/23 5h"); 
   tp7->AddText("CH0 is read 5 #mus after CH1");
   tp7->AddText("Cut on CH0>0 and 0<CH1<0.7");

*/
   auto report = df_new.Report(); // request cut report

   report->Print();           // print cut report on the terminal




}