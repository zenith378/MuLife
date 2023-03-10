#include "TTree.h"
#include "TH1.h"
#include <iostream>
#include "TCanvas.h"
#include "TGraph.h"
#include "TArrow.h"
#include "TText.h"
#include "TPaveText.h"

void sawstructure()
{
   //----------------------BLOCK 1------------------------//
   //------------------ Data Reading ---------------------//

   //---------- Define string for data handling----------//
   TString path_to_file="Dati/Sistematico/";

   TString fname = path_to_file+"risoluzione_28feb.dat";
   // auto df = ROOT::RDF::MakeCsvDataFrame(fname,false,'\t');
   TString date = "01/03/2023";
   TString authors = "G. Cordova, A. Giani";
   //---------- Define Tree ---------//
   TTree *tree = new TTree("tree", "tree");
   tree->ReadFile(fname, "x/D:y");

   Int_t N = tree->GetEntries();

   //---------- Tree Branches -------//
   Double_t time;

   //tree->SetBranchAddress("x", &channel);
   tree->SetBranchAddress("y", &time);

   Double_t t_vec[N];
   Double_t t_real[N];





   //-- Fill Histogram with stop-start signal--//
   auto tmax=85.8485;
   auto tmaxx=687.144;
   for (Int_t i = 0; i < N; i++)
   {
      tree->GetEntry(i);
      t_vec[i]=time;
      t_real[i]=i*0.050;
   }
      
   

   auto c = new TCanvas("c", "sawtooth", 950, 800);
   auto tp = new TPaveText(100,600,500,720);
   tp->AddText("MuLife");
   tp->AddText(authors);
   tp->AddText(date);
   auto L = new TLine(1100,tmaxx,1450,tmaxx);
   L->SetLineColor(kBlue);
   L->SetLineWidth(1);
   auto tmt= new TText(1180,300,"Restart at 687.144 s");
   tmt->SetTextSize(21);
   tmt->SetTextFont(43);
   auto L1 = new TLine (760,345,760,430);
   L1->SetLineColor(6);
   L1->SetLineWidth(1);
   TArrow *ar = new TArrow(760,150,760,300,0.02,"|>");
   TArrow *ar1 = new TArrow(1350,350,1350,670,0.02,"|>");
   TArrow *ar2 = new TArrow(1350,30,1350,270,0.02,"<|");

   auto tm = new TText(600,100,"Saw jump of 85.8485 s");
   tm->SetTextSize(21);
   tm->SetTextFont(43);
   auto g = new TGraph(N,t_real,t_vec);
   g->SetTitle("Saw-tooth structure of FPGA DAQ for 50 ms pulse;real time [s];read time [s]");
   g->Draw("AP");
   tp->Draw();
   L->Draw();
   tmt->Draw();
   L1->Draw();
   ar->Draw();
   ar1->Draw();
   ar2->Draw();
   tm->Draw();

   return;
}