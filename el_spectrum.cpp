//script per la genrazione dello spettro degli elettroni con endpoint a 52.8 MeV e a 61.6 MeV

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"

#include <iostream>
#include <vector>

void el_spectrum()
{
    TFile * f0 = new TFile("./el_spectrum0.root");
    TFile * f1 = new TFile("./el_spectrum1.root");
    TFile * f2 = new TFile("./el_spectrum2.root");
    TFile * f3 = new TFile("./el_spectrum3.root");
    TFile * f4 = new TFile("./el_spectrum4.root");
    TFile * f5 = new TFile("./el_spectrum5.root");
    TFile * f6 = new TFile("./el_spectrum6.root");

    //endpoint = 52.8 MeV
    TH1F * h1 = new TH1F("h1", "h1", 100, 0, 100);
    TH1F * h1Corr = new TH1F("h1Corr", "h1Corr", 100, 0, 100);
    TH1F * h2 = new TH1F("h2", "h2", 200, 0, 3000);

    //endpoint = 61.6 MeV
    TH1F * h3 = new TH1F("h3", "h3", 100, 0, 100);
    TH1F * h4 = new TH1F("h4", "h4", 200, 0, 3000);

    TTree * t0 = (TTree *)f0->Get("In the target");
    TTree * t1 = (TTree *)f1->Get("In the target");
    TTree * t2 = (TTree *)f2->Get("In the target");
    TTree * t3 = (TTree *)f3->Get("In the target");
    TTree * t4 = (TTree *)f4->Get("In the target");
    TTree * t5 = (TTree *)f5->Get("In the target");
    TTree * t6 = (TTree *)f6->Get("In the target");
    

    std::vector<double> vEdep, vPath, vEdepCorr;//vettori per riempire gli istogrammi a 52.8 MeV
    double Edep0, Edep1, Edep2, Edep3, Edep4, Edep5, Edep6, path0, path1, path2, path3, path4, path5, path6;

    t0->SetBranchAddress("EdepTarget", &Edep0);
    t0->SetBranchAddress("TargetPath", &path0);

    t1->SetBranchAddress("EdepTarget", &Edep1);
    t1->SetBranchAddress("TargetPath", &path1);

    t2->SetBranchAddress("EdepTarget", &Edep2);
    t2->SetBranchAddress("TargetPath", &path2);

    t3->SetBranchAddress("EdepTarget", &Edep3);
    t3->SetBranchAddress("TargetPath", &path3);
    
    t4->SetBranchAddress("EdepTarget", &Edep4);
    t4->SetBranchAddress("TargetPath", &path4);
    
    t5->SetBranchAddress("EdepTarget", &Edep5);
    t5->SetBranchAddress("TargetPath", &path5);

    t6->SetBranchAddress("EdepTarget", &Edep6);
    t6->SetBranchAddress("TargetPath", &path6);


    for (int i = 0; i < t0->GetEntries(); i++){
        t0->GetEntry(i);
        vEdep.push_back(Edep0);
        vPath.push_back(path0);
        
    }

    for (int i = 0; i < t1->GetEntries(); i++){
        t1->GetEntry(i);
        vEdep.push_back(Edep1);
        vPath.push_back(path1);
        
    }

    for (int i = 0; i < t2->GetEntries(); i++){
        t2->GetEntry(i);
        vEdep.push_back(Edep2);
        vPath.push_back(path2);
        
    }


    for (int i = 0; i < t3->GetEntries(); i++){
        t3->GetEntry(i);
        vEdep.push_back(Edep3);
        vPath.push_back(path3);
        
    }

    for (int i = 0; i < t4->GetEntries(); i++){
        t4->GetEntry(i);
        vEdep.push_back(Edep4);
        vPath.push_back(path4);
        
    }

    for (int i = 0; i < t5->GetEntries(); i++){
        t5->GetEntry(i);
        vEdep.push_back(Edep5);
        vPath.push_back(path5);
        
    }

    //istogramma per endpoint a 52.8 MeV
    for (double val : vEdep){
        
        h1->Fill(val);
        
    }

    for (double tmp : vPath){
        h2->Fill(tmp);
        
        //std::cout << val << std::endl;
    }

    
    


    for (int i = 0; i < t6->GetEntries(); i++){
        t6->GetEntry(i);
        vEdep.push_back(Edep6);
        vPath.push_back(path6);
        
    }

    //istogramma per endpoint a 61.6 MeV
    for (double val : vEdep){
        
        h3->Fill(val);
        
    }

    for (double tmp : vPath){
        h4->Fill(tmp);
        
    }
    

    TCanvas * c1  = new TCanvas("c1", "c1", 1);
    //h1->SetTitle("Energy loss through the target - 100000 generated muons; E [MeV]; Entries");
    h1->SetTitle("Energy loss through the target - end point: 52.8 MeV; E [MeV]; Entries");
    h1->Draw();

    TCanvas * c2  = new TCanvas("c2", "c2", 1);
    //h2->SetTitle("Total path inside the target - 100000 generated muons; Path [mm]; Entries");
    h2->SetTitle("Total path inside the target - end point: 52.8 MeV; Path [mm]; Entries");
    h2->Draw();

    TCanvas * c3  = new TCanvas("c3", "c3", 1);
    //h1->SetTitle("Energy loss through the target - 100000 generated muons; E [MeV]; Entries");
    h3->SetTitle("Energy loss through the target - end point: 61.6 MeV; E [MeV]; Entries");
    h3->Draw();

    TCanvas * c4  = new TCanvas("c4", "c4", 1);
    //h2->SetTitle("Total path inside the target - 100000 generated muons; Path [mm]; Entries");
    h4->SetTitle("Total path inside the target - end point: 61.6 MeV; Path [mm]; Entries");
    h4->Draw();
}