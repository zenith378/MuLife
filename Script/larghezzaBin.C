#include "TGraph.h"
#include "TF1.h"
#include <vector>
#include "TMath.h"
#include "TCanvas.h"
#include <iostream>

Double_t fInfo_loss(Double_t x, Double_t  pars)
{
    return 100*(1 - (  ( (pars/x)*(pars/x)*TMath::Exp(-pars/x) ) / ((1 - TMath::Exp(-pars/x))*(1 - TMath::Exp(-pars/x)))  ));
}

void larghezzaBin()
{
    Double_t fit_range = 20.-0.45; //in microsecondi
    Double_t tau = 2.;//approssimazione della vita media del muone nella palstica (in microsecondi)
    Double_t c = fit_range/tau;
    
    std::vector<Double_t> n_bins = {1., 5., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100.};
    std::vector<Double_t> res;

    //scorro su tutto l'array n_bins e salvo nella variabile val (sto fornendo l'indirizzo a cui salvare)
    for (auto  &val : n_bins)
    {
        res.push_back(fInfo_loss(val, c));
    }

    
    TGraph * graph = new TGraph(n_bins.size(), &n_bins[0], &res[0]);

    TCanvas * c1 = new TCanvas();

    graph->SetMarkerStyle(20);
    graph->SetMarkerColor(kRed);
    graph->SetMarkerSize(0.7);
    graph->SetTitle("Information loss vs number of bins; Number of bins; #deltaI [%]");
    //graph->GetXaxis()->SetLimits(0.,110.);
    //graph->GetYaxis()->SetRangeUser(0.,100.);

    graph->Draw("AP");

    c1->SaveAs("Plots/larghezzaBin.pdf");






}