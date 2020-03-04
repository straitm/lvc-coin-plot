#include <iostream>
#include <vector>
#include <stdio.h>
#include <math.h>
#include "TString.h"
#include "TH2.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TStyle.h"

using std::vector;

static char * gwname = NULL, * mass = NULL;

static int dig(const double n)
{
  return (Form("%.1g", n)[0] - '0' <= 2 &&
          Form("%.2g", n)[0] - '0' <= 2)? 2: 1;
}

static double il(const vector<double> & pr)
{
  double totprob = 0, accprob = 0;

  for(unsigned int i = 0; i < pr.size(); i++) totprob += pr[i];

  for(unsigned int i = 0; i < pr.size(); i++){
    accprob += pr[i];
    if(accprob/totprob > 0.9) return i;
  }

  fprintf(stderr, "Error: only found %f, not 90%%\n", accprob/totprob);
  return 0;
}

void draw_flux_limit(const vector<double> & norm,
                     const vector<double> & prob,
                     const double limit)
{
  const double size = 0.067;

  TCanvas * c1 = new TCanvas();
  c1->SetRightMargin(0.04);
  c1->SetTopMargin(0.03);
  c1->SetBottomMargin(0.15);
  c1->SetLeftMargin(0.14);
  c1->SetCanvasSize(600, 400);
  gStyle->SetOptStat(0);
  c1->SetTickx();
  c1->SetTicky();

  const double xmin = 0, xmax = 2,
               ymin = 0, ymax = 1.099;

  const int stepstep = 10;

  TH2D * dum = new TH2D("dum", "", 1, xmin, xmax, 1, ymin, ymax);

  dum->GetXaxis()->SetTickLength(0.02);
  dum->GetYaxis()->SetTickLength(0.02);
  dum->GetXaxis()->SetTitle(Form("%s signal strength relative to 10#kern[-0.5]{ }kpc",
                                 gwname));
  dum->GetYaxis()->SetTitle("Prob density (arb units)");
  dum->GetXaxis()->CenterTitle();
  dum->GetYaxis()->CenterTitle();
  dum->GetXaxis()->SetTitleSize(size);
  dum->GetYaxis()->SetTitleSize(size);
  dum->GetXaxis()->SetLabelSize(size);
  dum->GetYaxis()->SetLabelSize(size);
  dum->GetYaxis()->SetTitleOffset(1.0);
  dum->GetXaxis()->SetNdivisions(8);
  dum->GetYaxis()->SetDecimals();

  dum->Draw();

  double maxprob = 0;
  for(unsigned int i = 0; i < prob.size(); i++)
    if(prob[i] > maxprob)
      maxprob = prob[i];

  TGraph g, gfill;
  for(unsigned int i = 0; i < norm.size(); i+=stepstep){
    if(norm[i] > xmax) continue;

    g.SetPoint(g.GetN(), norm[i], prob[i]/maxprob);
    if(norm[i] < limit)
      gfill.SetPoint(g.GetN(), norm[i], prob[i]/maxprob);
  }
 
  gfill.SetPoint(gfill.GetN(), gfill.GetX()[gfill.GetN()-1], 0);
  gfill.SetPoint(gfill.GetN(), 0, 0);

  g.SetLineWidth(2);
  g.SetLineColor(kBlack);
  g.SetMarkerColor(kBlack);

  gfill.SetFillStyle(1001);
  gfill.SetFillColor(kGray);

  gfill.Draw("lf");
  g.Draw("l");

  TLegend leg(0.31, 0.65, 0.80, 0.93);
  leg.SetTextSize(size);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetTextFont(42);

  leg.AddEntry((TH1D*)NULL,
    Form("%s solar mass Garching", mass), "");
  leg.AddEntry((TH1D*)NULL, "Combined ND+FD", "");

  leg.AddEntry((TH1D*)NULL,
    Form("< %.1f of 10#kern[-0.5]{ }kpc at 90%% CL", limit), "");
  leg.AddEntry((TH1D*)NULL,
    Form("#Rightarrow  > %.0f#kern[-0.5]{ }kpc at 90%% CL", 10/sqrt(limit)), "");

  leg.Draw();

  c1->RedrawAxis();

  c1->SaveAs(Form("fluxlimit-%s-combined-%ssolarmass.pdf", gwname, mass));
}

int main(int argc, char ** argv)
{
  double n = 0, p1 = 0, p2 = 0;

  if(argc > 1) gwname = argv[1];
  if(argc > 1) mass = argv[2];

  vector<double> prob1, prob2, prob, norm;
  while(std::cin >> n >> p1 >> p2){
    prob.push_back(p1*p2);
    prob1.push_back(p1);
    prob2.push_back(p2);
    norm.push_back(n);
  }

  const int itot = il(prob);
  const int i1 = il(prob1);
  const int i2 = il(prob2);

  const double n1 = norm[itot];
  const double n2 = 10/sqrt(norm[itot]);
  const double n3 = 10/sqrt(norm[i1]);
  const double n4 = 10/sqrt(norm[i2]);
  printf("90%% CL at %.*g of 10kpc SN-like flux, %.*g kpc (%.*g, %.*g)",
         dig(n1), n1,
         dig(n2), n2,
         dig(n3), n3,
         dig(n4), n4);

  draw_flux_limit(norm, prob, n1);

  return 0;
}
