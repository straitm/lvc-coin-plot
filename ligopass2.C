#include "TH1.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TAxis.h"
#include "TStyle.h"

const double textsize = 0.05;

bool dividebylivetime = false;
const char * outbase = NULL;

TCanvas USN;

TH1D * getordie(const char * const histname, TDirectory * const ligofile)
{
  TH1D * h = dynamic_cast<TH1D *>(ligofile->Get(histname));
  if(!h){
    fprintf(stderr, "No histogram \"%s\" in your file\n", histname);
    exit(1);
  }
  return h;
}

// Count of something that helps me decide whether to put markers on the plot
double effectivecount(TH1D * h)
{
  double c = 0;
  for(int i = 1; i <= h->GetNbinsX(); i++){
    if(h->GetBinError(i) == 0) continue;
    c += h->GetBinContent(i)/h->GetBinError(i);
  }
  return c;
}

void stylehist(TH1D * h)
{
  TAxis* y = h->GetYaxis();
  TAxis* x = h->GetXaxis();

  y->SetTitleOffset(1.45);

  h->SetLineColor(kBlack);
  h->SetMarkerColor(kBlack);

  y->CenterTitle();
  x->CenterTitle();

  x->SetLabelSize(textsize);
  y->SetLabelSize(textsize);
  x->SetTitleSize(textsize);
  y->SetTitleSize(textsize);

  x->SetTitle("Time since GW event (s)");

  if(effectivecount(h) < h->GetNbinsX()){
    h->SetMarkerStyle(kOpenCircle);
  }
}

void stylecanvas(TCanvas * c)
{
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  const double scale = 5.0;
  c->SetCanvasSize(128*scale, 96*scale);

  c->  SetLeftMargin(0.14);
  c-> SetRightMargin(0.03);
  c->   SetTopMargin(0.06);
  c->SetBottomMargin(0.11);
}

void divide_livetime(TH1D * d, TH1D * live)
{
  if(!dividebylivetime) return;

  for(int i = 1; i <= live->GetNbinsX(); i++){
    const double bc = d->GetBinContent(i);
    const double be = d->GetBinError(i);
    if(live->GetBinContent(i) == 0){
      d->SetBinContent(i, 0);
      d->SetBinError  (i, 0);
    }
    else{
      d->SetBinContent(i, bc/live->GetBinContent(i));
      d->SetBinError  (i, be/live->GetBinContent(i));
    }
  }
}

void print2(const char * const name)
{
  USN.Print(Form("%s-%s.pdf", outbase, name));
  USN.Print(Form("%s.pdf", outbase));
}

void process_fullytracks(TH1D * tracks, TH1D * trackslive)
{
  divide_livetime(tracks, trackslive);

  stylehist(tracks);
  tracks->GetYaxis()->SetTitle("Number of slices with contained tracks/s");

  tracks->Draw("e");
  print2("containedtracks");
}

void process_halftracks(TH1D * tracks, TH1D * trackslive)
{
  divide_livetime(tracks, trackslive);

  stylehist(tracks);
  tracks->GetYaxis()->SetTitle("Number of slices with stopping tracks/s");

  tracks->Draw("e");
  print2("stoppingtracks");
}

void process_tracks(TH1D * tracks, TH1D * trackslive)
{
  divide_livetime(tracks, trackslive);

  stylehist(tracks);
  tracks->GetYaxis()->SetTitle("Number of tracks/s");

  tracks->Draw("e");
  print2("tracks");
}

void process_unslicedhits(TH1D * unsliced, TH1D * unslicedlive)
{
  divide_livetime(unsliced, unslicedlive);

  stylehist(unsliced);
  unsliced->GetYaxis()->SetTitle("Number of unsliced hits/s");

  unsliced->Draw("e");
  print2("unslicedhits");
}

void process_rawhits(TH1D * rawhits, TH1D * rawhitslive)
{
  divide_livetime(rawhits, rawhitslive);

  stylehist(rawhits);
  rawhits->GetYaxis()->SetTitle("Number of raw hits/s");

  rawhits->Draw("e");
  print2("rawhits");
}

int ligopass2(const char * const infilename, const char * outbase_,
              const bool dividebylivetime_)
{
  outbase = outbase_;
  dividebylivetime = dividebylivetime_;

  TFile * ligofile = new TFile(infilename, "read");
  if(!ligofile || ligofile->IsZombie()){
    fprintf(stderr, "Could not open your histogram file %s\n", infilename);
    exit(1);
  }

  TDirectory * ligodir = dynamic_cast<TDirectory*>(ligofile->Get("ligoanalysis"));

  if(!ligodir){
    fprintf(stderr, "No \"ligoanalysis\" directory in this file\n");
    exit(1);
  }

  #define GETORDIE(x) TH1D * x = getordie(#x, ligodir); \
                      TH1D * x##live = getordie(#x"live", ligodir)

  GETORDIE(rawhits);
  GETORDIE(unslice4ddhits);
  GETORDIE(tracks);
  GETORDIE(halfcontained_tracks);
  GETORDIE(fullycontained_tracks);

  stylecanvas(&USN);
  USN.Print(Form("%s.pdf[", outbase));
  process_rawhits(rawhits, rawhitslive);
  process_unslicedhits(unslice4ddhits, unslice4ddhitslive);
  process_tracks(tracks, trackslive);
  process_halftracks(halfcontained_tracks, halfcontained_trackslive);
  process_fullytracks(fullycontained_tracks, fullycontained_trackslive);
  USN.Print(Form("%s.pdf]", outbase));
}
