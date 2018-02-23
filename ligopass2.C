#include "TH1.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TAxis.h"
#include "TStyle.h"

const double textsize = 0.06;

bool dividebylivetime = false;
const char * outbase = NULL;

TCanvas USN;

TH1D * getordont(const char * const histname, TDirectory * const ligofile)
{
  TH1D * h = dynamic_cast<TH1D *>(ligofile->Get(histname));
  if(h == NULL){
    fprintf(stderr, "No histogram \"%s\" in your file\n", histname);
    return NULL;
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
  x->SetNdivisions(508);

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

  c->  SetLeftMargin(0.20);
  c-> SetRightMargin(0.01);
  c->   SetTopMargin(0.06);
  c->SetBottomMargin(0.13);
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
}

void process(TH1D * tracks, TH1D * trackslive, const char * const title)
{
  const int rebin = 10;
  
  tracks->Rebin(rebin);

  if(trackslive != NULL){
    trackslive->Rebin(rebin);
    divide_livetime(tracks, trackslive);
  }

  stylehist(tracks);
  tracks->GetYaxis()->SetTitle(Form("%s/%ds", title, rebin));

  tracks->Draw("e");
  print2(tracks->GetName());
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
  if(ligodir == NULL) ligodir = dynamic_cast<TDirectory*>(ligofile->Get("ligo"));

  if(ligodir == NULL){
    fprintf(stderr, "No \"ligoanalysis\" directory in this file\n");
    exit(1);
  }

  stylecanvas(&USN);

  #define GETORDONT(x, y) TH1D * x = getordont(#x, ligodir); \
                          TH1D * x##live = getordont(#x"live", ligodir); \
                          if(x != NULL) process(x, x##live, y)

  GETORDONT(rawtrigger,            "Raw triggers");
  GETORDONT(rawhits,               "Raw hits");
  GETORDONT(unslice4ddhits,        "Hits in the noise slice");
  GETORDONT(tracks,                "Slices with tracks");
  GETORDONT(halfcontained_tracks,  "Slices with stopping tracks");
  GETORDONT(fullycontained_tracks, "Slices with contained tracks");
  GETORDONT(contained_slices,      "Contained slices");
  GETORDONT(unslicedbighits,       "Unsliced big hits");
  GETORDONT(unslicedhitpairs,      "Supernova-like events");
  GETORDONT(upmu_tracks,           "Upward going muons");
  GETORDONT(energy_low_cut,        "Very high energy triggers");
  GETORDONT(energy_high_cut,       "Very very high energy triggers");
  GETORDONT(energy_low_cut_pertime,"Very high power events");
  GETORDONT(energy_high_cut_pertime,"Very very high power events");
}
