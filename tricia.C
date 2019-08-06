/*
 * "I think I wanted to see 1d histograms of the #events found per
 * timestamp for each of the trigger classes. One histogram from the 8am
 * triggers (from which your background thresholds were established),
 * then superimposed, the # of events found in each of the GW-1hour
 * timestamps, to show they are consistent."
 */

#include "TGaxis.h"
#include "TError.h"
#include "TLatex.h"
#include "TMinuit.h"
#include "Math/Math.h"
#include "Math/ProbFuncMathCore.h"
#include "Math/QuantFuncMathCore.h"
#include "TH1.h"
#include "TF1.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TAxis.h"
#include "TStyle.h"
#include "TEllipse.h"
#include <fstream>

#define SEPARATEPDFS

const double textsize = 0.06;

TCanvas bitwash;
const char * outbase = NULL;

void printstart()
{
  static bool started = false;
  if(!started)
    bitwash.Print(Form("%s.pdf[", outbase));
  started = true;
}

void printend()
{
  bitwash.Print(Form("%s.pdf]", outbase));

  printf("------------- Finished processing -------------\n");
}

void print2(const char * const name)
{
#ifdef SEPARATEPDFS
  bitwash.Print(Form("%s-%s.pdf", outbase, name));
#endif
  bitwash.Print(Form("%s.pdf", outbase));
}

void stylecanvas()
{
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  const double scale = 5.0;
  bitwash.SetCanvasSize(150*scale, 100*scale);

  gStyle->SetGridColor(kGray);
  gStyle->SetGridStyle(kSolid);

  bitwash.SetBottomMargin(0.143);
  bitwash.SetLeftMargin(0.155);
  bitwash.SetRightMargin(0.085);
  bitwash.SetTopMargin(0.02);
  bitwash.SetTickx(1);
  bitwash.SetTicky(1);
}

void stylehist(TH1D * h)
{
  TAxis* y = h->GetYaxis();
  TAxis* x = h->GetXaxis();

  h->SetLineColor(kBlack);
  h->SetMarkerColor(kBlack);

  h->SetLineWidth(2);

  y->CenterTitle();
  x->CenterTitle();

  x->SetLabelSize(textsize);
  y->SetLabelSize(textsize);
  x->SetTitleSize(textsize);
  y->SetTitleSize(textsize);

  y->SetDecimals();

  x->SetTickSize(0.03);
  y->SetTickSize(0.015);

  x->SetNdivisions(508);
  y->SetNdivisions(510);
}

void Add(TH1D & hist, const char * const filename,
         const char * const histname, const bool livetimediv)
{
  TFile * ligofile = new TFile(filename, "read");
  if(!ligofile || ligofile->IsZombie()){
    fprintf(stderr, "Could not open your histogram file %s\n", filename);
    exit(1);
  }

  TDirectory * ligodir =
    dynamic_cast<TDirectory*>(ligofile->Get("ligoanalysis"));

  if(ligodir == NULL){
    fprintf(stderr, "No \"ligoanalysis\" directory in %s\n", filename);
    exit(1);
  }

  TH1D * inhist = dynamic_cast<TH1D*>(ligodir->Get(histname));
  inhist->Sumw2();
  if(inhist == NULL){
    fprintf(stderr, "No histogram named %s in %s\n", histname, filename);
    exit(1);
  }

  if(livetimediv){
    TH1D * live = dynamic_cast<TH1D*>(ligodir->Get("trackslive"));
    if(live == NULL){
      fprintf(stderr, "No histogram named trackslive in %s\n", filename);
      exit(1);
    }
    for(int i = 1; i <= inhist->GetNbinsX(); i++)
      if(live->GetBinContent(i) > 0)
        hist.Fill(inhist->GetBinContent(i)/live->GetBinContent(i));
  }
  else{
    for(int i = 1; i <= inhist->GetNbinsX(); i++)
      hist.Fill(inhist->GetBinContent(i));
  }

  delete inhist;
  delete ligodir;
  delete ligofile;
}

struct opts{
  string histname;
  int nbins;
  double low, high;
  string xtitle;
};

opts mkopts(string histname_, int nbins_, double low_, double high_, string xtitle_)
{
  opts o;
  o.histname = histname_;
  o.nbins = nbins_;
  o.low = low_;
  o.high = high_;
  o.xtitle = xtitle_;
  return o;
}

double minz(TH1D * h)
{
  double ans = 1e300;
  for(int i = 1; i < h->GetNbinsX(); i++){
    const double b = h->GetBinContent(i);
    if(b > 0 && b < ans) ans = b;
  }
  return ans;
}

void tricia(const std::string & bgfilelistfile,
            const std::string & infilename,
            const std::string & trigname,
            const std::string & pdfbase, const bool livetimediv)
{
  outbase = pdfbase.c_str();

  TH1::SetDefaultSumw2();

  stylecanvas();

  vector<opts> todolist;
  if(trigname == "FD 10Hz trigger"){
    todolist.push_back(mkopts("halfcontained_tracks_point_0", 100, 0, 100,
                              "FD 10Hz, Stopping tracks, 1.3#circ, raw/s"));
    todolist.push_back(mkopts("fullycontained_tracks", 20, 0, 20,
                              "FD 10Hz, Contained tracks, raw/s"));
    todolist.push_back(mkopts("fullycontained_tracks_point_0", 20, 0, 20,
                              "FD 10Hz, Contained tracks, 1.3#circ, raw/s"));
    todolist.push_back(mkopts("fullycontained_tracks_point_1", 20, 0, 20,
                              "FD 10Hz, Contained tracks, 16#circ, raw/s"));
    todolist.push_back(mkopts("contained_slices", 20, 0, 20,
                              "FD 10Hz, Contained slices, raw/s"));
    todolist.push_back(mkopts("upmu_tracks", 20, 0, 20,
                              "FD 10Hz, Upward-going muons, raw/s"));
    todolist.push_back(mkopts("upmu_tracks_point_0", 20, 0, 20,
                              "FD 10Hz, Upward-going muons, 1.3#circ, raw/s"));
    todolist.push_back(mkopts("upmu_tracks_point_1", 20, 0, 20,
                              "FD 10Hz, Upward-going muons, 16#circ, raw/s"));
  }
  else if(trigname == "ND energy"){
    todolist.push_back(mkopts("tracks_point_0", 20, 0, 20,
                              "ND energy, Tracks, 1.3#circ, raw/s"));
    todolist.push_back(mkopts("halfcontained_tracks", 20, 0, 20,
                              "ND energy, Stopping tracks, raw/s"));
    todolist.push_back(mkopts("halfcontained_tracks_point_1", 20, 0, 20,
                              "ND energy, Stopping tracks, 16#circ, raw/s"));
    todolist.push_back(mkopts("halfcontained_tracks_point_0", 20, 0, 20,
                              "ND energy, Stopping tracks, 1.3#circ, raw/s"));
    todolist.push_back(mkopts("fullycontained_tracks", 20, 0, 20,
                              "ND energy, Contained tracks, raw/s"));
    todolist.push_back(mkopts("fullycontained_tracks_point_0", 20, 0, 20,
                              "ND energy, Contained tracks, 1.3#circ, raw/s"));
    todolist.push_back(mkopts("fullycontained_tracks_point_1", 20, 0, 20,
                              "ND energy, Contained tracks, 16#circ, raw/s"));
    todolist.push_back(mkopts("contained_slices", 20, 0, 20,
                              "ND energy, Contained slices, raw/s"));
  }
  else if(trigname == "FD energy"){
    todolist.push_back(mkopts("energy_high_cut", 20, 0, 20,
                              "FD energy, >2.5 ADC total/50#mus/s"));
    todolist.push_back(mkopts("energy_high_cut_pertime", 20, 0, 20,
                              "FD energy, >25 ADC total/50#mus/s"));
    todolist.push_back(mkopts("energy_vhigh_cut", 20, 0, 20,
                              "FD energy, >500 ADC total/s"));
    todolist.push_back(mkopts("energy_vhigh_cut_pertime", 20, 0, 20,
                              "FD energy, >250 ADC total/50#mus/s"));
  }

  bitwash.SetLogy();
  
  for(unsigned int i = 0; i < todolist.size(); i++){
    opts & o = todolist[i];

    TH1D * bg  = new TH1D("bg",  "", o.nbins, o.low, o.high);
    TH1D * sig = new TH1D("sig", "", o.nbins, o.low, o.high);

    Add(*sig, infilename.c_str(), o.histname.c_str(), livetimediv);
   
    ifstream bglist(bgfilelistfile.c_str());
    if(!bglist.is_open()){
      fprintf(stderr, "Couldn't open %s\n", bgfilelistfile.c_str());
      exit(1);
    }

    int nbgfiles = 0;
    std::string bgfilename;
    while(bglist >> bgfilename)
      Add(*bg, bgfilename.c_str(), o.histname.c_str(), livetimediv);

    printstart();

    sig->SetLineColor(kBlack);
    bg->SetLineColor(kRed);

    stylehist(sig);

    sig->GetXaxis()->SetTitle(o.xtitle.c_str());
    sig->GetYaxis()->SetTitle("Fraction of bins");

    if(bg->Integral()) bg->Scale(1/bg->Integral());
    if(sig->Integral()) sig->Scale(1/sig->Integral());

    const double ks = 100*sig->KolmogorovTest(bg);
    const double chi2 = 100*sig->Chi2Test(bg, "uw");

    sig->Draw("e");
    bg->Draw("sameehist");

    const double mn = minz(bg)*0.7,
                 mx = sig->GetMaximum()*1.3;
    if(mx > mn) sig->GetYaxis()->SetRangeUser(mn, mx);

    TLatex * l = new TLatex(0, 0, "");
    l->SetTextSize(textsize);
    l->SetTextFont(42);
    if(sig->GetEntries()){
      if(bg->Integral() == bg->GetBinContent(1) &&
        sig->Integral() == sig->GetBinContent(1)){
        l->DrawLatexNDC(0.3, 0.8, "No signal or background");
      }
      else{
        l->SetTextColor(ks < 5?kRed:kBlack);
        l->DrawLatexNDC(0.6, 0.8, Form("KS %.1f%%", ks));
        //l->SetTextColor(chi2 < 5?kRed:kBlack);
        //l->DrawLatexNDC(0.6, 0.7, Form("#chi^{2} %.1f%%", chi2));
      }
    }
   
    TLatex * l2 = new TLatex(0, 0, "");
    l2->SetTextSize(textsize);
    l2->SetTextFont(42);
    if(sig->GetEntries())
      l2->DrawLatexNDC(0.2, 0.9,
        Form("%s", pdfbase.substr(pdfbase.find('-')+1).c_str()));
   
    print2(o.histname.c_str());

    delete bg;
    delete sig;
  }
  printend();
}
