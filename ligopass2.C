#include "TLatex.h"
#include "TMinuit.h"
#include "Math/Math.h"
#include "Math/ProbFuncMathCore.h"
#include "TH1.h"
#include "TF1.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TAxis.h"
#include "TStyle.h"

const double textsize = 0.0666;

bool dividebylivetime = false;
const char * outbase = NULL;

TCanvas USN;
TPad * botpad, * midpad, * toppad;
const double divheight1 = 0.26, divheight2 = 0.47;

// Thanks ROOT for making this hard
// Ratio of top pane to middle pane
const double textratiomid = (1-divheight2)/(divheight2 - divheight1);
// Ratio of top pane to bottom pane
const double textratiobot = (1-divheight2)/divheight1;
// Ratio of top pane to whole canvas
const double textratiofull = (1-divheight2);

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

void stylehist(TH1D * h, const int can /* 0: full, 1: top, 2: mid, 3: bot */)
{
  TAxis* y = h->GetYaxis();
  TAxis* x = h->GetXaxis();

  y->SetTitleOffset(textsize/0.04);

  h->SetLineColor(kBlack);
  h->SetMarkerColor(kBlack);

  h->SetLineWidth(2);

  y->CenterTitle();
  x->CenterTitle();

  const double siz = can == 0? textsize:
    can == 1? textsize:
    can == 2? textsize*textratiomid:
              textsize*textratiobot;

  x->SetLabelSize(siz);
  y->SetLabelSize(siz);
  x->SetTitleSize(siz);
  y->SetTitleSize(siz);

  x->SetTickSize(0.03);
  y->SetTickSize(0.015);

  x->SetTitle("Time since GW event (s)");
  x->SetNdivisions(508);
  y->SetNdivisions(can == 0?510: can == 1? 508: 504);

  if(effectivecount(h) < h->GetNbinsX()){
    h->SetMarkerStyle(kOpenCircle);
  }
}

const double leftmargin = 0.12 * textsize/0.06;
const double rightmargin = 0.03;

void stylepad(TPad * pad)
{
  pad->SetBorderMode(0);
  pad->SetBorderSize(2);
  pad->SetTickx(1);
  pad->SetTicky(1);
  pad->SetFrameLineWidth(2);
  pad->SetRightMargin(rightmargin);
  pad->SetLeftMargin(leftmargin);
}

void stylecanvas(TCanvas * c)
{
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  const double scale = 5.0;
  c->SetCanvasSize(100*scale, 100*scale);

  if(dividebylivetime){
    toppad = new TPad("div", "div", 0, divheight2, 1,          1);
    midpad = new TPad("raw", "raw", 0, divheight1, 1, divheight2);
    botpad = new TPad("liv", "liv", 0,          0, 1, divheight1);
    stylepad(toppad);
    stylepad(midpad);
    stylepad(botpad);

    const double topmargin = 0.09 * textsize/0.07;

    toppad->SetTopMargin(topmargin);
    toppad->SetBottomMargin(0);

    midpad->SetBottomMargin(0);
    midpad->SetTopMargin(0);

    botpad->SetTopMargin(0);
    botpad->SetBottomMargin(0.078/divheight1);

    toppad->Draw();
    botpad->Draw();
    midpad->Draw();
  }
}

TH1D * divide_livetime(TH1D * d, TH1D * live)
{
  TH1D * answer = (TH1D*)d->Clone(Form("%sdiv", d->GetName()));

  for(int i = 1; i <= live->GetNbinsX(); i++){
    const double bc = d->GetBinContent(i);
    const double be = d->GetBinError(i);
    if(live->GetBinContent(i) == 0){
      answer->SetBinContent(i, 0);
      answer->SetBinError  (i, 0);
    }
    else{
      answer->SetBinContent(i, bc/live->GetBinContent(i));
      answer->SetBinError  (i, be/live->GetBinContent(i));
    }
  }
  return answer;
}

void print2(const char * const name)
{
  USN.Print(Form("%s-%s.pdf", outbase, name));
}


TH1D * fithist = NULL, * fithistlive = NULL;
static void fcn(__attribute__((unused)) int & np,
  __attribute__((unused)) double * gin, double & like, double *par,
  __attribute__((unused)) int flag)
{
  like = 0;

  const double flat = par[0];

  for(int i = 1; i <= fithist->GetNbinsX(); i++){

    const double model = flat * fithistlive->GetBinContent(i);
    const double data = fithist->GetBinContent(i);

    if(model > 0){
      like += model - data;
      if(data > 0) like += data * log(data/model);
    }
  }

  like *= 2;
}

double getpar(TMinuit & mn, int i) // 0-indexed!
{
  double answer, dum;
  mn.GetParameter(i, answer, dum);
  return answer;
}

double geterr(TMinuit & mn, int i) // 0-indexed!
{
  double val, err;
  mn.GetParameter(i, val, err);
  return err;
}

void bumphunt(TH1D * hist, TH1D * histlive, const bool verbose)
{
  TMinuit mn(1);
  mn.SetPrintLevel(-1);
  mn.SetFCN(fcn);
  int ierr;
  mn.mnparm(0, "flat", hist->Integral()/hist->GetNbinsX(), 10, 0, 0, ierr);

  fithist = hist;
  fithistlive = histlive;
  mn.Command("SET STRATEGY 2");
  mn.Command("MIGRAD");
  const double flat = getpar(mn, 0);

  toppad->cd();

  static TF1 * flatf = new TF1("flatf", "[0]", hist->GetBinLowEdge(1),
                        hist->GetBinLowEdge(hist->GetNbinsX()+1));
  flatf->SetParameter(0, flat);
  flatf->SetLineColor(kRed);
  flatf->SetNpx(500);
  flatf->SetLineWidth(2);
  flatf->Draw("same");

  if(verbose){
    double minprob = 1;
    double minprob_bin = 0;
    for(int i = 1; i <= hist->GetNbinsX(); i++){
      if(histlive->GetBinContent(i) == 0) continue;

      const double expected = flat * histlive->GetBinContent(i);
      const int actual = (int)hist->GetBinContent(i);

      const double prob = ROOT::Math::poisson_cdf_c(actual, expected);
      if(i == hist->GetXaxis()->FindBin(0.5))
        printf("Events in 0-1s: %d, %.1f expected. Prob of this or more: %.2g\n",
               actual, expected, prob);
      if(prob < minprob){
        minprob = prob;
        minprob_bin = i;
      }
    }

    const double expected = flat * histlive->GetBinContent(minprob_bin);
    const int actual = (int)hist->GetBinContent(minprob_bin);

    // This isn't quite right, because it could be more than 1.  Really what
    // I mean is not to boost the probability, but to reduce the threshold for
    // being surprised.  But this is easier to write.
    const double corrprob_almostanywhere =
      ROOT::Math::poisson_cdf_c(actual, expected) * (hist->GetNbinsX() - 1);

    printf("Biggest excess: %d - %ds: %d obs, %.1f expected. "
           "Prob of this or more somewhere: %.2g\n",
           (int)hist->GetBinLowEdge(minprob_bin),
           (int)hist->GetBinLowEdge(minprob_bin+1),
           actual, expected, corrprob_almostanywhere);

    // We'll consider the first N seconds after the event to be special
    const int nearlybins = 10;
    if(minprob_bin >= hist->GetXaxis()->FindBin(0.5) &&
       minprob_bin <  hist->GetXaxis()->FindBin(0.5)+nearlybins){

      const double corrprob_early =
        ROOT::Math::poisson_cdf_c(actual, expected) * nearlybins;
      printf("In fact, this was in 0-%ds. Prob of this: %.2g\n",
             nearlybins, corrprob_early);
    }
  }
}

void novapreliminary()
{
  static TLatex * t = new TLatex(0, 0, "NOvA Preliminary");
  t->SetTextColor(kBlue);
  t->SetTextSize(textsize*textratiofull);
  t->SetTextFont(42);
  t->SetNDC();
  t->SetTextAlign(32);
  t->SetX(1-rightmargin-0.005);
  t->SetY(0.978);
  t->Draw();
}

void process_rebin(TH1D *hist, TH1D * histlive, const char * title,
                   const int rebin)
{
  TH1D * rebinned = (TH1D*)hist->Rebin(rebin, "rebinned");
  TH1D * rebinnedlive = histlive == NULL? NULL:
    (TH1D*)histlive->Rebin(rebin, "rebinnedlive");

  stylehist(rebinned, 0);

  const bool preliminary = true;

  // Some trickery (static) here needed to avoid getting the first label stuck
  // on all the output PDFs.
  static TLatex * ltitle = new TLatex;
  ltitle->SetText(0.5 + leftmargin/2 - rightmargin/2
                  - preliminary * 0.1, 0.978, title);
  ltitle->SetTextSize(textsize*textratiofull);
  ltitle->SetTextFont(42);
  ltitle->SetTextAlign(22);
  ltitle->SetNDC();
  ltitle->Draw();

  novapreliminary();

  if(dividebylivetime){
    TH1D * divided = divide_livetime(rebinned, rebinnedlive);

    stylehist(divided, 1);
    stylehist(rebinned, 2);
    stylehist(rebinnedlive, 3);

    divided ->GetYaxis()->SetTitle("Events/s");
    rebinned->GetYaxis()->SetTitle("Raw");

    const double ytitleoff = 1.1 / sqrt(textsize/0.05);
    divided     ->GetYaxis()->SetTitleOffset(ytitleoff);
    rebinned    ->GetYaxis()->SetTitleOffset(ytitleoff/textratiomid);
    rebinnedlive->GetYaxis()->SetTitleOffset(ytitleoff/textratiobot);

    if(rebinned->Integral() > 0){
      divided->GetYaxis()->SetRangeUser(0.001, divided->GetMaximum() +
        std::max(divided->GetMaximum()/10, divided->GetBinError(divided->GetMaximumBin())*1.5));

      rebinned->GetYaxis()->SetRangeUser(0.001, rebinned->GetMaximum() +
        std::max(rebinned->GetMaximum()/10, rebinned->GetBinError(rebinned->GetMaximumBin())*1.5));
    }

    rebinnedlive->GetYaxis()->SetTitle("Livetime (%)");
    rebinnedlive->Scale(100./rebin);

    toppad->cd();
    divided->Draw("e");

    midpad->cd();
    rebinned->Draw("e");

    botpad->cd();
    rebinnedlive->Draw("hist");

    bumphunt(hist, histlive, rebin == 1);
  }
  else{
    rebinned->Draw("e");
  }

  print2(Form("%s-%ds", hist->GetName(), rebin));
}

void process(TH1D * hist, TH1D * histlive, const char * const title)
{
  process_rebin(hist, histlive, title,  1);
  process_rebin(hist, histlive, title, 10);
}

void ligopass2(const char * const infilename, const char * outbase_,
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

  #define GETORDONT(hname, title) \
    TH1D * hname = getordont(#hname, ligodir); \
    TH1D * hname##live = getordont(#hname"live", ligodir); \
    if(hname != NULL) process(hname, hname##live, title)

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
