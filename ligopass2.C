#include "TGaxis.h"
#include "TError.h"
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

struct popts{
  const char * name;

  // What order polynominal to fit the background time series to
  unsigned int polyorder;

  // Whether to run statistical tests on deviations from background
  bool stattest;

  popts(const char * const name_, const unsigned int polyorder_, const bool stattest_)
  {
    name = name_;
    polyorder = polyorder_;
    stattest = stattest_;
  }
};

const double textsize = 0.0666;

bool dividebylivetime = false;
bool longreadout = false;
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

  x->SetTitle("Time since event (s)");
  x->SetNdivisions(508);
  y->SetNdivisions(can == 0?510: can == 1? 508: 504);

  if(effectivecount(h) < h->GetNbinsX()){
    h->SetMarkerStyle(kOpenCircle);
  }

  if(h->GetEntries() < 10) y->SetNdivisions(h->GetEntries());
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
  c->SetCanvasSize(100*scale + (!dividebylivetime)*50*scale, 100*scale);

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
  else{
   c->SetBottomMargin(0.14);
   c->SetLeftMargin(0.14);
   c->SetRightMargin(0.02);
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
static void fcn(int & np, __attribute__((unused)) double * gin, double & like,
                double *par, __attribute__((unused)) int flag)
{ 
  like = 0;

  for(int i = 1; i <= fithist->GetNbinsX(); i++){
    double model = 0;
    for(int p = 0; p < np; p++){
      // Avoid warning about the function value not depending on the parameters
      // when fitting to a constant.
      const double pr = np == 1 ? fabs(par[p]): par[p];

      model += pr * pow(fithist->GetBinCenter(i), p);
    }
    model *= fithistlive->GetBinContent(i);

    const double data = fithist->GetBinContent(i);

    if(model > 0){
      like += model - data;
      if(data > 0) like += data * log(data/model);
    }
  }

  like *= 2; // Make chi2-like.  Alternatively could SET UP or whatever
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

void bumphunt(TH1D * hist, TH1D * histlive, const bool verbose, const unsigned int polyorder)
{
  TMinuit mn(polyorder+1);
  mn.SetPrintLevel(-1);
  mn.SetFCN(fcn);
  int ierr;
  mn.mnparm(0, "flat", 1, 10, 0, 0, ierr);
  for(unsigned int i = 1; i <= polyorder; i++)
    mn.mnparm(i, Form("p%d", i), 0, 1, 0, 0, ierr);

  fithist = hist;
  fithistlive = histlive;
  mn.Command("SET STRATEGY 2");
  mn.Command("MINIMIZE");

  toppad->cd();

  string fs = "[0]";
  for(unsigned int i = 1; i <= polyorder; i++)
    fs += Form("+[%d]*pow(x, %d)", polyorder, polyorder);

  TF1 * poly =
    new TF1("poly", fs.c_str(), hist->GetBinLowEdge(1),
                                hist->GetBinLowEdge(hist->GetNbinsX()+1));

  for(unsigned int i = 0; i <= polyorder; i++)
    poly->SetParameter(i, getpar(mn, i));
  poly->SetLineColor(kRed);
  poly->SetNpx(500);
  poly->SetLineWidth(2);
  poly->Draw("same");

  if(verbose){
    double minprob = 1;
    double minprob_bin = 0;
    for(int i = 1; i <= hist->GetNbinsX(); i++){
      if(histlive->GetBinContent(i) == 0) continue;

      const double expected = poly->Eval(hist->GetBinCenter(i)) * histlive->GetBinContent(i);
      const int actual = (int)hist->GetBinContent(i);

      const double prob = ROOT::Math::poisson_cdf_c(actual, expected);
      if(i == hist->GetXaxis()->FindBin(0.5))
        printf("Events in 0-1s: %d, %.3f expected. Prob of this or more: %.2g\n",
               actual, expected, prob);
      if(prob < minprob){
        minprob = prob;
        minprob_bin = i;
      }
    }

    const double expected = poly->Eval(hist->GetBinCenter(minprob_bin))
      * histlive->GetBinContent(minprob_bin);
    const int actual = (int)hist->GetBinContent(minprob_bin);

    // This isn't quite right, because it could be more than 1.  Really what
    // I mean is not to boost the probability, but to reduce the threshold for
    // being surprised.  But this is easier to write.
    const double corrprob_almostanywhere =
      ROOT::Math::poisson_cdf_c(actual, expected) * (hist->GetNbinsX() - 1);

    printf("Biggest excess: {%d, %d}s: %d obs, %.1f expected. "
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

void randomtime()
{
  static TLatex * t = new TLatex(0, 0, "Around 8:30 SNEWS Test");
  t->SetTextColorAlpha(kGray, 0.5);
  t->SetTextSize(textsize*textratiofull*2.5);
  t->SetTextFont(42);
  t->SetNDC();
  t->SetTextAlign(22);
  t->SetX(0.5);
  t->SetY(0.5);
  t->SetTextAngle(55);
  t->Draw();
}

void novapreliminary()
{
  static TLatex * t = new TLatex(0, 0, "NOvA Preliminary");
  t->SetTextColor(kBlue);
  t->SetTextSize(dividebylivetime?textsize*textratiofull:textsize);
  t->SetTextFont(42);
  t->SetNDC();
  t->SetTextAlign(32);
  t->SetX(1-rightmargin-0.005);
  t->SetY(0.978 - (!dividebylivetime)*0.02);
  t->Draw();
}

// polyorder: order of the polynominal to fit for the no-signal hypothesis.
//            Naively, zero.  For anything with pointing, at least 1.
void process_rebin(TH1D *hist, TH1D * histlive,
                   const int rebin, const popts opts)
{
  TH1D * rebinned = (TH1D*)hist->Rebin(rebin, "rebinned");
  TH1D * rebinnedlive = histlive == NULL? NULL:
    (TH1D*)histlive->Rebin(rebin, "rebinnedlive");

  stylehist(rebinned, 0);

  const bool preliminary = true;

  // Some trickery (static) here needed to avoid getting the first label stuck
  // on all the output PDFs.
  static TLatex * ltitle = new TLatex;
  ltitle->SetText(
    0.53 + leftmargin/2 - rightmargin/2 - preliminary * 0.2,
    0.978 - (!dividebylivetime)*0.02, opts.name);
  ltitle->SetTextSize(dividebylivetime?textsize*textratiofull:textsize);
  ltitle->SetTextFont(42);
  ltitle->SetTextAlign(22);
  ltitle->SetNDC();
  ltitle->Draw();

  if(preliminary) novapreliminary();
  randomtime();

  double mint = 0, maxt = 0;

  if(longreadout && rebin == 1){
    mint = -10, maxt = 45;
    rebinned->GetXaxis()->SetRangeUser(mint, maxt);
    if(rebinnedlive != NULL) rebinnedlive->GetXaxis()->SetRangeUser(mint, maxt);
  }


  if(dividebylivetime){
    TH1D * divided = divide_livetime(rebinned, rebinnedlive);
    if(rebin == 1) divided->GetXaxis()->SetRangeUser(mint, maxt);

    stylehist(divided, 1);
    stylehist(rebinned, 2);
    stylehist(rebinnedlive, 3);

    divided ->GetYaxis()->SetTitle("Events/s");
    rebinned->GetYaxis()->SetTitle("Raw");

    const double ytitleoff = 1.1 / sqrt(textsize/0.05);
    divided     ->GetYaxis()->SetTitleOffset(ytitleoff);
    rebinned    ->GetYaxis()->SetTitleOffset(ytitleoff/textratiomid);
    rebinnedlive->GetYaxis()->SetTitleOffset(ytitleoff/textratiobot);

    rebinnedlive->GetYaxis()->SetTitle("Livetime (%)");
    rebinnedlive->Scale(100./rebin);
    rebinnedlive->GetYaxis()->SetRangeUser(0, rebinnedlive->GetMaximum()*1.2);

    toppad->cd();
    divided->Draw("e");
    if(toppad->GetUymin() == 0) divided->GetYaxis()->ChangeLabel(1, -1, 0);

    midpad->cd();
    rebinned->Draw("e");
    if(midpad->GetUymin() == 0) rebinned->GetYaxis()->ChangeLabel(1, -1, 0);

    botpad->cd();
    rebinnedlive->Draw("hist");

    if(hist->Integral() > 0 && opts.stattest)
      bumphunt(hist, histlive, rebin == 1,
               std::min((unsigned int)(hist->Integral()), opts.polyorder));
  }
  else{
    rebinned->GetYaxis()->SetTitleOffset(
      (rebinned->GetEntries() < 10? 0.8: 1.1) / sqrt(textsize/0.05));
    rebinned->GetYaxis()->SetTitle("Events/s");
    rebinned->Draw("e");
    if(preliminary) novapreliminary();
    randomtime();
    ltitle->Draw();
  }

  print2(Form("%s-%ds", hist->GetName(), rebin));
}

void process(TH1D * hist, TH1D * histlive, const popts opts)
{
  printf("\n%s:\n", opts.name);
  process_rebin(hist, histlive,  1, opts);
  process_rebin(hist, histlive, 10, opts);
}

void ligopass2(const char * const infilename, const char * outbase_,
               const bool dividebylivetime_, const bool longreadout_)
{
  outbase = outbase_;
  dividebylivetime = dividebylivetime_;
  longreadout = longreadout_;

  // Don't print about making PDFs
  gErrorIgnoreLevel = kWarning;

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

  #define DOIT(hname, opts) \
    TH1D * hname = getordont(#hname, ligodir); \
    TH1D * hname##live = getordont(#hname"live", ligodir); \
    if(hname != NULL && hname##live != NULL) process(hname, hname##live, opts)

  DOIT(rawtrigger,                   popts("Raw triggers", 0, true));
  DOIT(rawhits,                      popts("Raw hits", 0, false));
  DOIT(unslicedbighits,              popts("Unsliced big hits", 0, false));
  DOIT(unslice4ddhits,               popts("Hits in noise slice", 0, true));
  DOIT(tracks,                       popts("Slices with tracks", 0, true));
  DOIT(tracks_point_0,               popts("Slices w/ tracks, 1.3#circ", 1, true));
  DOIT(tracks_point_1,               popts("Slices w/ tracks, 16#circ", 1, true));
  DOIT(fullycontained_tracks,        popts("Contained tracks", 0, true));
  DOIT(fullycontained_tracks_point_0,popts("Contained tracks, 1.3#circ", 1, true));
  DOIT(fullycontained_tracks_point_1,popts("Contained tracks, 16#circ", 1, true));
  DOIT(contained_slices,             popts("Contained slices", 0, true));
  DOIT(halfcontained_tracks,         popts("Stopping tracks", 0, true));
  DOIT(halfcontained_tracks_point_0, popts("Stopping tracks, 1.3#circ", 1, true));
  DOIT(halfcontained_tracks_point_1, popts("Stopping tracks, 16#circ", 1, true));
  DOIT(unslicedhitpairs,             popts("Supernova-like events", 0, true));
  DOIT(upmu_tracks,                  popts("Upward going muons", 0, true));
  DOIT(upmu_tracks_point_0,          popts("Up-#mu, 1.3#circ", 1, true));
  DOIT(upmu_tracks_point_1,          popts("Up-#mu, 16#circ", 1, true));
  DOIT(energy_low_cut,               popts("High energy triggers", 0, true));
  DOIT(energy_high_cut,              popts("Very high energy triggers", 0, true));
  DOIT(energy_low_cut_pertime,       popts("High energy events", 0, true));
  DOIT(energy_high_cut_pertime,      popts("Very high energy events", 0, true));
}
