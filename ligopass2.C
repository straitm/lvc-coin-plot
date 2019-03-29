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
const double topmargin = 0.09 * textsize/0.07;

bool dividebylivetime = false;
bool longreadout = false;
const char * outbase = NULL;
const char * trigname = NULL;

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

  TH1D * lastlive = NULL;

TH1D * getordont(const char * const histname, TDirectory * const ligofile)
{
  TH1D * h = dynamic_cast<TH1D *>(ligofile->Get(histname));

  const bool livehist =
    !strcmp(histname + strlen(histname) - (sizeof("live") - 1), "live");

  if(livehist && h != NULL) lastlive = h;

  if(h == NULL){
    printf("\nNo histogram \"%s\" in your file\n", histname);
    if(livehist){
      if(lastlive != NULL){
        printf(
          "Warning: Substituting previous livetime.  As long as different\n"
          "selections don't somehow have different livetimes, this is fine.\n");
        return lastlive;
      }
      else{
        printf("***ERROR***: I have to skip this one!\n");
      }
    }
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

  if(longreadout || h->GetNbinsX() < 300) h->SetLineWidth(2);
  else                                    h->SetLineWidth(1);

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

  y->SetDecimals();

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

  gStyle->SetGridColor(kGray);
  gStyle->SetGridStyle(kSolid);

  if(dividebylivetime){
    toppad = new TPad("div", "div", 0, divheight2, 1,          1);
    midpad = new TPad("raw", "raw", 0, divheight1, 1, divheight2);
    botpad = new TPad("liv", "liv", 0,          0, 1, divheight1);
    stylepad(toppad);
    stylepad(midpad);
    stylepad(botpad);

    toppad->SetTopMargin(topmargin);
    toppad->SetBottomMargin(0);

    midpad->SetBottomMargin(0);
    midpad->SetTopMargin(0);

    botpad->SetTopMargin(0);
    botpad->SetBottomMargin(0.078/divheight1);

    toppad->SetGrid(1, 0);
    botpad->SetGrid(1, 0);
    midpad->SetGrid(1, 0);

    toppad->Draw();
    botpad->Draw();
    midpad->Draw();
  }
  else{
    c->SetBottomMargin(0.143);
    c->SetLeftMargin(0.125);
    c->SetRightMargin(0.02);
    c->SetTopMargin(0.08);
    c->SetTickx(1);
    c->SetTicky(1);
  }
}

TH1D * divide_livetime(TH1D * d, TH1D * live)
{
  TH1D * answer = (TH1D*)d->Clone(Form("%sdiv", d->GetName()));

  for(int i = 1; i <= live->GetNbinsX(); i++){
    const double bc = d->GetBinContent(i);
    const double be = d->GetBinError(i);
    const double livefrac = live->GetBinContent(i);
    if(livefrac == 0){
      answer->SetBinContent(i, 0);
      answer->SetBinError  (i, 0);
    }
    else{
      answer->SetBinContent(i, bc/livefrac);
      answer->SetBinError  (i, be/livefrac);
    }
  }
  return answer;
}

void printstart()
{
  USN.Print(Form("%s.pdf[", outbase));
}

void printend()
{
  USN.Print(Form("%s.pdf]", outbase));
}

void print2(const char * const name)
{
  //USN.Print(Form("%s-%s.pdf", outbase, name));
  USN.Print(Form("%s.pdf", outbase));
}


TH1D * fithist = NULL, * fithistlive = NULL;
static void fcn(int & np, __attribute__((unused)) double * gin, double & like,
                double *par, __attribute__((unused)) int flag)
{
  like = 0;

  for(int bin = 1; bin <= fithist->GetNbinsX(); bin++){
    double model = 0;
    for(int p = 0; p < np; p++){
      // Avoid warning about the function value not depending on the parameters
      // when fitting to a constant.
      const double pr = np == 1 ? fabs(par[p]): par[p];

      model += pr * pow(fithist->GetBinCenter(bin), p);
    }
    model *= fithistlive->GetBinContent(bin);

    const double data = fithist->GetBinContent(bin);

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

void stylefunction(TF1 * f)
{
  f->SetLineColor(kRed);
  f->SetNpx(500);
  f->SetLineWidth(2);
}

double prob_this_or_more(const int actual, const double expected)
{
  if(actual == 0) return 1;
  // poisson_cdf_c(n, mu) gives the probability that *more* than n events
  // are seen for a mean of mu.
  return ROOT::Math::poisson_cdf_c(actual-1, expected);
}

double lookelsewhere(const double localprob, const int trials)
{
  return 1 - pow(1-localprob, trials);
}

void sigmacheck(const double p)
{
  const unsigned int nlev = 2;
  const double lev[nlev] = { 5, 3 };

  for(unsigned int i = 0; i < nlev; i++)
    if(p < (1-ROOT::Math::gaussian_cdf(lev[i], 1))*2){
      printf("*** That's more than %.0f sigma ***\n", lev[i]);
      return;
    }
}

void fit(TMinuit & mn, TH1D * hist, TH1D * histlive, const unsigned int polyorder)
{
  fithist = hist;
  fithistlive = histlive;

  mn.Command("SET STRATEGY 2");
  for(unsigned int p = 0; p <= polyorder; p++) mn.Command(Form("FIX %d", p+1));

  for(unsigned int i = 0; i <= polyorder; i++){
    // Sloppily try to prevent fitting more than we can chew
    if(i > 0 && hist->Integral() < 20*i) continue;

    mn.Command(Form("REL %d", i+1));
    mn.Command("MINIMIZE");

    // Since we use fabs() in the FCN to avoid warnings, we have to fix it here...
    if((polyorder == 0 || i == 0) && getpar(mn, 0) < 0){
      mn.Command(Form("SET PAR 1 %f\n", fabs(getpar(mn, 0))));
      mn.Command("MINIMIZE");
    }
  }
}

TH1D * dummy_histlive(TH1D * hist, const int rebin)
{
  TH1D * dum = new TH1D("dum", "",
                        hist->GetNbinsX(),
                        hist->GetBinLowEdge(1),
                        hist->GetBinLowEdge(hist->GetNbinsX()+1));
  for(int i = 1; i <= dum->GetNbinsX(); i++)
    dum->SetBinContent(i, 1./rebin);
  return dum;
}


void setup_mn(TMinuit & mn, const unsigned int polyorder)
{
  mn.SetPrintLevel(polyorder < 2?-1:0);
  mn.SetFCN(fcn);
  int ierr;
  mn.mnparm(0, "flat", 1, 1, 0, 0, ierr);
  for(unsigned int i = 1; i <= polyorder; i++)
    mn.mnparm(i, Form("p%d", i), 0, 0.001, 0, 0, ierr);
}


void bumphunt(TH1D * hist, TH1D * histlive, const int rebin,
              const unsigned int polyorder)
{
  if(histlive == NULL || !dividebylivetime) histlive = dummy_histlive(hist, rebin);

  const bool verbose = rebin == 1;
  TMinuit mn(polyorder+1);
  setup_mn(mn, polyorder);

  fit(mn, hist, histlive, polyorder);

  if(dividebylivetime) toppad->cd();

  std::string fs = "[0]";
  for(unsigned int i = 1; i <= polyorder; i++)
    fs += Form("+[%d]*pow(x, %d)", polyorder, polyorder);

  TF1 * poly =
    new TF1("poly", fs.c_str(), hist->GetBinLowEdge(1),
                                hist->GetBinLowEdge(hist->GetNbinsX()+1));

  for(unsigned int i = 0; i <= polyorder; i++)
    poly->SetParameter(i, getpar(mn, i));
  stylefunction(poly);
  poly->Draw("same");

  if(dividebylivetime){
    TH1D * fitraw = new TH1D(Form("%sfit%d", hist->GetName(), rebin), "",
                             hist->GetNbinsX(), hist->GetBinLowEdge(1),
                             hist->GetBinLowEdge(hist->GetNbinsX()+1));
    fitraw->SetLineWidth(2);
    fitraw->SetLineColor(kRed);
    fitraw->SetMarkerColor(kRed);

    for(int i = 1; i <= fitraw->GetNbinsX(); i += rebin){
      double sum = 0;
      for(int j = 0; j < rebin; j++)
        sum += poly->Eval(fitraw->GetBinCenter(i+j)) * histlive->GetBinContent(i+j);
      for(int j = 0; j < rebin; j++)
        fitraw->SetBinContent(i+j, sum);
    }

    midpad->cd();
    fitraw->Draw(longreadout?"samehist":"samehist][");
  }

  if(verbose){
    double minprob = 1;
    double minprob_bin = 0;
    for(int i = 1; i <= hist->GetNbinsX(); i++){
      if(histlive-> GetBinContent(i) == 0) continue;

      const double expected = poly->Eval(hist->GetBinCenter(i)) * histlive->GetBinContent(i);
      const int actual = (int)hist->GetBinContent(i);

      const double prob = prob_this_or_more(actual, expected);

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

    const double localprob = prob_this_or_more(actual, expected);
    const double corrprob_almostanywhere =
      lookelsewhere(localprob, hist->GetNbinsX());

    printf("Biggest excess: {%d, %d}s: %d obs, %.3f expected. "
           "Prob of this or more somewhere: %.2g\n",
           (int)hist->GetBinLowEdge(minprob_bin),
           (int)hist->GetBinLowEdge(minprob_bin+1),
           actual, expected, corrprob_almostanywhere);
    sigmacheck(corrprob_almostanywhere);

    // We'll consider the first N seconds after the event to be special
    const int nearlybins = 10;
    if(minprob_bin >= hist->GetXaxis()->FindBin(0.5) &&
       minprob_bin <  hist->GetXaxis()->FindBin(0.5)+nearlybins){

      const double corrprob_early =
        lookelsewhere(prob_this_or_more(actual, expected), nearlybins);
      printf("This was in the first %ds. Prob of this: %.2g\n",
             nearlybins, corrprob_early);
      sigmacheck(corrprob_early);
    }
  }
}

void randomtime()
{
  static TLatex * t = new TLatex(0, 0, "8:30 SNEWS Test");
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
  t->SetTextColorAlpha(kBlue, 0.5);
  t->SetTextSize(dividebylivetime?textsize*textratiofull:textsize);
  t->SetTextFont(42);
  t->SetNDC();
  t->SetTextAlign(32);
  t->SetX(1 - rightmargin - 0.02);
  t->SetY(1 -  topmargin  + 0.01 - (!dividebylivetime)*0.06);
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
    0.16 + leftmargin/2,
    0.978 - (!dividebylivetime)*0.02, Form("%s: %s", trigname, opts.name));
  ltitle->SetTextSize(dividebylivetime?textsize*textratiofull:textsize);
  ltitle->SetTextFont(42);
  ltitle->SetTextAlign(12);
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
    if(rebin == 1){
       divided->GetXaxis()->SetRangeUser(mint, maxt);

      double miny = divided->GetMaximum();
      double maxy = divided->GetMaximum();
      double extramax = 0, extramin = 0;

      for(int i = 1; i <= divided->GetNbinsX(); i++){
        const double vmin = divided->GetBinContent(i) - divided->GetBinError(i);
        const double vmax = divided->GetBinContent(i) + divided->GetBinError(i);
        const double e = divided->GetBinError(i);
        if(vmax > maxy){
           maxy = vmax;
           extramax = e * 0.5;
        }
        if(vmin > 0 && vmin < miny){
           miny = vmin;
           extramin = e * 0.5;
        }
      }

      if(maxy < 2*miny)
        divided->GetYaxis()->SetRangeUser(miny-extramin, maxy+extramax);
    }

    stylehist(divided, 1);
    stylehist(rebinned, 2);
    stylehist(rebinnedlive, 3);

    divided ->GetYaxis()->SetTitle("Events/s");
    rebinned->GetYaxis()->SetTitle("Raw");

    const double ytitleoff = 1.2 / sqrt(textsize/0.05);
    divided ->GetYaxis()->SetTitleOffset(ytitleoff);
    rebinned->GetYaxis()->SetTitleOffset(ytitleoff/textratiomid);

    // So it aligns to bottom of letters, not bottom of parentheses...
    rebinnedlive->GetYaxis()->SetTitleOffset(ytitleoff/textratiobot * 0.955);

    rebinnedlive->GetYaxis()->SetTitle("Livetime (%)");
    rebinnedlive->Scale(100./rebin);
    rebinnedlive->GetYaxis()->SetRangeUser(0, rebinnedlive->GetMaximum()*1.2);

    toppad->cd();
    divided->Draw("e");
    toppad->Update(); // Necessary to get Uymin correctly!
    if(toppad->GetUymin() == 0) divided->GetYaxis()->ChangeLabel(1, -1, 0);

    midpad->cd();
    rebinned->Draw("e");
    midpad->Update(); // Necessary to get Uymin correctly!
    if(midpad->GetUymin() == 0) rebinned->GetYaxis()->ChangeLabel(1, -1, 0);

    botpad->cd();
    rebinnedlive->Draw(longreadout?"hist":"hist][");
  }
  else{
    rebinned->GetYaxis()->SetTitleOffset(
      (rebinned->GetEntries() < 10? 0.8: 1.1) / sqrt(textsize/0.05));
    if(rebin == 1)
      rebinned->GetYaxis()->SetTitle("Events/s");
    else
      rebinned->GetYaxis()->SetTitle(Form("Events/%d#kern[-0.5]{ }s", rebin));

    rebinned->Draw("e");
    if(preliminary) novapreliminary();
    randomtime();
    ltitle->Draw();
  }

  if(hist->Integral() == 0){
    if(rebin == 1) printf("No events, so there's no excess\n");
  }
  else if(!opts.stattest){
    if(rebin == 1) printf("This distribution is known to be non-Poissonian\n");
  }
  else{
    bumphunt(hist, histlive, rebin,
             std::min((unsigned int)(hist->Integral()), opts.polyorder));
  }

  print2(Form("%s-%ds", hist->GetName(), rebin));
}

void process(TH1D * hist, TH1D * histlive, const popts opts)
{
  printf("\n%s: ", opts.name);
  process_rebin(hist, histlive,  1, opts);
  if(!longreadout) process_rebin(hist, histlive, 10, opts);
}

void ligopass2(const char * const infilename, const char * trigname_,
               const char * outbase_, const bool dividebylivetime_,
               const bool longreadout_)
{
  outbase = outbase_;
  dividebylivetime = dividebylivetime_;
  longreadout = longreadout_;
  trigname = trigname_;

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
  printstart();

  #define DOIT(hname, opts) \
    TH1D * hname = getordont(#hname, ligodir); \
    TH1D * hname##live = hname == NULL? NULL: getordont(#hname"live", ligodir); \
    if(hname != NULL && (!dividebylivetime || hname##live != NULL)) \
      process(hname, hname##live, opts)

  DOIT(rawhits,                      popts("Raw hits",                  0, trigname[0] == 'N'));
  DOIT(unslice4ddhits,               popts("Unsliced hits",             0, trigname[0] == 'N'));
  DOIT(unslicedbighits,              popts("Unsliced big hits",         0, trigname[0] == 'N'));
  DOIT(unslicedhitpairs,             popts("Supernova-like events",     0, true));

  // Skip these for the ND long readout, since they are redundant with the 100%
  // efficienct ND ddactivity1 trigger.
  if( strcmp(trigname, "ND long readout") ){
    DOIT(contained_slices,             popts("Contained slices",          0, true));
    DOIT(tracks,                       popts("Slices with tracks",        0, true));
    DOIT(tracks_point_1,               popts("Slices w/ tracks, 16#circ", 1, true));
    DOIT(tracks_point_0,               popts("Slices w/ tracks, 1.3#circ",1, true));
    DOIT(halfcontained_tracks,         popts("Stopping tracks",           0, true));
    DOIT(halfcontained_tracks_point_1, popts("Stopping tracks, 16#circ",  1, true));
    DOIT(halfcontained_tracks_point_0, popts("Stopping tracks, 1.3#circ", 1, true));
    DOIT(fullycontained_tracks,        popts("Contained tracks",          0, true));
    DOIT(fullycontained_tracks_point_1,popts("Contained tracks, 16#circ", 1, true));
    DOIT(fullycontained_tracks_point_0,popts("Contained tracks, 1.3#circ",1, true));
    DOIT(upmu_tracks,                  popts("Upward going muons",        0, true));
    DOIT(upmu_tracks_point_1,          popts("Up-#mu, 16#circ",           1, true));
    DOIT(upmu_tracks_point_0,          popts("Up-#mu, 1.3#circ",          1, true));

    DOIT(rawtrigger,                   popts("Raw triggers",              0, true));
    DOIT(energy_low_cut,               popts("> 5M ADC total",            0, true));
    DOIT(energy_high_cut,              popts("> 50M ADC total",           0, true));
    DOIT(energy_low_cut_pertime,       popts("> 5M ADC per 50#mus",       0, true));
    DOIT(energy_high_cut_pertime,      popts("> 50M ADC per 50#mus",      0, true));
  }

  printend();
}
