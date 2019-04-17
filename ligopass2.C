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

//#define DEBUG

// The total number of histograms searched for
const unsigned int NHIST = 52;

// If true, allow using external expectations for determining excess
// If false, always run fit internally (since this job will determine what
// the background is for other jobs)
bool gextexp = true;

enum stream_t {
  fardet_t02, neardet_ddactivity1, fardet_ddenergy, neardet_long, fardet_long,
  UNDEFINED_STREAM
};

struct popts{
  const char * name;

  // What order polynominal to fit the background time series to
  unsigned int polyorder;

  // Whether to run statistical tests on deviations from background
  bool stattest;

  // Background per second if this histogram is too low-stats to determine
  // it in one window.  Otherwise zero.
  double extexp;

  // Slope of background if the background is time varying
  double extexp1;

  popts(const char * const name_, const unsigned int polyorder_,
        const bool stattest_, const double extexp_ = 0, const double extexp1_ = 0)
  {
    name = name_;
    polyorder = polyorder_;
    stattest = stattest_;

    // If we're evaluating the background, do not set the background
    if(gextexp){
      extexp = extexp_;
      extexp1 = extexp1_;
    }
    else{
      extexp = extexp1 = 0;
    }
  }
};

const double textsize = 0.07;
const double topmargin = 0.09 * textsize/0.07;

bool dividebylivetime = false;
bool longreadout = false;

// Number of readout windows that the input data represents.  It is 1 for
// looking at signal samples and more than when for background samples.
int nwindows = 1;

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
const double textratiofull = 1-divheight2;

TH1D * lastlive = NULL;

TH1D * getordont(const char * const histname, TDirectory * const ligodir)
{
  TH1D * h = dynamic_cast<TH1D *>(ligodir->Get(histname));

  const bool livehist =
    !strcmp(histname + strlen(histname) - (sizeof("live") - 1), "live");

  if(livehist && h != NULL) lastlive = h;

  if(h == NULL){
#ifdef DEBUG
    printf("\nNo histogram \"%s\" in your file\n", histname);
#endif
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

void styledrawellipse(TEllipse * a)
{
  a->SetFillStyle(0);
  a->SetLineColor(kRed);
  a->Draw();
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

  if(longreadout || h->Integral() < 100){
    h->SetMarkerStyle(kOpenCircle);
    h->SetMarkerSize(dividebylivetime?0.7:1.4);
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
  c->SetCanvasSize((100 + (!dividebylivetime)*50)*scale, 100*scale);

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
    c->SetLeftMargin(0.175);
    c->SetRightMargin(0.02);
    c->SetTopMargin(0.085);
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

  printf("------------- Finished processing -------------\n");
}

void print2(const char * const name)
{
  USN.Print(Form("%s-%s.pdf", outbase, name));
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

void stylefunction(TF1 * f, const bool externalexpectation)
{
  f->SetLineColor(externalexpectation?kBlue:kRed);
  f->SetLineStyle(externalexpectation?9:1);
  f->SetLineWidth(externalexpectation?4:2);
  f->SetNpx(500);
}

double prob_this_or_more(const int actual, const double expected)
{
  if(actual == 0) return 1;
  // poisson_cdf_c(n, mu) gives the probability that *more* than n events
  // are seen for a mean of mu.
  return ROOT::Math::poisson_cdf_c(actual-1, expected);
}

// Evaluate look-elsewhere
double lookelse(const double localprob, const int trials)
{
  return 1 - pow(1-localprob, trials);
}

// Checks whether the probability 'p' is exciting, prints a message if
// it is, and returns whether it is.
bool sigmacheck(const double p)
{
  const unsigned int nlev = 2;
  const double lev[nlev] = { 5, 3 };

  for(unsigned int i = 0; i < nlev; i++)
    if(p < (1-ROOT::Math::gaussian_cdf(lev[i], 1))*2){
      printf("*******************************************\n"
             "*******************************************\n"
             "*******  That's more than %.0f sigma   *******\n"
             "*******************************************\n"
             "*******************************************\n", lev[i]);
      return true;
    }
  return false;
}

void fit(TMinuit & mn, TH1D * hist, TH1D * histlive, const unsigned int polyorder)
{
  fithist = hist;
  fithistlive = histlive;

  mn.Command("SET STRATEGY 2");
  for(unsigned int p = 0; p <= polyorder; p++) mn.Command(Form("FIX %d", p+1));

  for(unsigned int i = 0; i <= polyorder; i++){
    // Sloppily try to prevent fitting more than we can chew
    if(i > 0 && hist->Integral() < 50*i) continue;

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
  if(gROOT->FindObject("dum")) delete gROOT->FindObject("dum");
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


void stylefitraw(TH1D * fitraw)
{
  fitraw->SetLineWidth(2);
  fitraw->SetLineColor(kRed);
  fitraw->SetMarkerColor(kRed);
}

// Look for an excess anywhere in the histogram.  If 'extexp'
// is positive, use it as the background per second.  Otherwise,
// fit the histogram using a polnominal of order 'polyorder' and
// use that as the background.  Return the bin with the biggest
// excess, or -1 if none has an interesting excess
int bumphunt(TH1D * hist, TH1D * histlive, const int rebin,
             const double extexp, const double extexp1,
             unsigned int polyorder)
{
  if(histlive == NULL || !dividebylivetime)
    histlive = dummy_histlive(hist, rebin);

  const bool verbose = rebin == 1;
  TMinuit mn(polyorder+1);
  setup_mn(mn, polyorder);

  if(extexp > 0){
    polyorder = 1; // Just assume 1 (enforced by definition of popts) and
                   // give the draw function enough paramters to do its work
  }
  else{
    fit(mn, hist, histlive, polyorder);
  }

  if(dividebylivetime) toppad->cd();

  std::string fs = "[0]";
  for(unsigned int i = 1; i <= polyorder; i++)
    fs += Form("+[%d]*pow(x, %d)", polyorder, polyorder);

  TF1 * poly = new TF1("poly", fs.c_str(), hist->GetBinLowEdge(1),
                       hist->GetBinLowEdge(hist->GetNbinsX()+1));

  if(extexp > 0){
    poly->SetParameter(0, rebin*extexp);
    poly->SetParameter(1, rebin*extexp1);
  }
  else{
    for(unsigned int i = 0; i <= polyorder; i++)
      poly->SetParameter(i, getpar(mn, i));
  }
  stylefunction(poly, extexp > 0);
  poly->Draw("same");

  if(dividebylivetime){
    TH1D * fitraw = new TH1D(Form("%sfit%d", hist->GetName(), rebin), "",
                             hist->GetNbinsX(), hist->GetBinLowEdge(1),
                             hist->GetBinLowEdge(hist->GetNbinsX()+1));
    stylefitraw(fitraw);

    for(int i = 1; i <= fitraw->GetNbinsX(); i += rebin){
      double sum = 0;
      for(int j = 0; j < rebin; j++)
        sum += poly->Eval(fitraw->GetBinCenter(i+j))
               * histlive->GetBinContent(i+j);
      for(int j = 0; j < rebin; j++)
        fitraw->SetBinContent(i+j, sum);
    }

    midpad->cd();
    fitraw->Draw(longreadout?"samehist":"samehist][");
  }

  if(!verbose) return -1;

  double minprob = 1;
  double minprob_bin = 0;
  for(int i = 1; i <= hist->GetNbinsX(); i++){
    if(histlive->GetBinContent(i) == 0) continue;

    const double expected = poly->Eval(hist->GetBinCenter(i))
      * histlive->GetBinContent(i);
    const int actual = (int)hist->GetBinContent(i);

    const double prob = prob_this_or_more(actual, expected);

    if(i == hist->GetXaxis()->FindBin(0.5))
      printf("Events 0-1s: %d, %.3f exp. P: %.2g\n", actual, expected, prob);
    if(prob < minprob){
      minprob = prob;
      minprob_bin = i;
    }
  }

  const double expected = poly->Eval(hist->GetBinCenter(minprob_bin))
    * histlive->GetBinContent(minprob_bin);
  const int actual = (int)hist->GetBinContent(minprob_bin);

  const double localprob = prob_this_or_more(actual, expected);
  const double histprob = lookelse(localprob, hist->GetNbinsX());
  const double globprob = lookelse(localprob, NHIST*hist->GetNbinsX());

  printf("Highest: {%d, %d}s: %d, %.3f exp. P hist: %.2g (%.2g global)\n",
         (int)hist->GetBinLowEdge(minprob_bin),
         (int)hist->GetBinLowEdge(minprob_bin+1),
         actual, expected, histprob, globprob);
  int ret = -1;
  if(sigmacheck(globprob)) ret = minprob_bin;

  if(extexp == 0){
    if(polyorder == 0){
      printf("Number per second: %.5g\n", hist->Integral()/
        histlive->Integral()/histlive->GetBinWidth(1)
          /(dividebylivetime?1:nwindows));
    }
    else if(polyorder == 1){
      printf("Fit: %.5f + %.5ft\n",
             poly->GetParameter(0)/(dividebylivetime?1:nwindows),
             poly->GetParameter(1)/(dividebylivetime?1:nwindows));
    }
  }

  // We'll consider the first N seconds after the event to be special
  const int specialbins = 10;
  const int specialbin1 = hist->GetXaxis()->FindBin(0.5);
  if(minprob_bin >= specialbin1 && minprob_bin < specialbin1+specialbins){
    const double slocalprob = prob_this_or_more(actual, expected);
    const double shistprob = lookelse(localprob, specialbins);
    const double sglobprob     = lookelse(localprob, NHIST*specialbins);
    printf("This was in first %ds. P: %.2g (%.2g global)\n",
           specialbins, shistprob, sglobprob);
    sigmacheck(sglobprob);
  }

  return ret;
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

const double np_longreadoutmin = 10.0;
const double np_notlongreadoutmax = -5.0;

void draw_nonpoisson_f(const double n, const int width, const int color)
{
  static int c = 0;
  TF1 * f = new TF1(Form("npf%d", c++), Form("%f", n),
    longreadout?np_longreadoutmin:-1000,
    longreadout?1000:np_notlongreadoutmax);
  f->SetLineWidth(width);
  f->SetLineColor(color);
  f->SetNpx(200);
  f->Draw("same");

  TF1 * fdot = new TF1(Form("npfdot%d", c++), Form("%f", n),
    longreadout?-1000:np_notlongreadoutmax,
    longreadout?np_longreadoutmin:1000);
  fdot->SetLineStyle(7);
  fdot->SetLineWidth(width);
  fdot->SetLineColor(color);
  fdot->Draw("same");
}

// Search for excursions from normal behavior for a series that isn't
// Poissonian, like the total number of hits, which obviously has a lot
// of correlated activity.
void bumphunt_nonpoisson(TH1D * h)
{
  TH1D side("side", "", 100, h->GetMinimum(), h->GetMaximum()*1.001);

  for(int i = 1; i <= h->GetNbinsX(); i++){
    if(longreadout && h->GetBinCenter(i) < np_longreadoutmin) continue;
    if(!longreadout && h->GetBinCenter(i) > np_notlongreadoutmax) continue;

    side.Fill(h->GetBinContent(i));
  }

  if(side.Integral() == 0){
    printf("No data in baseline region\n");
    return;
  }

  side.Fit("gaus", "l0q");
  if(side.GetFunction("gaus") == NULL){
    fprintf(stderr, "Fit failed in bumphunt_nonpoisson on %s\n", h->GetName());
    return;
  }
  const double mean = side.GetFunction("gaus")->GetParameter(1);
  const double sigma = side.GetFunction("gaus")->GetParameter(2);

  toppad->cd();

  draw_nonpoisson_f(mean, 3, kRed);
  draw_nonpoisson_f(mean+sigma, 2, kRed);
  draw_nonpoisson_f(mean-sigma, 2, kRed);
  draw_nonpoisson_f(mean+2*sigma, 2, kRed+1);
  draw_nonpoisson_f(mean-2*sigma, 2, kRed+1);

  for(int i = 1; i <= h->GetNbinsX(); i++){
    if(longreadout && h->GetBinCenter(i) > np_longreadoutmin) continue;
    if(!longreadout && h->GetBinCenter(i) < np_notlongreadoutmax) continue;

    // XXX janky, but not as janky as just taking the central value.
    // In most cases, the error is quite small, but this protects against
    // bins with low livetime.
    const double content = h->GetBinContent(i) - h->GetBinError(i);

    const double localsigma = (content - mean)/sigma;
    const int trials = longreadout?40:500; // XXX fragile
    const double localprob = ROOT::Math::gaussian_cdf_c(localsigma);
    const double histprob = lookelse(localprob, trials);
    const double globprob = lookelse(localprob, NHIST*trials);

    const double histsigma = ROOT::Math::gaussian_quantile_c(histprob, 1);
    const double globsigma = ROOT::Math::gaussian_quantile_c(globprob, 1);

    if(sigmacheck(globprob)){
      printf("{%d, %d} %.1fsigma local, %.1fsigma in hist, %.1fsigma global\n",
         (int)h->GetBinLowEdge(i), (int)h->GetBinLowEdge(i+1),
         localsigma, histsigma, globsigma);
      styledrawellipse(new TEllipse(h->GetBinCenter(i), h->GetBinContent(i),
        longreadout?0.55:10, (toppad->GetUymax()-toppad->GetUymin())/30));
    }
  }

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
    0.978 - (!dividebylivetime)*0.02,
    nwindows == 1? Form("%s: %s", trigname, opts.name)
                 : Form("%s: %s x %d", trigname, opts.name, nwindows));
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


  TH1D * divided=dividebylivetime?divide_livetime(rebinned, rebinnedlive):rebinned;
  if(dividebylivetime){
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

    rebinnedlive->GetYaxis()->SetTitle("Live (%)");
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
      (rebinned->GetEntries() < 10? 0.8: rebinned->GetMaximum() < 1000?1.1:
                                         rebinned->GetMaximum() < 10000?1.4:1.6)
       / sqrt(textsize/0.05));
    if(rebin == 1)
      rebinned->GetYaxis()->SetTitle("Events/s");
    else
      rebinned->GetYaxis()->SetTitle(Form("Events/%d#kern[-0.5]{ }s", rebin));

    rebinned->Draw("e");
    if(preliminary) novapreliminary();
    randomtime();
    ltitle->Draw();

    if(opts.extexp > 0){
      rebinned->GetYaxis()->SetTitleOffset(rebinned->GetYaxis()->GetTitleOffset()*1.4);
      rebinned->GetYaxis()->SetRangeUser(opts.extexp/2, rebinned->GetMaximum()*8);
      USN.SetLogy();
    }
  }

  if(hist->Integral() == 0){
    if(rebin == 1) printf("No events, so there's no excess\n");
  }
  else if(!opts.stattest){
    bumphunt_nonpoisson(divided);
  }
  else{
    const int bestbin = bumphunt(hist, histlive, rebin, opts.extexp, opts.extexp1,
             std::min((unsigned int)(hist->Integral()), opts.polyorder));
    TPad * p = dividebylivetime? toppad: &USN;
    if(bestbin > 0){
      p->cd();
      styledrawellipse(new TEllipse(divided->GetBinCenter(bestbin),
                                    divided->GetBinContent(bestbin),
        longreadout?0.55:10, (p->GetUymax() - p->GetUymin())/30));
    }
  }

  print2(Form("%s-%ds", hist->GetName(), rebin));
}

void process(TH1D * hist, TH1D * histlive, const popts opts)
{
  printf("\n%s: ", opts.name);
  process_rebin(hist, histlive,  1, opts);
#if 0
  if(!longreadout) process_rebin(hist, histlive, 10, opts);
#endif
}

TDirectory * ligodir = NULL;

void DOIT(const char * const hname, const popts & opts)
{
  TH1D * h = getordont(hname, ligodir);
  if(h == NULL) return;

  TH1D * hlive = !dividebylivetime ? NULL
    : getordont(Form("%slive", hname), ligodir);

  if(!dividebylivetime || hlive != NULL) process(h, hlive, opts);
}

void ligopass2(const char * const infilename, const char * trigname_,
               const char * outbase_, const bool dividebylivetime_,
               const bool longreadout_, const int nwindows_)
{
  outbase = outbase_;
  dividebylivetime = dividebylivetime_;
  longreadout = longreadout_;
  trigname = trigname_;
  if(nwindows_ > 1){
    gextexp = false;
    nwindows = nwindows_;
    printf("Evaluating background using a sample of %d windows\n", nwindows);
  }

  enum stream_t stream = UNDEFINED_STREAM;
  if     (!strcmp(trigname, "FD 10Hz trigger")) stream = fardet_t02;
  else if(!strcmp(trigname, "ND energy"      )) stream = neardet_ddactivity1;
  else if(!strcmp(trigname, "FD energy"      )) stream = fardet_ddenergy;
  else if(!strcmp(trigname, "ND long readout")) stream = neardet_long;
  else if(!strcmp(trigname, "FD long readout")) stream = fardet_long;

  if(stream == UNDEFINED_STREAM){
    fprintf(stderr, "I don't know what to do with \"%s\"\n", trigname);
    exit(1);
  }

  // Don't print about making PDFs
  gErrorIgnoreLevel = kWarning;

  TFile * ligofile = new TFile(infilename, "read");
  if(!ligofile || ligofile->IsZombie()){
    fprintf(stderr, "Could not open your histogram file %s\n", infilename);
    exit(1);
  }

  ligodir = dynamic_cast<TDirectory*>(ligofile->Get("ligoanalysis"));
  if(ligodir == NULL) ligodir = dynamic_cast<TDirectory*>(ligofile->Get("ligo"));

  if(ligodir == NULL){
    fprintf(stderr, "No \"ligoanalysis\" directory in this file\n");
    exit(1);
  }

  stylecanvas(&USN);
  printstart();

  DOIT("supernovalike",   popts("Supernova-like events", 0, true));

  DOIT("unslicedbighits", popts("Unsliced big hits",     0, trigname[0] == 'N'));
  DOIT("unslicedhits",    popts("Unsliced hits",         0, trigname[0] == 'N'));


  // Skip these for the ND long readout, since they are redundant with the 100%
  // efficient ND ddactivity1 trigger.
  if(stream != neardet_long){

    // XXX Need a better way of loading background numbers in than manually
    // placing in the code, since they will be different for different signal events
    DOIT("contained_slices",
         popts("Contained slices", 0, true,
               stream == neardet_ddactivity1? 0.00096429:0));

    DOIT("tracks",                       popts("Slices with tracks",        0, true));
    DOIT("tracks_point_1",               popts("Slices w/ tracks, 16#circ", 1, true));
    DOIT("tracks_point_0",
         popts("Slices w/ tracks, 1.3#circ",1, true,
               stream == neardet_ddactivity1? 0.046: 0,
               stream == neardet_ddactivity1? 0.000: 0));


    DOIT("halfcontained_tracks",         popts("Stopping tracks",           0, true));
    DOIT("halfcontained_tracks_point_1",
         popts("Stopping tracks, 16#circ",  1, true,
               stream == neardet_ddactivity1? 0.05743: 0,
               stream == neardet_ddactivity1? 0.00000: 0));

    DOIT("halfcontained_tracks_point_0",
         popts("Stopping tracks, 1.3#circ", 1, true,
               stream == neardet_ddactivity1? 0.00087: 0,
               stream == neardet_ddactivity1? 0.00000: 0));

    DOIT("fullycontained_tracks",
         popts("Contained tracks",          0, true,
               stream == neardet_ddactivity1? 1.7857e-05: 0));

    DOIT("fullycontained_tracks_point_1",popts("Contained tracks, 16#circ", 1, true));

    DOIT("fullycontained_tracks_point_0",
         popts("Contained tracks, 1.3#circ",1, true));

    DOIT("upmu_tracks",                  popts("Upward going muons",        0, true));
    DOIT("upmu_tracks_point_1",          popts("Up-#mu, 16#circ",           1, true));
    DOIT("upmu_tracks_point_0",          popts("Up-#mu, 1.3#circ",          1, true));

    DOIT("rawtrigger",                   popts("Raw triggers",         0, true));
    DOIT("energy_low_cut",               popts(">5M ADC total",        0, true));
    DOIT("energy_low_cut_pertime",       popts(">2.5M ADC per 50#mus", 0, true));

    DOIT("energy_high_cut",              popts(">50M ADC total",       0, true, 0.0016));
    DOIT("energy_high_cut_pertime",      popts(">25M ADC per 50#mus",  0, true, 0.0027));

    DOIT("energy_vhigh_cut",             popts(">500M ADC total",      0, true, 4e-5));
    DOIT("energy_vhigh_cut_pertime",     popts(">250M ADC per 50#mus", 0, true, 4e-5));
  }

  printend();
}
