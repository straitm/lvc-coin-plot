#include "TGaxis.h"
#include "TError.h"
#include "TLatex.h"
#include "TMinuit.h"
#include "Math/Math.h"
#include "Math/ProbFuncMathCore.h"
#include "Math/QuantFuncMathCore.h"
#include "TH1.h"
#include "TH2.h"
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

//#define DEBUG
#define SEPARATEPDFS

// The total number of histograms searched for
const unsigned int NHIST = 52;

// If true, allow using external expectations for determining excess
// If false, always run fit internally (since this job will determine what
// the background is for other jobs)
bool gextexp = true;

bool is_ddactivity1 = false,
     is_fdminbias = false,
     is_ndminbias = false;

enum stream_t {
  fardet_t02 = 0x01, neardet_ddactivity1 = 0x02, fardet_ddenergy = 0x04,
  neardet_long = 0x08, fardet_long = 0x10,
  fardet_minbias = fardet_t02 | fardet_long, UNDEFINED_STREAM
};

static stream_t stream = UNDEFINED_STREAM;

struct popts{
  // Name to display on histograms and print in log
  const char * name;

  // Name to print in log and use for background lookup
  const char * codename;

  // What order polynomial to fit the background time series to
  unsigned int polyorder;

  // Whether to run statistical tests on deviations from background
  bool stattest;

  // Background per second if this histogram is too low-stats to determine
  // it in one window.  Otherwise zero.
  double extexp;

  // Slope of background if the background is time varying
  double extexp1;

  popts(const char * const name_, const char * const codename_,
        const unsigned int polyorder_, const bool stattest_,
        const double extexp_ = 0, const double extexp1_ = 0)
  {
    name = name_;
    codename = codename_;
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

const double textsize = 0.073;
const double topmargin = 0.085 * textsize/0.07;

bool dividebylivetime = false;
bool longreadout = false;

// Number of readout windows that the input data represents.  It is 1 for
// looking at signal samples and more than when for background samples.
int nwindows = 1;

const char * outbase = NULL;
const char * trigname = NULL;
const char * infilename = NULL;
const char * gwname = NULL;

bool yaxishack1000 = false;

TCanvas USN;
TPad * botpad, * midpad, * toppad;
const double divheight1 = 0.25, divheight2 = 0.47;

// Thanks ROOT for making this hard
// Ratio of top pane to middle pane
const double textratiomid = (1-divheight2)/(divheight2 - divheight1);
// Ratio of top pane to bottom pane
const double textratiobot = (1-divheight2)/divheight1;
// Ratio of top pane to whole canvas
const double textratiofull = 1-divheight2;

const double canxsize_nodiv = 150.;
const double canxsize_div = 100.;

const double canysize_nodiv = 80.;
const double canysize_div = 100.;

const double textratiodivnodiv = textratiofull*canxsize_nodiv/canxsize_div
 * canysize_div/canysize_nodiv;

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

  const int color = can == 3?kBlue+1: kBlack;

  h->SetLineColor(color);
  h->SetMarkerColor(color);

  if(longreadout || h->GetNbinsX() < 300) h->SetLineWidth(2);
  else                                    h->SetLineWidth(1);

  y->CenterTitle();
  x->CenterTitle();

  const double siz = can == 0? textsize*(dividebylivetime?1:textratiodivnodiv):
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
const double rightmargin = 0.01;

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
  c->SetCanvasSize((dividebylivetime?canxsize_div:canxsize_nodiv)*scale,
    (dividebylivetime?canysize_div:canysize_nodiv)*scale);

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
    botpad->SetBottomMargin(0.085/divheight1);

    toppad->Draw();
    botpad->Draw();
    midpad->Draw();
  }
  else{
    c->SetBottomMargin(0.155);
    c->SetLeftMargin(leftmargin);
    c->SetRightMargin(rightmargin);
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

  printf("------------- Finished processing -------------\n");
}

void print2(const char * const name)
{
#ifdef SEPARATEPDFS
  USN.Print(Form("%s-%s.eps", outbase, name));
#endif
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
// it is, and returns whether it is.  Check what the expected number
// of events was, though, because a very small number means that we
// don't really know what sigma is.
bool sigmacheck(const double p, const double expected)
{
  const unsigned int nlev = 4;
  const double lev[nlev] = { 5, 3, 2, 1 };

  for(unsigned int i = 0; i < nlev; i++)
    if(p < (1-ROOT::Math::gaussian_cdf(lev[i], 1))*2){
      if(expected > 1e-50){
        printf("*******************************************\n"
               "*******************************************\n"
               "******* That's more than %.1f sigma ********\n"
               "*******************************************\n"
               "*******************************************\n", lev[i]);
        return true;
      }
      else{
        printf("********************************************\n"
               "********************************************\n"
               "* inf sigma. Need a better background est! *\n"
               "********************************************\n"
               "********************************************\n");
        return false;
      }
    }
  return false;
}

void fit(TMinuit & mn, TH1D * hist, TH1D * histlive,
         const unsigned int polyorder)
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

    // Since we use fabs() in the FCN to avoid warnings, have to fix it here...
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

void stylefunction(TF1 * f, const bool externalexpectation)
{
  f->SetLineColor(externalexpectation?kBlue:kRed);
  f->SetLineStyle(1);
  f->SetLineWidth(3);
  if(!externalexpectation) f->SetNpx(500);
}

void stylefitraw(TH1D * fitraw, const bool externalexpectation)
{
  fitraw->SetLineWidth(2);
  fitraw->SetLineColor(externalexpectation?kBlue:kRed);
  fitraw->SetMarkerColor(externalexpectation?kBlue:kRed);
  fitraw->SetLineStyle(1);
}

enum stream_t decodestream(const char * const trigname)
{
  enum stream_t strm = UNDEFINED_STREAM;
  if     (!strcmp(trigname, "FD 10Hz trigger"  )) strm = fardet_t02;
  else if(!strcmp(trigname, "ND energy trigger")) strm = neardet_ddactivity1;
  else if(!strcmp(trigname, "FD energy trigger")) strm = fardet_ddenergy;
  else if(!strcmp(trigname, "ND long readout"  )) strm = neardet_long;
  else if(!strcmp(trigname, "FD long readout"  )) strm = fardet_long;

  if(strm == UNDEFINED_STREAM){
    fprintf(stderr, "I don't know what to do with \"%s\"\n", trigname);
    exit(1);
  }
  return strm;
}

// Look for an excess anywhere in the histogram.  If 'extexp'
// is positive, use it as the background per second.  Otherwise,
// fit the histogram using a polnominal of order 'polyorder' and
// use that as the background.  Return the bin with the biggest
// excess, or -1 if none has an interesting excess
int bumphunt(TH1D * hist, TH1D * histlive, const int rebin,
             const double extexp, const double extexp1,
             unsigned int polyorder, const char * codename)
{
  if(histlive == NULL || !dividebylivetime)
    histlive = dummy_histlive(hist, rebin);

  const bool verbose = rebin == 1;
  TMinuit mn(polyorder+1);
  setup_mn(mn, polyorder);

  if(extexp > 0){
    polyorder = 1; // Just assume 1 (enforced by definition of popts) and
                   // give the draw function enough parameters to do its work
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

  const double scale = yaxishack1000? 1e-3: 1;

  if(extexp > 0){
    poly->SetParameter(0, rebin*extexp*scale);
    poly->SetParameter(1, rebin*extexp1*scale);
  }
  else{
    for(unsigned int i = 0; i <= polyorder; i++){
      poly->SetParameter(i, getpar(mn, i)*scale);
      poly->SetParError(i, geterr(mn, i)*scale);
    }
  }
  stylefunction(poly, extexp > 0);
  poly->Draw("same");

  if(dividebylivetime){
    TH1D * fitraw = new TH1D(Form("%sfit%d", hist->GetName(), rebin), "",
                             hist->GetNbinsX(), hist->GetBinLowEdge(1),
                             hist->GetBinLowEdge(hist->GetNbinsX()+1));
    stylefitraw(fitraw, extexp > 0);

    for(int i = 1; i <= fitraw->GetNbinsX(); i += rebin){
      double sum = 0;
      for(int j = 0; j < rebin; j++)
        sum += poly->Eval(fitraw->GetBinCenter(i+j))
               * histlive->GetBinContent(i+j)/scale;
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
      * histlive->GetBinContent(i)/scale;
    const int actual = (int)hist->GetBinContent(i);

    const double prob = prob_this_or_more(actual, expected);

    if(prob < minprob){
      minprob = prob;
      minprob_bin = i;
    }
  }

  const double expected = poly->Eval(hist->GetBinCenter(minprob_bin))
    * histlive->GetBinContent(minprob_bin)/scale;
  const int actual = (int)hist->GetBinContent(minprob_bin);

  const int nbins = (stream == neardet_long || stream == fardet_long)? 45: 1000;

  const double localprob = prob_this_or_more(actual, expected);
  const double histprob = lookelse(localprob, nbins);
  const double globprob = lookelse(localprob, NHIST*nbins);

  //if(globprob < 0.5)
    printf("Highest: {%d, %d}s: %d, %.4f exp. P hist: %.2g (%.3g global)\n",
           (int)hist->GetBinLowEdge(minprob_bin),
           (int)hist->GetBinLowEdge(minprob_bin+1),
           actual, expected, histprob, globprob);
  int ret = -1;
  if(sigmacheck(globprob, expected)) ret = minprob_bin;

  if(extexp == 0){
    if(polyorder == 0){
      const double dem = histlive->Integral()*histlive->GetBinWidth(1)
          *(dividebylivetime?1:nwindows);
      const double mean = hist->Integral()/dem;
      const double error = sqrt(hist->Integral())/dem;

      printf("Number per second: %-6.3g\n"
             "                +- %-6.3g\n", mean, error);
    }
    else if(polyorder == 1){
      printf("Fit: (\n%s_%s_P0 %-6.3g\n"
             "   +- %-6.3g) + (\n%s_%s_P1 %-6.3g\n"
             "                 +- %-6.3g)t\n",
             codename, is_ddactivity1?"DDACTIVITY1":is_fdminbias?"FDMINBIAS":"?",
             poly->GetParameter(0)/scale/(dividebylivetime?1:nwindows),
             poly->GetParError(0)/scale/(dividebylivetime?1:nwindows),
             codename,is_ddactivity1?"DDACTIVITY1":is_fdminbias?"FDMINBIAS":"?",
             poly->GetParameter(1)/scale/(dividebylivetime?1:nwindows),
             poly->GetParError(1)/scale/(dividebylivetime?1:nwindows));
    }
  }

  // We'll consider the first N seconds after the event to be special
  const int specialbins = 10;
  const int specialbin1 = hist->GetXaxis()->FindBin(0.5);
  if(minprob_bin >= specialbin1 && minprob_bin < specialbin1+specialbins){
    const double slocalprob = prob_this_or_more(actual, expected);
    const double shistprob = lookelse(slocalprob, specialbins);
    const double sglobprob = lookelse(slocalprob, NHIST*specialbins);
    //if(sglobprob < 0.5)
      printf("This was in first %ds. P: %.2g (%.2g hist) (%.3g global)\n",
             specialbins, slocalprob, shistprob, sglobprob);
    sigmacheck(sglobprob, expected);
  }

  return ret;
}

void novapreliminary()
{
  static TLatex * t = new TLatex(0, 0, "NOvA Preliminary");
  t->SetTextColorAlpha(kBlue, 0.5);
  t->SetTextSize(dividebylivetime?textsize*textratiodivnodiv:textsize);
  t->SetTextFont(42);
  t->SetNDC();
  t->SetTextAlign(32);
  const double raise = dividebylivetime? 0.08:-0.085;
  const double shift = dividebylivetime? 0.02:0.00;
  t->SetX(1 - rightmargin - 0.02 + shift);
  t->SetY(1 -  topmargin  + 0.01 + raise);
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

__attribute__((unused)) void printhist(TH1D & h)
{
  for(int i = 1; i <= h.GetNbinsX(); i++){
    const int n = h.GetBinContent(i)*80/h.GetMaximum();
    for(int j = 0; j < n; j++) printf("#");
    printf("\n");
  }
  printf("\n");
}

static double minnonzero(TH1D * h)
{
  double min = 1e300;
  for(int i = 1; i <= h->GetNbinsX(); i++)
    if(h->GetBinContent(i) && h->GetBinContent(i) < min)
      min = h->GetBinContent(i);
  return min;
}

// Search for excursions from normal behavior for a series that isn't
// Poissonian, like the total number of hits, which obviously has a lot
// of correlated activity.
void bumphunt_nonpoisson(TH1D * h)
{
  TH1D side("side", "", 100, minnonzero(h)*0.999, h->GetMaximum()*1.001);

  for(int i = 1; i <= h->GetNbinsX(); i++){
    if(longreadout && h->GetBinCenter(i) < np_longreadoutmin) continue;
    if(!longreadout && h->GetBinCenter(i) > np_notlongreadoutmax) continue;

    side.Fill(h->GetBinContent(i));
  }

  if(side.Integral() == 0){
    printf("No data in baseline region\n");
    return;
  }

  TF1 Gaus("g", "[0]*exp(-pow(x-[1],2)/(2*[2]*[2]))",
    minnonzero(h)*0.999, h->GetMaximum()*1.001);
  Gaus.SetParameter(0, 10);
  Gaus.SetParameter(1, (h->GetMaximum()+minnonzero(h))/2);
  Gaus.SetParameter(2, (h->GetMaximum()-minnonzero(h))/10);

  // Important to prevent Fit from going crazy if there are outliers
  Gaus.SetParLimits(1, minnonzero(h), h->GetMaximum());

  side.Fit("g", "l0");

  if(side.GetFunction("g") == NULL){
    fprintf(stderr, "Fit failed in bumphunt_nonpoisson on %s\n",h->GetName());
    return;
  }
  const double mean = side.GetFunction("g")->GetParameter(1);
  const double sigma = side.GetFunction("g")->GetParameter(2);

  toppad->cd();

  draw_nonpoisson_f(mean, 3, kRed);
  draw_nonpoisson_f(mean+sigma, 2, kRed);
  draw_nonpoisson_f(mean-sigma, 2, kRed);
  draw_nonpoisson_f(mean+2*sigma, 2, kRed+1);
  draw_nonpoisson_f(mean-2*sigma, 2, kRed+1);

  for(int i = 1; i <= h->GetNbinsX(); i++){
    if(longreadout && h->GetBinCenter(i) > np_longreadoutmin) continue;
    if(!longreadout && h->GetBinCenter(i) < np_notlongreadoutmax) continue;

    const double content = h->GetBinContent(i);

    // Important to take the bin error into account, even though it is
    // a bit of double-counting, or else low livetime bins look like
    // big excesses.  The double-counting is only significant if the
    // errors are large compared to the spread, in which case, we probably
    // shouldn't be using this sort of search.
    const double localsigma = (content - mean)/
      sqrt(pow(sigma, 2) + pow(h->GetBinError(i), 2));

    const int trials = longreadout?40:500; // XXX fragile
    const double localprob = ROOT::Math::gaussian_cdf_c(localsigma);
    const double histprob = lookelse(localprob, trials);
    const double globprob = lookelse(localprob, NHIST*trials);

    const double histsigma = ROOT::Math::gaussian_quantile_c(histprob, 1);
    const double globsigma = ROOT::Math::gaussian_quantile_c(globprob, 1);

    if(sigmacheck(globprob, mean)){
      printf("{%d, %d} %f with %f expected\n"
             "P = %.1g local, %.1g in hist, %.3g global\n",
         (int)h->GetBinLowEdge(i), (int)h->GetBinLowEdge(i+1),
         content, mean,
         localprob, histprob, globprob);
      styledrawellipse(new TEllipse(h->GetBinCenter(i), h->GetBinContent(i),
        longreadout?0.55:10, (toppad->GetUymax()-toppad->GetUymin())/30));
    }
  }
}

// Given a number x, return a number close to x that isn't too near a multiple
// of 5.  This is for Y axes ranges so that the bottom label isn't cut off.
double notnear5(const double x)
{
  // Axes labels definitely won't be in multiples of 5 anyway
  if(x < 10) return x;
  if(x > 100) return x;

  if(x - int(x/5)*5 < 1) return int(x)/5 * 5 - 1;
  if(x - int(x/5)*5 > 4) return int(x)/5 * 5 + 4;

  return x;
}

// Given a number x, return a number close to x that isn't too near a multiple
// of 200.  This is for Y axes ranges so that the top label isn't cut off.
double notnear200(const double x)
{
  // Axes labels definitely won't be in multiples of 200 anyway
  if(x < 500) return x;
  if(x > 2000) return x;

  if(((int)x)%200 <  50) return int(x)/200 * 200 - 0.1;
  if(((int)x)%200 > 150) return int(x)/200 * 200 + 199.9;

  return x;
}

static TH1D * snhist = NULL, * snhistlive = NULL, * snmc = NULL;

// The objective function for MINUIT.
static void snlike(__attribute__((unused)) int & np,
                   __attribute__((unused)) double * gin,
                   double & llike, double * par, __attribute__((unused)) int flag)
{
  llike = 0;

  const double flat = par[0];
  const double snnorm = par[1];

  for(int i = 1; i <= snhist->GetNbinsX(); i++){
    const double data = snhist->GetBinContent(i);

    const double snmcexpect = snmc->GetBinContent(i) * snnorm;
    const double expect = snhistlive->GetBinContent(i)*(flat + snmcexpect);

    if(expect <= 0) continue;

    llike += expect;

    if(data > 0)
      llike += - data + data * log(data/expect);
  }

  llike *= 2;
}

static int dig(const double n)
{
  return (Form("%.1e", n)[0] - '0' <= 2 &&
          Form("%.2e", n)[0] - '0' <= 2)? 2: 1;
}

void draw_flux_limit(const double step, const double * const prob,
                     const unsigned int N, const double limit, const bool is27)
{
  const double size = 0.067;

  const double markersize = 1.5;

  TCanvas * c1 = new TCanvas();
  c1->SetRightMargin(0.04);
  c1->SetTopMargin(0.03);
  c1->SetBottomMargin(0.15);
  c1->SetLeftMargin(0.14);
  c1->SetCanvasSize(600, 400);
  gStyle->SetOptStat(0);
  c1->SetTickx();
  c1->SetTicky();

  const double xmin = 0, xmax = stream == fardet_long || stream == neardet_long?2: 20,
               ymin = 0, ymax = 1.099;

  const int stepstep = xmax < 2?1: xmax < 20?10:30;

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
  dum->GetXaxis()->SetNdivisions(stream == fardet_long || stream == neardet_long? 8:516);
  dum->GetYaxis()->SetDecimals();

  dum->Draw();

  TGraph g, gfill;
  for(unsigned int i = 0; i < N; i+=stepstep){
    if(i*step > xmax) continue;

    g.SetPoint(g.GetN(), i*step, prob[i]);
    if(i*step < limit)
      gfill.SetPoint(g.GetN(), i*step, prob[i]);
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
    Form("%s solar mass Garching", is27?"27":"9.6"), "");
  leg.AddEntry((TH1D*)NULL, 
    stream == fardet_t02? "FD 10Hz"
   :stream == fardet_long? "FD long readout"
   :stream == neardet_long?"ND long readout"
   :"???", "");

  leg.AddEntry((TH1D*)NULL,
    Form("< %.1f of 10#kern[-0.5]{ }kpc at 90%% CL", limit), "");
  leg.AddEntry((TH1D*)NULL,
    Form("#Rightarrow  > %.0f#kern[-0.5]{ }kpc at 90%% CL", 10/sqrt(limit)), "");

  leg.Draw();

  c1->RedrawAxis();

  c1->SaveAs(Form("fluxlimit-%s-%s-%ssolarmass.pdf", gwname,
    stream == fardet_t02? "FD10Hz"
   :stream == fardet_long? "FDlongreadout"
   :stream == neardet_long?"NDlongreadout"
   :"somethingelse",
    is27?"27":"9.6"));
}

// Given a distance limit in kpc, return a fluence limit in neutrinos/cm2
static double fllimit(const double distlimit, const bool is27)
{
  // We could do this in erg/cm2, but let's do it in neutrinos/cm2
  // to match the KamLAND paper.
#if 0
  // Found from
  // awk '!/^#/{sum += $2*($1-oldt); oldt=$1}END{print sum}'
  //   neutrino_signal_nu_{x,xbar,e,ebar}-LS220-z9.6co.txt
  //
  // same for neutrino_signal_nu_{x,xbar,e,ebar}-LS220-s27.0co.txt
  const double totE = is27? 224.504e51: 125.637e51;
#endif

  // awk '!/^#/{if($3 != 0) 
  //   sum += $2*($1-oldt) * 1e51/($3 / 624150.91); oldt=$1}
  //   END{printf("%g\n", sum)}'
  // neutrino_signal_nu_{x,xbar,e,ebar}-LS220-s27.0co.txt
  const double totN = is27? 1.11556e+58: 6.78803e+57;

  const double benchmarkdist = 10.; // kpc

  const double cmperkpc = 3.0856776e+21;

  const double cm2at10kpc = 4*M_PI*pow(benchmarkdist * cmperkpc, 2);

  const double fluencerat = pow(benchmarkdist/distlimit, 2);

  const double benchmarkfluence = totN/cm2at10kpc;

  return fluencerat * benchmarkfluence;
}

static TH1D * sn_10kpc_27 = NULL, * sn_10kpc_96 = NULL, * sn_10kpc = NULL;

std::pair<double, double>
supernova_flux_limit(TH1D *hist, TH1D * histlive, const double bg, const bool is27)
{
  sn_10kpc_27 = (TH1D*)hist->Clone("sn");
  sn_10kpc_96 = (TH1D*)hist->Clone("sn");
  sn_10kpc_27->Reset();
  sn_10kpc_96->Reset();
  sn_10kpc = is27? sn_10kpc_27: sn_10kpc_96;

  if(is_fdminbias){
    // For 9.6 solar mass Garching flux, as found in Andrey's files.  I just had
    // to do a little post-processing since somehow all the timestamps are zero
    // in his file.  However, as long as I assume that it's a progression of
    // contiguous 5ms events, it's fine.
    //
    // The 6th bin and onward for 27 solar masses are taken from an analytic
    // Garching curve, matched to the shape of the MC.
    sn_10kpc_27->SetBinContent(sn_10kpc->FindBin(0), 239                  );
    sn_10kpc_27->SetBinContent(sn_10kpc->FindBin(1),  77                  );
    sn_10kpc_27->SetBinContent(sn_10kpc->FindBin(2),  41                  );
    sn_10kpc_27->SetBinContent(sn_10kpc->FindBin(3),  26                  );
    sn_10kpc_27->SetBinContent(sn_10kpc->FindBin(4),  15                  );
    sn_10kpc_27->SetBinContent(sn_10kpc->FindBin(5),  10                  );
    sn_10kpc_27->SetBinContent(sn_10kpc->FindBin(6), 52.7209*pow(2./10, 2));
    sn_10kpc_27->SetBinContent(sn_10kpc->FindBin(7), 35.4769*pow(2./10, 2));
    sn_10kpc_27->SetBinContent(sn_10kpc->FindBin(8), 24.0291*pow(2./10, 2));
    sn_10kpc_27->SetBinContent(sn_10kpc->FindBin(9), 15.4561*pow(2./10, 2));
    sn_10kpc_27->SetBinContent(sn_10kpc->FindBin(10), 8.54  *pow(2./10, 2));
    sn_10kpc_27->SetBinContent(sn_10kpc->FindBin(11), 4.70  *pow(2./10, 2));

    sn_10kpc_96->SetBinContent(sn_10kpc->FindBin(0), 47);
    sn_10kpc_96->SetBinContent(sn_10kpc->FindBin(1), 42);
    sn_10kpc_96->SetBinContent(sn_10kpc->FindBin(2), 16);
    sn_10kpc_96->SetBinContent(sn_10kpc->FindBin(3),  9);
    sn_10kpc_96->SetBinContent(sn_10kpc->FindBin(4),  7);
    sn_10kpc_96->SetBinContent(sn_10kpc->FindBin(5),  1);
    sn_10kpc_96->SetBinContent(sn_10kpc->FindBin(6),  0);
    sn_10kpc_96->SetBinContent(sn_10kpc->FindBin(7),  0);
    sn_10kpc_96->SetBinContent(sn_10kpc->FindBin(8),  0);
    sn_10kpc_96->SetBinContent(sn_10kpc->FindBin(9),  0);
    sn_10kpc_96->SetBinContent(sn_10kpc->FindBin(10), 0);
    sn_10kpc_96->SetBinContent(sn_10kpc->FindBin(11), 0);
  }
  else{
    // For 9.6 solar mass Garching flux, as found in Andrey's files.  His 2kpc
    // file has 89 events total, which makes 3.6 at 10kpc, or 18 true events at
    // 10kpc assuming 20% efficiency.  I think that efficiency means fiducial
    // volume, so that's really 27 true events or so.
    //
    // This is close to the 33 IBD events predicted by GVKM, as found on
    // Justin's 2017 APS poster.  It's still close if elastic scattering is
    // included in the Garching MC, which looks to be true.
    //
    // XXX the stats are really too low here

    sn_10kpc_27->SetBinContent(sn_10kpc->FindBin(0), 5.22378  );
    sn_10kpc_27->SetBinContent(sn_10kpc->FindBin(1), 1.73077  );
    sn_10kpc_27->SetBinContent(sn_10kpc->FindBin(2), 0.314685 );
    sn_10kpc_27->SetBinContent(sn_10kpc->FindBin(3), 0.251748 );
    sn_10kpc_27->SetBinContent(sn_10kpc->FindBin(4), 0.0314685);
    sn_10kpc_27->SetBinContent(sn_10kpc->FindBin(5), 0.0314685);

    sn_10kpc_96->SetBinContent(sn_10kpc->FindBin(0), 1.1014   );
    sn_10kpc_96->SetBinContent(sn_10kpc->FindBin(1), 0.818182 );
    sn_10kpc_96->SetBinContent(sn_10kpc->FindBin(2), 0.660839 );
    sn_10kpc_96->SetBinContent(sn_10kpc->FindBin(3), 0.22028  );
    sn_10kpc_96->SetBinContent(sn_10kpc->FindBin(4), 0.0629371);
    sn_10kpc_96->SetBinContent(sn_10kpc->FindBin(5), 0.0314685);
  }

  const double replive = histlive->GetBinContent(histlive->FindBin(1.5));
  const double mostlylive = replive > 0.1;

  TMinuit mn(2);
  mn.fGraphicsMode = false;
  int ierr;
  mn.mnparm(0, "const", is_fdminbias?400:0.4, is_fdminbias?30.:0.05, 0, 0, ierr);
  mn.mnparm(1, "sn",    0.1 , 0.1, 0, mostlylive?2:400, ierr);

  mn.SetFCN(snlike);

  if(bg != 0){
    mn.Command(Form("SET PAR 1 %f", bg));
    mn.Command("FIX 1");
  }

  snhist = hist;
  snhistlive = histlive;
  snmc = sn_10kpc;

  mn.Command("SCAN");
  mn.Command("MIGRAD");
  if(bg == 0) mn.Command("MNCONT 1 2 100");

  const std::pair<double, double> ans(getpar(mn, 0), getpar(mn, 1));

  const double gmin = mn.fAmin;

  const double step = 0.001;
  const int N = 100000;
  double prob[N];
  memset(prob, 0, sizeof(double)*N);

  mn.Command("FIX 2");

  mn.Command("SET PRINT -1");
  bool worthit = true;
  for(int i = 0; i < N; i++){
    if(worthit){
      mn.Command(Form("SET PAR 2 %f", i*step));
      if(bg == 0) mn.Command("MIGRAD");
      else        mn.Command("CALL");
      prob[i] = exp(-(mn.fAmin-gmin));
      printf("snprob%s %10f %10g\n", is27?"27":"9.6", i*step, prob[i]);
      if(i > 0 && prob[i] < prob[i-1] && prob[i] < 1e-30) worthit = false;
    }
    else{
      printf("snprob%s %10f 0\n", is27?"27":"9.6", i*step);
    }
  }
  mn.Command("SET PRINT 1");

  double totprob = 0;
  for(int i = 0; i < N; i++) totprob += prob[i];

  double accprob = 0;
  double limit = 0;
  for(int i = 0; i < N; i++){
    accprob += prob[i];
    if(accprob/totprob > 0.9){
      limit = i*step;
      const double distlimit = 10/sqrt(limit);
      const double fluencelimit = fllimit(distlimit, is27);
      const int ndigf = dig(limit);
      const int ndigd = dig(distlimit);
      const int ndigfl = dig(fluencelimit);
      printf("90%% CL at %.*g of 10kpc %sSN-like flux, %.*g kpc, %#8.*g\n",
             ndigf, limit,
             is27?"27":"9.6",
             ndigd, distlimit,
             ndigfl, fluencelimit/1e12);
      break;
    }
  }

  draw_flux_limit(step, prob, N, limit, is27);

#if 0
  exit(0);
#endif

  return ans;
}

// Process, possibly rebinning.  (If rebin=1, don't rebin.)
// polyorder: order of the polynomial to fit for the no-signal hypothesis.
//            Naively, zero.  For anything with pointing, at least 1.
void process_rebin(TH1D *hist, TH1D * histlive,
                   const int rebin, const popts opts)
{
  TH1D * rebinned = (TH1D*)hist->Rebin(rebin, "rebinned");
  TH1D * rebinnedlive = histlive == NULL? NULL:
    (TH1D*)histlive->Rebin(rebin, "rebinnedlive");

  stylehist(rebinned, 0);

  std::pair<double, double> snfit27;
  if(!strcmp(opts.name, "Supernova-like"))
    snfit27 = supernova_flux_limit(hist, histlive, opts.extexp, true);

  const std::pair<double, double> snfit = !strcmp(opts.name, "Supernova-like")?
     supernova_flux_limit(hist, histlive, opts.extexp, false)
    : std::pair<double, double>(0, 0);

  const bool setlogy = opts.extexp > 1 || !dividebylivetime;

  const bool preliminary = false;
  const bool usegwname = true;

  // Some trickery (static) here needed to avoid getting the first label stuck
  // on all the output PDFs.
  static TLatex * ltitle = new TLatex;

  char title_event[1024];
  strncpy(title_event, infilename, 1024);
  if(NULL != strstr(title_event, "Z"))
    *(strstr(title_event, "Z")+1) = '\0';

  if(usegwname){
     strncpy(title_event, gwname, std::min((size_t)1024, strlen(gwname)+1));
     title_event[1023] = '\0';
   }

  ltitle->SetText(
    0.5 + leftmargin/2 - rightmargin/2,
    dividebylivetime?0.975:0.955,
    nwindows == 1? Form("%s #minus %s: %s", title_event, trigname, opts.name)
                 : Form("#splitline{%s}{%s: %s x %d}", title_event, trigname,
                        opts.name, nwindows));
  ltitle->SetTextSize(dividebylivetime?textsize*textratiofull:
                                       textsize*textratiodivnodiv);
  ltitle->SetTextFont(42);
  ltitle->SetTextAlign(22);
  ltitle->SetNDC();

  USN.cd();
  ltitle->Draw();

  if(preliminary) novapreliminary();

  double mint = -500, maxt = 500;

  if(longreadout && rebin == 1)
    mint = -10, maxt = 45;

  rebinned->GetXaxis()->SetRangeUser(mint, maxt);
  if(rebinnedlive != NULL)
    rebinnedlive->GetXaxis()->SetRangeUser(mint, maxt);

  TH1D * divided
    = dividebylivetime?divide_livetime(rebinned, rebinnedlive):rebinned;
  if(dividebylivetime){
    if(rebin == 1){
      divided->GetXaxis()->SetRangeUser(mint, maxt);

      double miny = divided->GetMaximum();
      double maxy = -1e100;
      double extramax = 0, extramin = 0;

      for(int i = 1; i <= divided->GetNbinsX(); i++){
        const double vmin = divided->GetBinContent(i)-divided->GetBinError(i);
        const double vmax = divided->GetBinContent(i)+divided->GetBinError(i);
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

      // Special cases to make particular plots work
      if(miny - extramin < 350 && miny - extramin > 345) miny += 10;
      if(miny - extramin < 360 && miny - extramin > 357) miny += 3;
      if(miny - extramin < 520e3 && miny - extramin > 515e3) miny += 4e3;

      if(!setlogy)
        divided->GetYaxis()->SetRangeUser(0.001, maxy+extramax);
      if(maxy < 2*miny)
        divided->GetYaxis()->SetRangeUser(miny-extramin, maxy+extramax);
    }

    stylehist(divided, 1);
    stylehist(rebinned, 2);
    stylehist(rebinnedlive, 3);

    yaxishack1000 = divided->GetMaximum() > 10e3;
    if(yaxishack1000) divided->Scale(1e-3);
    divided->GetYaxis()->SetTitle(Form("%sCorrected events/s",
      yaxishack1000?"10^{3} ":""));
    rebinned->GetYaxis()->SetTitle("Raw events/s");

    const double ytitleoff = 1.2 / sqrt(textsize/0.05);
    divided->GetYaxis()->SetTitleOffset(ytitleoff
      + 0.08 * (divided->GetMaximum() > 10000 && divided->GetMaximum() < 100000));
    if(!strcmp(hist->GetName(), "supernovalike")){
      if(!strcmp(trigname, "ND long readout"))
        divided->GetYaxis()->SetRangeUser(0.01, 7);
      else if(!strcmp(trigname, "FD long readout"))
        divided->GetYaxis()->SetRangeUser(380.1, 730);
    }
    else{
      const double hi = divided->GetMaximum(), lo = minnonzero(divided);
      divided->GetYaxis()->SetRangeUser(notnear5(lo - (hi-lo)*0.22),
                                                 hi + (hi-lo)*0.12);
    }

    rebinned->GetYaxis()->SetTitleOffset(ytitleoff/textratiomid);

    // So it aligns to bottom of letters, not bottom of parentheses...
    rebinnedlive->GetYaxis()->SetTitleOffset(ytitleoff/textratiobot * 0.955);

    rebinnedlive->GetYaxis()->SetTitle("Live (%)");
    rebinnedlive->Scale(100./rebin);
    rebinnedlive->GetYaxis()->SetRangeUser(0, rebinnedlive->GetMaximum()*1.23);

    toppad->cd();
    divided->Draw("e");
    toppad->Update(); // Necessary to get Uymin correctly!
    if(setlogy && toppad->GetUymin() == 0) divided->GetYaxis()->ChangeLabel(1, -1, 0);

    midpad->cd();

    const double midpadrange = rebinned->GetMaximum() - rebinned->GetMinimum();
    rebinned->GetYaxis()->SetRangeUser(
      std::max(0., rebinned->GetMinimum() - midpadrange*0.2),
      notnear200(rebinned->GetMaximum() + midpadrange*0.2)
    );

    rebinned->Draw("hist");
    midpad->Update(); // Necessary to get Uymin correctly!
    if(midpad->GetUymin() == 0) rebinned->GetYaxis()->ChangeLabel(1, -1, 0);

    botpad->cd();
    rebinnedlive->Draw(longreadout?"hist":"hist][");
  }
  else{
    rebinned->GetYaxis()->SetTitleOffset(
      (rebinned->GetEntries() < 10? 0.8: rebinned->GetMaximum()<1000?1.1:
                                         rebinned->GetMaximum()<10000?1.4:1.6)
       / sqrt(textsize/0.05) * canysize_nodiv/canysize_div);
    if(rebin == 1)
      rebinned->GetYaxis()->SetTitle("Events/s");
    else
      rebinned->GetYaxis()->SetTitle(Form("Events/%d#kern[-0.5]{ }s", rebin));

    rebinned->Draw("e");
    if(preliminary) novapreliminary();
    ltitle->Draw();
  }

  // Values like 1e-100 mean we have no background measurement
  if(setlogy && opts.extexp > 1e-50){
    TH1 * h = dividebylivetime? divided: rebinned;
    h->GetYaxis()->SetTitleOffset(h->GetYaxis()->GetTitleOffset()
                                  *(dividebylivetime? 1: 1.4));
    h->GetYaxis()->SetRangeUser(opts.extexp/2, h->GetMaximum()*8);
  }

  (dividebylivetime?toppad:&USN)->SetLogy(setlogy);

  if(hist->Integral() == 0){
    if(rebin == 1){
      printf("No events, so there's no excess\n");
      if(opts.polyorder == 1){
        printf("%s_%s_P0 1e-100\n"
               "%s_%s_P1 1e-100\n",
          opts.codename, is_ddactivity1?"DDACTIVITY1":is_fdminbias?"FDMINBIAS":"?",
          opts.codename, is_ddactivity1?"DDACTIVITY1":is_fdminbias?"FDMINBIAS":"?");
      }
    }
  }
  else if(!opts.stattest){
    bumphunt_nonpoisson(divided);
  }
  else{
    const int bestbin =
      bumphunt(hist, histlive, rebin, opts.extexp, opts.extexp1,
               std::min((unsigned int)(hist->Integral()), opts.polyorder),
               opts.codename);
    TPad * p = dividebylivetime? midpad: &USN;
    if(bestbin > 0){
      p->cd();
      styledrawellipse(new TEllipse(hist->GetBinCenter(bestbin),
                                    hist->GetBinContent(bestbin),
        (maxt-mint)/50, (p->GetUymax() - p->GetUymin())/20));
    }
  }

  if(!strcmp(opts.name, "Supernova-like")){
    toppad->cd();
    for(int i = 1; i <= sn_10kpc_27->GetNbinsX(); i++)
      sn_10kpc_27->SetBinContent(i, sn_10kpc_27->GetBinContent(i) + snfit27.first);
    sn_10kpc_27->SetLineWidth(3);
    sn_10kpc_27->SetLineStyle(7);
    sn_10kpc_27->SetLineColor(is_fdminbias?kRed:kBlue);
    sn_10kpc_27->Draw("same][");

    for(int i = 1; i <= sn_10kpc_96->GetNbinsX(); i++)
      sn_10kpc_96->SetBinContent(i, sn_10kpc_96->GetBinContent(i) + snfit.first);
    sn_10kpc_96->SetLineWidth(2);
    sn_10kpc_96->SetLineStyle(1);
    sn_10kpc_96->SetLineColor(is_fdminbias?kRed:kBlue);
    sn_10kpc_96->Draw("same][");
  }

  print2(Form("%s-%ds", hist->GetName(), rebin));
}

void process(TH1D * hist, TH1D * histlive, const popts & opts)
{
  if(strlen(opts.codename) > 0)
    printf("\n%s (%s_%s):\n", opts.name, opts.codename,
           is_ddactivity1?"DDACTIVITY1":
           is_fdminbias?  "FDMINBIAS":
           is_ndminbias?  "": "?");
  else
    printf("\n%s:\n", opts.name);
  process_rebin(hist, histlive,  1, opts);
#if 0
  if(!longreadout) process_rebin(hist, histlive, 10, opts);
#endif
}

TDirectory * ligodir = NULL;

// Return true iff we did something
bool DOIT(const char * const hname, const popts & opts)
{
  TH1D * h = getordont(hname, ligodir);
  if(h == NULL) return false;

  TH1D * hlive = !dividebylivetime ? NULL
    : getordont(Form("%slive", hname), ligodir);

  if(!dividebylivetime || hlive != NULL){
    process(h, hlive, opts);
    return true;
  }

  return false;
}

std::map<std::string, double> readbg(const char * const bgfile)
{
  std::ifstream infile(bgfile);
  if(!infile.is_open()){
    printf("Could not open %s\n", bgfile);
    exit(1);
  }
  std::string key;
  double value;
  std::map<std::string, double> extbg;
  while(infile >> key >> value) extbg[key] = value;
  printf("Read in %u external background numbers from %s\n",
    (unsigned int)extbg.size(), bgfile);
  return extbg;
}

void ligopass2(const char * const infilename_, const char * trigname_,
               const char * outbase_, const bool dividebylivetime_,
               const bool longreadout_, const int nwindows_,
               const char * const bgfile, const char * const gwname_)
{
  infilename = infilename_;
  outbase = outbase_;
  dividebylivetime = dividebylivetime_;
  longreadout = longreadout_;
  trigname = trigname_;
  gwname = gwname_;

  TFile * ligofile = new TFile(infilename, "read");
  if(!ligofile || ligofile->IsZombie()){
    fprintf(stderr, "Could not open your histogram file %s\n", infilename);
    exit(1);
  }

  if(nwindows_ > 1){
    gextexp = false;
    nwindows = nwindows_;
    printf("Evaluating background using a sample of %d windows\n", nwindows);
  }

  stream = decodestream(trigname);

  gErrorIgnoreLevel = kWarning; // Don't print about making PDFs

  ligodir = dynamic_cast<TDirectory*>(ligofile->Get("ligoanalysis"));

  if(ligodir == NULL){
    fprintf(stderr, "No \"ligoanalysis\" directory in this file\n");
    exit(1);
  }

  stylecanvas(&USN);
  printstart();

  is_ddactivity1 = stream == neardet_ddactivity1,
  is_fdminbias = (stream & fardet_minbias),
  is_ndminbias = stream == neardet_long;

  if(DOIT("blind", popts("Livetime report only", "", 0, true))){
    printend();
    return; // Make double-sure we stay blind if this is a blind report
  }

  DOIT("supernovalike",
       popts("Supernova-like", "", 0, true, is_ndminbias? 0.53: 0));

  DOIT("unslicedbighits", popts("Sub-supernova", "", 0, trigname[0] == 'N'));
  DOIT("unslicedhits",    popts("Sub-supernova, no threshold", "",     0, trigname[0] == 'N'));

  // Skip rest, since it's redundant w/ 100% efficient ddactivity1.
  if(stream == neardet_long){
    printend();
    return;
  }

  DOIT("tracks", popts("All tracks", "", 0, true));

  DOIT("tracks_point_1", popts("All tracks, 16#circ", "", 1, true));

  std::map<std::string, double> extbg = readbg(bgfile);

  DOIT("tracks_point_0",
    popts("All tracks, 1.3#circ", "TRACKS_POINT_0", 1, true,
          is_ddactivity1? extbg.at("TRACKS_POINT_0_DDACTIVITY1_P0"): 0,
          is_ddactivity1? extbg.at("TRACKS_POINT_0_DDACTIVITY1_P1"): 0));

  DOIT("halfcontained_tracks",
       popts("Stopping tracks", "", 0, true,
             is_ddactivity1? 0.401: 0));

  DOIT("halfcontained_tracks_point_1",
    popts("Stopping tracks, 16#circ", "HALFCONTAINED_TRACKS_POINT_1", 1, true,
    is_ddactivity1?extbg.at("HALFCONTAINED_TRACKS_POINT_1_DDACTIVITY1_P0"): 0,
    is_ddactivity1?extbg.at("HALFCONTAINED_TRACKS_POINT_1_DDACTIVITY1_P1"): 0)
    );

  DOIT("halfcontained_tracks_point_0",
    popts("Stopping tracks, 1.3#circ", "HALFCONTAINED_TRACKS_POINT_0", 1, true,
    is_ddactivity1? extbg.at("HALFCONTAINED_TRACKS_POINT_0_DDACTIVITY1_P0"):
    is_fdminbias?   extbg.at("HALFCONTAINED_TRACKS_POINT_0_FDMINBIAS_P0"): 0,
    is_ddactivity1? extbg.at("HALFCONTAINED_TRACKS_POINT_0_DDACTIVITY1_P1"):
    is_fdminbias?   extbg.at("HALFCONTAINED_TRACKS_POINT_0_FDMINBIAS_P1"): 0)
    );

  DOIT("fullycontained_tracks",
       popts("Contained tracks", "", 0, true,
             is_ddactivity1? 2.8e-5:
             is_fdminbias? 0.45: 0));

  DOIT("fullycontained_tracks_point_1",
    popts("Contained tracks, 16#circ", "FULLYCONTAINED_TRACKS_POINT_1", 1, true,
    is_ddactivity1? extbg.at("FULLYCONTAINED_TRACKS_POINT_1_DDACTIVITY1_P0"):
    is_fdminbias?   extbg.at("FULLYCONTAINED_TRACKS_POINT_1_FDMINBIAS_P0"): 0,
    is_ddactivity1? extbg.at("FULLYCONTAINED_TRACKS_POINT_1_DDACTIVITY1_P1"):
    is_fdminbias?   extbg.at("FULLYCONTAINED_TRACKS_POINT_1_FDMINBIAS_P1"): 0)
    );

  DOIT("fullycontained_tracks_point_0",
    popts("Contained tracks, 1.3#circ", "FULLYCONTAINED_TRACKS_POINT_0",1, true,
    is_ddactivity1? extbg.at("FULLYCONTAINED_TRACKS_POINT_0_DDACTIVITY1_P0"):
    is_fdminbias?   extbg.at("FULLYCONTAINED_TRACKS_POINT_0_FDMINBIAS_P0"): 0,
    is_ddactivity1? extbg.at("FULLYCONTAINED_TRACKS_POINT_0_DDACTIVITY1_P1"):
    is_fdminbias?   extbg.at("FULLYCONTAINED_TRACKS_POINT_0_FDMINBIAS_P1"): 0)
    );

  DOIT("contained_slices",
       popts("Contained events", "", 0, true,
             is_ddactivity1? 0.00093:
             is_fdminbias? 24.085: 0));

  DOIT("upmu_tracks",
       popts("Upward going muons", "", 0, true,
             is_fdminbias? 1.98: 0));

  DOIT("upmu_tracks_point_1",
       popts("Upward going muons, 16#circ", "UPMU_TRACKS_POINT_1", 1, true,
             is_fdminbias? extbg.at("UPMU_TRACKS_POINT_1_FDMINBIAS_P0"): 0,
             is_fdminbias? extbg.at("UPMU_TRACKS_POINT_1_FDMINBIAS_P1"): 0));

  DOIT("upmu_tracks_point_0",
       popts("Upward going muons, 1.3#circ", "UPMU_TRACKS_POINT_0", 1, true,
             is_fdminbias? extbg.at("UPMU_TRACKS_POINT_0_FDMINBIAS_P0"): 0,
             is_fdminbias? extbg.at("UPMU_TRACKS_POINT_0_FDMINBIAS_P1"): 0));

  DOIT("rawtrigger",
       popts("Raw triggers", "",         0, true));
  DOIT("energy_low_cut",
       popts(">5M ADC total", "",        0, true));
  DOIT("energy_low_cut_pertime",
       popts(">2.5M ADC per 50#mus", "", 0, true));
  DOIT("energy_high_cut",
       popts(">50M ADC total", "",       0, true, 0.0010));
  DOIT("energy_high_cut_pertime",
       popts(">25M ADC per 50#mus", "",  0, true, 0.0021));
  DOIT("energy_vhigh_cut",
      popts(">500M ADC total", "",      0, true, 1e-100));
  DOIT("energy_vhigh_cut_pertime",
       popts(">250M ADC per 50#mus", "", 0, true, 1e-100));

  printend();
}
