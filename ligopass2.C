#include "TH1.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TFile.h"

int ligopass2(const char * const infilename)
{
  TFile * ligofile = new TFile(infilename, "read");
  if(!ligofile || ligofile.IsZombie()){
    fprintf(stderr, "Could not open your histogram file %s\n", infilename);
    exit(1);
  }
}
