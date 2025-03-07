//
// dltdc4_module.cxx
//
// K.Olchanski
//

#include <stdio.h>
#include <assert.h>
#include <unistd.h>

#include "manalyzer.h"
#include "midasio.h"

#include "dltdc.h"

#include "DlFlow.h"

#include "ncfm.h"

#include <deque>

#include <TStyle.h>

#ifdef HAVE_ROOT
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#endif

#if 0
static double amin(double a, double b)
{
   if (a < b)
      return a;
   else
      return b;
}

static std::string toString(int i)
{
   char buf[256];
   sprintf(buf, "%d", i);
   return buf;
}

static std::string toString(const char* fmt, double v)
{
   char buf[256];
   sprintf(buf, fmt, v);
   return buf;
}
#endif

class DlTdcFlags
{
public:
   bool fEnabled = false;
   bool fTriggered = false;
   bool fDebug = false;
   bool fPrint = false;
   bool fTWC = false;
};

#define NUM_TDC_CHAN (32+3)
#define MAX_TDC_CHAN (NUM_TDC_CHAN-1)

//#define CHANA 0
//#define CHANB 1
//#define CHAN1 2
//#define CHAN2 3
//#define CHAN3 4
//#define CHAN4 5
//#define CHAN5 6
//#define CHAN6 7
//#define CHAN7 8
//#define CHAN8 9
//#define CHANT 10

#define CHAN1 (fMap4->fChan1)
#define CHAN2 (fMap4->fChan2)
#define CHAN3 (fMap4->fChan3)
#define CHAN4 (fMap4->fChan4)
#define CHAN5 (fMap4->fChan5)
#define CHAN6 (fMap4->fChan6)
#define CHAN7 (fMap4->fChan7)
#define CHAN8 (fMap4->fChan8)


class DlTdcMap4
{
public:
   int fChan1 = 2;
   int fChan2 = 3;
   int fChan3 = 4;
   int fChan4 = 5;
   int fChan5 = 6;
   int fChan6 = 7;
   int fChan7 = 8;
   int fChan8 = 9;
   int fChanA = 0;
   int fChanB = 1;
   int fChanT = 10;

public:
   bool Init(int runno);
};

bool DlTdcMap4::Init(int runno)
{
   if (runno >= 17 && runno < 900000) {
      printf("DlTdcMap4 for run %d!\n", runno);
      fChan1 = 16;
      fChan2 = 17;
      fChan3 = 22;
      fChan4 = 23;
      fChan5 = 26;
      fChan6 = 27;
      fChan7 = 30;
      fChan8 = 31;
      fChanA = 32;
      fChanB = 33;
      fChanT = 34;
      return true;
   }

   return false;
};

class DlTdc4Module: public TARunObject
{
public:
   DlTdcFlags* fFlags = NULL;
   
#ifdef HAVE_ROOT
#if 0
   TH1D* fHphaseLe[MAX_TDC_CHAN+1];
   TH1D* fHphaseTe[MAX_TDC_CHAN+1];

   TH1D* fHfineLe[MAX_TDC_CHAN+1];
   TH1D* fHfineTe[MAX_TDC_CHAN+1];

   TH1D* fHwidth[MAX_TDC_CHAN+1];

   TH1D* fHpulserLeAll = NULL;
   TH1D* fHpulserTeAll = NULL;
   TH1D* fHpulserWiAll = NULL;
   TH1D* fHpulserLe[MAX_TDC_CHAN+1];
   TH1D* fHpulserTe[MAX_TDC_CHAN+1];
   TH1D* fHpulserWi[MAX_TDC_CHAN+1];

   TH1D* fHhitdt1ns = NULL;
   TH1D* fHhitdt2ns = NULL;
   TH1D* fHhitdt3ns = NULL;
   TH1D* fHhitdt4ns = NULL;

   TH1D* fHeventdt1ns = NULL;
   TH1D* fHeventdt2ns = NULL;
   TH1D* fHeventdt3ns = NULL;
   TH1D* fHeventdt4ns = NULL;

   TH1D* fHdt1Le[MAX_TDC_CHAN+1];
   TH1D* fHdt2Le[MAX_TDC_CHAN+1];
   TH1D* fHdt3Le[MAX_TDC_CHAN+1];
   TH1D* fHdt4Le[MAX_TDC_CHAN+1];

   std::vector<TH1D*> fHunphysical_pair_ns;
#endif

   TH1D* fHtABns = NULL;

   TH1D* fHt12ns = NULL;
   TH1D* fHt34ns = NULL;
   TH1D* fHt56ns = NULL;
   TH1D* fHt78ns = NULL;

   TH1D* fHt14ns = NULL;
   TH1D* fHt23ns = NULL;
   TH1D* fHt58ns = NULL;
   TH1D* fHt67ns = NULL;

   TH1D* fHt12ns_cut = NULL;
   TH1D* fHt34ns_cut = NULL;
   TH1D* fHt56ns_cut = NULL;
   TH1D* fHt78ns_cut = NULL;

   TH1D* fHt14ns_cut = NULL;
   TH1D* fHt23ns_cut = NULL;
   TH1D* fHt58ns_cut = NULL;
   TH1D* fHt67ns_cut = NULL;

   TH1D* fHt14ns_twc = NULL;
   TH1D* fHt23ns_twc = NULL;
   TH1D* fHt58ns_twc = NULL;
   TH1D* fHt67ns_twc = NULL;

   TH1D* fHw1ns = NULL;
   TH1D* fHw2ns = NULL;
   TH1D* fHw3ns = NULL;
   TH1D* fHw4ns = NULL;
   TH1D* fHw5ns = NULL;
   TH1D* fHw6ns = NULL;
   TH1D* fHw7ns = NULL;
   TH1D* fHw8ns = NULL;

   TH1D* fHw1ns_cut = NULL;
   TH1D* fHw2ns_cut = NULL;
   TH1D* fHw3ns_cut = NULL;
   TH1D* fHw4ns_cut = NULL;
   TH1D* fHw5ns_cut = NULL;
   TH1D* fHw6ns_cut = NULL;
   TH1D* fHw7ns_cut = NULL;
   TH1D* fHw8ns_cut = NULL;

   TH1D* fHwAns = NULL;
   TH1D* fHwBns = NULL;
   TH1D* fHwTns = NULL;

   TH2D* fHw14ns = NULL;
   TH2D* fHw23ns = NULL;
   TH2D* fHw58ns = NULL;
   TH2D* fHw67ns = NULL;

   TH2D* fHt14ns_w1ns = NULL;
   TH2D* fHt14ns_w4ns = NULL;
   TH2D* fHt23ns_w2ns = NULL;
   TH2D* fHt23ns_w3ns = NULL;
   TH2D* fHt58ns_w5ns = NULL;
   TH2D* fHt58ns_w8ns = NULL;
   TH2D* fHt67ns_w6ns = NULL;
   TH2D* fHt67ns_w7ns = NULL;

   TH2D* fHt14ns_w1ns_cutw4 = NULL;
   TH2D* fHt14ns_w4ns_cutw1 = NULL;
   TH2D* fHt23ns_w2ns_cutw3 = NULL;
   TH2D* fHt23ns_w3ns_cutw2 = NULL;
   TH2D* fHt58ns_w5ns_cutw8 = NULL;
   TH2D* fHt58ns_w8ns_cutw5 = NULL;
   TH2D* fHt67ns_w6ns_cutw7 = NULL;
   TH2D* fHt67ns_w7ns_cutw6 = NULL;

   TH2D* fHt14ns_w1ns_cutw4_twc = NULL;
   TH2D* fHt14ns_w4ns_cutw1_twc = NULL;
   TH2D* fHt23ns_w2ns_cutw3_twc = NULL;
   TH2D* fHt23ns_w3ns_cutw2_twc = NULL;
   TH2D* fHt58ns_w5ns_cutw8_twc = NULL;
   TH2D* fHt58ns_w8ns_cutw5_twc = NULL;
   TH2D* fHt67ns_w6ns_cutw7_twc = NULL;
   TH2D* fHt67ns_w7ns_cutw6_twc = NULL;
   
   // Gabby
   // Initialize historgrams
   // Time difference between t6 (t7 next) and average of t1 and t4 vs width 6 (width 7 next)
   TH2D* fHt146ns_w6ns_cut = NULL;
   TH2D* fHt147ns_w7ns_cut= NULL;

   // QUAD 14*58

   TH1D* fHtof_1458   = NULL;
   TH1D* fHtof_1458_cut   = NULL;
   TH1D* fHtof_1458_cut_twc = NULL;

   TH1D* fHt14ns_1458 = NULL;
   TH1D* fHt58ns_1458 = NULL;
   TH2D* fHt14t58ns_1458 = NULL;
   TH2D* fHt14t58ns_1458_cut = NULL;
   TH2D* fHt14t58ns_1458_cut_twc = NULL;

   TH1D* fHt158ns_1458_twc = NULL;
   TH1D* fHt458ns_1458_twc = NULL;
   TH1D* fHt145ns_1458_twc = NULL;
   TH1D* fHt148ns_1458_twc = NULL;

   TH1D* fHt158ns_1458_twc_cut = NULL;
   TH1D* fHt458ns_1458_twc_cut = NULL;
   TH1D* fHt145ns_1458_twc_cut = NULL;
   TH1D* fHt148ns_1458_twc_cut = NULL;

   TH2D* fHw14ns_1458 = NULL;
   TH2D* fHw58ns_1458 = NULL;

   TH2D* fHtof_w1ns_1458 = NULL;
   TH2D* fHtof_w4ns_1458 = NULL;
   TH2D* fHtof_w5ns_1458 = NULL;
   TH2D* fHtof_w8ns_1458 = NULL;

   TH2D* fHtof_w1ns_1458_cut = NULL;
   TH2D* fHtof_w4ns_1458_cut = NULL;
   TH2D* fHtof_w5ns_1458_cut = NULL;
   TH2D* fHtof_w8ns_1458_cut = NULL;

   TH2D* fHtof_w1ns_1458_twc = NULL;
   TH2D* fHtof_w4ns_1458_twc = NULL;
   TH2D* fHtof_w5ns_1458_twc = NULL;
   TH2D* fHtof_w8ns_1458_twc = NULL;

   TH2D* fHt158ns_w1ns_1458 = NULL;
   TH2D* fHt458ns_w4ns_1458 = NULL;
   TH2D* fHt145ns_w5ns_1458 = NULL;
   TH2D* fHt148ns_w8ns_1458 = NULL;

   TH2D* fHt158ns_w1ns_1458_twc = NULL;
   TH2D* fHt458ns_w4ns_1458_twc = NULL;
   TH2D* fHt145ns_w5ns_1458_twc = NULL;
   TH2D* fHt148ns_w8ns_1458_twc = NULL;
   
   // QUAD 23*58

   TH1D* fHtof_2358   = NULL;
   TH1D* fHtof_2358_cut   = NULL;
   TH1D* fHtof_2358_cut_twc = NULL;
   
   TH1D* fHt23ns_2358 = NULL;
   TH1D* fHt58ns_2358 = NULL;
   TH2D* fHt23t58ns_2358 = NULL;
   TH2D* fHt23t58ns_2358_cut = NULL;
   TH2D* fHt23t58ns_2358_cut_twc = NULL;

   TH2D* fHw23ns_2358 = NULL;
   TH2D* fHw58ns_2358 = NULL;

   TH2D* fHtof_w2ns_2358_twc = NULL;
   TH2D* fHtof_w3ns_2358_twc = NULL;
   TH2D* fHtof_w5ns_2358_twc = NULL;
   TH2D* fHtof_w8ns_2358_twc = NULL;

   // QUAD 14*67

   TH1D* fHtof_1467   = NULL;
   TH1D* fHtof_1467_cut   = NULL;
   TH1D* fHtof_1467_cut_twc = NULL;
   
   TH1D* fHt14ns_1467 = NULL;
   TH1D* fHt67ns_1467 = NULL;
   TH2D* fHt14t67ns_1467 = NULL;
   TH2D* fHt14t67ns_1467_cut = NULL;
   TH2D* fHt14t67ns_1467_cut_twc = NULL;

   TH2D* fHw14ns_1467 = NULL;
   TH2D* fHw67ns_1467 = NULL;

   TH2D* fHt67ns_w6ns_1467 = NULL;
   TH2D* fHt67ns_w7ns_1467 = NULL;

   TH2D* fHtof_w1ns_1467_twc = NULL;
   TH2D* fHtof_w4ns_1467_twc = NULL;
   TH2D* fHtof_w6ns_1467_twc = NULL;
   TH2D* fHtof_w7ns_1467_twc = NULL;
   
   // QUAD 23*67

   TH1D* fHtof_2367   = NULL;
   TH1D* fHtof_2367_cut   = NULL;
   TH1D* fHtof_2367_cut_twc = NULL;
   
   TH1D* fHt23ns_2367 = NULL;
   TH1D* fHt67ns_2367 = NULL;
   TH2D* fHt23t67ns_2367 = NULL;
   TH2D* fHt23t67ns_2367_cut = NULL;
   TH2D* fHt23t67ns_2367_cut_twc = NULL;

   TH1D* fHt267ns_2367_twc = NULL;
   TH1D* fHt367ns_2367_twc = NULL;
   TH1D* fHt236ns_2367_twc = NULL;
   TH1D* fHt237ns_2367_twc = NULL;

   TH1D* fHt267ns_2367_twc_cut = NULL;
   TH1D* fHt367ns_2367_twc_cut = NULL;
   TH1D* fHt236ns_2367_twc_cut = NULL;
   TH1D* fHt237ns_2367_twc_cut = NULL;

   TH2D* fHw23ns_2367 = NULL;
   TH2D* fHw67ns_2367 = NULL;

   TH2D* fHtof_w2ns_2367_twc = NULL;
   TH2D* fHtof_w3ns_2367_twc = NULL;
   TH2D* fHtof_w6ns_2367_twc = NULL;
   TH2D* fHtof_w7ns_2367_twc = NULL;

   TH2D* fHt267ns_w2ns_2367 = NULL;
   TH2D* fHt367ns_w3ns_2367 = NULL;
   TH2D* fHt236ns_w6ns_2367 = NULL;
   TH2D* fHt237ns_w7ns_2367 = NULL;

   TH2D* fHt267ns_w2ns_2367_twc = NULL;
   TH2D* fHt367ns_w3ns_2367_twc = NULL;
   TH2D* fHt236ns_w6ns_2367_twc = NULL;
   TH2D* fHt237ns_w7ns_2367_twc = NULL;

   int fCount14 = 0;
   int fCount23 = 0;
   int fCount58 = 0;
   int fCount67 = 0;

   int fCount1458 = 0;
   int fCount2358 = 0;
   int fCount1467 = 0;
   int fCount2367 = 0;

   int fCountMyA = 0;
   int fCountMyB = 0;
   int fCountMyT = 0;

   int fCountA = 0;
   int fCountB = 0;
   int fCountT = 0;

#endif

   int fCountCut2367 = 0;
   int fCountCut1 = 0;
   int fCountCut4 = 0;
   int fCountCut5 = 0;
   int fCountCut8 = 0;

   double fPrevEventTimeSec = 0;

   DlTdcMap4 *fMap4 = NULL;

   bool fTrace = false;

   Ncfm* fCfm = NULL;
   
   DlTdc4Module(TARunInfo* runinfo, DlTdcFlags* flags)
      : TARunObject(runinfo)
   {
      if (fTrace)
         printf("DlTdc4Module::ctor!\n");

      fModuleName = "dltdc4_module";
      fFlags   = flags;

      fCfm = new Ncfm("dlcfmdb");

      fMap4 = new DlTdcMap4();
   }

   ~DlTdc4Module()
   {
      if (fTrace)
         printf("DlTdc4Module::dtor!\n");

      if (fMap4) {
         delete fMap4;
         fMap4 = NULL;
      }

      if (fCfm) {
         delete fCfm;
         fCfm = NULL;
      }
   }

   void BeginRun(TARunInfo* runinfo)
   {
      if (fTrace)
         printf("DlTdcModule::BeginRun, run %d, file %s\n", runinfo->fRunNo, runinfo->fFileName.c_str());

      if (!fFlags->fEnabled)
         return;

      bool conf_ok = fMap4->Init(runinfo->fRunNo);
      if (!conf_ok) {
         printf("Cannot load TDC map for run %d\n", runinfo->fRunNo);
         exit(123);
      }

#ifdef HAVE_ROOT
      runinfo->fRoot->fOutputFile->cd(); // select correct ROOT directory
      TDirectory* dir = gDirectory->mkdir("dltdc4");
      dir->cd(); // select correct ROOT directory

      fHtABns = new TH1D("tABns", "time between 1234 trigger and 4567 trigger, A-B, ns", 200, -10, 10);

      fHt12ns = new TH1D("t12ns", "sipm board 1, t2-t1, ns", 200, -10, 10);
      fHt34ns = new TH1D("t34ns", "sipm board 2, t4-t3, ns", 200, -10, 10);
      fHt56ns = new TH1D("t56ns", "sipm board 3, t6-t5, ns", 200, -10, 10);
      fHt78ns = new TH1D("t78ns", "sipm board 4, t8-t7, ns", 200, -10, 10);

      fHt14ns = new TH1D("t14ns", "Paddle 1 time difference, t4-t1 (ns)", 200, -10, 10);
      fHt23ns = new TH1D("t23ns", "Paddle 2 time difference, t3-t2 (ns)", 200, -10, 10);
      fHt58ns = new TH1D("t58ns", "Paddle 3 time difference, t8-t5 (ns)", 200, -10, 10);
      fHt67ns = new TH1D("t67ns", "Paddle 4 time difference, t7-t6 (ns)", 200, -10, 10);

      fHt12ns_cut = new TH1D("t12ns_cut", "sipm board 1, t2-t1, ns with cut", 200, -10, 10);
      fHt34ns_cut = new TH1D("t34ns_cut", "sipm board 2, t4-t3, ns with cut", 200, -10, 10);
      fHt56ns_cut = new TH1D("t56ns_cut", "sipm board 3, t6-t5, ns with cut", 200, -10, 10);
      fHt78ns_cut = new TH1D("t78ns_cut", "sipm board 4, t8-t7, ns with cut", 200, -10, 10);

      fHt14ns_cut = new TH1D("t14ns_cut", "Paddle 1 time difference, t4-t1 (ns) with cut", 200, -10, 10);
      fHt23ns_cut = new TH1D("t23ns_cut", "Paddle 2 time difference, t3-t2 (ns) with cut", 200, -10, 10);
      fHt58ns_cut = new TH1D("t58ns_cut", "Paddle 3 time difference, t8-t5 (ns) with cut", 200, -10, 10);
      fHt67ns_cut = new TH1D("t67ns_cut", "Paddle 4 time difference, t7-t6 (ns) with cut", 200, -10, 10);

      fHt14ns_twc = new TH1D("t14ns_twc", "Paddle 1 time difference, t4-t1 (ns) with twc", 200, -10, 10);
      fHt23ns_twc = new TH1D("t23ns_twc", "Paddle 2 time difference, t3-t2 (ns) with twc", 200, -10, 10);
      fHt58ns_twc = new TH1D("t58ns_twc", "Paddle 3 time difference, t8-t5 (ns) with twc", 200, -10, 10);
      fHt67ns_twc = new TH1D("t67ns_twc", "Paddle 4 time difference, t7-t6 (ns) with twc", 200, -10, 10);

      fHw1ns = new TH1D("w1ns", "w1ns", 100, 0, 100);
      fHw2ns = new TH1D("w2ns", "w2ns", 100, 0, 100);
      fHw3ns = new TH1D("w3ns", "w3ns", 100, 0, 100);
      fHw4ns = new TH1D("w4ns", "w4ns", 100, 0, 100);
      fHw5ns = new TH1D("w5ns", "w5ns", 100, 0, 100);
      fHw6ns = new TH1D("w6ns", "w6ns", 100, 0, 100);
      fHw7ns = new TH1D("w7ns", "w7ns", 100, 0, 100);
      fHw8ns = new TH1D("w8ns", "w8ns", 100, 0, 100);

      fHw1ns_cut = new TH1D("w1ns_cut", "w1ns with cut", 100, 0, 100);
      fHw2ns_cut = new TH1D("w2ns_cut", "w2ns with cut", 100, 0, 100);
      fHw3ns_cut = new TH1D("w3ns_cut", "w3ns with cut", 100, 0, 100);
      fHw4ns_cut = new TH1D("w4ns_cut", "w4ns with cut", 100, 0, 100);
      fHw5ns_cut = new TH1D("w5ns_cut", "w5ns with cut", 100, 0, 100);
      fHw6ns_cut = new TH1D("w6ns_cut", "w6ns with cut", 100, 0, 100);
      fHw7ns_cut = new TH1D("w7ns_cut", "w7ns with cut", 100, 0, 100);
      fHw8ns_cut = new TH1D("w8ns_cut", "w8ns with cut", 100, 0, 100);

      fHwAns = new TH1D("wAns", "wAns", 100, 0, 100);
      fHwBns = new TH1D("wBns", "wBns", 100, 0, 100);
      fHwTns = new TH1D("wTns", "wTns", 100, 0, 100);

      fHw14ns = new TH2D("w14ns", "w4ns vs w1ns", 100, 0, 100, 100, 0, 100);
      fHw23ns = new TH2D("w23ns", "w3ns vs w2ns", 100, 0, 100, 100, 0, 100);
      fHw58ns = new TH2D("w58ns", "w8ns vs w5ns", 100, 0, 100, 100, 0, 100);
      fHw67ns = new TH2D("w67ns", "w7ns vs w6ns", 100, 0, 100, 100, 0, 100);

      fHt14ns_w1ns = new TH2D("t14ns_w1ns", "t4-t1 (ns) vs w1 (ns)", 100, 0, 100, 200, -10, 10);
      fHt14ns_w4ns = new TH2D("t14ns_w4ns", "t4-t1 (ns) vs w4 (ns)", 100, 0, 100, 200, -10, 10);
      fHt23ns_w2ns = new TH2D("t23ns_w2ns", "t3-t2 (ns) vs w2 (ns)", 100, 0, 100, 200, -10, 10);
      fHt23ns_w3ns = new TH2D("t23ns_w3ns", "t3-t2 (ns) vs w3 (ns)", 100, 0, 100, 200, -10, 10);
      fHt58ns_w5ns = new TH2D("t58ns_w5ns", "t8-t5 (ns) vs w5 (ns)", 100, 0, 100, 200, -10, 10);
      fHt58ns_w8ns = new TH2D("t58ns_w8ns", "t8-t5 (ns) vs w8 (ns)", 100, 0, 100, 200, -10, 10);
      fHt67ns_w6ns = new TH2D("t67ns_w6ns", "t7-t6 (ns) vs w6 (ns)", 100, 0, 100, 200, -10, 10);
      fHt67ns_w7ns = new TH2D("t67ns_w7ns", "t7-t6 (ns) vs w7 (ns)", 100, 0, 100, 200, -10, 10);

      fHt14ns_w1ns_cutw4 = new TH2D("t14ns_w1ns_cutw4", "t4-t1 (ns) vs w1 (ns) cut on w4", 100, 0, 100, 200, -10, 10);
      fHt14ns_w4ns_cutw1 = new TH2D("t14ns_w4ns_cutw1", "t4-t1 (ns) vs w4 (ns) cut on w1", 100, 0, 100, 200, -10, 10);
      fHt23ns_w2ns_cutw3 = new TH2D("t23ns_w2ns_cutw3", "t3-t2 (ns) vs w2 (ns) cut on w3", 100, 0, 100, 200, -10, 10);
      fHt23ns_w3ns_cutw2 = new TH2D("t23ns_w3ns_cutw2", "t3-t2 (ns) vs w3 (ns) cut on w2", 100, 0, 100, 200, -10, 10);
      fHt58ns_w5ns_cutw8 = new TH2D("t58ns_w5ns_cutw8", "t8-t5 (ns) vs w5 (ns) cut on w8", 100, 0, 100, 200, -10, 10);
      fHt58ns_w8ns_cutw5 = new TH2D("t58ns_w8ns_cutw5", "t8-t5 (ns) vs w8 (ns) cut on w5", 100, 0, 100, 200, -10, 10);
      fHt67ns_w6ns_cutw7 = new TH2D("t67ns_w6ns_cutw7", "t7-t6 (ns) vs w6 (ns) cut on w7", 100, 0, 100, 200, -10, 10);
      fHt67ns_w7ns_cutw6 = new TH2D("t67ns_w7ns_cutw6", "t7-t6 (ns) vs w7 (ns) cut on w6", 100, 0, 100, 200, -10, 10);

      fHt14ns_w1ns_cutw4_twc = new TH2D("t14ns_w1ns_cutw4_twc", "t4-t1 (ns) vs w1 (ns) cut on w4 w/ twc", 100, 0, 100, 200, -10, 10);
      fHt14ns_w4ns_cutw1_twc = new TH2D("t14ns_w4ns_cutw1_twc", "t4-t1 (ns) vs w4 (ns) cut on w1 w/ twc", 100, 0, 100, 200, -10, 10);
      fHt23ns_w2ns_cutw3_twc = new TH2D("t23ns_w2ns_cutw3_twc", "t3-t2 (ns) vs w2 (ns) cut on w3 w/ twc", 100, 0, 100, 200, -10, 10);
      fHt23ns_w3ns_cutw2_twc = new TH2D("t23ns_w3ns_cutw2_twc", "t3-t2 (ns) vs w3 (ns) cut on w2 w/ twc", 100, 0, 100, 200, -10, 10);
      fHt58ns_w5ns_cutw8_twc = new TH2D("t58ns_w5ns_cutw8_twc", "t8-t5 (ns) vs w5 (ns) cut on w8 w/ twc", 100, 0, 100, 200, -10, 10);
      fHt58ns_w8ns_cutw5_twc = new TH2D("t58ns_w8ns_cutw5_twc", "t8-t5 (ns) vs w8 (ns) cut on w5 w/ twc", 100, 0, 100, 200, -10, 10);
      fHt67ns_w6ns_cutw7_twc = new TH2D("t67ns_w6ns_cutw7_twc", "t7-t6 (ns) vs w6 (ns) cut on w7 w/ twc", 100, 0, 100, 200, -10, 10);
      fHt67ns_w7ns_cutw6_twc = new TH2D("t67ns_w7ns_cutw6_twc", "t7-t6 (ns) vs w7 (ns) cut on w6 w/ twc", 100, 0, 100, 200, -10, 10);

      // Gabby
      // Make the form of the histogram to store the data in (file name, graph title, axis and bin information)
      fHt146ns_w6ns_cut = new TH2D("t146ns_w6ns_cut", "t6-0.5(t1+t4) (ns) vs w6 (ns) w/  cut", 100, 0, 100, 200, -10, 10);
      fHt147ns_w7ns_cut = new TH2D("t147ns_w7ns_cut", "t7-0.5(t1+t4) (ns) vs w7 (ns) w/ cut", 100, 0, 100, 200, -10, 10);

      // QUAD 14*58

      fHtof_1458 = new TH1D("tof_1458", "TOF 14 vs 58 (ns)", 200, -10, 10);
      fHtof_1458_cut = new TH1D("tof_1458_cut", "TOF 14 vs 58 (ns) with cut", 200, -10, 10);
      fHtof_1458_cut_twc = new TH1D("tof_1458_cut_twc", "TOF 14 vs 58 (ns) with cut and twc", 200, -10, 10);

      fHt14ns_1458 = new TH1D("t14ns_1458", "Paddle 1 time difference, t4-t1 (ns), 1*4*5*8", 200, -10, 10);
      fHt58ns_1458 = new TH1D("t58ns_1458", "Paddle X time difference, t8-t5 (ns), 1*4*5*8", 200, -10, 10);
      fHt14t58ns_1458 = new TH2D("t14t58ns_1458", "t14 vs t58 (ns), 1*4*5*8", 200, -10, 10, 200, -10, 10);
      fHt14t58ns_1458_cut = new TH2D("t14t58ns_1458_cut", "t14 vs t58 (ns), 1*4*5*8 with cut", 200, -10, 10, 200, -10, 10);
      fHt14t58ns_1458_cut_twc = new TH2D("t14t58ns_1458_cut_twc", "t14 vs t58 (ns), 1*4*5*8 with cut and twc", 200, -10, 10, 200, -10, 10);
      //fHt14ns_1458_cut = new TH1D("t14ns_1458_cut", "Paddle 1 time difference, t4-t1 (ns), 1*4*5*8 with cut", 200, -10, 10);
      //fHt58ns_1458_cut = new TH1D("t58ns_1458_cut", "Paddle X time difference, t8-t5 (ns), 1*4*5*8 with cut", 200, -10, 10);

      fHt158ns_1458_twc = new TH1D("t158ns_1458_twc", "t1 - 0.5*(t5 + t8) (ns), CH1 minus lower paddle w/ twc", 200, -10, 10);
      fHt458ns_1458_twc = new TH1D("t458ns_1458_twc", "t4 - 0.5*(t5 + t8) (ns), CH4 minus lower paddle w/ twc", 200, -10, 10);
      fHt145ns_1458_twc = new TH1D("t145ns_1458_twc", "t5 - 0.5*(t1 + t4) (ns), CH5 minus upper paddle w/ twc", 200, -10, 10);
      fHt148ns_1458_twc = new TH1D("t148ns_1458_twc", "t8 - 0.5*(t1 + t4) (ns), CH8 minus upper paddle w/ twc", 200, -10, 10);

      fHt158ns_1458_twc_cut = new TH1D("t158ns_1458_twc_cut", "t1 - 0.5*(t5 + t8) (ns), CH1 minus lower paddle w/ twc and cut", 200, -10, 10);
      fHt458ns_1458_twc_cut = new TH1D("t458ns_1458_twc_cut", "t4 - 0.5*(t5 + t8) (ns), CH4 minus lower paddle w/ twc and cut", 200, -10, 10);
      fHt145ns_1458_twc_cut = new TH1D("t145ns_1458_twc_cut", "t5 - 0.5*(t1 + t4) (ns), CH5 minus upper paddle w/ twc and cut", 200, -10, 10);
      fHt148ns_1458_twc_cut = new TH1D("t148ns_1458_twc_cut", "t8 - 0.5*(t1 + t4) (ns), CH8 minus upper paddle w/ twc and cut", 200, -10, 10); 

      fHw14ns_1458 = new TH2D("w14ns_1458", "w4ns vs w1ns, 1*4*5*8", 100, 0, 100, 100, 0, 100);
      fHw58ns_1458 = new TH2D("w58ns_1458", "w8ns vs w5ns, 1*4*5*8", 100, 0, 100, 100, 0, 100);

      fHtof_w1ns_1458 = new TH2D("tof_w1ns_1458", "tof vs w1 (ns), 1*4*5*8", 100, 0, 100, 100, -10, 10);
      fHtof_w4ns_1458 = new TH2D("tof_w4ns_1458", "tof vs w4 (ns), 1*4*5*8", 100, 0, 100, 100, -10, 10);
      fHtof_w5ns_1458 = new TH2D("tof_w5ns_1458", "tof vs w5 (ns), 1*4*5*8", 100, 0, 100, 100, -10, 10);
      fHtof_w8ns_1458 = new TH2D("tof_w8ns_1458", "tof vs w8 (ns), 1*4*5*8", 100, 0, 100, 100, -10, 10);

      fHtof_w1ns_1458_cut = new TH2D("tof_w1ns_1458_cut", "tof vs w1 (ns), 1*4*5*8 with cut", 100, 0, 100, 100, -10, 10);
      fHtof_w4ns_1458_cut = new TH2D("tof_w4ns_1458_cut", "tof vs w4 (ns), 1*4*5*8 with cut", 100, 0, 100, 100, -10, 10);
      fHtof_w5ns_1458_cut = new TH2D("tof_w5ns_1458_cut", "tof vs w5 (ns), 1*4*5*8 with cut", 100, 0, 100, 100, -10, 10);
      fHtof_w8ns_1458_cut = new TH2D("tof_w8ns_1458_cut", "tof vs w8 (ns), 1*4*5*8 with cut", 100, 0, 100, 100, -10, 10);

      fHtof_w1ns_1458_twc = new TH2D("tof_w1ns_1458_twc", "tof vs w1 (ns), 1*4*5*8 with twc", 100, 0, 100, 100, -10, 10);
      fHtof_w4ns_1458_twc = new TH2D("tof_w4ns_1458_twc", "tof vs w4 (ns), 1*4*5*8 with twc", 100, 0, 100, 100, -10, 10);
      fHtof_w5ns_1458_twc = new TH2D("tof_w5ns_1458_twc", "tof vs w5 (ns), 1*4*5*8 with twc", 100, 0, 100, 100, -10, 10);
      fHtof_w8ns_1458_twc = new TH2D("tof_w8ns_1458_twc", "tof vs w8 (ns), 1*4*5*8 with twc", 100, 0, 100, 100, -10, 10);

      fHt158ns_w1ns_1458 = new TH2D("t158ns_w1ns_1458", "t1 - 0.5*(t5 + t8) (ns) vs w1 (ns)", 100, 0, 100, 100, -10, 10);
      fHt458ns_w4ns_1458 = new TH2D("t458ns_w4ns_1458", "t4 - 0.5*(t5 + t8) (ns) vs w4 (ns)", 100, 0, 100, 100, -10, 10);
      fHt145ns_w5ns_1458 = new TH2D("t145ns_w5ns_1458", "t5 - 0.5*(t1 + t4) (ns) vs w5 (ns)", 100, 0, 100, 100, -10, 10);
      fHt148ns_w8ns_1458 = new TH2D("t148ns_w8ns_1458", "t8 - 0.5*(t1 + t4) (ns) vs w8 (ns)", 100, 0, 100, 100, -10, 10);

      fHt158ns_w1ns_1458_twc = new TH2D("t158ns_w1ns_1458_twc", "t1 - 0.5*(t5 + t8) (ns) vs w1 (ns) w/ twc", 100, 0, 100, 100, -10, 10);
      fHt458ns_w4ns_1458_twc = new TH2D("t458ns_w4ns_1458_twc", "t4 - 0.5*(t5 + t8) (ns) vs w4 (ns) w/ twc", 100, 0, 100, 100, -10, 10);
      fHt145ns_w5ns_1458_twc = new TH2D("t145ns_w5ns_1458_twc", "t5 - 0.5*(t1 + t4) (ns) vs w5 (ns) w/ twc", 100, 0, 100, 100, -10, 10);
      fHt148ns_w8ns_1458_twc = new TH2D("t148ns_w8ns_1458_twc", "t8 - 0.5*(t1 + t4) (ns) vs w8 (ns) w/ twc", 100, 0, 100, 100, -10, 10);

      // QUAD 23*58

      fHtof_2358 = new TH1D("tof_2358", "TOF 23 vs 58 (ns)", 200, -10, 10);
      fHtof_2358_cut = new TH1D("tof_2358_cut", "TOF 23 vs 58 (ns) with cut", 200, -10, 10);
      fHtof_2358_cut_twc = new TH1D("tof_2358_cut_twc", "TOF 23 vs 58 (ns) with cut and twc", 200, -10, 10);

      fHt23ns_2358 = new TH1D("t23ns_2358", "t3-t2 (ns), 2*3*5*8", 200, -10, 10);
      fHt58ns_2358 = new TH1D("t58ns_2358", "t8-t5 (ns), 2*3*5*8", 200, -10, 10);
      fHt23t58ns_2358 = new TH2D("t23t58ns_2358", "t23 vs t58 (ns), 2*3*5*8", 200, -10, 10, 200, -10, 10);
      fHt23t58ns_2358_cut = new TH2D("t23t58ns_2358_cut", "t23 vs t58 (ns), 2*3*5*8 with cut", 200, -10, 10, 200, -10, 10);
      fHt23t58ns_2358_cut_twc = new TH2D("t23t58ns_2358_cut_twc", "t23 vs t58 (ns), 2*3*5*8 with cut and twc", 200, -10, 10, 200, -10, 10);
      //fHt23ns_2358_cut = new TH1D("t23ns_2358_cut", "t3-t2 (ns), 2*3*5*8 with cut", 200, -10, 10);
      //fHt58ns_2358_cut = new TH1D("t58ns_2358_cut", "t8-t5 (ns), 2*3*5*8 with cut", 200, -10, 10);

      fHw23ns_2358 = new TH2D("w23ns_2358", "w3ns vs w2ns, 2*3*5*8", 100, 0, 100, 100, 0, 100);
      fHw58ns_2358 = new TH2D("w58ns_2358", "w8ns vs w5ns, 2*3*5*8", 100, 0, 100, 100, 0, 100);

      fHtof_w2ns_2358_twc = new TH2D("tof_w2ns_2358_twc", "tof vs w2 (ns), 2*3*5*8 with twc", 100, 0, 100, 100, -10, 10);
      fHtof_w3ns_2358_twc = new TH2D("tof_w3ns_2358_twc", "tof vs w3 (ns), 2*3*5*8 with twc", 100, 0, 100, 100, -10, 10);
      fHtof_w5ns_2358_twc = new TH2D("tof_w5ns_2358_twc", "tof vs w5 (ns), 2*3*5*8 with twc", 100, 0, 100, 100, -10, 10);
      fHtof_w8ns_2358_twc = new TH2D("tof_w8ns_2358_twc", "tof vs w8 (ns), 2*3*5*8 with twc", 100, 0, 100, 100, -10, 10);

      // QUAD 14*67

      fHtof_1467 = new TH1D("tof_1467", "TOF 14 vs 67 (ns)", 200, -10, 10);
      fHtof_1467_cut = new TH1D("tof_1467_cut", "TOF 14 vs 67 (ns) with cut", 200, -10, 10);
      fHtof_1467_cut_twc = new TH1D("tof_1467_cut_twc", "TOF 14 vs 67 (ns) with cut and twc", 200, -10, 10);

      fHt14ns_1467 = new TH1D("t14ns_1467", "Paddle 1 time difference, t4-t1 (ns), 1*4*6*7", 200, -10, 10);
      fHt67ns_1467 = new TH1D("t67ns_1467", "Paddle X time difference, t7-t6 (ns), 1*4*6*7", 200, -10, 10);
      fHt14t67ns_1467 = new TH2D("t14t67ns_1467", "t14 vs t67 (ns), 1*4*6*7", 200, -10, 10, 200, -10, 10);
      fHt14t67ns_1467_cut = new TH2D("t14t67ns_1467_cut", "t14 vs t67 (ns), 1*4*6*7 with cut", 200, -10, 10, 200, -10, 10);
      fHt14t67ns_1467_cut_twc = new TH2D("t14t67ns_1467_cut_twc", "t14 vs t67 (ns), 1*4*6*7 with cut and twc", 200, -10, 10, 200, -10, 10);
      //fHt14ns_1467_cut = new TH1D("t14ns_1467_cut", "Paddle 1 time difference, t4-t1 (ns), 1*4*6*7 with cut", 200, -10, 10);
      //fHt67ns_1467_cut = new TH1D("t67ns_1467_cut", "Paddle X time difference, t7-t6 (ns), 1*4*6*7 with cut", 200, -10, 10);

      fHw14ns_1467 = new TH2D("w14ns_1467", "w4ns vs w1ns, 1*4*6*7", 100, 0, 100, 100, 0, 100);
      fHw67ns_1467 = new TH2D("w67ns_1467", "w7ns vs w6ns, 1*4*6*7", 100, 0, 100, 100, 0, 100);

      fHt67ns_w6ns_1467 = new TH2D("t67ns_w6ns_1467", "t7-t6 (ns) vs w6 (ns), 1*4*6*7", 100, 0, 100, 200, -10, 10);
      fHt67ns_w7ns_1467 = new TH2D("t67ns_w7ns_1467", "t7-t6 (ns) vs w7 (ns), 1*4*6*7", 100, 0, 100, 200, -10, 10);

      fHtof_w1ns_1467_twc = new TH2D("tof_w1ns_1467_twc", "tof vs w1 (ns), 1*4*6*7 with twc", 100, 0, 100, 100, -10, 10);
      fHtof_w4ns_1467_twc = new TH2D("tof_w4ns_1467_twc", "tof vs w4 (ns), 1*4*6*7 with twc", 100, 0, 100, 100, -10, 10);
      fHtof_w6ns_1467_twc = new TH2D("tof_w6ns_1467_twc", "tof vs w6 (ns), 1*4*6*7 with twc", 100, 0, 100, 100, -10, 10);
      fHtof_w7ns_1467_twc = new TH2D("tof_w7ns_1467_twc", "tof vs w7 (ns), 1*4*6*7 with twc", 100, 0, 100, 100, -10, 10);

      // QUAD 23*67

      fHtof_2367 = new TH1D("tof_2367", "TOF 23 vs 67 (ns)", 200, -10, 10);
      fHtof_2367_cut = new TH1D("tof_2367_cut", "TOF 23 vs 67 (ns) with cut", 200, -10, 10);
      fHtof_2367_cut_twc = new TH1D("tof_2367_cut_twc", "TOF 23 vs 67 (ns) with cut and twc", 200, -10, 10);

      fHt23ns_2367 = new TH1D("t23ns_2367", "Paddle 1 time difference, t3-t2 (ns), 2*3*6*7", 200, -10, 10);
      fHt67ns_2367 = new TH1D("t67ns_2367", "Paddle X time difference, t7-t6 (ns), 2*3*6*7", 200, -10, 10);
      fHt23t67ns_2367 = new TH2D("t23t67ns_2367", "t23 vs t67 (ns), 2*3*6*7", 200, -10, 10, 200, -10, 10);
      fHt23t67ns_2367_cut = new TH2D("t23t67ns_2367_cut", "t23 vs t67 (ns), 2*3*6*7 with cut", 200, -10, 10, 200, -10, 10);
      fHt23t67ns_2367_cut_twc = new TH2D("t23t67ns_2367_cut_twc", "t23 vs t67 (ns), 2*3*6*7 with cut and twc", 200, -10, 10, 200, -10, 10);
      //fHt23ns_2367_cut = new TH1D("t23ns_2367_cut", "Paddle 1 time difference, t3-t2 (ns), 2*3*6*7 with cut", 200, -10, 10);
      //fHt67ns_2367_cut = new TH1D("t67ns_2367_cut", "Paddle X time difference, t7-t6 (ns), 2*3*6*7 with cut", 200, -10, 10);

      fHt267ns_2367_twc = new TH1D("t267ns_2367_twc", "t2 - 0.5*(t6 + t7) (ns), CH2 minus lower paddle w/ twc", 200, -10, 10);
      fHt367ns_2367_twc = new TH1D("t367ns_2367_twc", "t3 - 0.5*(t6 + t7) (ns), CH3 minus lower paddle w/ twc", 200, -10, 10);
      fHt236ns_2367_twc = new TH1D("t236ns_2367_twc", "t6 - 0.5*(t2 + t3) (ns), CH6 minus upper paddle w/ twc", 200, -10, 10);
      fHt237ns_2367_twc = new TH1D("t237ns_2367_twc", "t7 - 0.5*(t2 + t3) (ns), CH7 minus upper paddle w/ twc", 200, -10, 10);

      fHt267ns_2367_twc_cut = new TH1D("t267ns_2367_twc_cut", "t2 - 0.5*(t6 + t7) (ns), CH2 minus lower paddle w/ twc and cut", 200, -10, 10);
      fHt367ns_2367_twc_cut = new TH1D("t367ns_2367_twc_cut", "t3 - 0.5*(t6 + t7) (ns), CH3 minus lower paddle w/ twc and cut", 200, -10, 10);
      fHt236ns_2367_twc_cut = new TH1D("t236ns_2367_twc_cut", "t6 - 0.5*(t2 + t3) (ns), CH6 minus upper paddle w/ twc and cut", 200, -10, 10);
      fHt237ns_2367_twc_cut = new TH1D("t237ns_2367_twc_cut", "t7 - 0.5*(t2 + t3) (ns), CH7 minus upper paddle w/ twc and cut", 200, -10, 10);

      fHw23ns_2367 = new TH2D("w23ns_2367", "w3ns vs w2ns, 2*3*6*7", 100, 0, 100, 100, 0, 100);
      fHw67ns_2367 = new TH2D("w67ns_2367", "w7ns vs w6ns, 2*3*6*7", 100, 0, 100, 100, 0, 100);

      fHtof_w2ns_2367_twc = new TH2D("tof_w2ns_2367_twc", "tof vs w2 (ns), 2*3*6*7 with twc", 100, 0, 100, 100, -10, 10);
      fHtof_w3ns_2367_twc = new TH2D("tof_w3ns_2367_twc", "tof vs w3 (ns), 2*3*6*7 with twc", 100, 0, 100, 100, -10, 10);
      fHtof_w6ns_2367_twc = new TH2D("tof_w6ns_2367_twc", "tof vs w6 (ns), 2*3*6*7 with twc", 100, 0, 100, 100, -10, 10);
      fHtof_w7ns_2367_twc = new TH2D("tof_w7ns_2367_twc", "tof vs w7 (ns), 2*3*6*7 with twc", 100, 0, 100, 100, -10, 10);

      fHt267ns_w2ns_2367 = new TH2D("t267ns_w2ns_2367", "t2 - 0.5*(t6 + t7) (ns) vs w2 (ns)", 100, 0, 100, 100, -10, 10);
      fHt367ns_w3ns_2367 = new TH2D("t367ns_w3ns_2367", "t3 - 0.5*(t6 + t7) (ns) vs w3 (ns)", 100, 0, 100, 100, -10, 10);
      fHt236ns_w6ns_2367 = new TH2D("t236ns_w6ns_2367", "t6 - 0.5*(t2 + t3) (ns) vs w6 (ns)", 100, 0, 100, 100, -10, 10);
      fHt237ns_w7ns_2367 = new TH2D("t237ns_w7ns_2367", "t7 - 0.5*(t2 + t3) (ns) vs w7 (ns)", 100, 0, 100, 100, -10, 10);

      fHt267ns_w2ns_2367_twc = new TH2D("t267ns_w2ns_2367_twc", "t2 - 0.5*(t6 + t7) (ns) vs w2 (ns) w/ twc", 100, 0, 100, 100, -10, 10);
      fHt367ns_w3ns_2367_twc = new TH2D("t367ns_w3ns_2367_twc", "t3 - 0.5*(t6 + t7) (ns) vs w3 (ns) w/ twc", 100, 0, 100, 100, -10, 10);
      fHt236ns_w6ns_2367_twc = new TH2D("t236ns_w6ns_2367_twc", "t6 - 0.5*(t2 + t3) (ns) vs w6 (ns) w/ twc", 100, 0, 100, 100, -10, 10);
      fHt237ns_w7ns_2367_twc = new TH2D("t237ns_w7ns_2367_twc", "t7 - 0.5*(t2 + t3) (ns) vs w7 (ns) w/ twc", 100, 0, 100, 100, -10, 10);
#endif
   }

   void PreEndRun(TARunInfo* runinfo)
   {
      if (fTrace)
         printf("DlTdcModule::PreEndRun, run %d\n", runinfo->fRunNo);
   }
   
   void EndRun(TARunInfo* runinfo)
   {
      if (fTrace)
         printf("DlTdcModule::EndRun, run %d\n", runinfo->fRunNo);

      if (!fFlags->fEnabled)
         return;

      printf("Run %d coincidences: 1-4: %d, 2-3: %d, 5-8: %d, 6-7: %d, 14-58: %d, 23-58: %d, 14-67: %d, 23-67: %d, A: %d/%d, B: %d/%d, T: %d/%d\n",
             runinfo->fRunNo,
             fCount14,
             fCount23,
             fCount58,
             fCount67,
             fCount1458,
             fCount2358,
             fCount1467,
             fCount2367,
             fCountMyA, fCountA,
             fCountMyB, fCountB,
             fCountMyT, fCountT);

      printf("Run %d calibration: t12ns: %8.3f, t34ns: %8.3f, t56ns: %8.3f, t78ns: %8.3f, t14ns: %8.3f, t58ns: %8.3f, t1458ns: %8.3f\n",
             runinfo->fRunNo,
             fHt12ns->GetMean(),
             fHt34ns->GetMean(),
             fHt56ns->GetMean(),
             fHt78ns->GetMean(),
             fHt14ns->GetMean(),
             fHt58ns->GetMean(),
             fHtof_1458->GetMean()
             );

      printf("cuts: 2367: %d, w1 %d, w4 %d, w5 %d, w8 %d\n", fCountCut2367, fCountCut1, fCountCut4, fCountCut5, fCountCut8);
      
      if (fFlags->fTWC) {
         std::vector<TH2D *> fit_histograms (8);
         
         fit_histograms[0] = fHt158ns_w1ns_1458;
         fit_histograms[1] = fHt267ns_w2ns_2367;
         fit_histograms[2] = fHt367ns_w3ns_2367;
         fit_histograms[3] = fHt458ns_w4ns_1458;
         fit_histograms[4] = fHt145ns_w5ns_1458;
         fit_histograms[5] = fHt236ns_w6ns_2367;
         fit_histograms[6] = fHt237ns_w7ns_2367;
         fit_histograms[7] = fHt148ns_w8ns_1458;
         
         
         TF1 *f_twc = new TF1("f_twc", "[0] + [1]/sqrt(x)",0,100);
         f_twc->SetParameters(-1,10);
         f_twc->SetParNames("Offset", "W");
         
         printf("----------\n\nTWC W PARAMETERS:\n");
         for (int i=0; i<8; i++) {
            fit_histograms[i]->Fit(f_twc);
            TF1 *fit  = fit_histograms[i]->GetFunction("f_twc");
            Double_t chi2 = fit->GetChisquare();
            Double_t W = fit->GetParameter(1);
            Double_t W_err = fit->GetParError(1);
            
            printf("FOR CH%d: Chi squared %f, W = %f +/- %f\n",i+1, chi2, W, W_err);
         }
         printf("\n----------\n");
      }
   }
   
   void PauseRun(TARunInfo* runinfo)
   {
      if (fTrace)
         printf("DlTdcModule::PauseRun, run %d\n", runinfo->fRunNo);
   }

   void ResumeRun(TARunInfo* runinfo)
   {
      if (fTrace)
         printf("DlTdcModule::ResumeRun, run %d\n", runinfo->fRunNo);
   }

   void AnalyzeTdcEvent(const DlTdcEvent& t)
   {
      if (fFlags->fDebug) {
         printf("EVENT %d %d %d %d %d %d %d %d, ABT %d%d%d\n", t.HaveCh(CHAN1), t.HaveCh(CHAN2), t.HaveCh(CHAN3), t.HaveCh(CHAN4), t.HaveCh(CHAN5), t.HaveCh(CHAN6), t.HaveCh(CHAN7), t.HaveCh(CHAN8), t.HaveCh(fMap4->fChanA), t.HaveCh(fMap4->fChanB), t.HaveCh(fMap4->fChanT));
      }
      
      ///////// check for triggered event ///////////

      if (fFlags->fTriggered && !t.HaveCh(fMap4->fChanT)) {
         return;
      }

      ///////// plot trigger A-B time ///////////
      
      if (t.HaveCh(fMap4->fChanA) && t.HaveCh(fMap4->fChanB)) {
         double tAB_ns = subtract_ns(t.GetCh(fMap4->fChanA).fLe, t.GetCh(fMap4->fChanB).fLe);
         fHtABns->Fill(tAB_ns);
      }

      //if (!t.HaveCh(fMap4->fChanB)) return;

      // cut on t6-t7
      
      //if (t.HaveCh(CHAN6) && t.HaveCh(CHAN7)) {
      //   double t67_ns = subtract_ns(t.GetCh(CHAN7).fLe, t.GetCh(CHAN6).fLe);
      //
      //   //if (t67_ns > -0.6) return;
      //   if (t67_ns < -0.4) return;
      //} else {
      //   return;
      //}

      ///////// SET WIDTH CUTOFFS (NS)  /////////
      
      double w1_cut = 10;
      double w2_cut = 10;
      double w3_cut = 10;
      double w4_cut = 10;
      double w5_cut = 10;
      double w6_cut = 10;
      double w7_cut = 10;
      double w8_cut = 10;

      ///////// COMPUTE WIDTH AND PULSE HEIGHT ///////////

      double w1_ns = -9999;

      if (t.HaveCh(CHAN1)) {
         w1_ns = t.GetCh(CHAN1).fWidthNs;
         if (w1_ns < 0.01) {
            printf("WWW: BAD WIDTH chan1 %f!\n", w1_ns);
            w1_ns = -9999;
         } else {
         }
      }

      double w2_ns = -9999;

      if (t.fHits[CHAN2].fDown) {
         w2_ns = t.fHits[CHAN2].fWidthNs;
         if (w2_ns < 0.01) {
            printf("WWW: BAD WIDTH chan2 %f!\n", w2_ns);
            w2_ns = -9999;
         } else {
         }
      }

      double w3_ns = -9999;

      if (t.fHits[CHAN3].fDown) {
         w3_ns = t.fHits[CHAN3].fWidthNs;
         if (w3_ns < 0.01) {
            printf("WWW: BAD WIDTH chan3 %f!\n", w3_ns);
            w3_ns = -9999;
         } else {
         }
      }

      double w4_ns = -9999;

      if (t.fHits[CHAN4].fDown) {
         w4_ns = t.fHits[CHAN4].fWidthNs;
         if (w4_ns < 0.01) {
            printf("WWW: BAD WIDTH chan4 %f!\n", w4_ns);
            w4_ns = -9999;
         } else {
         }
      }

      double w5_ns = -9999;

      if (t.fHits[CHAN5].fDown) {
         w5_ns = t.fHits[CHAN5].fWidthNs;
         if (w5_ns < 0.01) {
            printf("WWW: BAD WIDTH chan5 %f!\n", w5_ns);
            w5_ns = -9999;
         } else {
         }

      }

      double w6_ns = -9999;

      if (t.fHits[CHAN6].fDown) {
         w6_ns = t.fHits[CHAN6].fWidthNs;
         if (w6_ns < 0.01) {
            printf("WWW: BAD WIDTH chan6 %f!\n", w6_ns);
            w6_ns = -9999;
         } else {
         }
      }

      double w7_ns = -9999;

      if (t.fHits[CHAN7].fDown) {
         w7_ns = t.fHits[CHAN7].fWidthNs;
         if (w7_ns < 0.01) {
            printf("WWW: BAD WIDTH chan7 %f!\n", w7_ns);
            w7_ns = -9999;
         } else {
         }
      }

      double w8_ns = -9999;

      if (t.fHits[CHAN8].fDown) {
         w8_ns = t.fHits[CHAN8].fWidthNs;
         if (w8_ns < 0.01) {
            printf("WWW: BAD WIDTH chan8 %f!\n", w8_ns);
            w8_ns = -9999;
         } else {
         }
      }

      // SET UP TIME WALK CORRECT VALUES //
      //double W_twc = 8.0;

      // OPTIMIZED WITH CUT
      double W1_twc = 9.32; // 8.41;
      double W2_twc = 9.89; // 8.32;
      double W3_twc = 11.07; // 8.67;
      double W4_twc = 10.54; // 9.26;
      double W5_twc = 8.36; // 7.75;
      double W6_twc = 9.39; // 8.73;
      double W7_twc = 9.01; // 8.15;
      double W8_twc = 9.11; // 8.05;

      // UNIFORM FOR GENERAL USE
      //double W1_twc = W_twc;
      //double W2_twc = W_twc;
      //double W3_twc = W_twc;
      //double W4_twc = W_twc;
      //double W5_twc = W_twc;
      //double W6_twc = W_twc;
      //double W7_twc = W_twc;
      //double W8_twc = W_twc;

      ///////// TRIGGER CHANNALS A, B and T not used yet ///////////

      if (t.fHits[fMap4->fChanA].fDown) {
         double wA_ns = t.fHits[fMap4->fChanA].fWidthNs;

         if (wA_ns < 0.01) {
            printf("WWW: BAD WIDTH chanA %f!\n", wA_ns);
         }

         if (fFlags->fPrint) {
            printf("new dlsc event A, le %.9f, w %.0f!\n", t.GetCh(fMap4->fChanA).fLe.time_sec, wA_ns);
         }

         fHwAns->Fill(wA_ns);

         fCountA++;
      }

      if (t.fHits[fMap4->fChanB].fDown) {
         double wB_ns = t.fHits[fMap4->fChanB].fWidthNs;

         if (wB_ns < 0.01) {
            printf("WWW: BAD WIDTH chanB %f!\n", wB_ns);
         }

         if (fFlags->fPrint) {
            printf("new dlsc event B, le %.9f, w %.0f!\n", t.GetCh(fMap4->fChanB).fLe.time_sec, wB_ns);
         }

         fHwBns->Fill(wB_ns);

         fCountB++;
      }

      if (t.fHits[fMap4->fChanT].fDown) {
         double wT_ns = t.fHits[fMap4->fChanT].fWidthNs;

         if (wT_ns < 0.01) {
            printf("WWW: BAD WIDTH chanT %f!\n", wT_ns);
         }

         if (fFlags->fPrint) {
            printf("new dlsc event T, le %.9f, w %.0f!\n", t.GetCh(fMap4->fChanT).fLe.time_sec, wT_ns);
         }

         fHwTns->Fill(wT_ns);
         
         fCountT++;
      }

#if 0
      bool cut = false;

      if (w2_ns>0&&w3_ns>0&&w6_ns>0&&w7_ns>0) {
         bool flag = false;
         if (w1_ns>0) { fCountCut1++; flag = true; }
         if (w4_ns>0) { fCountCut4++; flag = true; }
         if (w5_ns>0) { fCountCut5++; flag = true; }
         if (w8_ns>0) { fCountCut8++; flag = true; }
         if (flag)
            cut = true;
      } else if (w1_ns>0&&w4_ns>0&&w6_ns>0&&w7_ns>0) {
         bool flag = false;
         if (w2_ns>0) { flag = true; }
         if (w3_ns>0) { flag = true; }
         if (w5_ns>0) { flag = true; }
         if (w8_ns>0) { flag = true; }
         if (flag)
            cut = true;
      } else {
         //fCountCut2367++;
         cut = true;
      }

      if (cut) return;
      //if (!cut) return;
#endif

      ///////// SINGLES ///////////

      if (w1_ns > 0) {
         if (fFlags->fPrint) {
            printf("new dlsc event 1, le %.9f sec, w1 %.0f ns!\n", t.GetCh(CHAN1).fLe.time_sec, w1_ns);
         }

         fHw1ns->Fill(w1_ns);

         if (w1_ns > w1_cut) { 
            fHw1ns_cut->Fill(w1_ns);
         }
      }

      if (w2_ns > 0) {
         if (fFlags->fPrint) {
            printf("new dlsc event 2, le %.9f sec, w2 %.0f ns!\n", t.GetCh(CHAN2).fLe.time_sec, w2_ns);
         }

         fHw2ns->Fill(w2_ns);

         if (w2_ns > w2_cut) { 
            fHw2ns_cut->Fill(w2_ns);
         }
      }

      if (w3_ns > 0) {
         if (fFlags->fPrint) {
            printf("new dlsc event 3, le %.9f sec, w3 %.0f ns!\n", t.GetCh(CHAN3).fLe.time_sec, w3_ns);
         }

         fHw3ns->Fill(w3_ns);

         if (w3_ns > w3_cut) { 
            fHw3ns_cut->Fill(w3_ns);
         }
      }

      if (w4_ns > 0) {
         if (fFlags->fPrint) {
            printf("new dlsc event 4, le %.9f sec, w4 %.0f ns!\n", t.GetCh(CHAN4).fLe.time_sec, w4_ns);
         }

         fHw4ns->Fill(w4_ns);

         if (w4_ns > w4_cut) { 
            fHw4ns_cut->Fill(w4_ns);
         }
      }

      if (w5_ns > 0) {
         if (fFlags->fPrint) {
            printf("new dlsc event 5, le %.9f sec, w5 %.0f ns!\n", t.GetCh(CHAN5).fLe.time_sec, w5_ns);
         }

         fHw5ns->Fill(w5_ns);

         if (w5_ns > w5_cut) { 
            fHw5ns_cut->Fill(w5_ns);
         }
      }

      if (w6_ns > 0) {
         if (fFlags->fPrint) {
            printf("new dlsc event 6, le %.9f sec, w6 %.0f ns!\n", t.GetCh(CHAN6).fLe.time_sec, w6_ns);
         }

         fHw6ns->Fill(w6_ns);

         if (w6_ns > w6_cut) { 
            fHw6ns_cut->Fill(w6_ns);
         }
      }

      if (w7_ns > 0) {
         if (fFlags->fPrint) {
            printf("new dlsc event 7, le %.9f sec, w7 %.0f ns!\n", t.GetCh(CHAN7).fLe.time_sec, w7_ns);
         }

         fHw7ns->Fill(w7_ns);

         if (w7_ns > w7_cut) { 
            fHw7ns_cut->Fill(w7_ns);
         }
      }

      if (w8_ns > 0) {
         if (fFlags->fPrint) {
            printf("new dlsc event 8, le %.9f sec, w8 %.0f ns!\n", t.GetCh(CHAN8).fLe.time_sec, w8_ns);
         }

         fHw8ns->Fill(w8_ns);

         if (w8_ns > w8_cut) { 
            fHw8ns_cut->Fill(w8_ns);
         }
      }

      ///////// CROSS BOARD COINCIDENCES ///////////
      
      if (w1_ns > 0 && w2_ns > 0) {
         double t12_ns = subtract_ns(t.GetCh(CHAN2).fLe, t.GetCh(CHAN1).fLe);
         fHt12ns->Fill(t12_ns);
         
         if (w1_ns > w1_cut && w2_ns > w1_cut) {
            fHt12ns_cut->Fill(t12_ns);
         }
      }
      
      if (w3_ns > 0 && w4_ns > 0) {
         double t34_ns = subtract_ns(t.GetCh(CHAN4).fLe, t.GetCh(CHAN3).fLe);
         fHt34ns->Fill(t34_ns);

         if (w3_ns > w3_cut && w4_ns > w4_cut) {
            fHt34ns_cut->Fill(t34_ns);
         }
      }
      
      if (w5_ns > 0 && w6_ns > 0) {
         double t56_ns = subtract_ns(t.GetCh(CHAN6).fLe, t.GetCh(CHAN5).fLe);
         fHt56ns->Fill(t56_ns);

         if (w5_ns > w5_cut && w6_ns > w6_cut) {
            fHt56ns_cut->Fill(t56_ns);
         }
      }

      if (w7_ns > 0 && w8_ns > 0) {
         double t78_ns = subtract_ns(t.GetCh(CHAN8).fLe, t.GetCh(CHAN7).fLe);
         fHt78ns->Fill(t78_ns);

         if (w7_ns > w7_cut && w8_ns > w8_cut) {
            fHt78ns_cut->Fill(t78_ns);
         }
      }

      ///////// PAIR COINCIDENCES ///////////

      double t14_ns = -9999;
      double t14_ns_twc = -9999;

      if (w1_ns > 0 && w4_ns > 0) {
         t14_ns = subtract_ns(t.GetCh(CHAN4).fLe, t.GetCh(CHAN1).fLe);
         t14_ns_twc = t14_ns + (W1_twc / sqrt(w1_ns) - W4_twc / sqrt(w4_ns));

         if (fFlags->fPrint) {
            printf("new dlsc event 1*4, le %.9f %.9f sec, diff14 %.0f ns, w1 %.0f, w4 %.0f ns!\n", t.GetCh(CHAN1).fLe.time_sec, t.GetCh(CHAN4).fLe.time_sec, t14_ns, w1_ns, w4_ns);
         }

         fCount14++;
         fCountMyA++;
         
         fHt14ns->Fill(t14_ns);
         fHw14ns->Fill(w1_ns, w4_ns);

         fHt14ns_w1ns->Fill(w1_ns, t14_ns);
         fHt14ns_w4ns->Fill(w4_ns, t14_ns);

         if (w4_ns > w4_cut) {
            fHt14ns_w1ns_cutw4->Fill(w1_ns, t14_ns);
            fHt14ns_w1ns_cutw4_twc->Fill(w1_ns, t14_ns_twc);
         }

         if (w1_ns > w1_cut) {
            fHt14ns_w4ns_cutw1->Fill(w4_ns, t14_ns);
            fHt14ns_w4ns_cutw1_twc->Fill(w4_ns, t14_ns_twc);
         }

         if (w1_ns > w1_cut && w4_ns > w4_cut) {
            fHt14ns_cut->Fill(t14_ns);
         }
         
         fHt14ns_twc->Fill(t14_ns_twc);
      }

      double t23_ns = -9999;
      double t23_ns_twc = -9999;

      if (w2_ns > 0 && w3_ns > 0) {
         t23_ns = subtract_ns(t.GetCh(CHAN3).fLe, t.GetCh(CHAN2).fLe);
         t23_ns_twc = t23_ns + (W2_twc / sqrt(w2_ns) - W3_twc / sqrt(w3_ns));
         
         if (fFlags->fPrint) {
            printf("new dlsc event 2*3, le %.9f %.9f sec, diff23 %.0f ns, w2 %.0f, w3 %.0f ns!\n", t.GetCh(CHAN2).fLe.time_sec, t.GetCh(CHAN3).fLe.time_sec, t23_ns, w2_ns, w3_ns);
         }
         
         fCount23++;
         fCountMyA++;

         fHt23ns->Fill(t23_ns);
         fHw23ns->Fill(w2_ns, w3_ns);

         fHt23ns_w2ns->Fill(w2_ns, t23_ns);
         fHt23ns_w3ns->Fill(w3_ns, t23_ns);

         if (w3_ns > w3_cut) {
            fHt23ns_w2ns_cutw3->Fill(w2_ns, t23_ns);
            fHt23ns_w2ns_cutw3_twc->Fill(w2_ns, t23_ns_twc);
         }

         if (w2_ns > w2_cut) {
            fHt23ns_w3ns_cutw2->Fill(w3_ns, t23_ns);
            fHt23ns_w3ns_cutw2_twc->Fill(w3_ns, t23_ns_twc);
         }

         if (w2_ns > w2_cut && w3_ns > w3_cut) {
            fHt23ns_cut->Fill(t23_ns);
         }

         fHt23ns_twc->Fill(t23_ns_twc);
      }

      double t58_ns = -9999;
      double t58_ns_twc = -9999;

      if (w5_ns > 0 && w8_ns > 0) {
         t58_ns = subtract_ns(t.GetCh(CHAN8).fLe, t.GetCh(CHAN5).fLe);
         t58_ns_twc = t58_ns + (W5_twc / sqrt(w5_ns) - W8_twc / sqrt(w8_ns));

         if (fFlags->fPrint) {
            printf("new dlsc event 5*8, le %.9f %.9f sec, diff58 %.0f ns, w5 %.0f, w8 %.0f ns!\n", t.GetCh(CHAN5).fLe.time_sec, t.GetCh(CHAN8).fLe.time_sec, t58_ns, w5_ns, w8_ns);
         }

         fCount58++;
         fCountMyB++;

         fHt58ns->Fill(t58_ns);
         fHw58ns->Fill(w5_ns, w8_ns);

         fHt58ns_w5ns->Fill(w5_ns, t58_ns);
         fHt58ns_w8ns->Fill(w8_ns, t58_ns);

         if (w8_ns > w8_cut) {
            fHt58ns_w5ns_cutw8->Fill(w5_ns, t58_ns);
            fHt58ns_w5ns_cutw8_twc->Fill(w5_ns, t58_ns_twc);
         }

         if (w5_ns > w5_cut) {
            fHt58ns_w8ns_cutw5->Fill(w8_ns, t58_ns);
            fHt58ns_w8ns_cutw5_twc->Fill(w8_ns, t58_ns_twc);
         }

         if (w5_ns > w5_cut && w8_ns > w8_cut) {
            fHt58ns_cut->Fill(t58_ns);
         }

         fHt58ns_twc->Fill(t58_ns_twc);
      }

      double t67_ns = -9999;
      double t67_ns_twc = -9999;
   
      if (w6_ns > 0 && w7_ns > 0) {
         t67_ns = subtract_ns(t.GetCh(CHAN7).fLe, t.GetCh(CHAN6).fLe);
         t67_ns_twc = t67_ns + (W6_twc / sqrt(w6_ns) - W7_twc / sqrt(w7_ns));

         if (fFlags->fPrint) {
            printf("new dlsc event 6*7, le %.9f %.9f sec, diff67 %.0f ns, w6 %.0f, w7 %.0f ns!\n", t.GetCh(CHAN6).fLe.time_sec, t.GetCh(CHAN7).fLe.time_sec, t67_ns, w6_ns, w7_ns);
         }

         fCount67++;
         fCountMyB++;

         fHt67ns->Fill(t67_ns);
         fHw67ns->Fill(w6_ns, w7_ns);

         fHt67ns_w6ns->Fill(w6_ns, t67_ns);
         fHt67ns_w7ns->Fill(w7_ns, t67_ns);

         if (w7_ns > w7_cut) {
            fHt67ns_w6ns_cutw7->Fill(w6_ns, t67_ns);
            fHt67ns_w6ns_cutw7_twc->Fill(w6_ns, t67_ns_twc);
         }

         if (w6_ns > w6_cut) {
            fHt67ns_w7ns_cutw6->Fill(w7_ns, t67_ns);
            fHt67ns_w7ns_cutw6_twc->Fill(w7_ns, t67_ns_twc);
         }

         if (w6_ns > w6_cut && w7_ns > w7_cut) {
            fHt67ns_cut->Fill(t67_ns);
         }

         fHt67ns_twc->Fill(t67_ns_twc);
      }
      
      ///////// QUAD COINCIDENCES ///////////

      if (t14_ns > -9999 && t58_ns > -9999) {
         if (fFlags->fPrint) {
            printf("new dlsc event 1*4*5*8, le %.9f %.9f sec, diff58 %.0f ns, w5 %.0f, w8 %.0f ns!\n", t.GetCh(CHAN5).fLe.time_sec, t.GetCh(CHAN8).fLe.time_sec, t58_ns, w5_ns, w8_ns);
         }

         fCount1458++;
         fCountMyT++;
         
         //double avg14 = 0.5*(t.chan1le.time_sec + t.chan4le.time_sec);
         //double avg58 = 0.5*(t.GetCh(CHAN5).fLe.time_sec + t.GetCh(CHAN8).fLe.time_sec);
         //
         //double tof1458 = avg58 - avg14;
         
         double t15_ns = subtract_ns(t.GetCh(CHAN5).fLe, t.GetCh(CHAN1).fLe);
         double t18_ns = subtract_ns(t.GetCh(CHAN8).fLe, t.GetCh(CHAN1).fLe);
         double t45_ns = subtract_ns(t.GetCh(CHAN5).fLe, t.GetCh(CHAN4).fLe);
         double t48_ns = subtract_ns(t.GetCh(CHAN8).fLe, t.GetCh(CHAN4).fLe);

         
         double tof1458 = 0.5*(t15_ns + t48_ns);
         double tof1458_twc = tof1458 + 0.5*(W1_twc / sqrt(w1_ns) + W4_twc / sqrt(w4_ns) - W5_twc / sqrt(w5_ns) - W8_twc / sqrt(w8_ns));

         double t145_ns = 0.5*(t15_ns + t45_ns);
         double t148_ns = 0.5*(t18_ns + t48_ns);
         double t158_ns = -0.5*(t15_ns + t18_ns);
         double t458_ns = -0.5*(t45_ns + t48_ns);
         double t145_ns_twc = t145_ns + 0.5*(W1_twc / sqrt(w1_ns) + W4_twc / sqrt(w4_ns) - 2*W5_twc / sqrt(w5_ns));
         double t148_ns_twc = t148_ns + 0.5*(W1_twc / sqrt(w1_ns) + W4_twc / sqrt(w4_ns) - 2*W8_twc / sqrt(w8_ns));
         double t158_ns_twc = t158_ns + 0.5*(W5_twc / sqrt(w5_ns) + W8_twc / sqrt(w8_ns) - 2*W1_twc / sqrt(w1_ns));
         double t458_ns_twc = t458_ns + 0.5*(W5_twc / sqrt(w5_ns) + W8_twc / sqrt(w8_ns) - 2*W4_twc / sqrt(w4_ns));
         
         fHtof_1458->Fill(tof1458);

         fHt14ns_1458->Fill(t14_ns);
         fHt58ns_1458->Fill(t58_ns);
         fHt14t58ns_1458->Fill(t14_ns, t58_ns);

         fHt158ns_1458_twc->Fill(t158_ns_twc);
         fHt458ns_1458_twc->Fill(t458_ns_twc);
         fHt145ns_1458_twc->Fill(t145_ns_twc);
         fHt148ns_1458_twc->Fill(t148_ns_twc);
         
         fHw14ns_1458->Fill(w1_ns, w4_ns);
         fHw58ns_1458->Fill(w5_ns, w8_ns);
         
         fHtof_w1ns_1458->Fill(w1_ns, tof1458);
         fHtof_w4ns_1458->Fill(w4_ns, tof1458);

         fHtof_w5ns_1458->Fill(w5_ns, tof1458);
         fHtof_w8ns_1458->Fill(w8_ns, tof1458);

         if (w1_ns > w1_cut && w4_ns > w4_cut && w5_ns > w5_cut && w8_ns > w8_cut) {
            fHtof_1458_cut->Fill(tof1458);
            fHtof_1458_cut_twc->Fill(tof1458_twc);

            fHtof_w1ns_1458_cut->Fill(w1_ns, tof1458);
            fHtof_w4ns_1458_cut->Fill(w4_ns, tof1458);

            fHtof_w5ns_1458_cut->Fill(w5_ns, tof1458);
            fHtof_w8ns_1458_cut->Fill(w8_ns, tof1458);

            //fHt14ns_1458_cut->Fill(t14_ns);
            //fHt58ns_1458_cut->Fill(t58_ns);
            fHt14t58ns_1458_cut->Fill(t14_ns, t58_ns);
            fHt14t58ns_1458_cut_twc->Fill(t14_ns_twc, t58_ns_twc);

            fHt158ns_1458_twc_cut->Fill(t158_ns_twc);
            fHt458ns_1458_twc_cut->Fill(t458_ns_twc);
            fHt145ns_1458_twc_cut->Fill(t145_ns_twc);
            fHt148ns_1458_twc_cut->Fill(t148_ns_twc);
         }

         fHtof_w1ns_1458_twc->Fill(w1_ns, tof1458_twc);
         fHtof_w4ns_1458_twc->Fill(w4_ns, tof1458_twc);

         fHtof_w5ns_1458_twc->Fill(w5_ns, tof1458_twc);
         fHtof_w8ns_1458_twc->Fill(w8_ns, tof1458_twc);

         fHt158ns_w1ns_1458->Fill(w1_ns, t158_ns);
         fHt458ns_w4ns_1458->Fill(w4_ns, t458_ns);
         fHt145ns_w5ns_1458->Fill(w5_ns, t145_ns);
         fHt148ns_w8ns_1458->Fill(w8_ns, t148_ns);

         fHt158ns_w1ns_1458_twc->Fill(w1_ns, t158_ns_twc);
         fHt458ns_w4ns_1458_twc->Fill(w4_ns, t458_ns_twc);
         fHt145ns_w5ns_1458_twc->Fill(w5_ns, t145_ns_twc);
         fHt148ns_w8ns_1458_twc->Fill(w8_ns, t148_ns_twc);

      }

      if (t23_ns > -9999 && t58_ns > -9999) {
         t23_ns = subtract_ns(t.GetCh(CHAN3).fLe, t.GetCh(CHAN2).fLe);

         if (fFlags->fPrint) {
            printf("new dlsc event 2*3*5*8, le %.9f %.9f sec, diff58 %.0f ns, w5 %.0f, w8 %.0f ns!\n", t.GetCh(CHAN5).fLe.time_sec, t.GetCh(CHAN8).fLe.time_sec, t58_ns, w5_ns, w8_ns);
         }
         
         fCount2358++;
         fCountMyT++;

         double t25_ns = subtract_ns(t.GetCh(CHAN5).fLe, t.GetCh(CHAN2).fLe);
         double t38_ns = subtract_ns(t.GetCh(CHAN8).fLe, t.GetCh(CHAN3).fLe);
         
         double tof2358 = 0.5*(t25_ns + t38_ns);
         double tof2358_twc = tof2358 + 0.5*(W2_twc / sqrt(w2_ns) + W3_twc / sqrt(w3_ns) - W5_twc / sqrt(w5_ns) - W8_twc / sqrt(w8_ns));         

         fHtof_2358->Fill(tof2358);
         
         fHt23ns_2358->Fill(t23_ns);
         fHt58ns_2358->Fill(t58_ns);
         fHt23t58ns_2358->Fill(t23_ns, t58_ns);
         
         fHw23ns_2358->Fill(w2_ns, w3_ns);
         fHw58ns_2358->Fill(w5_ns, w8_ns);

         if (w2_ns > w2_cut && w3_ns > w3_cut && w5_ns > w5_cut && w8_ns > w8_cut) {
            fHtof_2358_cut->Fill(tof2358);
            fHtof_2358_cut_twc->Fill(tof2358_twc);

            //fHt23ns_2358_cut->Fill(t23_ns);
            //fHt58ns_2358_cut->Fill(t58_ns);
            fHt23t58ns_2358_cut->Fill(t23_ns, t58_ns);
            fHt23t58ns_2358_cut_twc->Fill(t23_ns_twc, t58_ns_twc);
         }

         fHtof_w2ns_2358_twc->Fill(w2_ns, tof2358_twc);
         fHtof_w3ns_2358_twc->Fill(w3_ns, tof2358_twc);

         fHtof_w5ns_2358_twc->Fill(w5_ns, tof2358_twc);
         fHtof_w8ns_2358_twc->Fill(w8_ns, tof2358_twc);

      }

      if (t14_ns > -9999 && t67_ns > -9999) {

         if (fFlags->fPrint) {
            printf("new dlsc event 1*4*6*7, le %.9f %.9f sec, diff58 %.0f ns, w5 %.0f, w8 %.0f ns!\n", t.GetCh(CHAN5).fLe.time_sec, t.GetCh(CHAN8).fLe.time_sec, t67_ns, w5_ns, w8_ns);
         }
         
         fCount1467++;
         fCountMyT++;

         double t16_ns = subtract_ns(t.GetCh(CHAN6).fLe, t.GetCh(CHAN1).fLe);
	 double t17_ns = subtract_ns(t.GetCh(CHAN7).fLe, t.GetCh(CHAN1).fLe);
         double t46_ns = subtract_ns(t.GetCh(CHAN6).fLe, t.GetCh(CHAN4).fLe);
         double t47_ns = subtract_ns(t.GetCh(CHAN7).fLe, t.GetCh(CHAN4).fLe);
         
         double tof1467 = 0.5*(t16_ns + t47_ns);
         double tof1467_twc = tof1467 + 0.5*(W1_twc / sqrt(w1_ns) + W4_twc / sqrt(w4_ns) - W6_twc / sqrt(w6_ns) - W7_twc / sqrt(w7_ns));         
         
         // Gabby
         // t6 - 0.5*(t1+t4)
         double t146_ns = 0.5*(t16_ns + t46_ns);

         // t7 - 0.5*(t1+t4)
         double t147_ns = 0.5*(t17_ns + t47_ns);

         fHtof_1467->Fill(tof1467);
         
         fHt14ns_1467->Fill(t14_ns);
         fHt67ns_1467->Fill(t67_ns);
         fHt14t67ns_1467->Fill(t14_ns, t67_ns);
         
         fHw14ns_1467->Fill(w1_ns, w4_ns);
         fHw67ns_1467->Fill(w6_ns, w7_ns);

         fHt67ns_w6ns_1467->Fill(w6_ns, t67_ns);
         fHt67ns_w7ns_1467->Fill(w7_ns, t67_ns);

         if (w1_ns > w1_cut && w4_ns > w4_cut && w6_ns > w6_cut && w7_ns > w7_cut) {
            fHtof_1467_cut->Fill(tof1467);
            fHtof_1467_cut_twc->Fill(tof1467_twc);

            //fHt14ns_1467_cut->Fill(t14_ns);
            //fHt67ns_1467_cut->Fill(t67_ns);
            fHt14t67ns_1467_cut->Fill(t14_ns, t67_ns);
            fHt14t67ns_1467_cut_twc->Fill(t14_ns_twc, t67_ns_twc);

            // Gabby
            fHt146ns_w6ns_cut->Fill(w6_ns,t146_ns);
            fHt147ns_w7ns_cut->Fill(w7_ns,t147_ns);
         }

         fHtof_w1ns_1467_twc->Fill(w1_ns, tof1467_twc);
         fHtof_w4ns_1467_twc->Fill(w4_ns, tof1467_twc);

         fHtof_w6ns_1467_twc->Fill(w6_ns, tof1467_twc);
         fHtof_w7ns_1467_twc->Fill(w7_ns, tof1467_twc);
      }

      if (t23_ns > -9999 && t67_ns > -9999) {
         
         if (fFlags->fPrint) {
            printf("new dlsc event 2*3*6*7, le %.9f %.9f sec, diff58 %.0f ns, w5 %.0f, w8 %.0f ns!\n", t.GetCh(CHAN5).fLe.time_sec, t.GetCh(CHAN8).fLe.time_sec, t67_ns, w5_ns, w8_ns);
         }

         fCount2367++;
         fCountMyT++;
         
         double t26_ns = subtract_ns(t.GetCh(CHAN6).fLe, t.GetCh(CHAN2).fLe);
         double t27_ns = subtract_ns(t.GetCh(CHAN7).fLe, t.GetCh(CHAN2).fLe);
         double t36_ns = subtract_ns(t.GetCh(CHAN6).fLe, t.GetCh(CHAN3).fLe);
         double t37_ns = subtract_ns(t.GetCh(CHAN7).fLe, t.GetCh(CHAN3).fLe);
         
         double tof2367 = 0.5*(t26_ns + t37_ns);
         double tof2367_twc = tof2367 + 0.5*(W2_twc / sqrt(w2_ns) + W3_twc / sqrt(w3_ns) - W6_twc / sqrt(w6_ns) - W7_twc / sqrt(w7_ns));        

         double t236_ns = 0.5*(t26_ns + t36_ns);
         double t237_ns = 0.5*(t27_ns + t37_ns);
         double t267_ns = -0.5*(t26_ns + t27_ns);
         double t367_ns = -0.5*(t36_ns + t37_ns);
         double t236_ns_twc = t236_ns + 0.5*(W2_twc / sqrt(w2_ns) + W3_twc / sqrt(w3_ns) - 2*W6_twc / sqrt(w6_ns));
         double t237_ns_twc = t237_ns + 0.5*(W2_twc / sqrt(w2_ns) + W3_twc / sqrt(w3_ns) - 2*W7_twc / sqrt(w7_ns));
         double t267_ns_twc = t267_ns + 0.5*(W6_twc / sqrt(w6_ns) + W7_twc / sqrt(w7_ns) - 2*W2_twc / sqrt(w2_ns));
         double t367_ns_twc = t367_ns + 0.5*(W6_twc / sqrt(w6_ns) + W7_twc / sqrt(w7_ns) - 2*W3_twc / sqrt(w3_ns));
         
         fHtof_2367->Fill(tof2367);
         
         fHt23ns_2367->Fill(t23_ns);
         fHt67ns_2367->Fill(t67_ns);
         fHt23t67ns_2367->Fill(t23_ns, t67_ns);

         fHt267ns_2367_twc->Fill(t267_ns_twc);
         fHt367ns_2367_twc->Fill(t367_ns_twc);
         fHt236ns_2367_twc->Fill(t236_ns_twc);
         fHt237ns_2367_twc->Fill(t237_ns_twc);
         
         fHw23ns_2367->Fill(w2_ns, w3_ns);
         fHw67ns_2367->Fill(w6_ns, w7_ns);

         if (w2_ns > w2_cut && w3_ns > w3_cut && w6_ns > w6_cut && w7_ns > w7_cut) {
            fHtof_2367_cut->Fill(tof2367);
            fHtof_2367_cut_twc->Fill(tof2367_twc);

            //fHt23ns_2367_cut->Fill(t23_ns);
            //fHt67ns_2367_cut->Fill(t67_ns);
            fHt23t67ns_2367_cut->Fill(t23_ns, t67_ns);
            fHt23t67ns_2367_cut_twc->Fill(t23_ns_twc, t67_ns_twc);

            fHt267ns_2367_twc_cut->Fill(t267_ns_twc);
            fHt367ns_2367_twc_cut->Fill(t367_ns_twc);
            fHt236ns_2367_twc_cut->Fill(t236_ns_twc);
            fHt237ns_2367_twc_cut->Fill(t237_ns_twc);
         }

         fHtof_w2ns_2367_twc->Fill(w2_ns, tof2367_twc);
         fHtof_w3ns_2367_twc->Fill(w3_ns, tof2367_twc);

         fHtof_w6ns_2367_twc->Fill(w6_ns, tof2367_twc);
         fHtof_w7ns_2367_twc->Fill(w7_ns, tof2367_twc);

         fHt267ns_w2ns_2367->Fill(w2_ns, t267_ns);
         fHt367ns_w3ns_2367->Fill(w3_ns, t367_ns);
         fHt236ns_w6ns_2367->Fill(w6_ns, t236_ns);
         fHt237ns_w7ns_2367->Fill(w7_ns, t237_ns);

         fHt267ns_w2ns_2367_twc->Fill(w2_ns, t267_ns_twc);
         fHt367ns_w3ns_2367_twc->Fill(w3_ns, t367_ns_twc);
         fHt236ns_w6ns_2367_twc->Fill(w6_ns, t236_ns_twc);
         fHt237ns_w7ns_2367_twc->Fill(w7_ns, t237_ns_twc);

      }
   }

   TAFlowEvent* Analyze(TARunInfo* runinfo, TMEvent* event, TAFlags* flags, TAFlowEvent* flow)
   {
      //printf("DlTdcModule::Analyze, run %d, event serno %d, id 0x%04x, data size %d\n", runinfo->fRunNo, event->serial_number, (int)event->event_id, event->data_size);
      return flow;
   }

   TAFlowEvent* AnalyzeFlowEvent(TARunInfo* runinfo, TAFlags* flags, TAFlowEvent* flow)
   {
      if (!fFlags->fEnabled)
         return flow;

      DlTdcEventFlow* f = flow->Find<DlTdcEventFlow>();
      if (f) {
         AnalyzeTdcEvent(*f->fDlTdcEvent);
      }
      return flow;
   }

   void AnalyzeSpecialEvent(TARunInfo* runinfo, TMEvent* event)
   {
      if (fTrace)
         printf("DlTdcModule::AnalyzeSpecialEvent, run %d, event serno %d, id 0x%04x, data size %d\n", runinfo->fRunNo, event->serial_number, (int)event->event_id, event->data_size);
   }
};

class DlTdc4ModuleFactory: public TAFactory
{
public:
   DlTdcFlags fFlags;

public:
   void Usage()
   {
      printf("DlTdc4ModuleFactory flags:\n");
      printf("--dltdc4 -- enable analysis of 4 paddle data\n");
      printf("--dltdc4-triggered -- analyze only events with hit in channel T (coincidence of A and B)\n");
      printf("--dltdc4-debug -- print detailed information\n");
      printf("--dltdc4-print -- print events\n");
      printf("--dltdc4-twc -- get twc W parameter\n");
   }

   void Init(const std::vector<std::string> &args)
   {
      printf("DlTdc4ModuleFactory::Init!\n");

      for (unsigned i=0; i<args.size(); i++) {
         if (args[i] == "--dltdc4") {
            fFlags.fEnabled = true;
         }
         if (args[i] == "--dltdc4-triggered") {
            fFlags.fTriggered = true;
         }
         if (args[i] == "--dltdc4-debug") {
            fFlags.fDebug = true;
         }
         if (args[i] == "--dltdc4-print") {
            fFlags.fPrint = true;
         }
         if (args[i] == "--dltdc4-twc") {
            fFlags.fTWC = true;
         }
      }
   }

   void Finish()
   {
      printf("DlTdc4ModuleFactory::Finish!\n");
   }
   
   TARunObject* NewRunObject(TARunInfo* runinfo)
   {
      printf("DlTdc4ModuleFactory::NewRunObject, run %d, file %s\n", runinfo->fRunNo, runinfo->fFileName.c_str());
      return new DlTdc4Module(runinfo, &fFlags);
   }
};

static TARegister tar(new DlTdc4ModuleFactory);

/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */
