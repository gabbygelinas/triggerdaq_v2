//
// dltdc8_module.cxx
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

class DlTdcMap8
{
public:
   int fChanA = 32;
   int fChanB = 33;
   int fChanT = 34;

   std::vector<int> fMap;

   std::vector<int> fPair1;
   std::vector<int> fPair2;

public:
   bool Init(int runno);
};

bool DlTdcMap8::Init(int runno)
{
   if (runno < 156) {
      printf("DlTdcMap8 old 1458 map for run %d!\n", runno);

      fMap.resize(16+1);
      
      fMap[1] = 16; // chan1
      fMap[2] = 17; // chan2
      fMap[3] = 22; // chan3
      fMap[4] = 23; // chan4
      fMap[5] = 26; // chan5
      fMap[6] = 27; // chan6
      fMap[7] = 30; // chan7
      fMap[8] = 31; // chan8
      
      fMap[9]  = 18;
      fMap[10] = 19;
      fMap[11] = 20;
      fMap[12] = 21;
      fMap[13] = 24;
      fMap[14] = 25;
      fMap[15] = 28;
      fMap[16] = 29;
      
      fPair1.resize(8+1);
      fPair2.resize(8+1);
      
      fPair1[1] =  1; fPair2[1] =  4; // chan14
      fPair1[2] =  2; fPair2[2] =  3; // chan23
      fPair1[3] =  5; fPair2[3] =  8; // chan58
      fPair1[4] =  6; fPair2[4] =  7; // chan67
      fPair1[5] =  9; fPair2[5] = 12;
      fPair1[6] = 10; fPair2[6] = 11;
      fPair1[7] = 13; fPair2[7] = 16;
      fPair1[8] = 14; fPair2[8] = 15;
   } else {
      printf("DlTdcMap8 for run %d!\n", runno);

      fMap.resize(16+1);
      
      fMap[1]  =  0; // chan1
      fMap[2]  =  1; // chan2
      fMap[3]  = 10; // chan3
      fMap[4]  = 11; // chan4
      fMap[5]  =  2 + 16; // chan5
      fMap[6]  =  3 + 16; // chan6
      fMap[7]  =  8 + 16; // chan7
      fMap[8]  =  9 + 16; // chan8
      
      fMap[9]  = 15;
      fMap[10] = 14;
      fMap[11] =  7;
      fMap[12] =  6;
      fMap[13] = 13 + 16;
      fMap[14] = 12 + 16;
      fMap[15] =  5 + 16;
      fMap[16] =  4 + 16;
      
      fPair1.resize(8+1);
      fPair2.resize(8+1);
      
      fPair1[1] =  1; fPair2[1] =  9; // chan14
      fPair1[2] =  2; fPair2[2] = 10; // chan23
      fPair1[3] =  3; fPair2[3] = 11; // chan58
      fPair1[4] =  4; fPair2[4] = 12; // chan67
      fPair1[5] =  5; fPair2[5] = 13;
      fPair1[6] =  6; fPair2[6] = 14;
      fPair1[7] =  7; fPair2[7] = 15;
      fPair1[8] =  8; fPair2[8] = 16;
   }

   return true;
};

class DlTdc8Module: public TARunObject
{
public:
   DlTdcFlags* fFlags = NULL;
   
#ifdef HAVE_ROOT
   TH1D* fHtABns = NULL;

   std::vector<TH1D*> fHw_ns;
   std::vector<TH1D*> fHw_ns_cut;

   TH2D* fHw_ns_chanmap = NULL;
   TH2D* fHw_ns_chanmap_cut = NULL;

   TH1D* fHwAns = NULL;
   TH1D* fHwBns = NULL;
   TH1D* fHwTns = NULL;

   std::vector<TH2D*> fHwpair_ns;

   std::vector<TH1D*> fHtpair_ns;
   std::vector<TH2D*> fHtpair_w1_ns;
   std::vector<TH2D*> fHtpair_w2_ns;
   std::vector<TH2D*> fHtpair_w1_ns_cut_w2;
   std::vector<TH2D*> fHtpair_w2_ns_cut_w1;
   std::vector<TH1D*> fHtpair_ns_cut;
   std::vector<TH2D*> fHtpair_w1_ns_cut_w2_twc;
   std::vector<TH2D*> fHtpair_w2_ns_cut_w1_twc;
   std::vector<TH1D*> fHtpair_ns_cut_twc;
   //std::vector<TH2D*> fHpairPP;
   std::vector<TH2D*> fHpairFF;

   TH1D* fHquad11 = NULL;
   TH1D* fHquad15 = NULL;
   TH1D* fHquad16 = NULL;
   TH1D* fHquad17 = NULL;
   TH1D* fHquad18 = NULL;

   TH1D* fHquad11AB = NULL;
   TH1D* fHquad15AB = NULL;
   TH1D* fHquad16AB = NULL;
   TH1D* fHquad17AB = NULL;
   TH1D* fHquad18AB = NULL;

#if 0
   TH1D* fHt12ns = NULL;
   TH1D* fHt34ns = NULL;
   TH1D* fHt56ns = NULL;
   TH1D* fHt78ns = NULL;

   TH1D* fHt12ns_cut = NULL;
   TH1D* fHt34ns_cut = NULL;
   TH1D* fHt56ns_cut = NULL;
   TH1D* fHt78ns_cut = NULL;

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
#endif

   int fCountA = 0;
   int fCountB = 0;
   int fCountT = 0;

#endif

#if 0
   int fCountCut2367 = 0;
   int fCountCut1 = 0;
   int fCountCut4 = 0;
   int fCountCut5 = 0;
   int fCountCut8 = 0;
#endif

   double fPrevEventTimeSec = 0;

   DlTdcMap8 *fMap8 = NULL;

   bool fTrace = false;

   Ncfm* fCfm = NULL;
   
   DlTdc8Module(TARunInfo* runinfo, DlTdcFlags* flags)
      : TARunObject(runinfo)
   {
      if (fTrace)
         printf("DlTdc8Module::ctor!\n");

      fModuleName = "dltdc8_module";
      fFlags   = flags;

      fCfm = new Ncfm("dlcfmdb");

      fMap8 = new DlTdcMap8();
   }

   ~DlTdc8Module()
   {
      if (fTrace)
         printf("DlTdc8Module::dtor!\n");

      if (fMap8) {
         delete fMap8;
         fMap8 = NULL;
      }

      if (fCfm) {
         delete fCfm;
         fCfm = NULL;
      }
   }

   void BeginRun(TARunInfo* runinfo)
   {
      if (fTrace)
         printf("DlTdc8Module::BeginRun, run %d, file %s\n", runinfo->fRunNo, runinfo->fFileName.c_str());

      if (!fFlags->fEnabled)
         return;

      bool conf_ok = fMap8->Init(runinfo->fRunNo);
      if (!conf_ok) {
         printf("Cannot load TDC map for run %d\n", runinfo->fRunNo);
         exit(123);
      }

#ifdef HAVE_ROOT
      runinfo->fRoot->fOutputFile->cd(); // select correct ROOT directory
      TDirectory* dir = gDirectory->mkdir("dltdc8");
      dir->cd(); // select correct ROOT directory

      fHtABns = new TH1D("dltdc8_tABns", "time between 1234 trigger and 4567 trigger, A-B, ns", 200, -10, 10);

      fHw_ns.resize(fMap8->fMap.size() + 1);

      for (size_t i=1; i<fMap8->fMap.size(); i++) {
         char name[256];
         char title[256];

         sprintf(name, "dltdc8_w%02d_ns", (int)i);
         sprintf(title, "pulse width chan %2d, ns", (int)i);
         fHw_ns[i] = new TH1D(name, title, 100, 0, 100);
      }
         
      fHw_ns_cut.resize(fMap8->fMap.size() + 1);

      for (size_t i=1; i<fMap8->fMap.size(); i++) {
         char name[256];
         char title[256];

         sprintf(name, "dltdc8_w%02d_ns_cut", (int)i);
         sprintf(title, "pulse width chan %2d with pulse height cut, ns", (int)i);
         fHw_ns_cut[i] = new TH1D(name, title, 100, 0, 100);
      }

      fHw_ns_chanmap = new TH2D("dltdc8_w_ns_chanmap", "Pulse width, ns, for each channel", fMap8->fMap.size(), 0.5, fMap8->fMap.size()+0.5, 100, 0, 100);

      fHw_ns_chanmap_cut = new TH2D("dltdc8_w_ns_chanmap_cut", "Pulse width, ns, for each channel", fMap8->fMap.size(), 0.5, fMap8->fMap.size()+0.5, 100, 0, 100);

      fHwAns = new TH1D("dltdc8_wA_ns", "width of A, ns", 100, 0, 100);
      fHwBns = new TH1D("dltdc8_wB_ns", "width of B, ns", 100, 0, 100);
      fHwTns = new TH1D("dltdc8_wT_ns", "width of T, ns", 100, 0, 100);

      fHwpair_ns.resize(fMap8->fPair1.size() + 1);

      for (size_t p=1; p<fMap8->fPair1.size(); p++) {
         int c1 = fMap8->fPair1[p];
         int c2 = fMap8->fPair2[p];

         char name[256];
         char title[256];

         sprintf(name, "dltdc8_pair%d_w_%02d_%02d_ns", (int)p, c1, c2);
         sprintf(title, "pair %d pulse width chan %2d vs chan %2d, ns", (int)p, c2, c1);
         fHwpair_ns[p] = new TH2D(name, title, 100, 0, 100, 100, 0, 100);
      }

      //fHpairPP.resize(fMap8->fPair1.size() + 1);
      fHpairFF.resize(fMap8->fPair1.size() + 1);

      for (size_t p=1; p<fMap8->fPair1.size(); p++) {
         int c1 = fMap8->fPair1[p];
         int c2 = fMap8->fPair2[p];

         int tdc1 = fMap8->fMap[c1];
         int tdc2 = fMap8->fMap[c2];

         char name[256];
         char title[256];

         //sprintf(name, "dltdc8_pair%d_t_%02d_%02d_PP", (int)p, c1, c2);
         //sprintf(title, "pair %d time bin chan %2d vs chan %2d, tdc %2d vs %2d", (int)p, c2, c1, tdc2, tdc1);
         //fHpairPP[p] = new TH2D(name, title, 101, -50.5, 50.5, 101, -50.5, 50.5);

         sprintf(name, "dltdc8_pair%d_t_%02d_%02d_FF", (int)p, c1, c2);
         sprintf(title, "pair %d fine time chan %2d vs %2d, tdc %2d vs %2d, ns vs ns", (int)p, c2, c1, tdc2, tdc1);
         fHpairFF[p] = new TH2D(name, title, 200, -5, 15, 200, -5, 15);
      }

      fHtpair_ns.resize(fMap8->fPair1.size() + 1);

      for (size_t p=1; p<fMap8->fPair1.size(); p++) {
         int c1 = fMap8->fPair1[p];
         int c2 = fMap8->fPair2[p];

         char name[256];
         char title[256];

         sprintf(name, "dltdc8_pair%d_t_%02d_%02d_ns", (int)p, c1, c2);
         sprintf(title, "pair %d time difference chan %2d minus %2d, ns", (int)p, c2, c1);
         fHtpair_ns[p] = new TH1D(name, title, 200, -10, 10);
      }

      fHtpair_w1_ns.resize(fMap8->fPair1.size() + 1);
      fHtpair_w2_ns.resize(fMap8->fPair1.size() + 1);

      for (size_t p=1; p<fMap8->fPair1.size(); p++) {
         int c1 = fMap8->fPair1[p];
         int c2 = fMap8->fPair2[p];

         char name[256];
         char title[256];

         sprintf(name, "dltdc8_pair%d_t_w1_%02d_%02d_ns", (int)p, c1, c2);
         sprintf(title, "pair %d time difference chan %2d minus %2d, ns vs pulse width chan %2d", (int)p, c2, c1, c1);
         fHtpair_w1_ns[p] = new TH2D(name, title, 100, 0, 100, 200, -10, 10);

         sprintf(name, "dltdc8_pair%d_t_w2_%02d_%02d_ns", (int)p, c1, c2);
         sprintf(title, "pair %d time difference chan %2d minus %2d, ns vs pulse width chan %2d", (int)p, c2, c1, c2);
         fHtpair_w2_ns[p] = new TH2D(name, title, 100, 0, 100, 200, -10, 10);
      }

      fHtpair_w1_ns_cut_w2.resize(fMap8->fPair1.size() + 1);
      fHtpair_w2_ns_cut_w1.resize(fMap8->fPair1.size() + 1);

      for (size_t p=1; p<fMap8->fPair1.size(); p++) {
         int c1 = fMap8->fPair1[p];
         int c2 = fMap8->fPair2[p];

         char name[256];
         char title[256];

         sprintf(name, "dltdc8_pair%d_t_w1_%02d_%02d_ns_cut_w2", (int)p, c1, c2);
         sprintf(title, "pair %d time difference chan %2d minus %2d, ns vs pulse width chan %2d, cut on width of %2d", (int)p, c2, c1, c1, c2);
         fHtpair_w1_ns_cut_w2[p] = new TH2D(name, title, 100, 0, 100, 200, -10, 10);

         sprintf(name, "dltdc8_pair%d_t_w2_%02d_%02d_ns_cut_w1", (int)p, c1, c2);
         sprintf(title, "pair %d time difference chan %2d minus %2d, ns vs pulse width chan %2d, cut on width of %2d", (int)p, c2, c1, c2, c1);
         fHtpair_w2_ns_cut_w1[p] = new TH2D(name, title, 100, 0, 100, 200, -10, 10);
      }

      fHtpair_ns_cut.resize(fMap8->fPair1.size() + 1);

      for (size_t p=1; p<fMap8->fPair1.size(); p++) {
         int c1 = fMap8->fPair1[p];
         int c2 = fMap8->fPair2[p];

         char name[256];
         char title[256];

         sprintf(name, "dltdc8_pair%d_t_%02d_%02d_ns_cut", (int)p, c1, c2);
         sprintf(title, "pair %d time difference chan %2d minus %2d, ns, with cut on width", (int)p, c2, c1);
         fHtpair_ns_cut[p] = new TH1D(name, title, 200, -10, 10);
      }

      fHtpair_w1_ns_cut_w2_twc.resize(fMap8->fPair1.size() + 1);
      fHtpair_w2_ns_cut_w1_twc.resize(fMap8->fPair1.size() + 1);

      for (size_t p=1; p<fMap8->fPair1.size(); p++) {
         int c1 = fMap8->fPair1[p];
         int c2 = fMap8->fPair2[p];

         char name[256];
         char title[256];

         sprintf(name, "dltdc8_pair%d_t_w1_%02d_%02d_ns_twc", (int)p, c1, c2);
         sprintf(title, "pair %d time difference chan %2d minus %2d, ns vs pulse width chan %2d, cut on width of %2d, with time walk correction", (int)p, c2, c1, c1, c2);
         fHtpair_w1_ns_cut_w2_twc[p] = new TH2D(name, title, 100, 0, 100, 200, -10, 10);

         sprintf(name, "dltdc8_pair%d_t_w2_%02d_%02d_ns_twc", (int)p, c1, c2);
         sprintf(title, "pair %d time difference chan %2d minus %2d, ns vs pulse width chan %2d, cut on width of %2d, with time walk correction", (int)p, c2, c1, c2, c1);
         fHtpair_w2_ns_cut_w1_twc[p] = new TH2D(name, title, 100, 0, 100, 200, -10, 10);
      }

      fHtpair_ns_cut_twc.resize(fMap8->fPair1.size() + 1);

      for (size_t p=1; p<fMap8->fPair1.size(); p++) {
         int c1 = fMap8->fPair1[p];
         int c2 = fMap8->fPair2[p];

         char name[256];
         char title[256];

         sprintf(name, "dltdc8_pair%d_t_%02d_%02d_ns_cut_twc", (int)p, c1, c2);
         sprintf(title, "pair %d time difference chan %2d minus %2d, ns, with cut on width and time walk correction", (int)p, c2, c1);
         fHtpair_ns_cut_twc[p] = new TH1D(name, title, 200, -10, 10);
      }

      // quad coincidences

      fHquad11 = new TH1D("quad11", "quad11", 200, -10, 10);
      fHquad15 = new TH1D("quad15", "quad15", 200, -10, 10);
      fHquad16 = new TH1D("quad16", "quad16", 200, -10, 10);
      fHquad17 = new TH1D("quad17", "quad17", 200, -10, 10);
      fHquad18 = new TH1D("quad18", "quad18", 200, -10, 10);

      fHquad11AB = new TH1D("quad11AB", "quad11 A-B, ns", 200, -10, 10);
      fHquad15AB = new TH1D("quad15AB", "quad15 A-B, ns", 200, -10, 10);
      fHquad16AB = new TH1D("quad16AB", "quad16 A-B, ns", 200, -10, 10);
      fHquad17AB = new TH1D("quad17AB", "quad17 A-B, ns", 200, -10, 10);
      fHquad18AB = new TH1D("quad18AB", "quad18 A-B, ns", 200, -10, 10);

#if 0
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
#endif
   }

   void PreEndRun(TARunInfo* runinfo)
   {
      if (fTrace)
         printf("DlTdc8Module::PreEndRun, run %d\n", runinfo->fRunNo);
   }
   
   void EndRun(TARunInfo* runinfo)
   {
      if (fTrace)
         printf("DlTdc8Module::EndRun, run %d\n", runinfo->fRunNo);

      if (!fFlags->fEnabled)
         return;

#if 0
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
#endif
   }
   
   void PauseRun(TARunInfo* runinfo)
   {
      if (fTrace)
         printf("DlTdc8Module::PauseRun, run %d\n", runinfo->fRunNo);
   }

   void ResumeRun(TARunInfo* runinfo)
   {
      if (fTrace)
         printf("DlTdc8Module::ResumeRun, run %d\n", runinfo->fRunNo);
   }

   void AnalyzeTdcEvent(const DlTdcEvent& t)
   {
      // compute number of fired channels

      if (0) {
         int nhits = 0;
         double first_le_sec = 0;
         double last_le_sec = 0;
         double last_te_sec = 0;
         
         for (size_t i=1; i<fMap8->fMap.size(); i++) {
            int tdc_ch = fMap8->fMap[i];
            if (t.HaveCh(tdc_ch)) {
               double le_sec = t.GetCh(tdc_ch).fLe.time_sec;
               double te_sec = t.GetCh(tdc_ch).fTe.time_sec;
               if (first_le_sec==0 || le_sec < first_le_sec) {
                  first_le_sec = le_sec;
               }
               if (last_le_sec==0 || le_sec > last_le_sec) {
                  last_le_sec = le_sec;
               }
               if (last_te_sec==0 || te_sec > last_te_sec) {
                  last_te_sec = te_sec;
               }
               nhits++;
            } else {
            }
         }
         
         double le_le_ns = sec_to_ns(last_le_sec-first_le_sec);
         double le_te_ns = sec_to_ns(last_te_sec-first_le_sec);
         
         if (le_le_ns > 20.0)
            printf("nhits %d, le-le %.3f ns, le-te %.3f ns\n", nhits, le_le_ns, le_te_ns);
      }

      int nhits = 0;
      for (size_t i=1; i<fMap8->fMap.size(); i++) {
         int tdc_ch = fMap8->fMap[i];
         if (t.HaveCh(tdc_ch)) {
            nhits++;
         }
      }

      // debug coincidences

      //if (nhits != 4)
      //   return;

      if (fFlags->fDebug) {
         std::string s = "";
         for (size_t i=1; i<fMap8->fMap.size(); i++) {
            int tdc_ch = fMap8->fMap[i];
            if (t.HaveCh(tdc_ch)) {
               s += "H";
            } else {
               s += ".";
            }
         }

         printf("EVENT %s, ABT %d%d%d, %d hits\n", s.c_str(), t.HaveCh(fMap8->fChanA), t.HaveCh(fMap8->fChanB), t.HaveCh(fMap8->fChanT), nhits);
      }

      ///////// check for triggered event ///////////
         
      if (fFlags->fTriggered && !t.HaveCh(fMap8->fChanT)) {
         return;
      }

      ///////// plot trigger A-B time ///////////

      double tAB_ns = -9999;
      
      if (t.HaveCh(fMap8->fChanA) && t.HaveCh(fMap8->fChanB)) {
         tAB_ns = subtract_ns(t.GetCh(fMap8->fChanA).fLe, t.GetCh(fMap8->fChanB).fLe);
         fHtABns->Fill(tAB_ns);
      }

      //if (!t.HaveCh(fMap8->fChanB)) return;

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
      
      double ww_cut_ns = 20;
      //double ww_cut_ns = 2;

      ///////// COMPUTE WIDTH AND PULSE HEIGHT ///////////

      std::vector<double> ww_ns(fMap8->fMap.size() + 1);

      for (size_t i=1; i<fMap8->fMap.size(); i++) {
         int tdc_ch = fMap8->fMap[i];

         double w_ns = -9999;

         if (t.HaveCh(tdc_ch)) {
            w_ns = t.GetCh(tdc_ch).fWidthNs;

            if (w_ns < 0.01) {
               printf("WWW: BAD WIDTH chan %02d, tdc_ch %02d: %f ns!\n", (int)i, (int)tdc_ch, w_ns);
               w_ns = -9999;
            }
         }

         ww_ns[i] = w_ns;
      }

      // SET UP TIME WALK CORRECT VALUES //
      double ww_twc = 9.0;

#if 0
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
#endif

      ///////// TRIGGER CHANNALS A, B and T not used yet ///////////

      if (t.fHits[fMap8->fChanA].fDown) {
         double wA_ns = t.fHits[fMap8->fChanA].fWidthNs;

         if (wA_ns < -2.0) {
            printf("WWW: BAD WIDTH chanA %f!\n", wA_ns);
         }

         if (fFlags->fPrint) {
            printf("new dlsc event A, le %.9f, w %.0f!\n", t.GetCh(fMap8->fChanA).fLe.time_sec, wA_ns);
         }

         fHwAns->Fill(wA_ns);

         fCountA++;
      }

      if (t.fHits[fMap8->fChanB].fDown) {
         double wB_ns = t.fHits[fMap8->fChanB].fWidthNs;

         if (wB_ns < -2.0) {
            printf("WWW: BAD WIDTH chanB %f!\n", wB_ns);
         }

         if (fFlags->fPrint) {
            printf("new dlsc event B, le %.9f, w %.0f!\n", t.GetCh(fMap8->fChanB).fLe.time_sec, wB_ns);
         }

         fHwBns->Fill(wB_ns);

         fCountB++;
      }

      if (t.fHits[fMap8->fChanT].fDown) {
         double wT_ns = t.fHits[fMap8->fChanT].fWidthNs;

         if (wT_ns < -2.0) {
            printf("WWW: BAD WIDTH chanT %f!\n", wT_ns);
         }

         if (fFlags->fPrint) {
            printf("new dlsc event T, le %.9f, w %.0f!\n", t.GetCh(fMap8->fChanT).fLe.time_sec, wT_ns);
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

      for (size_t i=1; i<fMap8->fMap.size(); i++) {
         if (ww_ns[i] > 0) {
            fHw_ns[i]->Fill(ww_ns[i]);
            fHw_ns_chanmap->Fill(i, ww_ns[i]);
            
            if (ww_ns[i] > ww_cut_ns) { 
               fHw_ns_cut[i]->Fill(ww_ns[i]);
               fHw_ns_chanmap_cut->Fill(i, ww_ns[i]);
            }
         }
      }

#if 0
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
#endif

      ///////// PAIR COINCIDENCES ///////////

      std::vector<bool>   ttpair_hit;
      std::vector<bool>   ttpair_hit_cut;
      std::vector<double> ttpair_ns;
      std::vector<double> ttpair_ns_twc;

      ttpair_hit.resize(fMap8->fPair1.size() + 1);
      ttpair_hit_cut.resize(fMap8->fPair1.size() + 1);
      ttpair_ns.resize(fMap8->fPair1.size() + 1);
      ttpair_ns_twc.resize(fMap8->fPair1.size() + 1);

      for (size_t p=1; p<fMap8->fPair1.size(); p++) {
         int c1 = fMap8->fPair1[p];
         int c2 = fMap8->fPair2[p];
         if (ww_ns[c1]>0 && ww_ns[c2]>0) {
            ttpair_hit[p] = true;
         }
      }

      int npairhits = 0;

      for (size_t p=1; p<fMap8->fPair1.size(); p++) {
         if (ttpair_hit[p]) {
            npairhits++;
         }
      }

      if (0) {
         if (npairhits != 2)
            return;

         //if (ttpair_hit[8] && (ttpair_hit[1] || ttpair_hit[2] || ttpair_hit[3] || ttpair_hit[4])) {
         if (ttpair_hit[1] && (ttpair_hit[5] || ttpair_hit[6] || ttpair_hit[7] || ttpair_hit[8])) {
            //if (ttpair_hit[1] && (ttpair_hit[5] || ttpair_hit[8])) {
            
         } else {
            return;
         }
      }

      for (size_t p=1; p<fMap8->fPair1.size(); p++) {
         if (ttpair_hit[p]) {
            int c1 = fMap8->fPair1[p];
            int c2 = fMap8->fPair2[p];

            int tdc1 = fMap8->fMap[c1];
            int tdc2 = fMap8->fMap[c2];

            double tpair_ns = subtract_ns(t.GetCh(tdc2).fLe, t.GetCh(tdc1).fLe);
            double tpair_ns_twc = tpair_ns + (ww_twc / sqrt(ww_ns[c1]) - ww_twc / sqrt(ww_ns[c2]));

            //if (fFlags->fPrint) {
            //   printf("new dlsc event 1*4, le %.9f %.9f sec, diff14 %.0f ns, w1 %.0f, w4 %.0f ns!\n", t.GetCh(CHAN1).fLe.time_sec, t.GetCh(CHAN4).fLe.time_sec, t14_ns, w1_ns, w4_ns);
            //}

            ttpair_hit[p] = true;
            ttpair_ns[p] = tpair_ns;
            ttpair_ns_twc[p] = tpair_ns_twc;
            
            //fCount14++;
            //fCountMyA++;
            
            fHtpair_ns[p]->Fill(tpair_ns);

#if 0
            if (t.fHits[tdc1].fLe.phase == -10) continue;
            if (t.fHits[tdc1].fLe.phase == -20) continue;
            if (t.fHits[tdc1].fLe.phase == -30) continue;
            if (t.fHits[tdc1].fLe.phase == -40) continue;

            if (t.fHits[tdc2].fLe.phase == -10) continue;
            if (t.fHits[tdc2].fLe.phase == -20) continue;
            if (t.fHits[tdc2].fLe.phase == -30) continue;
            if (t.fHits[tdc2].fLe.phase == -40) continue;

            if (t.fHits[tdc1].fLe.phase ==  10) continue;
            if (t.fHits[tdc1].fLe.phase ==  20) continue;
            if (t.fHits[tdc1].fLe.phase ==  30) continue;
            if (t.fHits[tdc1].fLe.phase ==  40) continue;

            if (t.fHits[tdc2].fLe.phase ==  10) continue;
            if (t.fHits[tdc2].fLe.phase ==  20) continue;
            if (t.fHits[tdc2].fLe.phase ==  30) continue;
            if (t.fHits[tdc2].fLe.phase ==  40) continue;
#endif

            //fHpairPP[p]->Fill(t.fHits[tdc1].fLe.phase, t.fHits[tdc2].fLe.phase);
            fHpairFF[p]->Fill(t.fHits[tdc1].fLe.fine_ns, t.fHits[tdc2].fLe.fine_ns);

            if (ww_ns[c1] < ww_cut_ns && ww_ns[c2] < ww_cut_ns) {
            } else {
               fHwpair_ns[p]->Fill(ww_ns[c1], ww_ns[c2]);
            }
            
            fHtpair_w1_ns[p]->Fill(ww_ns[c1], tpair_ns);
            fHtpair_w2_ns[p]->Fill(ww_ns[c2], tpair_ns);
            
            if (ww_ns[c2] > ww_cut_ns) {
               fHtpair_w1_ns_cut_w2[p]->Fill(ww_ns[c1], tpair_ns);
               fHtpair_w1_ns_cut_w2_twc[p]->Fill(ww_ns[c1], tpair_ns_twc);
            }
            
            if (ww_ns[c1] > ww_cut_ns) {
               fHtpair_w2_ns_cut_w1[p]->Fill(ww_ns[c2], tpair_ns);
               fHtpair_w2_ns_cut_w1_twc[p]->Fill(ww_ns[c2], tpair_ns_twc);
            }
            
            if (ww_ns[c1] > ww_cut_ns && ww_ns[c2] > ww_cut_ns) {
               fHtpair_ns_cut[p]->Fill(tpair_ns);
               fHtpair_ns_cut_twc[p]->Fill(tpair_ns_twc);
               ttpair_hit_cut[p] = true;
            }
         }
      }

      if (fFlags->fDebug) {
         for (size_t p=1; p<fMap8->fPair1.size(); p++) {
            if (ttpair_hit[p]) {
               int c1 = fMap8->fPair1[p];
               int c2 = fMap8->fPair2[p];
               printf("pair %02zu chan %02d-%02d ttpair: hit %d, ns_twc %6.3f ns\n", p, c1, c2, (int)ttpair_hit[p], ttpair_ns_twc[p]);
            }
         }
      }

      ///////// QUAD COINCIDENCES ///////////

      if (0) {
         if (ttpair_hit_cut[5]) {
            if (drand48() > 38418.0/39118.0)
               ttpair_hit_cut[5] = false;
         }
         if (ttpair_hit_cut[6]) {
            if (drand48() > 38418.0/41452.0)
               ttpair_hit_cut[6] = false;
         }
         if (ttpair_hit_cut[7]) {
            if (drand48() > 38418.0/41306.0)
               ttpair_hit_cut[7] = false;
         }
         if (ttpair_hit_cut[8]) {
            if (drand48() > 38418.0/38418.0)
               ttpair_hit_cut[8] = false;
         }
      }

      if (0) {
         if (ttpair_hit_cut[5]) {
            if (drand48() > 0.85)
               ttpair_hit_cut[5] = false;
         }
         if (ttpair_hit_cut[6]) {
            if (drand48() > 0.95)
               ttpair_hit_cut[6] = false;
         }
         if (ttpair_hit_cut[7]) {
            if (drand48() > 0.90)
               ttpair_hit_cut[7] = false;
         }
         if (ttpair_hit_cut[8]) {
            if (drand48() > 0.90)
               ttpair_hit_cut[8] = false;
         }
      }

      if (1) {
         int p=3;
         if (ttpair_hit_cut[p]) {
            if (ttpair_hit_cut[5]) {
               fHquad11->Fill(ttpair_ns_twc[p]);
               fHquad15->Fill(ttpair_ns_twc[p]);
               fHquad11AB->Fill(tAB_ns);
               fHquad15AB->Fill(tAB_ns);
            }
            if (ttpair_hit_cut[6]) {
               fHquad11->Fill(ttpair_ns_twc[p]);
               fHquad16->Fill(ttpair_ns_twc[p]);
               fHquad11AB->Fill(tAB_ns);
               fHquad16AB->Fill(tAB_ns);
            }
            if (ttpair_hit_cut[7]) {
               fHquad11->Fill(ttpair_ns_twc[p]);
               fHquad17->Fill(ttpair_ns_twc[p]);
               fHquad11AB->Fill(tAB_ns);
               fHquad17AB->Fill(tAB_ns);
            }
            if (ttpair_hit_cut[8]) {
               fHquad11->Fill(ttpair_ns_twc[p]);
               fHquad18->Fill(ttpair_ns_twc[p]);
               fHquad11AB->Fill(tAB_ns);
               fHquad18AB->Fill(tAB_ns);
            }
         }
      }

#if 0      
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
#endif
   }

   TAFlowEvent* Analyze(TARunInfo* runinfo, TMEvent* event, TAFlags* flags, TAFlowEvent* flow)
   {
      //printf("DlTdc8Module::Analyze, run %d, event serno %d, id 0x%04x, data size %d\n", runinfo->fRunNo, event->serial_number, (int)event->event_id, event->data_size);
      return flow;
   }

   TAFlowEvent* AnalyzeFlowEvent(TARunInfo* runinfo, TAFlags* flags, TAFlowEvent* flow)
   {
      //printf("DlTdc8Module::AnalyzeFlowEvent, run %d\n", runinfo->fRunNo);

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
         printf("DlTdc8Module::AnalyzeSpecialEvent, run %d, event serno %d, id 0x%04x, data size %d\n", runinfo->fRunNo, event->serial_number, (int)event->event_id, event->data_size);
   }
};

class DlTdc8ModuleFactory: public TAFactory
{
public:
   DlTdcFlags fFlags;

public:
   void Usage()
   {
      printf("DlTdc8ModuleFactory flags:\n");
      printf("--dltdc8 -- enable analysis of 4 paddle data\n");
      printf("--dltdc8-triggered -- analyze only events with hit in channel T (coincidence of A and B)\n");
      printf("--dltdc8-debug -- print detailed information\n");
      printf("--dltdc8-print -- print events\n");
      printf("--dltdc8-twc -- get twc W parameter\n");
   }

   void Init(const std::vector<std::string> &args)
   {
      printf("DlTdc8ModuleFactory::Init!\n");

      //for (unsigned i=0; i<args.size(); i++) {
      //   printf("arg[%d] is [%s]\n", i, args[i].c_str());
      //}

      //printf("test1 %d\n", (std::string("--dltdc8") == "--dltdc8"));
      //printf("test2 %d\n", (args[1] == "--dltdc8"));
      //printf("test3 length %zu\n", args[1].length());
      //for (size_t i=0; i<args[1].length(); i++) {
      //   printf("char %zu is %d [%c]\n", i, args[1][i], args[1][i]);
      //}

      for (unsigned i=0; i<args.size(); i++) {
         if (args[i] == "--dltdc8") {
            fFlags.fEnabled = true;
         }
         if (args[i] == "--dltdc8-triggered") {
            fFlags.fTriggered = true;
         }
         if (args[i] == "--dltdc8-debug") {
            fFlags.fDebug = true;
         }
         if (args[i] == "--dltdc8-print") {
            fFlags.fPrint = true;
         }
         if (args[i] == "--dltdc8-twc") {
            fFlags.fTWC = true;
         }
      }
   }

   void Finish()
   {
      printf("DlTdc8ModuleFactory::Finish!\n");
   }
   
   TARunObject* NewRunObject(TARunInfo* runinfo)
   {
      printf("DlTdc8ModuleFactory::NewRunObject, run %d, file %s\n", runinfo->fRunNo, runinfo->fFileName.c_str());
      return new DlTdc8Module(runinfo, &fFlags);
   }
};

static TARegister tar(new DlTdc8ModuleFactory);

/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */
