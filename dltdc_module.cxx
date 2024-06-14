//
// dltdc_module.cxx
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
#endif

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

static double sec_to_ns(double t)
{
   return t*1e9;
}

static double subtract_ns(const DlTdcHit& h1, const DlTdcHit& h2)
{
   return (h1.coarse_sec- h2.coarse_sec)*1e9 + ((h1.fine_ns + h1.offset_ns) - (h2.fine_ns + h2.offset_ns));
}

class DlTdcFlags
{
public:
   bool fEnabled = true;
   bool fTriggered = false;
   bool fCalibFineTime = false;
   bool fCalibOffsets  = false;
   bool fDebug = false;
   bool fPrint = false;
   bool fTWC = false;
};

class DlTdcHit2
{
public:
   DlTdcHit fLe;
   DlTdcHit fTe;
   bool     fUp   = false;
   bool     fDown = false;
   int      fCount = 0;

public:
   double fTimeSec  = 0;
   double fWidthNs = 0;

public:
   void Clear()
   {
      fLe.Clear();
      fTe.Clear();
      fUp   = false;
      fDown = false;
      fCount = 0;
      fTimeSec  = 0;
      fWidthNs = 0;
   }

   void AddHit(const DlTdcHit& h)
   {
      bool complain = true; // false;

      if (h.le) {
         if (!fUp) {
            fUp = true;
            fDown = false;
            if (fCount == 0) {
               fLe = h;
               fTimeSec = h.time_sec;
            } else {
               double te_to_le_ns = sec_to_ns(h.time_sec - fTe.time_sec);
               if (te_to_le_ns < 5.0) {
                  if (complain)
                     printf("TTT: ch %d: LE-TE-LE, MERGE TE to LE %.3f ns\n", h.ch, te_to_le_ns);
                  fCount = 0;
               } else {
                  if (complain)
                     printf("TTT: ch %d: MULTIPLE HIT, time %.9f -> %.9f sec, dt %.3f ns, TE to LE %.3f ns, count %d\n", h.ch, fTimeSec, h.time_sec, sec_to_ns(h.time_sec - fTimeSec), sec_to_ns(h.time_sec - fTe.time_sec), fCount);
                  //fLe.Print(); printf("\n");
                  //fTe.Print(); printf("\n");
                  //h.Print(); printf("\n");
               }
            }
         } else {
            double le_le_dt_ns = sec_to_ns(h.time_sec - fLe.time_sec);
            if (le_le_dt_ns < 80.0) {
               if (complain)
                  printf("TTT: ch %d: LE-LE, LE-LE %.3f ns\n", h.ch, le_le_dt_ns);
            } else {
               if (complain)
                  printf("TTT: ch %d: LE-LE, LE FROM WRONG EVENT, LE-LE %.3f ns, count %d\n", h.ch, le_le_dt_ns, fCount);
            }
         }
      } else if (h.te) {
         if (fUp) {
            fUp = false;
            fDown = true;
            if (fCount == 0) {
               fTe = h;
               fWidthNs = subtract_ns(fTe, fLe);
            } else {
               if (complain)
                  printf("TTT: ch %d: TE multiple hit, count %d\n", h.ch, fCount);
            }
            fCount++;
         } else {
            if (fCount == 0) {
               if (complain)
                  printf("TTT: ch %d: TE without LE\n", h.ch);
            } else {
               double le_te_dt_ns = sec_to_ns(fTe.time_sec - fLe.time_sec);
               double te_te_dt_ns = sec_to_ns(h.time_sec - fTe.time_sec);
               if (te_te_dt_ns < 80.0) {
                  if (complain)
                     printf("TTT: ch %d: LE-TE-TE, MISSING LE, LE-TE %.3f, TE-TE %.3f ns, count %d\n", h.ch, le_te_dt_ns, te_te_dt_ns, fCount);
               } else {
                  if (complain)
                     printf("TTT: ch %d: LE-TE-TE, TE FROM WRONG EVENT, LE-TE %.3f, TE-TE %.3f ns, count %d\n", h.ch, le_te_dt_ns, te_te_dt_ns, fCount);
               }
               //fLe.Print(); printf("\n");
               //fTe.Print(); printf("\n");
               //h.Print(); printf("\n");
            }
         }
      }
   };

   void Print() const
   {
      printf("DlTdcHit2: u/d %d%d, time_sec %.9f, w_ns %.3f\n", fUp, fDown, fTimeSec, fWidthNs);
      fLe.Print();
      printf("\n");
      fTe.Print();
      printf("\n");
   }
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

#define CHAN1 (fConf->fChan1)
#define CHAN2 (fConf->fChan2)
#define CHAN3 (fConf->fChan3)
#define CHAN4 (fConf->fChan4)
#define CHAN5 (fConf->fChan5)
#define CHAN6 (fConf->fChan6)
#define CHAN7 (fConf->fChan7)
#define CHAN8 (fConf->fChan8)


class DlTdcConfig
{
public:
   std::vector<double> fLeOffset;
   std::vector<double> fTeOffset;

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
   DlTdcConfig(int num);
   bool ReadJson(int runno);
};

DlTdcConfig::DlTdcConfig(int num)
{
   fLeOffset.resize(num);
   fTeOffset.resize(num);
}

bool DlTdcConfig::ReadJson(int runno)
{
   if (runno >= 17 && runno < 900000) {
      printf("NEW MAP!\n");
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
   }

   return true;
};

class DlTdcEvent
{
public: // matching with ADC data
   double time_sec = 0;
   double dt = 0;

public: // time of first and last TDC hits
   double first_time_sec = 0;
   double last_time_sec = 0;

   double min_time_sec = 0;
   double max_time_sec = 0;

public: // TDC hits
   DlTdcHit2 fHits[MAX_TDC_CHAN+1];

public:
   void Clear()
   {
      first_time_sec = 0;
      last_time_sec = 0;

      min_time_sec = 0;
      max_time_sec = 0;

      for (int i=0; i<=MAX_TDC_CHAN; i++) {
         fHits[i].Clear();
      }
   }

   void AddHit8(const DlTdcHit& h)
   {
      if (first_time_sec == 0) {
         first_time_sec = h.time_sec;
         min_time_sec = h.time_sec;
         max_time_sec = h.time_sec;
      }

      last_time_sec = h.time_sec;

      if (h.time_sec < min_time_sec)
         min_time_sec = h.time_sec;
      
      if (h.time_sec > max_time_sec)
         max_time_sec = h.time_sec;

      assert(h.ch >= 0);
      assert(h.ch <= MAX_TDC_CHAN);

      fHits[h.ch].AddHit(h);
   }

   bool HaveCh(int ch) const
   {
      return fHits[ch].fDown;
   }

   const DlTdcHit2& GetCh(int ch) const
   {
      return fHits[ch];
   }
};

class DlTdcModule: public TARunObject
{
public:
   DlTdcFlags* fFlags = NULL;
   DlTdcUnpack* fU = NULL;
   
#ifdef HAVE_ROOT
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

   DlTdcEvent *fCt = NULL;

   DlTdcConfig *fConf = NULL;

   bool fTrace = false;

   Ncfm* fCfm = NULL;
   
   DlTdcModule(TARunInfo* runinfo, DlTdcFlags* flags)
      : TARunObject(runinfo)
   {
      if (fTrace)
         printf("DlTdcModule::ctor!\n");

      fModuleName = "dltdc_module";
      fFlags   = flags;

      fCfm = new Ncfm("dlcfmdb");

      fU = new DlTdcUnpack(NUM_TDC_CHAN);

      fConf = new DlTdcConfig(NUM_TDC_CHAN);
   }

   ~DlTdcModule()
   {
      if (fTrace)
         printf("DlTdcModule::dtor!\n");

      if (fU) {
         delete fU;
         fU = NULL;
      }

      if (fConf) {
         delete fConf;
         fConf = NULL;
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

      bool conf_ok = fConf->ReadJson(runinfo->fRunNo);
      if (!conf_ok) {
         printf("Cannot load TDC configuration for run %d\n", runinfo->fRunNo);
         exit(123);
      }

      if (!fFlags->fCalibFineTime) {
         std::string fine_time_json = fCfm->GetFilename("dltdc", "finetime", runinfo->fRunNo, "json");
         bool load_ok = fU->LoadCalib(fine_time_json.c_str());

         printf("json file %s load_ok %d\n", fine_time_json.c_str(), load_ok);

         //bool load_ok = fU->LoadCalib(runinfo->fRunNo);
         if (!load_ok) {
            printf("Cannot load TDC fine time calibration for run %d\n", runinfo->fRunNo);
            exit(123);
         }
      }

      if (!fFlags->fCalibOffsets) {

         if (runinfo->fRunNo > 906000) {

            // mean of corresponding tNM_ns plot with same sign
            double xt12 = 3; // pulser run 906008.
            double xt34 = 2.159;
            double xt56 = 5.333;
            double xt78 = -3.452;
            double xt14 = -0.774;
            double xt58 = -6.328;
            double xt1458 = -2.9;

            if (runinfo->fRunNo >= 906028) {
               // pass1
               xt12 = 2.792; // pulser run 906028.
               xt34 = 3.209;
               xt56 = 2.690;
               xt78 = -2.563;
               // pass2
               xt14 = -0.069;
               xt58 = -2.337;
               // pass3
               xt1458 = -2.614;
            }

            if (runinfo->fRunNo >= 906041) {
               // pass1
               xt12 = 2.792-0.556224+0.042; // pulser run 906042.
               xt34 = 3.209+0.029908-0.074;
               xt56 = 2.690-1.881289+0.737;
               xt78 = -2.563-0.028077+0.193;
               // pass2
               xt14 = -0.069-0.332+0.748-0.042;
               xt58 = -2.337-0.737+1.400-0.737;
               // pass3
               xt1458 = -2.614-0.509-0.138+0.674-0.348;
            }

            fU->fCalib[0].lepos.fOffsetNs = fU->fCalib[0].leneg.fOffsetNs = -3.865; // A
            fU->fCalib[1].lepos.fOffsetNs = fU->fCalib[1].leneg.fOffsetNs = -1.578; // B
            fU->fCalib[2].lepos.fOffsetNs = fU->fCalib[2].leneg.fOffsetNs = xt12 + xt14 + xt1458; // chan1
            fU->fCalib[3].lepos.fOffsetNs = fU->fCalib[3].leneg.fOffsetNs = xt14 + xt1458; // chan2
            fU->fCalib[4].lepos.fOffsetNs = fU->fCalib[4].leneg.fOffsetNs = xt34 + xt1458; // chan3
            fU->fCalib[5].lepos.fOffsetNs = fU->fCalib[5].leneg.fOffsetNs = xt1458; // chan4
            fU->fCalib[6].lepos.fOffsetNs = fU->fCalib[6].leneg.fOffsetNs = xt56 + xt58; // chan5
            fU->fCalib[7].lepos.fOffsetNs = fU->fCalib[7].leneg.fOffsetNs = xt58; // chan6
            fU->fCalib[8].lepos.fOffsetNs = fU->fCalib[8].leneg.fOffsetNs = xt78; // chan7
            fU->fCalib[9].lepos.fOffsetNs = fU->fCalib[9].leneg.fOffsetNs = 0; // chan8
            fU->fCalib[10].lepos.fOffsetNs = fU->fCalib[10].leneg.fOffsetNs = -4.739; // T
            fU->fCalib[11].lepos.fOffsetNs = fU->fCalib[11].leneg.fOffsetNs = 0; // nc

            // number from width plot with opposite sign
            // cosmic run 906005
            fU->fCalib[0].tepos.fOffsetNs  = fU->fCalib[0].teneg.fOffsetNs  = fU->fCalib[0].lepos.fOffsetNs  -0.100 +3.905; // A
            fU->fCalib[1].tepos.fOffsetNs  = fU->fCalib[1].teneg.fOffsetNs  = fU->fCalib[1].lepos.fOffsetNs  +1.000 +0.400 -0.120; // B
            fU->fCalib[2].tepos.fOffsetNs  = fU->fCalib[2].teneg.fOffsetNs  = fU->fCalib[2].lepos.fOffsetNs  -0.600 +1.100; // chan1
            fU->fCalib[3].tepos.fOffsetNs  = fU->fCalib[3].teneg.fOffsetNs  = fU->fCalib[3].lepos.fOffsetNs  +0.100 +0.300 +3.358; // chan2
            fU->fCalib[4].tepos.fOffsetNs  = fU->fCalib[4].teneg.fOffsetNs  = fU->fCalib[4].lepos.fOffsetNs  +0.300 +0.600 + 0.400 -1.361; // chan3
            fU->fCalib[5].tepos.fOffsetNs  = fU->fCalib[5].teneg.fOffsetNs  = fU->fCalib[5].lepos.fOffsetNs  -0.300 +0.100 +3.923; // chan4
            fU->fCalib[6].tepos.fOffsetNs  = fU->fCalib[6].teneg.fOffsetNs  = fU->fCalib[6].lepos.fOffsetNs  -0.700 +0.900 -0.867; // chan5
            fU->fCalib[7].tepos.fOffsetNs  = fU->fCalib[7].teneg.fOffsetNs  = fU->fCalib[7].lepos.fOffsetNs  +0.500 +0.300 +3.373; // chan6
            fU->fCalib[8].tepos.fOffsetNs  = fU->fCalib[8].teneg.fOffsetNs  = fU->fCalib[8].lepos.fOffsetNs  +0.400 +0.100 +3.412; // chan7
            fU->fCalib[9].tepos.fOffsetNs  = fU->fCalib[9].teneg.fOffsetNs  = fU->fCalib[9].lepos.fOffsetNs  -0.100 +0.500 -0.583; // chan8
            fU->fCalib[10].tepos.fOffsetNs = fU->fCalib[10].teneg.fOffsetNs = fU->fCalib[10].lepos.fOffsetNs +0.500 +3.951; // T
            fU->fCalib[11].tepos.fOffsetNs = fU->fCalib[11].teneg.fOffsetNs = fU->fCalib[11].lepos.fOffsetNs +0; // nc

            if (runinfo->fRunNo >= 906095) {
               // number from pulser tdc chanNN LE table:
               fU->fCalib[0].lepos.fOffsetNs  = fU->fCalib[0].leneg.fOffsetNs  -= -0.650; // A
               fU->fCalib[1].lepos.fOffsetNs  = fU->fCalib[1].leneg.fOffsetNs  -=  0.315; // B
               fU->fCalib[2].lepos.fOffsetNs  = fU->fCalib[2].leneg.fOffsetNs  -=  0.000; // chan1
               fU->fCalib[3].lepos.fOffsetNs  = fU->fCalib[3].leneg.fOffsetNs  -= -0.424; // chan2
               fU->fCalib[4].lepos.fOffsetNs  = fU->fCalib[4].leneg.fOffsetNs  -=  0.559; // chan3
               fU->fCalib[5].lepos.fOffsetNs  = fU->fCalib[5].leneg.fOffsetNs  -=  0.025; // chan4
               fU->fCalib[6].lepos.fOffsetNs  = fU->fCalib[6].leneg.fOffsetNs  -=  0.234; // chan5
               fU->fCalib[7].lepos.fOffsetNs  = fU->fCalib[7].leneg.fOffsetNs  -= -0.140; // chan6
               fU->fCalib[8].lepos.fOffsetNs  = fU->fCalib[8].leneg.fOffsetNs  -= -0.417; // chan7
               fU->fCalib[9].lepos.fOffsetNs  = fU->fCalib[9].leneg.fOffsetNs  -= -1.040; // chan8
               fU->fCalib[10].lepos.fOffsetNs = fU->fCalib[10].leneg.fOffsetNs -= -1.330; // T
               fU->fCalib[11].lepos.fOffsetNs = fU->fCalib[11].leneg.fOffsetNs = 0; // nc

               // number from pulser tdc chanNN TE table:
               fU->fCalib[0].tepos.fOffsetNs  = fU->fCalib[0].teneg.fOffsetNs  -= -0.057 + 0.442; // A
               fU->fCalib[1].tepos.fOffsetNs  = fU->fCalib[1].teneg.fOffsetNs  -= -0.313  -0.761; // B
               fU->fCalib[2].tepos.fOffsetNs  = fU->fCalib[2].teneg.fOffsetNs  -=  0.000; // chan1
               fU->fCalib[3].tepos.fOffsetNs  = fU->fCalib[3].teneg.fOffsetNs  -= -1.352 + 1.405; // chan2
               fU->fCalib[4].tepos.fOffsetNs  = fU->fCalib[4].teneg.fOffsetNs  -=  1.250  -0.538; // chan3
               fU->fCalib[5].tepos.fOffsetNs  = fU->fCalib[5].teneg.fOffsetNs  -= -0.460 + 0.616; // chan4
               fU->fCalib[6].tepos.fOffsetNs  = fU->fCalib[6].teneg.fOffsetNs  -=  2.412  -1.914; // chan5
               fU->fCalib[7].tepos.fOffsetNs  = fU->fCalib[7].teneg.fOffsetNs  -= -1.652 + 1.754; // chan6
               fU->fCalib[8].tepos.fOffsetNs  = fU->fCalib[8].teneg.fOffsetNs  -= -2.276 + 1.783; // chan7
               fU->fCalib[9].tepos.fOffsetNs  = fU->fCalib[9].teneg.fOffsetNs  -= -2.919 + 0.857 -0.148; // chan8
               fU->fCalib[10].tepos.fOffsetNs = fU->fCalib[10].teneg.fOffsetNs -= -1.253 + 0.148 -0.148; // T
               fU->fCalib[11].tepos.fOffsetNs = fU->fCalib[11].teneg.fOffsetNs = 0; // nc
            }
         } else {
#if 0
            if (runinfo->fRunNo >= 22) {
               // number from pulser tdc chanNN LE table:
               fU->fCalib[16].lepos.fOffsetNs = fU->fCalib[16].leneg.fOffsetNs  =  0;     // chan1
               fU->fCalib[17].lepos.fOffsetNs = fU->fCalib[17].leneg.fOffsetNs -= -0.359; // chan2
               fU->fCalib[22].lepos.fOffsetNs = fU->fCalib[22].leneg.fOffsetNs -= -1.495; // chan3
               fU->fCalib[23].lepos.fOffsetNs = fU->fCalib[23].leneg.fOffsetNs -=  0.123; // chan4
               fU->fCalib[26].lepos.fOffsetNs = fU->fCalib[26].leneg.fOffsetNs -=  0.435; // chan5
               fU->fCalib[27].lepos.fOffsetNs = fU->fCalib[27].leneg.fOffsetNs -= -1.175; // chan6
               fU->fCalib[30].lepos.fOffsetNs = fU->fCalib[30].leneg.fOffsetNs -=  0.614; // chan7
               fU->fCalib[31].lepos.fOffsetNs = fU->fCalib[31].leneg.fOffsetNs -= -0.731; // chan8
               fU->fCalib[32].lepos.fOffsetNs = fU->fCalib[32].leneg.fOffsetNs -=  2.059; // chanA
               fU->fCalib[33].lepos.fOffsetNs = fU->fCalib[33].leneg.fOffsetNs -=  0.893; // chanB
               fU->fCalib[34].lepos.fOffsetNs = fU->fCalib[34].leneg.fOffsetNs -=  2.238; // chanT

               // number from pulser tdc chanNN TE table:
               fU->fCalib[16].tepos.fOffsetNs = fU->fCalib[16].teneg.fOffsetNs  =  0;     // chan1
               fU->fCalib[17].tepos.fOffsetNs = fU->fCalib[17].teneg.fOffsetNs -= -1.918; // chan2
               fU->fCalib[22].tepos.fOffsetNs = fU->fCalib[22].teneg.fOffsetNs -= -1.714; // chan3
               fU->fCalib[23].tepos.fOffsetNs = fU->fCalib[23].teneg.fOffsetNs -= -1.888; // chan4
               fU->fCalib[26].tepos.fOffsetNs = fU->fCalib[26].teneg.fOffsetNs -= -1.589; // chan5
               fU->fCalib[27].tepos.fOffsetNs = fU->fCalib[27].teneg.fOffsetNs -= -2.000; // chan6
               fU->fCalib[30].tepos.fOffsetNs = fU->fCalib[30].teneg.fOffsetNs -= -1.211; // chan7
               fU->fCalib[31].tepos.fOffsetNs = fU->fCalib[31].teneg.fOffsetNs -= -1.239; // chan8
               fU->fCalib[32].tepos.fOffsetNs = fU->fCalib[32].teneg.fOffsetNs -= -0.729; // chanA
               fU->fCalib[33].tepos.fOffsetNs = fU->fCalib[33].teneg.fOffsetNs -=  0.373; // chanB
               fU->fCalib[34].tepos.fOffsetNs = fU->fCalib[34].teneg.fOffsetNs -= -0.168; // chanT
            }
#endif

            std::string offset_json = fCfm->GetFilename("dltdc", "offsets", runinfo->fRunNo, "json");

            bool load_ok = fU->LoadOffsets(offset_json.c_str());
            
            printf("json file %s load_ok %d\n", offset_json.c_str(), load_ok);
            
            if (!load_ok) {
               printf("Cannot load TDC offset calibration for run %d\n", runinfo->fRunNo);
               exit(123);
            }
         }
      }

#ifdef HAVE_ROOT
      runinfo->fRoot->fOutputFile->cd(); // select correct ROOT directory
      TDirectory* dir = gDirectory->mkdir("dltdc");
      dir->cd(); // select correct ROOT directory

      dir->mkdir("fine_time")->cd();

      for (int i=0; i<=MAX_TDC_CHAN; i++) {
         char name[256];
         char title[256];

         sprintf(name,  "tdc%02d_phase_le", i);
         sprintf(title, "tdc%02d_phase_le", i);
         fHphaseLe[i] = new TH1D(name, title, 101, -50, 50);

         sprintf(name,  "tdc%02d_phase_te", i);
         sprintf(title, "tdc%02d_phase_te", i);
         fHphaseTe[i] = new TH1D(name, title, 101, -50, 50);

         sprintf(name,  "tdc%02d_fine_le", i);
         sprintf(title, "tdc%02d_fine_le, ns", i);
         fHfineLe[i] = new TH1D(name, title, 200, -5, 15);

         sprintf(name,  "tdc%02d_fine_te", i);
         sprintf(title, "tdc%02d_fine_te, ns", i);
         fHfineTe[i] = new TH1D(name, title, 200, -5, 15);

         sprintf(name,  "tdc%02d_width_ns", i);
         sprintf(title, "tdc%02d_width_ns, ns", i);
         fHwidth[i] = new TH1D(name, title, 100, -5, 5);
      }

      dir->mkdir("pulser")->cd();

      fHpulserLeAll = new TH1D("tdc_pulser_le_all", "tdc pulser le all, ns", 400, -40, 40);
      fHpulserTeAll = new TH1D("tdc_pulser_te_all", "tdc_pulser_te_all, ns", 400, -40, 40);
      fHpulserWiAll = new TH1D("tdc_pulser_width_all", "tdc_pulser_width_all, ns", 400, 0, 80);

      for (int i=0; i<=MAX_TDC_CHAN; i++) {
         char name[256];
         char title[256];

         sprintf(name,  "tdc%02d_pulser_le", i);
         sprintf(title, "tdc%02d_pulser_le, ns", i);
         fHpulserLe[i] = new TH1D(name, title, 400, -40, 40);

         sprintf(name,  "tdc%02d_pulser_te", i);
         sprintf(title, "tdc%02d_pulser_te, ns", i);
         fHpulserTe[i] = new TH1D(name, title, 400, -40, 40);

         sprintf(name,  "tdc%02d_pulser_width", i);
         sprintf(title, "tdc%02d_pulser_width, ns", i);
         fHpulserWi[i] = new TH1D(name, title, 400, 0, 80);
      }

      dir->mkdir("poisson")->cd();

      fHhitdt1ns = new TH1D("hitdt1ns", "hit dt 100 ns", 100, 0, 100); // 100 ns
      fHhitdt2ns = new TH1D("hitdt2ns", "hit dt 1000 ns", 100, 0, 1000); // 1 usec
      fHhitdt3ns = new TH1D("hitdt3ns", "hit dt 1000 us", 100, 0, 1000000); // 1 msec
      fHhitdt4ns = new TH1D("hitdt4ns", "hit dt 1000 ms", 100, 0, 1000000000); // 1 sec

      fHeventdt1ns = new TH1D("eventdt1ns", "event dt 100 ns", 100, 0, 100); // 100 ns
      fHeventdt2ns = new TH1D("eventdt2ns", "event dt 1000 ns", 100, 0, 1000); // 1 usec
      fHeventdt3ns = new TH1D("eventdt3ns", "event dt 1000 us", 100, 0, 1000000); // 1 msec
      fHeventdt4ns = new TH1D("eventdt4ns", "event dt 1000 ms", 100, 0, 1000000000); // 1 sec

      for (int i=0; i<=MAX_TDC_CHAN; i++) {
         char name[256];
         char title[256];

         sprintf(name,  "tdc%02d_dt1_le", i);
         sprintf(title, "tdc%02d_dt1_le, ns", i);
         fHdt1Le[i] = new TH1D(name, title, 100, 0, 100); // 10 ns

         sprintf(name,  "tdc%02d_dt2_le", i);
         sprintf(title, "tdc%02d_dt2_le, ns", i);
         fHdt2Le[i] = new TH1D(name, title, 100, 0, 1000); // 1 usec

         sprintf(name,  "tdc%02d_dt3_le", i);
         sprintf(title, "tdc%02d_dt3_le, ns", i);
         fHdt3Le[i] = new TH1D(name, title, 100, 0, 1000000); // 1 msec

         sprintf(name,  "tdc%02d_dt4_le", i);
         sprintf(title, "tdc%02d_dt4_le, ns", i);
         fHdt4Le[i] = new TH1D(name, title, 100, 0, 1000000000); // 1 sec
      }

      dir->mkdir("unphysical")->cd();

      fHunphysical_pair_ns.resize((MAX_TDC_CHAN+1)*(MAX_TDC_CHAN+1));
      for (int i=0; i<=MAX_TDC_CHAN; i++) {
         for (int j=i+1; j<=MAX_TDC_CHAN; j++) {
            char name[256];
            char title[256];

            sprintf(name,  "tdc%02d_tdc%02d_le_ns", i, j);
            sprintf(title, "tdc%02d_le - tdc%02d_le, ns", i, j);
            fHunphysical_pair_ns[i*(MAX_TDC_CHAN+1) + j] = new TH1D(name, title, 200, -50, 50);
            fHunphysical_pair_ns[i*(MAX_TDC_CHAN+1) + j]->SetMinimum(0);
         }
      }

      dir->cd();

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

      if (fFlags->fCalibFineTime) {
         printf("DlTdcModule::EndRun: Saving TDC fine time calibrations for run %d\n", runinfo->fRunNo);
         fU->UpdateCalib();
         fU->SaveCalib(runinfo->fRunNo);
      }

      if (fFlags->fCalibOffsets) {
         printf("DlTdcModule::EndRun: Saving TDC offsets calibrations for run %d\n", runinfo->fRunNo);
         //fU->SaveOffsets(runinfo->fRunNo);
      }

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

      printf("Run %d pulser:\n", runinfo->fRunNo);

      std::string s;
      s += "{\n";
      s += "\"dltdc_num_channels\":";
      s += toString(MAX_TDC_CHAN+1);
      s += ",\n";
      s += "\"offset_ns\":[\n";

      for (size_t ch=0; ch<=MAX_TDC_CHAN; ch++) {
         printf("tdc chan%02zu: LE %8.3f, TE %8.3f, Width: %8.3f, RMS %.3f %.3f %.3f\n",
                ch,
                fHpulserLe[ch]->GetMean(),
                fHpulserTe[ch]->GetMean(),
                fHpulserWi[ch]->GetMean(),
                fHpulserLe[ch]->GetRMS(),
                fHpulserTe[ch]->GetRMS(),
                fHpulserWi[ch]->GetRMS());

         if (ch > 0)
            s += ",\n";
         s += "[";
         s += toString("%8.3f, ", fHpulserLe[ch]->GetMean());
         s += toString("%8.3f, ", fHpulserTe[ch]->GetMean());
         s += toString("%8.3f, ", fHpulserWi[ch]->GetMean());
         s += "   ";
         s += toString("%.3f, ",  fHpulserLe[ch]->GetRMS());
         s += toString("%.3f, ",  fHpulserTe[ch]->GetRMS());
         s += toString("%.3f",    fHpulserWi[ch]->GetRMS());
         s += "]";
      }

      s += "\n";
      s += "],\n";

         printf("tdc sumary: LE %8.3f, LE %8.3f, Width: %8.3f, RMS %.3f %.3f %.3f\n",
                fHpulserLeAll->GetMean(),
                fHpulserTeAll->GetMean(),
                fHpulserWiAll->GetMean(),
                fHpulserLeAll->GetRMS(),
                fHpulserTeAll->GetRMS(),
                fHpulserWiAll->GetRMS());

      s += "\"offset_all_ns\":[\n";
      s += " ";
      s += toString("%8.3f, ", fHpulserLeAll->GetMean());
      s += toString("%8.3f, ", fHpulserTeAll->GetMean());
      s += toString("%8.3f, ", fHpulserWiAll->GetMean());
      s += "   ";
      s += toString("%.3f, ",  fHpulserLeAll->GetRMS());
      s += toString("%.3f, ",  fHpulserTeAll->GetRMS());
      s += toString("%.3f",    fHpulserWiAll->GetRMS());
      s += "]\n";
      s += "},\n";

      printf("offsets json:\n%s\n", s.c_str());

      if (fFlags->fCalibOffsets) {
         char fname[256];
         sprintf(fname, "dlcfmdb/dltdc_offsets_%06d.json", runinfo->fRunNo);
         FILE *fp = fopen(fname, "w");
         if (fp) {
            fprintf(fp, "%s", s.c_str());
            fclose(fp);
         }
      }

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

   void FinishEventT(double prev_event_time_sec, const DlTdcEvent& t)
   {
      if (fFlags->fDebug) {
         printf("EVENT %d %d %d %d %d %d %d %d, ABT %d%d%d\n", t.HaveCh(CHAN1), t.HaveCh(CHAN2), t.HaveCh(CHAN3), t.HaveCh(CHAN4), t.HaveCh(CHAN5), t.HaveCh(CHAN6), t.HaveCh(CHAN7), t.HaveCh(CHAN8), t.HaveCh(fConf->fChanA), t.HaveCh(fConf->fChanB), t.HaveCh(fConf->fChanT));
      }
      
      ///////// check for triggered event ///////////

      if (fFlags->fTriggered && !t.HaveCh(fConf->fChanT)) {
         return;
      }

      ///////// pulser calibration histograms ///////////
         
      size_t pulserMaster = CHAN1;
      bool   pulserMasterHave = t.HaveCh(pulserMaster);

      if (pulserMasterHave) {
         for (size_t ch=0; ch<=MAX_TDC_CHAN; ch++) {
            //printf("ch %d, up down %d %d\n", ch, t.fHits[ch].fUp, t.fHits[ch].fDown);
            if (t.HaveCh(ch)) {
               if (ch != pulserMaster) {
                  fHpulserLeAll->Fill(subtract_ns(t.GetCh(ch).fLe, t.GetCh(pulserMaster).fLe));
                  fHpulserTeAll->Fill(subtract_ns(t.GetCh(ch).fTe, t.GetCh(pulserMaster).fTe));

                  fHpulserLe[ch]->Fill(subtract_ns(t.GetCh(ch).fLe, t.GetCh(pulserMaster).fLe));
                  fHpulserTe[ch]->Fill(subtract_ns(t.GetCh(ch).fTe, t.GetCh(pulserMaster).fTe));
               }

               fHpulserWiAll->Fill(t.fHits[ch].fWidthNs);
               fHpulserWi[ch]->Fill(t.fHits[ch].fWidthNs);
            }
         }
      }

      ///////// create width calibration histograms ///////////

      for (size_t ch=0; ch<=MAX_TDC_CHAN; ch++) {
         //printf("ch %d, up down %d %d\n", ch, t.fHits[ch].fUp, t.fHits[ch].fDown);
         if (!t.fHits[ch].fUp && t.fHits[ch].fDown) {
            fHwidth[ch]->Fill(t.fHits[ch].fWidthNs);
         }
      }

      ///////// plot trigger A-B time ///////////
      
      if (t.HaveCh(fConf->fChanA) && t.HaveCh(fConf->fChanB)) {
         double tAB_ns = subtract_ns(t.GetCh(fConf->fChanA).fLe, t.GetCh(fConf->fChanB).fLe);
         fHtABns->Fill(tAB_ns);
      }

      //if (!t.HaveCh(fConf->fChanB)) return;

      // cut on t6-t7
      
      //if (t.HaveCh(CHAN6) && t.HaveCh(CHAN7)) {
      //   double t67_ns = subtract_ns(t.GetCh(CHAN7).fLe, t.GetCh(CHAN6).fLe);
      //
      //   //if (t67_ns > -0.6) return;
      //   if (t67_ns < -0.4) return;
      //} else {
      //   return;
      //}

      ///////// plot unphysical time pairs /////////

      for (int i=0; i<=MAX_TDC_CHAN; i++) {
         for (int j=i+1; j<=MAX_TDC_CHAN; j++) {
            if (t.HaveCh(i) && t.HaveCh(j)) {
               fHunphysical_pair_ns[i*(MAX_TDC_CHAN+1) + j]->Fill(subtract_ns(t.GetCh(i).fLe, t.GetCh(j).fLe));
            }
         }
      }

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
      double W_twc = 8.0;

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

      if (t.fHits[fConf->fChanA].fDown) {
         double wA_ns = t.fHits[fConf->fChanA].fWidthNs;

         if (wA_ns < 0.01) {
            printf("WWW: BAD WIDTH chanA %f!\n", wA_ns);
         }

         if (fFlags->fPrint) {
            printf("new dlsc event A, le %.9f, w %.0f!\n", t.GetCh(fConf->fChanA).fLe.time_sec, wA_ns);
         }

         fHwAns->Fill(wA_ns);

         fCountA++;
      }

      if (t.fHits[fConf->fChanB].fDown) {
         double wB_ns = t.fHits[fConf->fChanB].fWidthNs;

         if (wB_ns < 0.01) {
            printf("WWW: BAD WIDTH chanB %f!\n", wB_ns);
         }

         if (fFlags->fPrint) {
            printf("new dlsc event B, le %.9f, w %.0f!\n", t.GetCh(fConf->fChanB).fLe.time_sec, wB_ns);
         }

         fHwBns->Fill(wB_ns);

         fCountB++;
      }

      if (t.fHits[fConf->fChanT].fDown) {
         double wT_ns = t.fHits[fConf->fChanT].fWidthNs;

         if (wT_ns < 0.01) {
            printf("WWW: BAD WIDTH chanT %f!\n", wT_ns);
         }

         if (fFlags->fPrint) {
            printf("new dlsc event T, le %.9f, w %.0f!\n", t.GetCh(fConf->fChanT).fLe.time_sec, wT_ns);
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

   double fLastXtime = 0;

   bool fFirstTrig = true;

   std::vector<DlTdcHit> fPrevLe;

   TAFlowEvent* Analyze(TARunInfo* runinfo, TMEvent* event, TAFlags* flags, TAFlowEvent* flow)
   {
      //printf("DlTdcModule::Analyze, run %d, event serno %d, id 0x%04x, data size %d\n", runinfo->fRunNo, event->serial_number, (int)event->event_id, event->data_size);

      if (event->event_id != 4)
         return flow;

      if (!fFlags->fEnabled)
         return flow;

      if (fPrevLe.size() != MAX_TDC_CHAN+1) {
         //printf("RESIZE!\n");
         fPrevLe.resize(MAX_TDC_CHAN+1);
      }
      
      TMBank* tdcbank = event->FindBank("CBT2");

      bool calib = fFlags->fCalibFineTime;

      if (tdcbank) {
         int tdc_nw64 = tdcbank->data_size/8;
         uint32_t* tdc_data = (uint32_t*)event->GetBankData(tdcbank);
         //printf("nw: %d\n", tdc_nw64);

         if (fCt == NULL) {
            fCt = new DlTdcEvent;
         }
     
	 for (int i=0; i<tdc_nw64; i++) {
            uint32_t wlo = tdc_data[i*2+0];
            uint32_t whi = tdc_data[i*2+1];

            //printf("%2d: 0x%08x 0x%08x\n", i, whi, wlo);

            DlTdcHit h;
            fU->Unpack(&h, wlo, whi);

            if (fFlags->fDebug) {
               h.Print(); printf("\n");
            }

            //static double prev_time_sec = 0;
            //
            //if (h.le) {
            //   if (fFlags->fDebug) {h.Print(); printf(", dt %4.0f ns\n", sec_to_ns(h.time_sec - prev_time_sec));}
            //   prev_time_sec = h.time_sec;
            //}

            if (calib) {
               fU->fCalib[h.ch].AddHit(h);
            }

            if (h.le) {
               double dt = subtract_ns(h, fPrevLe[h.ch]);
               //printf("ch %d: dt %.3f ns\n", h.ch, dt);
               fHdt1Le[h.ch]->Fill(dt);
               fHdt2Le[h.ch]->Fill(dt);
               fHdt3Le[h.ch]->Fill(dt);
               fHdt4Le[h.ch]->Fill(dt);
               fPrevLe[h.ch] = h;
            }

            double hit_dt_ns = sec_to_ns(h.time_sec - fCt->max_time_sec);

            fHhitdt1ns->Fill(hit_dt_ns);
            fHhitdt2ns->Fill(hit_dt_ns);
            fHhitdt3ns->Fill(hit_dt_ns);
            fHhitdt4ns->Fill(hit_dt_ns);

            //printf("TTX %.9f %.9f sec, first %.9f, last %.9f sec\n", fCt->min_time_sec, fCt->max_time_sec, fCt->first_time_sec, fCt->last_time_sec);

            if (h.le && hit_dt_ns > 80.0) {
               double event_dt_ns = sec_to_ns(h.time_sec - fCt->min_time_sec);


               if (0) {
                  printf("finish event ch %d, lete %d%d, time %.9f sec, max %.9f sec, dt %.3f ns\n", h.ch, h.le, h.te, h.time_sec, fCt->max_time_sec, event_dt_ns);
                  //h.Print();
                  //printf("===\n");
               }

               //for (size_t ch=0; ch<MAX_TDC_CHAN; ch++) {
               //   if (fCt->fHits[ch].fUp) {
               //      printf("TTT: ch %zu no TE\n", ch);
               //   }
               //}

               if (1) {
                  //printf("TTT %.9f -> %.9f sec, dt %.3f ns\n", fCt->min_time_sec, h.time_sec, event_dt_ns);
                  fHeventdt1ns->Fill(event_dt_ns);
                  fHeventdt2ns->Fill(event_dt_ns);
                  fHeventdt3ns->Fill(event_dt_ns);
                  fHeventdt4ns->Fill(event_dt_ns);
                  
                  FinishEventT(fPrevEventTimeSec, *fCt);
                  
                  fPrevEventTimeSec = fCt->min_time_sec;
                  
                  fCt->Clear();
               }
            }

#ifdef HAVE_ROOT
            if (h.ch >= 0 && h.ch <= MAX_TDC_CHAN) {
               if (h.le) {
                  fHphaseLe[h.ch]->Fill(h.phase);
                  fHfineLe[h.ch]->Fill(h.fine_ns);
               }
               if (h.te) {
                  fHphaseTe[h.ch]->Fill(h.phase);
                  fHfineTe[h.ch]->Fill(h.fine_ns);
               }
            }
#endif
            fCt->AddHit8(h);

         } // loop over data

         time_t now = time(NULL);
         static time_t last = 0;
         
         if (now - last > 5) {
            if (calib) {
               for (auto& c: fU->fCalib) {
                  c.Update();
                  if (fFlags->fDebug) c.Print();
               }
            }

            last = time(NULL);
         }
         
      }
      
      return flow;
   }

   double fPrevFlowTime = 0;

   TAFlowEvent* AnalyzeFlowEvent(TARunInfo* runinfo, TAFlags* flags, TAFlowEvent* flow)
   {
      return flow;
   }

   void AnalyzeSpecialEvent(TARunInfo* runinfo, TMEvent* event)
   {
      if (fTrace)
         printf("DlTdcModule::AnalyzeSpecialEvent, run %d, event serno %d, id 0x%04x, data size %d\n", runinfo->fRunNo, event->serial_number, (int)event->event_id, event->data_size);
   }
};

class DlTdcModuleFactory: public TAFactory
{
public:
   DlTdcFlags fFlags;

public:
   void Usage()
   {
      printf("DlTdcModuleFactory flags:\n");
      printf("--dltdc -- enable dltdc code (default)\n");
      printf("--no-dltdc -- disable dltdc code\n");
      printf("--dltdc-finetime -- calibrate dltdc fine time\n");
      printf("--dltdc-offsets  -- calibrate dltdc offsets\n");
      printf("--dltdc-debug -- print detailed information\n");
      printf("--dltdc-print -- print events\n");
      printf("--dltdc-twc -- get twc W parameter");
   }

   void Init(const std::vector<std::string> &args)
   {
      printf("DlTdcModuleFactory::Init!\n");

      for (unsigned i=0; i<args.size(); i++) {
         if (args[i] == "--no-dltdc") {
            fFlags.fEnabled = false;
         }
         if (args[i] == "--dltdc") {
            fFlags.fEnabled = true;
         }
         if (args[i] == "--dltdc-finetime") {
            fFlags.fCalibFineTime = true;
         }
         if (args[i] == "--dltdc-offsets") {
            fFlags.fCalibOffsets = true;
         }
         if (args[i] == "--dltdc-triggered") {
            fFlags.fTriggered = true;
         }
         if (args[i] == "--dltdc-debug") {
            fFlags.fDebug = true;
         }
         if (args[i] == "--dltdc-print") {
            fFlags.fPrint = true;
         }
         if (args[i] == "--dltdc-twc") {
            fFlags.fTWC = true;
         }
      }
   }

   void Finish()
   {
      printf("DlTdcModuleFactory::Finish!\n");
   }
   
   TARunObject* NewRunObject(TARunInfo* runinfo)
   {
      printf("DlTdcModuleFactory::NewRunObject, run %d, file %s\n", runinfo->fRunNo, runinfo->fFileName.c_str());
      return new DlTdcModule(runinfo, &fFlags);
   }
};

static TARegister tar(new DlTdcModuleFactory);

/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */
