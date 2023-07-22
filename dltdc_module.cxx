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

#include "AgFlow.h"

#include <deque>

#ifdef HAVE_ROOT
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#endif

static double amin(double a, double b)
{
   if (a < b)
      return a;
   else
      return b;
}

static double sec_to_ns(double t)
{
   return t*1e9;
}

class DlTdcFlags
{
public:
   bool fEnabled = false;
   bool fCalib   = false;
   bool fHaveAdc = false;
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
   DlTdcHit h0le;
   DlTdcHit h0te;
   
   DlTdcHit h1le;
   DlTdcHit h1te;
   
   //DlTdcHit h2le;
   //DlTdcHit h2te;
   
   //DlTdcHit h3le;
   //DlTdcHit h3te;
   
   DlTdcHit h4le;
   DlTdcHit h4te;
   
   DlTdcHit h5le;
   DlTdcHit h5te;

   DlTdcHit h8le;
   DlTdcHit h8te;

   DlTdcHit h9le;
   DlTdcHit h9te;

   DlTdcHit hxle;
   DlTdcHit hxte;

   bool have0le = false;
   bool have1le = false;
   bool have0te = false;
   bool have1te = false;
   
   bool have4le = false;
   bool have5le = false;
   bool have4te = false;
   bool have5te = false;
   
   bool have8le = false;
   bool have9le = false;
   bool have8te = false;
   bool have9te = false;

   bool havexle = false;
   bool havexte = false;

   DlTdcHit chan1le;
   DlTdcHit chan1te;

   DlTdcHit chan2le;
   DlTdcHit chan2te;

   DlTdcHit chan3le;
   DlTdcHit chan3te;

   DlTdcHit chan4le;
   DlTdcHit chan4te;

   bool havechan1le = false;
   bool havechan1te = false;

   bool havechan2le = false;
   bool havechan2te = false;

   bool havechan3le = false;
   bool havechan3te = false;

   bool havechan4le = false;
   bool havechan4te = false;

public:
   void Reset()
   {
      first_time_sec = 0;
      last_time_sec = 0;

      min_time_sec = 0;
      max_time_sec = 0;

      have0le = false;
      have1le = false;
      have0te = false;
      have1te = false;
      
      have4le = false;
      have5le = false;
      have4te = false;
      have5te = false;
      
      have8le = false;
      have9le = false;
      have8te = false;
      have9te = false;
      
      havexle = false;
      havexte = false;

      havechan1le = false;
      havechan1te = false;

      havechan2le = false;
      havechan2te = false;

      havechan3le = false;
      havechan3te = false;

      havechan4le = false;
      havechan4te = false;
   }

   void AddHit(const DlTdcHit& h)
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
      
      if (h.ch == 0 && h.le && !have0le) {
         h0le = h;
         have0le = true;
      }
      
      if (h.ch == 0 && h.te && !have0te) {
         h0te = h;
         have0te = true;
      }
      
      if (h.ch == 1 && h.le && !have1le) {
         h1le = h;
         have1le = true;
      }
      
      if (h.ch == 1 && h.te && !have1te) {
         h1te = h;
         have1te = true;
      }
      
      if (h.ch == 4 && h.le && !have4le) {
         h4le = h;
         have4le = true;
      }
      
      if (h.ch == 4 && h.te && !have4te) {
         h4te = h;
         have4te = true;
      }
      
      if (h.ch == 5 && h.le && !have5le) {
         h5le = h;
         have5le = true;
      }
      
      if (h.ch == 5 && h.te && !have5te) {
         h5te = h;
         have5te = true;
      }
      
      if (h.ch == 6 && h.le && !have8le) {
         h8le = h;
         have8le = true;
      }
      
      if (h.ch == 6 && h.te && !have8te) {
         h8te = h;
         have8te = true;
      }
      
      if (h.ch == 7 && h.le && !have9le) {
         h9le = h;
         have9le = true;
      }
      
      if (h.ch == 7 && h.te && !have9te) {
         h9te = h;
         have9te = true;
      }
      
      if (h.ch == 2 && h.le && !havexle) {
         hxle = h;
         havexle = true;
      }
      
      if (h.ch == 2 && h.te && !havexte) {
         hxte = h;
         havexte = true;
      }

      if (h.ch == 8 && h.le && !havechan1le) {
         chan1le = h;
         havechan1le = true;
      }
      
      if (h.ch == 8 && h.te && !havechan1te) {
         chan1te = h;
         havechan1te = true;
      }

      if (h.ch == 9 && h.le && !havechan2le) {
         chan2le = h;
         havechan2le = true;
      }
      
      if (h.ch == 9 && h.te && !havechan2te) {
         chan2te = h;
         havechan2te = true;
      }

      if (h.ch == 10 && h.le && !havechan3le) {
         chan3le = h;
         havechan3le = true;
      }
      
      if (h.ch == 10 && h.te && !havechan3te) {
         chan3te = h;
         havechan3te = true;
      }

      if (h.ch == 11 && h.le && !havechan4le) {
         chan4le = h;
         havechan4le = true;
      }
      
      if (h.ch == 11 && h.te && !havechan4te) {
         chan4te = h;
         havechan4te = true;
      }
   }
};

class AdcEvent
{
public:
   double time_sec = 0;
   double dt = 0;
   AgAwHit a7;
   AgAwHit a15;

   AgAwHit a0;
   AgAwHit a1;
   AgAwHit a10;
   AgAwHit a11;
};

static double ns_to_mv(double x)
{
#if 0
   if (x < 25)
      return 250;

   return 211 + 0.704*x + 0.0309*x*x;
#endif

   if (x < 0)
      return 0;

   if (x > 300)
      return 9999;

   return 191 + 1.45*x + 0.0259*x*x;
}

#define MAX_TDC_CHAN 11

class DlTdcModule: public TARunObject
{
public:
   DlTdcFlags* fFlags = NULL;
   DlTdcUnpack* fU = NULL;
   
#ifdef HAVE_ROOT

   TCanvas* gWindow = NULL;
   TH1D* hphasele[2];
   TH1D* hphasete[2];
   TH1D* hfinele[2];
   TH1D* hfinete[2];
   TH1D* hdiff01le = NULL;
   TH1D* hdiff01te = NULL;
   TH1D* hdiff01le_fine = NULL;
   TH1D* hdiff01te_fine = NULL;
   TH1D* hwid0 = NULL;
   TH1D* hwid0_fine = NULL;
   int icd = 0;

   TH1D* fHphaseLe[MAX_TDC_CHAN+1];
   TH1D* fHphaseTe[MAX_TDC_CHAN+1];

   TH1D* hcalle[MAX_TDC_CHAN+1];
   TH1D* hcalte[MAX_TDC_CHAN+1];
   TH1D* hcalle_fine[MAX_TDC_CHAN+1];
   TH1D* hcalte_fine[MAX_TDC_CHAN+1];

   TH1D* fHhitdt1ns = NULL;
   TH1D* fHhitdt2ns = NULL;
   TH1D* fHhitdt3ns = NULL;
   TH1D* fHhitdt4ns = NULL;

   TH1D* fHeventdt1ns = NULL;
   TH1D* fHeventdt2ns = NULL;
   TH1D* fHeventdt3ns = NULL;
   TH1D* fHeventdt4ns = NULL;

   TCanvas* fDL1 = NULL;
   TH1D* fHt0;
   TH1D* fHt1;
   TH1D* fHt4;
   TH1D* fHt5;
   TH1D* fHw0;
   TH1D* fHw1;
   TH1D* fHw4;
   TH1D* fHw5;
   int   fDL1icd = 0;

   TH2D* fHw0w1 = NULL;
   TH2D* fHw0w4 = NULL;
   TH2D* fHw0w5 = NULL;
   TH2D* fHw1w4 = NULL;
   TH2D* fHw1w5 = NULL;
   TH2D* fHw4w5 = NULL;

   TH2D* fHtw0 = NULL;
   TH2D* fHtw1 = NULL;
   TH2D* fHtw4 = NULL;
   TH2D* fHtw5 = NULL;

   TH2D* fHt0t1 = NULL;
   TH2D* fHt0t4 = NULL;
   TH2D* fHt0t5 = NULL;
   TH2D* fHt1t4 = NULL;
   TH2D* fHt1t5 = NULL;
   TH2D* fHt4t5 = NULL;

   TH1D* fHt0m4 = NULL;
   TH1D* fHt1m5 = NULL;

   TH2D* fHt0m4xt1m5 = NULL;

   TH2D* fHt0m4w0 = NULL;
   TH2D* fHt0m4w4 = NULL;

   TH2D* fHt1m5w1 = NULL;
   TH2D* fHt1m5w5 = NULL;

   TH1D* fHt8mt9 = NULL;
   TH1D* fHw8 = NULL;
   TH1D* fHw9 = NULL;

   TH1D* fHa7 = NULL;
   TH1D* fHa15 = NULL;

   TH2D* fHw8a15 = NULL;
   TH2D* fHw9a7  = NULL;

   TH1D* fHa0  = NULL;
   TH1D* fHa1  = NULL;
   TH1D* fHa10 = NULL;
   TH1D* fHa11 = NULL;

   TH2D* fHw0a11 = NULL;
   TH2D* fHw1a10 = NULL;
   TH2D* fHw4a0 = NULL;
   TH2D* fHw5a1 = NULL;

   //TH2D* fX1 = NULL;
   //TH2D* fX2 = NULL;
   //TH2D* fX3 = NULL;
   //TH2D* fX4 = NULL;

   TH1D* fHt12ns = NULL;
   TH1D* fHt34ns = NULL;

   TH1D* fHt14ns = NULL;
   TH1D* fHt23ns = NULL;

   TH1D* fHw1ns = NULL;
   TH1D* fHw2ns = NULL;
   TH1D* fHw3ns = NULL;
   TH1D* fHw4ns = NULL;

   TH2D* fHw14ns = NULL;
   TH2D* fHw23ns = NULL;

   TH2D* fHt14w1ns = NULL;
   TH2D* fHt14w4ns = NULL;

   TH2D* fHt23w2ns = NULL;
   TH2D* fHt23w3ns = NULL;

   TH1D* fHa1mv = NULL;
   TH1D* fHa2mv = NULL;
   TH1D* fHa3mv = NULL;
   TH1D* fHa4mv = NULL;

   TH2D* fH_a1mv_t14ns = NULL;
   TH2D* fH_a4mv_t14ns = NULL;

   TH2D* fH_a2mv_t23ns = NULL;
   TH2D* fH_a3mv_t23ns = NULL;

   TH1D* fHt14ns_cut = NULL;
   TH1D* fHt23ns_cut = NULL;

#endif

   double fPrevEventTimeSec = 0;

   DlTdcEvent *fCt = NULL;

   std::deque<DlTdcEvent*> fTq;
   std::deque<AdcEvent*> fAq;

   bool fTrace = false;
   
   DlTdcModule(TARunInfo* runinfo, DlTdcFlags* flags)
      : TARunObject(runinfo)
   {
      if (fTrace)
         printf("DlTdcModule::ctor!\n");

      fModuleName = "dltdc_module";
      fFlags   = flags;
      fU = new DlTdcUnpack(35);
   }

   ~DlTdcModule()
   {
      if (fTrace)
         printf("DlTdcModule::dtor!\n");
      if (fU) {
         delete fU;
         fU = NULL;
      }
      for (auto t: fTq) {
         delete t;
      }
      for (auto a: fAq) {
         delete a;
      }
   }

   void BeginRun(TARunInfo* runinfo)
   {
      if (fTrace)
         printf("DlTdcModule::BeginRun, run %d, file %s\n", runinfo->fRunNo, runinfo->fFileName.c_str());

      if  (!fFlags->fCalib) {
         bool load_ok = fU->Load(runinfo->fRunNo);
         if (!load_ok) {
            printf("Cannot load TDC calibration for run %d\n", runinfo->fRunNo);
            exit(123);
         }
      }

#ifdef HAVE_ROOT
      runinfo->fRoot->fOutputFile->cd(); // select correct ROOT directory
      TDirectory* dir = gDirectory->mkdir("dltdc");
      dir->cd(); // select correct ROOT directory

      gWindow = new TCanvas("dltdc_calib", "dltdc_calib", 800, 800);
      gWindow->Clear();
      gWindow->Divide(2, 7);
      
      icd = 1;

      for (int i=0; i<2; i++) {
         char name[256];
         char title[256];
         sprintf(name,  "hphase%dle", i);
         sprintf(title, "hphase%dle", i);
         hphasele[i] = new TH1D(name, title, 101, -50, 50);
         gWindow->cd(icd++);
         hphasele[i]->Draw();
      }
      
      for (int i=0; i<2; i++) {
         char name[256];
         char title[256];
         sprintf(name,  "hphase%dte", i);
         sprintf(title, "hphase%dte", i);
         hphasete[i] = new TH1D(name, title, 101, -50, 50);
         gWindow->cd(icd++);
         hphasete[i]->Draw();
      }
      
      for (int i=0; i<2; i++) {
         char name[256];
         char title[256];
         sprintf(name,  "hfine%dle", i);
         sprintf(title, "hfine%dle, ns", i);
         hfinele[i] = new TH1D(name, title, 100, -15, 15);
         gWindow->cd(icd++);
         hfinele[i]->Draw();
      }
      
      for (int i=0; i<2; i++) {
         char name[256];
         char title[256];
         sprintf(name,  "hfine%dte", i);
         sprintf(title, "hfine%dte, ns", i);
         hfinete[i] = new TH1D(name, title, 100, -15, 15);
         gWindow->cd(icd++);
         hfinete[i]->Draw();
      }
      
      //TH1D* hper0le = new TH1D("hper0le", "period ch0 le, ns", 100, -50, 50);
      //gWindow->cd(9);
      //hper0le->Draw();
      
      //TH1D* hper0te = new TH1D("hper0te", "period ch0 te, ns", 100, -50, 50);
      //gWindow->cd(10);
      //hper0te->Draw();
      
      hdiff01le = new TH1D("hdiff01le", "diff 1-0 le, ns", 101, -25, 25);
      gWindow->cd(icd++);
      hdiff01le->Draw();
      
      hdiff01te = new TH1D("hdiff01te", "diff 1-0 te, ns", 101, -25, 25);
      gWindow->cd(icd++);
      hdiff01te->Draw();
      
      hdiff01le_fine = new TH1D("hdiff01le_fine", "diff 1-0 le fine, ns", 101, -5, 5);
      gWindow->cd(icd++);
      hdiff01le_fine->Draw();
      
      hdiff01te_fine = new TH1D("hdiff01te_fine", "diff 1-0 te fine, ns", 101, -5, 5);
      gWindow->cd(icd++);
      hdiff01te_fine->Draw();
      
      hwid0 = new TH1D("hwid0", "width ch0, ns", 101, 0, 1000);
      gWindow->cd(icd++);
      hwid0->Draw();
      
      hwid0_fine = new TH1D("hwid0_fine", "width ch0 fine, ns", 101, 545.0-5, 545.0+5);
      gWindow->cd(icd++);
      hwid0_fine->Draw();
      
      gWindow->Modified();
      gWindow->Update();

      for (int i=0; i<=MAX_TDC_CHAN; i++) {
         char name[256];
         char title[256];

         sprintf(name,  "tdc%02d_phase_le", i);
         sprintf(title, "tdc%02d_phase_le", i);
         fHphaseLe[i] = new TH1D(name, title, 101, -50, 50);

         sprintf(name,  "tdc%02d_phase_te", i);
         sprintf(title, "tdc%02d_phase_te", i);
         fHphaseTe[i] = new TH1D(name, title, 101, -50, 50);
      }

      for (int i=0; i<=MAX_TDC_CHAN; i++) {
         char name[256];
         char title[256];
         sprintf(name,  "hcalle%d", i+1);
         sprintf(title, "hcalle%d", i+1);
         hcalle[i] = new TH1D(name, title, 101, -25, 25);
         sprintf(name,  "hcalte%d", i+1);
         sprintf(title, "hcalte%d", i+1);
         hcalte[i] = new TH1D(name, title, 101, -25, 25);
         sprintf(name,  "hcalle%d_fine", i+1);
         sprintf(title, "hcalle%d_fine", i+1);
         hcalle_fine[i] = new TH1D(name, title, 101, -5, 5);
         sprintf(name,  "hcalte%d_fine", i+1);
         sprintf(title, "hcalte%d_fine", i+1);
         hcalte_fine[i] = new TH1D(name, title, 101, -5, 5);
      }

      fHhitdt1ns = new TH1D("hitdt1ns", "hit dt 100 ns", 100, 0, 100); // 100 ns
      fHhitdt2ns = new TH1D("hitdt2ns", "hit dt 1000 ns", 100, 100, 1000); // 1 usec
      fHhitdt3ns = new TH1D("hitdt3ns", "hit dt 1000 us", 100, 1000, 1000000); // 1 msec
      fHhitdt4ns = new TH1D("hitdt4ns", "hit dt 1000 ms", 100, 1000000, 1000000000); // 1 sec

      fHeventdt1ns = new TH1D("eventdt1ns", "event dt 100 ns", 100, 0, 100); // 100 ns
      fHeventdt2ns = new TH1D("eventdt2ns", "event dt 1000 ns", 100, 100, 1000); // 1 usec
      fHeventdt3ns = new TH1D("eventdt3ns", "event dt 1000 us", 100, 1000, 1000000); // 1 msec
      fHeventdt4ns = new TH1D("eventdt4ns", "event dt 1000 ms", 100, 1000000, 1000000000); // 1 sec

      if (fFlags->fHaveAdc) {
         fDL1 = new TCanvas("dl1", "dl1", 1600, 800);
         fDL1->Clear();
         fDL1->Divide(4, 2);
         
         fDL1icd = 1;
         
         fHt0 = new TH1D("t0", "t0", 100, 0, 20);
         fDL1->cd(fDL1icd++);
         fHt0->Draw();
         fHt1 = new TH1D("t1", "t1", 100, 0, 20);
         fDL1->cd(fDL1icd++);
         fHt1->Draw();
         fHt4 = new TH1D("t4", "t4", 100, 0, 20);
         fDL1->cd(fDL1icd++);
         fHt4->Draw();
         fHt5 = new TH1D("t5", "t5", 100, 0, 20);
         fDL1->cd(fDL1icd++);
         fHt5->Draw();
         
         fHw0 = new TH1D("w0", "w0", 100, 0, 200);
         fDL1->cd(fDL1icd++);
         fHw0->Draw();
         fHw1 = new TH1D("w1", "w1", 100, 0, 200);
         fDL1->cd(fDL1icd++);
         fHw1->Draw();
         fHw4 = new TH1D("w4", "w4", 100, 0, 200);
         fDL1->cd(fDL1icd++);
         fHw4->Draw();
         fHw5 = new TH1D("w5", "w5", 100, 0, 200);
         fDL1->cd(fDL1icd++);
         fHw5->Draw();
         
         fDL1->Modified();
         fDL1->Update();
         
         fHw0w1 = new TH2D("w0w1", "w0w1", 50, 0, 200, 50, 0, 200);
         fHw0w4 = new TH2D("w0w4", "w0w4", 50, 0, 200, 50, 0, 200);
         fHw0w5 = new TH2D("w0w5", "w0w5", 50, 0, 200, 50, 0, 200);
         fHw1w4 = new TH2D("w1w4", "w1w4", 50, 0, 200, 50, 0, 200);
         fHw1w5 = new TH2D("w1w5", "w1w5", 50, 0, 200, 50, 0, 200);
         fHw4w5 = new TH2D("w4w5", "w4w5", 50, 0, 200, 50, 0, 200);
         
         fHtw0 = new TH2D("tw0", "tw0", 50, 0, 20, 50, 0, 200);
         fHtw1 = new TH2D("tw1", "tw1", 50, 0, 20, 50, 0, 200);
         fHtw4 = new TH2D("tw4", "tw4", 50, 0, 20, 50, 0, 200);
         fHtw5 = new TH2D("tw5", "tw5", 50, 0, 20, 50, 0, 200);
         
         fHt0t1 = new TH2D("t0t1", "t0t1", 50, 0, 20, 50, 0, 20);
         fHt0t4 = new TH2D("t0t4", "t0t4", 50, 0, 20, 50, 0, 20);
         fHt0t5 = new TH2D("t0t5", "t0t5", 50, 0, 20, 50, 0, 20);
         fHt1t4 = new TH2D("t1t4", "t1t4", 50, 0, 20, 50, 0, 20);
         fHt1t5 = new TH2D("t1t5", "t1t5", 50, 0, 20, 50, 0, 20);
         fHt4t5 = new TH2D("t4t5", "t4t5", 50, 0, 20, 50, 0, 20);
         
         fHt0m4 = new TH1D("t0m4", "t0m4", 100, -20, 20);
         fHt1m5 = new TH1D("t1m5", "t1m5", 100, -20, 20);
         
         fHt0m4xt1m5 = new TH2D("t0m4xt1m5", "t0m4xt1m5", 100, -20, 20, 100, -20, 20);
         
         fHt0m4w0 = new TH2D("t0m4w0", "t0m4w0", 100, -20, 20, 50, 0, 200);
         fHt0m4w4 = new TH2D("t0m4w4", "t0m4w4", 100, -20, 20, 50, 0, 200);
         
         fHt1m5w1 = new TH2D("t1m5w1", "t1m5w1", 100, -20, 20, 50, 0, 200);
         fHt1m5w5 = new TH2D("t1m5w5", "t1m5w5", 100, -20, 20, 50, 0, 200);
         
         fHt8mt9 = new TH1D("t8mt9", "t8mt9", 100, -20, 20);
         fHw8 = new TH1D("w8", "w8", 100, 0, 50);
         fHw9 = new TH1D("w9", "w9", 100, 0, 50);
         
         fHa7  = new TH1D("a7",  "a7",  100, 0, 20000);
         fHa15 = new TH1D("a15", "a15", 100, 0, 20000);
         
         fHw8a15 = new TH2D("w8a15", "w8a15", 50, 0, 50, 50, 0, 20000);
         fHw9a7  = new TH2D("w9a7",  "w9a7",  50, 0, 50, 50, 0, 20000);
         
         fHa0  = new TH1D("a0",  "a0",  100, 0, 20000);
         fHa1  = new TH1D("a1",  "a1",  100, 0, 20000);
         fHa10 = new TH1D("a10", "a10", 100, 0, 20000);
         fHa11 = new TH1D("a11", "a11", 100, 0, 20000);
         
         fHw0a11 = new TH2D("w0a11", "w0a11", 50, 0, 200, 50, 0, 20000);
         fHw1a10 = new TH2D("w1a10", "w1a10", 50, 0, 200, 50, 0, 20000);
         fHw4a0  = new TH2D("w4a0",  "w4a0",  50, 0, 200, 50, 0, 20000);
         fHw5a1  = new TH2D("w5a1",  "w5a1",  50, 0, 200, 50, 0, 20000);
         
         //fX1 = new TH2D("x1", "x1", 50, 0, 200, 50, 0, 20000);
         //fX2 = new TH2D("x2", "x2", 50, 0, 200, 50, 0, 20000);
         //fX3 = new TH2D("x3", "x3", 50, 0, 200, 50, 0, 20000);
         //fX4 = new TH2D("x4", "x4", 50, 0, 200, 50, 0, 20000);
      }
         
      fHt12ns = new TH1D("t12ns", "sipm board 1, t2-t1, ns", 500, -10, 10);
      fHt34ns = new TH1D("t34ns", "sipm board 2, t4-t3, ns", 500, -10, 10);

      fHt14ns = new TH1D("t14ns", "paddle 1, t4-t1, ns", 500, -10, 10);
      fHt23ns = new TH1D("t23ns", "paddle 2, t3-t2, ns", 500, -10, 10);

      fHw1ns = new TH1D("w1ns", "w1ns", 100, 0, 400);
      fHw2ns = new TH1D("w2ns", "w2ns", 100, 0, 400);
      fHw3ns = new TH1D("w3ns", "w3ns", 100, 0, 400);
      fHw4ns = new TH1D("w4ns", "w4ns", 100, 0, 400);

      fHw14ns = new TH2D("w14ns", "w4ns vs w1ns", 100, 0, 400, 100, 0, 400);
      fHw23ns = new TH2D("w23ns", "w3ns vs w2ns", 100, 0, 400, 100, 0, 400);

      fHt14w1ns = new TH2D("t14w1ns", "w1ns vs t14ns", 100, -5, 5, 100, 0, 400);
      fHt14w4ns = new TH2D("t14w4ns", "w4ns vs t14ns", 100, -5, 5, 100, 0, 400);

      fHt23w2ns = new TH2D("t23w2ns", "w2ns vs t23ns", 100, -5, 5, 100, 0, 400);
      fHt23w3ns = new TH2D("t23w3ns", "w3ns vs t23ns", 100, -5, 5, 100, 0, 400);

      fHa1mv = new TH1D("a1mv", "calculated amp 1, mV", 100, 0, 2000);
      fHa2mv = new TH1D("a2mv", "calculated amp 2, mV", 100, 0, 2000);
      fHa3mv = new TH1D("a3mv", "calculated amp 3, mV", 100, 0, 2000);
      fHa4mv = new TH1D("a4mv", "calculated amp 4, mV", 100, 0, 2000);

      fH_a1mv_t14ns = new TH2D("a1mv_t14ns", "t4-t1 (ns) vs a1 (mV)", 100, 0, 2000, 100, -5, 5);
      fH_a4mv_t14ns = new TH2D("a4mv_t14ns", "t4-t1 (ns) vs a4 (mV)", 100, 0, 2000, 100, -5, 5);

      fH_a2mv_t23ns = new TH2D("a2mv_t23ns", "t3-t2 (ns) vs a2 (mV)", 100, 0, 2000, 100, -5, 5);
      fH_a3mv_t23ns = new TH2D("a3mv_t23ns", "t3-t2 (ns) vs a3 (mV)", 100, 0, 2000, 100, -5, 5);

      fHt14ns_cut = new TH1D("t14ns_cut", "t14ns_cut", 500, -10, 10);
      fHt23ns_cut = new TH1D("t23ns_cut", "t23ns_cut", 500, -10, 10);

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

      if (fFlags->fCalib) {
         printf("DlTdcModule::EndRun: Saving TDC calibrations for run %d\n", runinfo->fRunNo);
         fU->Save(runinfo->fRunNo);
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
      double dt = t.min_time_sec - prev_event_time_sec;
      
      printf("dlsc %d %d %d %d, %.9f %.9f sec, dt %.9f sec\n", t.havechan1le, t.havechan2le, t.havechan3le, t.havechan4le, prev_event_time_sec, t.min_time_sec, dt);
      
      if (t.havechan1le && t.havechan2le) {
         double t12_ns = sec_to_ns(t.chan2le.time_sec - t.chan1le.time_sec);
         fHt12ns->Fill(t12_ns);
      }
      
      if (t.havechan3le && t.havechan4le) {
         double t34_ns = sec_to_ns(t.chan4le.time_sec - t.chan3le.time_sec);
         fHt34ns->Fill(t34_ns);
      }
      
      if (t.havechan1le && t.havechan4le) {
         double t14_ns = sec_to_ns(t.chan4le.time_sec - t.chan1le.time_sec);
         double w1_ns = sec_to_ns(t.chan1te.time_sec - t.chan1le.time_sec);
         double w4_ns = sec_to_ns(t.chan4te.time_sec - t.chan4le.time_sec);

         double a1_mv = ns_to_mv(w1_ns);
         double a4_mv = ns_to_mv(w4_ns);
         
         printf("new dlsc event 1*4, le %.9f %.9f sec, diff14 %.0f ns, w1 %.0f, w4 %.0f ns, a1 %.0f, a4 %.0f mV!\n", t.chan1le.time_sec, t.chan4le.time_sec, t14_ns, w1_ns, w4_ns, a1_mv, a4_mv);
         
         fHt14ns->Fill(t14_ns);
         fHw1ns->Fill(w1_ns);
         fHw4ns->Fill(w4_ns);
         fHw14ns->Fill(w1_ns, w4_ns);
         fHt14w1ns->Fill(t14_ns, w1_ns);
         fHt14w4ns->Fill(t14_ns, w4_ns);

         fHa1mv->Fill(a1_mv);
         fHa4mv->Fill(a4_mv);
         fH_a1mv_t14ns->Fill(a1_mv, t14_ns);
         fH_a4mv_t14ns->Fill(a4_mv, t14_ns);
         
         //if (w1_ns > 150 && w4_ns > 150) {
         //   fHt14ns_cut->Fill(t14_ns);
         //}
         
         if (w1_ns > 160 && w1_ns < 170) {
            if (w4_ns > 175 && w4_ns < 185) {
               fHt14ns_cut->Fill(t14_ns);
            }
         }
      }
      
      if (t.havechan2le && t.havechan3le) {
         double t23_ns = sec_to_ns(t.chan3le.time_sec - t.chan2le.time_sec);
         double w2_ns = sec_to_ns(t.chan2te.time_sec - t.chan2le.time_sec);
         double w3_ns = sec_to_ns(t.chan3te.time_sec - t.chan3le.time_sec);

         double a2_mv = ns_to_mv(w2_ns);
         double a3_mv = ns_to_mv(w3_ns);
         
         printf("new dlsc event 2*3, le %.9f %.9f sec, diff23 %.0f ns, w2 %.0f, w3 %.0f ns!\n", t.chan2le.time_sec, t.chan3le.time_sec, t23_ns, w2_ns, w3_ns);
         
         fHt23ns->Fill(t23_ns);
         fHw2ns->Fill(w2_ns);
         fHw3ns->Fill(w3_ns);
         fHw23ns->Fill(w2_ns, w3_ns);
         fHt23w2ns->Fill(t23_ns, w2_ns);
         fHt23w3ns->Fill(t23_ns, w3_ns);
         
         fHa2mv->Fill(a2_mv);
         fHa3mv->Fill(a3_mv);
         fH_a2mv_t23ns->Fill(a2_mv, t23_ns);
         fH_a3mv_t23ns->Fill(a3_mv, t23_ns);
         
         if (w2_ns > 150 && w3_ns > 150) {
            fHt23ns_cut->Fill(t23_ns);
         }
      }
   }

   double fPrevTdcTime = 0;

   void FinishEventTA(DlTdcEvent& t, AdcEvent& a)
   {
      double tdc_time = t.hxle.time_sec;

      double dt = tdc_time - fPrevTdcTime;

      double w = 0;

      if (t.havexle && t.havexte) {
         w = sec_to_ns(t.hxte.time_sec - t.hxle.time_sec);
      } else {
         return;
      }

      printf("time %.6f, dt %.6f, have %d %d - %d %d %d %d %d %d %d %d - %d %d %d %d - width %.0f\n", tdc_time, dt, t.havexle, t.havexte, t.have0le, t.have0te, t.have1le, t.have1te, t.have4le, t.have4te, t.have5le, t.have5te, t.have8le, t.have8te, t.have9le, t.have9te, w);

      if (w > 500) {
         fPrevTdcTime = tdc_time;
      }

      if (t.have0le && t.have0te) {
         double wid_ns = (t.h0te.time_sec - t.h0le.time_sec)*1e9;
         printf("xwidth 0: %.3f %.3f ns\n", (t.h0te.coarse_sec - t.h0le.coarse_sec)*1e9, wid_ns);
         //printf("WWW: "); h0le.Print(); printf("\n");
         //printf("WWW: "); h0te.Print(); printf("\n");
#ifdef HAVE_ROOT
         hwid0->Fill(wid_ns);
         hwid0_fine->Fill(wid_ns);
#endif
      }

      if (t.have0le && t.have0te) {
         double le0 = t.h0le.time_sec;
         double te0 = t.h0te.time_sec;

         int i=0;
         hcalle[i]->Fill(sec_to_ns(t.h1le.time_sec - le0));
         hcalte[i]->Fill(sec_to_ns(t.h1te.time_sec - te0));
         hcalle_fine[i]->Fill(sec_to_ns(t.h1le.time_sec - le0));
         hcalte_fine[i]->Fill(sec_to_ns(t.h1te.time_sec - te0));

         i++;
         hcalle[i]->Fill(sec_to_ns(t.hxle.time_sec - le0));
         hcalte[i]->Fill(sec_to_ns(t.hxte.time_sec - te0));
         hcalle_fine[i]->Fill(sec_to_ns(t.hxle.time_sec - le0));
         hcalte_fine[i]->Fill(sec_to_ns(t.hxte.time_sec - te0));

         i++;
         hcalle[i]->Fill(sec_to_ns(t.h4le.time_sec - le0));
         hcalte[i]->Fill(sec_to_ns(t.h4te.time_sec - te0));
         hcalle_fine[i]->Fill(sec_to_ns(t.h4le.time_sec - le0));
         hcalte_fine[i]->Fill(sec_to_ns(t.h4te.time_sec - te0));

         i++;
         hcalle[i]->Fill(sec_to_ns(t.h5le.time_sec - le0));
         hcalte[i]->Fill(sec_to_ns(t.h5te.time_sec - te0));
         hcalle_fine[i]->Fill(sec_to_ns(t.h5le.time_sec - le0));
         hcalte_fine[i]->Fill(sec_to_ns(t.h5te.time_sec - te0));

         i++;
         hcalle[i]->Fill(sec_to_ns(t.h8le.time_sec - le0));
         hcalte[i]->Fill(sec_to_ns(t.h8te.time_sec - te0));
         hcalle_fine[i]->Fill(sec_to_ns(t.h8le.time_sec - le0));
         hcalte_fine[i]->Fill(sec_to_ns(t.h8te.time_sec - te0));

         i++;
         hcalle[i]->Fill(sec_to_ns(t.h9le.time_sec - le0));
         hcalte[i]->Fill(sec_to_ns(t.h9te.time_sec - te0));
         hcalle_fine[i]->Fill(sec_to_ns(t.h9le.time_sec - le0));
         hcalte_fine[i]->Fill(sec_to_ns(t.h9te.time_sec - te0));
      }

      if (t.have0le && t.have0te && t.have1le && t.have1te) {
         double diff01le_ns = (t.h1le.time_sec - t.h0le.time_sec)*1e9;
         double diff01te_ns = (t.h1te.time_sec - t.h0te.time_sec)*1e9;
         printf("xdiff 01 le %.1f, te %.1f ns\n", diff01le_ns, diff01te_ns);
            
#ifdef HAVE_ROOT
         hdiff01le->Fill(diff01le_ns);
         hdiff01te->Fill(diff01te_ns);
         hdiff01le_fine->Fill(diff01le_ns);
         hdiff01te_fine->Fill(diff01te_ns);
#endif
      }

      if (t.have8le && t.have8te && t.have9le && t.have9te) {
         double t8mt9 = sec_to_ns(t.h8le.time_sec - t.h9le.time_sec);
         double w8 = sec_to_ns(t.h8te.time_sec - t.h8le.time_sec); 
         double w9 = sec_to_ns(t.h9te.time_sec - t.h9le.time_sec); 

         printf("t8mt9 %7.3f, w8 %7.3f w9 %7.3f, a7 %5.0f, a15 %5.0f\n", t8mt9, w8, w9, a.a7.amp, a.a15.amp);

         fHt8mt9->Fill(t8mt9);
         fHw8->Fill(w8);
         fHw9->Fill(w9);
         
         fHa7->Fill(a.a7.amp);
         fHa15->Fill(a.a15.amp);

         if (a.a15.amp > 0)
            fHw8a15->Fill(w8, a.a15.amp);
         if (a.a7.amp > 0)
            fHw9a7->Fill(w9, a.a7.amp);
      }

      if (t.have0le && t.have0te && t.have1le && t.have1te && t.have4le && t.have4te && t.have5le && t.have5te) {
         // DL analysis

         double tt = t.h0le.time_sec;
         tt = amin(tt, t.h1le.time_sec);
         tt = amin(tt, t.h4le.time_sec);
         tt = amin(tt, t.h5le.time_sec);

         double t0 = sec_to_ns(t.h0le.time_sec - tt);
         double t1 = sec_to_ns(t.h1le.time_sec - tt);
         double t4 = sec_to_ns(t.h4le.time_sec - tt);
         double t5 = sec_to_ns(t.h5le.time_sec - tt);

         double w0 = sec_to_ns(t.h0te.time_sec - t.h0le.time_sec);
         double w1 = sec_to_ns(t.h1te.time_sec - t.h1le.time_sec);
         double w4 = sec_to_ns(t.h4te.time_sec - t.h4le.time_sec);
         double w5 = sec_to_ns(t.h5te.time_sec - t.h5le.time_sec);

         printf("xtime %.6f %.6f %.6f %.6f tt %.6f, hit time %.3f %.3f %.3f %.3f, width %.3f %.3f %.3f %.3f\n", t.h0le.time_sec, t.h1le.time_sec, t.h4le.time_sec, t.h5le.time_sec, tt, t0, t1, t4, t5, w0, w1, w4, w5);

#ifdef HAVE_ROOT
         fHt0->Fill(t0);
         fHt1->Fill(t1);
         fHt4->Fill(t4);
         fHt5->Fill(t5);

         fHw0->Fill(w0);
         fHw1->Fill(w1);
         fHw4->Fill(w4);
         fHw5->Fill(w5);

         fHw0w1->Fill(w0, w1);
         fHw0w4->Fill(w0, w4);
         fHw0w5->Fill(w0, w5);
         fHw1w4->Fill(w1, w4);
         fHw1w5->Fill(w1, w5);
         fHw4w5->Fill(w4, w5);

         fHtw0->Fill(t0, w0);
         fHtw1->Fill(t1, w1);
         fHtw4->Fill(t4, w4);
         fHtw5->Fill(t5, w5);

         fHt0t1->Fill(t0, t1);
         fHt0t4->Fill(t0, t4);
         fHt0t5->Fill(t0, t5);
         fHt1t4->Fill(t1, t4);
         fHt1t5->Fill(t1, t5);
         fHt4t5->Fill(t4, t5);

         fHt0m4->Fill(t0-t4);
         fHt1m5->Fill(t1-t5);

         fHt0m4xt1m5->Fill(t0-t4, t1-t5);

         fHt0m4w0->Fill(t0-t4, w0);
         fHt0m4w4->Fill(t0-t4, w4);

         fHt1m5w1->Fill(t1-t5, w1);
         fHt1m5w5->Fill(t1-t5, w5);

         fHa0->Fill(a.a0.amp);
         fHa1->Fill(a.a1.amp);
         fHa10->Fill(a.a10.amp);
         fHa11->Fill(a.a11.amp);

         if (a.a11.amp > 0)
            fHw0a11->Fill(w0, a.a11.amp);
         if (a.a10.amp > 0)
            fHw1a10->Fill(w1, a.a10.amp);
         if (a.a0.amp > 0)
            fHw4a0->Fill(w4, a.a0.amp);
         if (a.a1.amp > 0)
            fHw5a1->Fill(w5, a.a1.amp);

         //fX1->Fill(w4, a.a0.amp);
         //fX2->Fill(w4, a.a1.amp);
         //fX3->Fill(w4, a.a10.amp);
         //fX4->Fill(w4, a.a11.amp);
#endif
      }
   }

   void A()
   {
      printf("A %zu %zu\n", fTq.size(), fAq.size());

      size_t n = fTq.size();
      if (n > fAq.size())
         n = fAq.size();

      for (size_t i=0; i<n; i++) {
         printf("slot %zu: time %.6f %.6f dt %.6f %.6f, ddt %.6f\n", i, fTq[i]->hxte.time_sec, fAq[i]->time_sec, fTq[i]->dt, fAq[i]->dt, fTq[i]->dt - fAq[i]->dt);
      }

      if (n < 5)
         return;

      double eps = 0.000050;

      if (fabs(fTq[0]->dt - fAq[0]->dt) < eps) {
         // good we are sycnhed
         printf("SYNC OK!\n");
         FinishEventTA(*fTq[0], *fAq[0]);
         fTq.pop_front();
         fAq.pop_front();
      } else if (fabs(fTq[1]->dt - fAq[0]->dt) < eps) {
         // off by one
         fTq.pop_front();
         printf("SYNC POP T!\n");
      } else if (fabs(fTq[2]->dt - fAq[0]->dt) < eps) {
         // off by one
         fTq.pop_front();
         fTq.pop_front();
         printf("SYNC POP 2 T!\n");
      } else if (fabs(fTq[0]->dt - fAq[1]->dt) < eps) {
         // off by one
         fAq.pop_front();
         printf("SYNC POP A!\n");
      } else if (fabs(fTq[0]->dt - fAq[2]->dt) < eps) {
         // off by one
         fAq.pop_front();
         fAq.pop_front();
         printf("SYNC POP 2 A!\n");
      } else if (fabs(fTq[1]->dt - fAq[1]->dt) < eps) {
         fTq.pop_front();
         fAq.pop_front();
         printf("SYNC POP T and A!\n");
      } else {
         fTq.pop_front();
         fAq.pop_front();
         printf("SYNC BAD POP both!\n");
      }

#if 0
      n = fTq.size();
      if (n > fAq.size())
         n = fAq.size();

      for (size_t i=0; i<n; i++) {
         printf("slot %zu: time %.6f %.6f dt %.6f %.6f\n", i, fTq[i]->hxte.time_sec, fAq[i]->time_sec, fTq[i]->dt, fAq[i]->dt);
      }
#endif
   }

   double fLastXtime = 0;

   bool fFirstTrig = true;

   TAFlowEvent* Analyze(TARunInfo* runinfo, TMEvent* event, TAFlags* flags, TAFlowEvent* flow)
   {
      //printf("DlTdcModule::Analyze, run %d, event serno %d, id 0x%04x, data size %d\n", runinfo->fRunNo, event->serial_number, (int)event->event_id, event->data_size);

      if (event->event_id != 4)
         return flow;

      if (!fFlags->fEnabled)
         return flow;

      TMBank* tdcbank = event->FindBank("CBT2");

      bool calib = fFlags->fCalib;

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
            
            h.Print(); printf("\n");

            if (calib) {
               fU->fCalib[h.ch].AddHit(h);
            }

            double hit_dt_ns = sec_to_ns(h.time_sec - fCt->max_time_sec);

            fHhitdt1ns->Fill(hit_dt_ns);
            fHhitdt2ns->Fill(hit_dt_ns);
            fHhitdt3ns->Fill(hit_dt_ns);
            fHhitdt4ns->Fill(hit_dt_ns);

            //printf("TTX %.9f %.9f sec, first %.9f, last %.9f sec\n", fCt->min_time_sec, fCt->max_time_sec, fCt->first_time_sec, fCt->last_time_sec);

            if (hit_dt_ns > 300.0) {
               double event_dt_ns = sec_to_ns(h.time_sec - fCt->min_time_sec);

               //printf("TTT %.9f -> %.9f sec, dt %.3f ns\n", fCt->min_time_sec, h.time_sec, event_dt_ns);

               fHeventdt1ns->Fill(event_dt_ns);
               fHeventdt2ns->Fill(event_dt_ns);
               fHeventdt3ns->Fill(event_dt_ns);
               fHeventdt4ns->Fill(event_dt_ns);

               if (fFlags->fHaveAdc) {
                  if (fCt->havexle && fCt->havexte) {
                     if (fFirstTrig) {
                        fU->Reset();
                        fFirstTrig = false;
                     } else {
                        fCt->time_sec = fCt->hxle.time_sec;
                        fCt->dt = fCt->hxle.time_sec - fLastXtime;
                        fLastXtime = fCt->hxle.time_sec;
                        fTq.push_back(fCt);
                        fCt = new DlTdcEvent;
                     }
                  } else {
                     fCt->Reset();
                  }
               } else {
                  FinishEventT(fPrevEventTimeSec, *fCt);
               }

               fPrevEventTimeSec = fCt->min_time_sec;

               fCt->Reset();
            }

#ifdef HAVE_ROOT
            if (h.ch == 0 || h.ch == 1) {
               if (h.le) {
                  hphasele[h.ch]->Fill(h.phase);
                  //if (h.phase > 0)
                  hfinele[h.ch]->Fill(h.fine_ns);
                  //else
                  // hfinele[h.ch]->Fill(-h.fine_ns);
               }
               
               if (h.te) {
                  hphasete[h.ch]->Fill(h.phase);
                  //if (h.phase > 0)
                  hfinete[h.ch]->Fill(h.fine_ns);
                  //else
                  //hfinete[h.ch]->Fill(-h.fine_ns);
               }
            }

            if (h.ch >= 0 && h.ch <= MAX_TDC_CHAN) {
               if (h.le) {
                  fHphaseLe[h.ch]->Fill(h.phase);
               }
               if (h.te) {
                  fHphaseTe[h.ch]->Fill(h.phase);
               }
            }
#endif

            fCt->AddHit(h);

         } // loop over data

         time_t now = time(NULL);
         static time_t last = 0;
         
         if (now - last > 5) {

#ifdef HAVE_ROOT
            for (int i=1; i<=icd; i++)
               gWindow->cd(i)->Modified();
            
            gWindow->Modified();
            gWindow->Update();
            gWindow->SaveAs("cal1a.root");
            gWindow->SaveAs("cal1a.pdf");

            if (fDL1) {
               for (int i=1; i<=fDL1icd; i++)
                  fDL1->cd(i)->Modified();
               
               fDL1->Modified();
               fDL1->Update();
               fDL1->SaveAs("dl1a.root");
               fDL1->SaveAs("dl1a.pdf");
            }
#endif

            if (calib) {
               for (auto& c: fU->fCalib) {
                  c.Update();
                  c.Print();
               }
            }

            last = time(NULL);
         }
         
      }

      if (fFlags->fHaveAdc) {
         A();
      }
      
      return flow;
   }

   double fPrevFlowTime = 0;

   TAFlowEvent* AnalyzeFlowEvent(TARunInfo* runinfo, TAFlags* flags, TAFlowEvent* flow)
   {
      if (!fFlags->fEnabled)
         return flow;
      
      //printf("DltdcModule::AnalyzeFlowEvent!\n");

      AgEventFlow *ef = flow->Find<AgEventFlow>();

      if (!ef || !ef->fEvent)
         return flow;

      AgAwHitsFlow* eawh = flow->Find<AgAwHitsFlow>();
      //AgPadHitsFlow* eph = flow->Find<AgPadHitsFlow>();
      //AgBscAdcHitsFlow* eba = flow->Find<AgBscAdcHitsFlow>();

      if (eawh) {
         AdcEvent* a = new AdcEvent;
         a->time_sec = ef->fEvent->time;
         a->a0.amp = 0;
         a->a1.amp = 0;
         a->a10.amp = 0;
         a->a11.amp = 0;
         a->a7.amp = 0;
         a->a15.amp = 0;

         for (const auto& ah: eawh->fAwHits) {
            if (ah.adc_module != 15)
               continue;
            //printf("ADC hit %d %d %f %f!\n", ah.adc_module, ah.adc_chan, ah.time, ah.amp);
            if (ah.adc_chan == 0)
               a->a0 = ah;
            if (ah.adc_chan == 1)
               a->a1 = ah;
            if (ah.adc_chan == 10)
               a->a10 = ah;
            if (ah.adc_chan == 11)
               a->a11 = ah;
            if (ah.adc_chan == 7)
               a->a7 = ah;
            if (ah.adc_chan == 15)
               a->a15 = ah;
         }

         a->dt = a->time_sec - fPrevFlowTime;

         printf("time %.6f, dt %.6f, adc %f %f ---------\n", a->time_sec, a->dt, a->a7.amp, a->a15.amp);

         fPrevFlowTime = a->time_sec;

         fAq.push_back(a);

         A();
      }

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
      printf("--dltdc -- enable dltdc code\n");
      printf("--dltdc-calib -- calibrate dltdc");
      printf("--dltdc-adc -- have ADC data");
   }

   void Init(const std::vector<std::string> &args)
   {
      printf("DlTdcModuleFactory::Init!\n");

      for (unsigned i=0; i<args.size(); i++) {
         if (args[i] == "--dltdc") {
            fFlags.fEnabled = true;
         }
         if (args[i] == "--dltdc-calib") {
            fFlags.fCalib = true;
         }
         if (args[i] == "--dltdc-adc") {
            fFlags.fHaveAdc = true;
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
