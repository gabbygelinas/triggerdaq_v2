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
#include <list>

#include "TStyle.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"

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

class DlTdcFlags
{
public:
   bool fEnabled = true;
   bool fTriggered = false;
   bool fCalibFineTime = false;
   bool fCalibOffsets  = false;
   bool fDebug = false;
   bool fPrintRawHits    = false;
   bool fPrintSortedHits = false;
   bool fPrintEventHits  = false;
   bool fAdd40 = false;
   bool fSub40 = false;
};

#define NUM_TDC_CHAN (32+3)
#define MAX_TDC_CHAN (NUM_TDC_CHAN-1)

class DlTdcConfig
{
public:
   int fChanA = 0;
   int fChanB = 1;
   int fChanT = 10;

public:
   DlTdcConfig(int num);
   bool ReadJson(int runno);
};

DlTdcConfig::DlTdcConfig(int num)
{
}

bool DlTdcConfig::ReadJson(int runno)
{
   if (runno >= 17 && runno < 900000) {
      printf("NEW MAP!\n");
      fChanA = 32;
      fChanB = 33;
      fChanT = 34;
   }

   return true;
};

class DlTdcModule: public TARunObject
{
public:
   DlTdcFlags* fFlags = NULL;
   DlTdcUnpack* fU = NULL;
   
#ifdef HAVE_ROOT
   TH2D* fHphaseLePerTdcChan = NULL;
   TH2D* fHphaseTePerTdcChan = NULL;

   TH2D* fHfineLePerTdcChan = NULL;
   TH2D* fHfineTePerTdcChan = NULL;

   TH1D* fHphaseLe[MAX_TDC_CHAN+1];
   TH1D* fHphaseTe[MAX_TDC_CHAN+1];

   std::vector<TProfile*> fHfineBinsLe;
   std::vector<TProfile*> fHfineBinsTe;

   TH1D* fHfineLe[MAX_TDC_CHAN+1];
   TH1D* fHfineTe[MAX_TDC_CHAN+1];

   TH1D* fHfineLeP[MAX_TDC_CHAN+1];
   TH1D* fHfineLeN[MAX_TDC_CHAN+1];
   TH1D* fHfineTeP[MAX_TDC_CHAN+1];
   TH1D* fHfineTeN[MAX_TDC_CHAN+1];

   TH2D* fHfineLeP0 = NULL;
   TH2D* fHfineTeN0 = NULL;
   TH2D* fHfineLeP1 = NULL;
   TH2D* fHfineTeN1 = NULL;

   TH1D* fHwidth[MAX_TDC_CHAN+1];

   TH2D* fHpulserLePerTdcChan = NULL;
   TH2D* fHpulserTePerTdcChan = NULL;
   TH2D* fHpulserWiPerTdcChan = NULL;
   TH1D* fHpulserLeAll = NULL;
   TH1D* fHpulserTeAll = NULL;
   TH1D* fHpulserWiAll = NULL;
   TH1D* fHpulserLe[MAX_TDC_CHAN+1];
   TH1D* fHpulserTe[MAX_TDC_CHAN+1];
   TH1D* fHpulserWi[MAX_TDC_CHAN+1];
   TH2D* fHpulserPP[MAX_TDC_CHAN+1];
   TH2D* fHpulserFF[MAX_TDC_CHAN+1];

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

   //std::vector<TH2D*> fHdt_PP;
   std::vector<TH2D*> fHdt_FF;

   TH1D* fHunphysical_map  = NULL;
   TH2D* fHunphysical_map2 = NULL;

   TH2D* fHunphysical_ns_01_14_w01 = NULL;
   TH2D* fHunphysical_ns_01_14_w14 = NULL;

   TH2D* fHunphysical_ns_01_14_c01 = NULL;
   TH2D* fHunphysical_ns_01_14_c14 = NULL;

   std::vector<TH1D*> fHunphysical_ns;
   //std::vector<TH2D*> fHunphysical_PP;
   std::vector<TH2D*> fHunphysical_FF;
#endif

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

      if (fCt) {
         delete fCt;
         fCt = NULL;
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
         std::string offset_json = fCfm->GetFilename("dltdc", "offsets", runinfo->fRunNo, "json");
         
         bool load_ok = fU->LoadOffsets(offset_json.c_str());
         
         printf("json file %s load_ok %d\n", offset_json.c_str(), load_ok);
         
         if (!load_ok) {
            printf("Cannot load TDC offset calibration for run %d\n", runinfo->fRunNo);
            exit(123);
         }
      }
      
#if 0
      if (fFlags->fAdd40) {
         // shift B-cable signals
         for (int i=16; i<32; i++) {
            fU->fCalib[i].lepos.fOffsetNs += 40.0;
            fU->fCalib[i].leneg.fOffsetNs += 40.0;
            fU->fCalib[i].tepos.fOffsetNs += 40.0;
            fU->fCalib[i].teneg.fOffsetNs += 40.0;
         }
         
         // shift B signal
         fU->fCalib[33].lepos.fOffsetNs += 40.0;
         fU->fCalib[33].leneg.fOffsetNs += 40.0;
         fU->fCalib[33].tepos.fOffsetNs += 40.0;
         fU->fCalib[33].teneg.fOffsetNs += 40.0;
      }
      
      if (fFlags->fSub40) {
         // shift B-cable signals
         for (int i=16; i<32; i++) {
            fU->fCalib[i].lepos.fOffsetNs -= 40.0;
            fU->fCalib[i].leneg.fOffsetNs -= 40.0;
            fU->fCalib[i].tepos.fOffsetNs -= 40.0;
            fU->fCalib[i].teneg.fOffsetNs -= 40.0;
         }
         
         // shift B signal
         fU->fCalib[33].lepos.fOffsetNs -= 40.0;
         fU->fCalib[33].leneg.fOffsetNs -= 40.0;
         fU->fCalib[33].tepos.fOffsetNs -= 40.0;
         fU->fCalib[33].teneg.fOffsetNs -= 40.0;
      }
#endif
      
#ifdef HAVE_ROOT

#define TDC_RANGE MAX_TDC_CHAN+2, 0-0.5, MAX_TDC_CHAN+1-0.5+1

      runinfo->fRoot->fOutputFile->cd(); // select correct ROOT directory
      TDirectory* dir = gDirectory->mkdir("dltdc");
      dir->cd(); // select correct ROOT directory

      dir->mkdir("fine_time")->cd();

      double max_time_hour = 3600;
      double max_time_day  = 24*max_time_hour;
      double max_time = 6.0*max_time_day;

      if (1) {
         int i=0;
         char name[256];
         char title[256];
         sprintf(name,  "aaa_tdc%02d_fine_time_le_pos", i);
         sprintf(title, "aaa_tdc%02d_LE_pos fine time distribution, ns", i);
         fHfineLeP0 = new TH2D(name, title, 100, -5, 15, 50, 0, max_time);
      }

      if (1) {
         int i=0;
         char name[256];
         char title[256];
         sprintf(name,  "aaa_tdc%02d_fine_time_te_neg", i);
         sprintf(title, "aaa_tdc%02d_TE_neg fine time distribution, ns", i);
         fHfineTeN0 = new TH2D(name, title, 100, -5, 15, 50, 0, max_time);
      }

      if (1) {
         int i=1;
         char name[256];
         char title[256];
         sprintf(name,  "aaa_tdc%02d_fine_time_le_pos", i);
         sprintf(title, "aaa_tdc%02d_LE_pos fine time distribution, ns", i);
         fHfineLeP1 = new TH2D(name, title, 100, -5, 15, 50, 0, max_time);
      }

      if (1) {
         int i=1;
         char name[256];
         char title[256];
         sprintf(name,  "aaa_tdc%02d_fine_time_te_neg", i);
         sprintf(title, "aaa_tdc%02d_TE_neg fine time distribution, ns", i);
         fHfineTeN1 = new TH2D(name, title, 100, -5, 15, 50, 0, max_time);
      }

      fHphaseLePerTdcChan = new TH2D("tdc_fine_bins_le_vs_tdc_chan", "TDC LE fine bin occupancy;tdc chan;fine bin number", TDC_RANGE, 100, 0-0.5, 100-0.5);
      fHphaseTePerTdcChan = new TH2D("tdc_fine_bins_te_vs_tdc_chan", "TDC TE fine bin occupancy;tdc chan;fine bin number", TDC_RANGE, 100, 0-0.5, 100-0.5);

      fHfineLePerTdcChan = new TH2D("tdc_fine_time_le_vs_tdc_chan", "TDC LE fine time occupancy;tdc chan;fine time, ns", TDC_RANGE, 200, -5, 15);
      fHfineTePerTdcChan = new TH2D("tdc_fine_time_te_vs_tdc_chan", "TDC TE fine time occupancy;tdc chan;fine time, ns", TDC_RANGE, 200, -5, 15);

      for (int i=0; i<=MAX_TDC_CHAN; i++) {
         char name[256];
         char title[256];

         sprintf(name,  "tdc%02d_fine_bins_le", i);
         sprintf(title, "tdc%02d_LE fine time bin occupancy", i);
         fHphaseLe[i] = new TH1D(name, title, 100, 0-0.5, 100-0.5);

         sprintf(name,  "tdc%02d_fine_bins_te", i);
         sprintf(title, "tdc%02d_TE fine time bin occupancy", i);
         fHphaseTe[i] = new TH1D(name, title, 100, 0-0.5, 100-0.5);

         BookFineBins(i);

         sprintf(name,  "tdc%02d_fine_time_le", i);
         sprintf(title, "tdc%02d_LE fine time distribution, ns", i);
         fHfineLe[i] = new TH1D(name, title, 200, -5, 15);

         sprintf(name,  "tdc%02d_fine_time_le_pos", i);
         sprintf(title, "tdc%02d_LE_pos fine time distribution, ns", i);
         fHfineLeP[i] = new TH1D(name, title, 200, -5, 15);

         sprintf(name,  "tdc%02d_fine_time_le_neg", i);
         sprintf(title, "tdc%02d_LE_neg fine time distribution, ns", i);
         fHfineLeN[i] = new TH1D(name, title, 200, -5, 15);

         sprintf(name,  "tdc%02d_fine_time_te", i);
         sprintf(title, "tdc%02d_TE fine time distribution, ns", i);
         fHfineTe[i] = new TH1D(name, title, 200, -5, 15);

         sprintf(name,  "tdc%02d_fine_time_te_pos", i);
         sprintf(title, "tdc%02d_TE_pos fine time distribution, ns", i);
         fHfineTeP[i] = new TH1D(name, title, 200, -5, 15);

         sprintf(name,  "tdc%02d_fine_time_te_neg", i);
         sprintf(title, "tdc%02d_TE_neg fine time distribution, ns", i);
         fHfineTeN[i] = new TH1D(name, title, 200, -5, 15);

         sprintf(name,  "tdc%02d_width_ns", i);
         sprintf(title, "tdc%02d pulse width, ns", i);
         fHwidth[i] = new TH1D(name, title, 100, -5, 5);
      }

      dir->mkdir("pulser")->cd();

      fHpulserLePerTdcChan = new TH2D("tdc_pulser_le_vs_tdc_chan", "Pulser LE time vs TDC channel;tdc channel;pulse le time, ns", TDC_RANGE, 400, -40, 40);
      fHpulserTePerTdcChan = new TH2D("tdc_pulser_te_vs_tdc_chan", "Pulser TE time vs TDC channel;tdc channel;pulse te time, ns", TDC_RANGE, 400, -40, 40);
      fHpulserWiPerTdcChan = new TH2D("tdc_pulser_width_vs_tdc_chan", "Pulser width vs TDC channel;tdc channel;pulse width, ns", TDC_RANGE, 400, 0, 80);

      fHpulserLeAll = new TH1D("tdc_pulser_le_all", "tdc pulser le all, ns", 400, -40, 40);
      fHpulserTeAll = new TH1D("tdc_pulser_te_all", "tdc_pulser_te_all, ns", 400, -40, 40);
      fHpulserWiAll = new TH1D("tdc_pulser_width_all", "tdc_pulser_width_all, ns", 400, 0, 80);

      for (int i=0; i<=MAX_TDC_CHAN; i++) {
         char name[256];
         char title[256];

         double shift40 = 0;

         if ((i>=16&&i<32) || (i==33)) { // B cable or B signal
            if (fFlags->fAdd40)
               shift40 = +40;
            if (fFlags->fSub40)
               shift40 = -40;
         }

         sprintf(name,  "tdc%02d_pulser_le", i);
         sprintf(title, "tdc%02d_pulser_le, ns", i);
         fHpulserLe[i] = new TH1D(name, title, 400, -40+shift40, 40+shift40);

         sprintf(name,  "tdc%02d_pulser_te", i);
         sprintf(title, "tdc%02d_pulser_te, ns", i);
         fHpulserTe[i] = new TH1D(name, title, 400, -40+shift40, 40+shift40);

         sprintf(name,  "tdc%02d_pulser_width", i);
         sprintf(title, "tdc%02d_pulser_width, ns", i);
         fHpulserWi[i] = new TH1D(name, title, 400, 0, 80);

         sprintf(name,  "tdc%02d_pulser_pp", i);
         sprintf(title, "tdc%02d_pulser_pp, phase TE vs phase LE;LE fine bin;TE fine bin", i);
         fHpulserPP[i] = new TH2D(name, title, 101, -50.5, 50.5, 101, -50.5, 50.5);

         sprintf(name,  "tdc%02d_pulser_tt", i);
         sprintf(title, "tdc%02d_pulser_tt, TE fine time vs LE fine time;LE fine time, ns;TE fine time, ns", i);
         fHpulserFF[i] = new TH2D(name, title, 200, -5, 15, 200, -5, 15);
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

      //fHdt_PP.resize(MAX_TDC_CHAN+1);
      fHdt_FF.resize(MAX_TDC_CHAN+1);

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

         //sprintf(name,  "tdc%02d_pp", i);
         //sprintf(title, "tdc%02d phase vs phase", i);
         //fHdt_PP[i] = new TH2D(name, title, 101, -50.5, 50.5, 101, -50.5, 50.5);

         sprintf(name,  "tdc%02d_ff", i);
         sprintf(title, "tdc%02d fine time vs fine time", i);
         fHdt_FF[i] = new TH2D(name, title, 200, -5, 15, 200, -5, 15);
      }

      dir->mkdir("unphysical")->cd();

      fHunphysical_map  = new TH1D("tdc_map", " TDC hit map", TDC_RANGE);
      fHunphysical_map2 = new TH2D("tdc_map2", "TDC pairs of hits map", TDC_RANGE, TDC_RANGE);

      fHunphysical_ns_01_14_w01 = new TH2D("tdc01_tdc14_vs_w01", "tdc01_le - tdc14_le, ns vs tdc01_width, ns", 200, -50, 50, 100, 0, 100);
      fHunphysical_ns_01_14_w14 = new TH2D("tdc01_tdc14_vs_w14", "tdc01_le - tdc14_le, ns vs tdc14_width, ns", 200, -50, 50, 100, 0, 100);

      fHunphysical_ns_01_14_c01 = new TH2D("tdc01_tdc14_vs_c10", "tdc01_le - tdc14_le, ns vs tdc01_hit_count", 200, -50, 50, 20, 0, 20);
      fHunphysical_ns_01_14_c14 = new TH2D("tdc01_tdc14_vs_c14", "tdc01_le - tdc14_le, ns vs tdc14_hit_count", 200, -50, 50, 20, 0, 20);

      fHunphysical_ns.resize((MAX_TDC_CHAN+1)*(MAX_TDC_CHAN+1));
      for (int i=0; i<=MAX_TDC_CHAN; i++) {
         for (int j=i+1; j<=MAX_TDC_CHAN; j++) {
            char name[256];
            char title[256];

            sprintf(name,  "tdc%02d_tdc%02d_le_ns", i, j);
            sprintf(title, "tdc%02d_le - tdc%02d_le, ns", i, j);
            fHunphysical_ns[i*(MAX_TDC_CHAN+1) + j] = new TH1D(name, title, 200, -50, 50);
            fHunphysical_ns[i*(MAX_TDC_CHAN+1) + j]->SetMinimum(0);
         }
      }

      dir->mkdir("unphysical_FF")->cd();

      //fHunphysical_PP.resize((MAX_TDC_CHAN+1)*(MAX_TDC_CHAN+1));
      fHunphysical_FF.resize((MAX_TDC_CHAN+1)*(MAX_TDC_CHAN+1));

      for (int i=0; i<=MAX_TDC_CHAN; i++) {
         for (int j=i+1; j<=MAX_TDC_CHAN; j++) {
            char name[256];
            char title[256];

            //sprintf(name,  "tdc%02d_tdc%02d_le_PP", i, j);
            //sprintf(title, "TDC fine bin tdc%02d_le vs tdc%02d_le, ns", i, j);
            //fHunphysical_PP[i*(MAX_TDC_CHAN+1) + j] = new TH2D(name, title, 101, -50.5, 50.5, 101, -50.5, 50.5);

            sprintf(name,  "tdc%02d_tdc%02d_le_FF", i, j);
            sprintf(title, "TDC fine time tdc%02d_le vs tdc%02d_le, ns", i, j);
            fHunphysical_FF[i*(MAX_TDC_CHAN+1) + j] = new TH2D(name, title, 200, -5, 15, 200, -5, 15);
         }
      }
#endif

      if (!fFlags->fCalibFineTime) {
         PlotFineBins();
      }
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

      ProcessSorted(runinfo, true); // flush accumulated data

      if (fFlags->fCalibFineTime) {
         printf("DlTdcModule::EndRun: Saving TDC fine time calibrations for run %d\n", runinfo->fRunNo);
         fU->UpdateCalib();
         fU->SaveCalib(runinfo->fRunNo);
         PlotFineBins();
      }

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

      if (fFlags->fCalibOffsets) {
         printf("offsets json:\n%s\n", s.c_str());

         printf("DlTdcModule::EndRun: Saving TDC offsets calibrations for run %d\n", runinfo->fRunNo);
         char fname[256];
         sprintf(fname, "dlcfmdb/dltdc_offsets_%06d.json", runinfo->fRunNo);
         FILE *fp = fopen(fname, "w");
         if (fp) {
            fprintf(fp, "%s", s.c_str());
            fclose(fp);
         }
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

   size_t fPulserMaster = 0;
   bool   fPulserMasterHave = false;

   void AnalyzeTdcEvent(const DlTdcEvent& t)
   {
      if (fFlags->fDebug) {
         printf("EVENT ABT %d%d%d\n", t.HaveCh(fConf->fChanA), t.HaveCh(fConf->fChanB), t.HaveCh(fConf->fChanT));
      }
      
      ///////// check for triggered event ///////////

      if (fFlags->fTriggered && !t.HaveCh(fConf->fChanT)) {
         return;
      }

      ///////// debug spikes in dark noise data ///////////

      if (0) {
         if (t.HaveCh(1) && t.HaveCh(14)) {
            double dt_ns = subtract_ns(t.GetCh(1).fLe, t.GetCh(14).fLe);
            //if (t.GetCh(1).fWidthNs > 10.0)
            //   return;
            //if (t.GetCh(14).fWidthNs > 10.0)
            //   return;
            //if (fabs(dt_ns) > 5.0)
            //   return;
            printf("dt 01-14: %6.3f ns, width %6.3f %6.3f ns\n", dt_ns, t.GetCh(1).fWidthNs, t.GetCh(14).fWidthNs);
         } else {
            return;
         }
      }
      
      int ntdchits = 0;

      for (int tdc_ch=0; tdc_ch<=MAX_TDC_CHAN; tdc_ch++) {
         if (t.HaveCh(tdc_ch)) {
            ntdchits++;
         }
      }

      if (fFlags->fDebug) {
         std::string s = "";
         for (int tdc_ch=0; tdc_ch<=MAX_TDC_CHAN; tdc_ch++) {
            if (t.HaveCh(tdc_ch)) {
               s += "H";
            } else {
               s += ".";
            }
         }

         printf("EVENT %s, ABT %d%d%d, %d hits\n", s.c_str(), t.HaveCh(fConf->fChanA), t.HaveCh(fConf->fChanB), t.HaveCh(fConf->fChanT), ntdchits);
      }

      //if (ntdchits != 7)
      //   return;

      if (0) {
         printf("have %zu hits:\n", t.fTdcHits.size());
         double t1 = -9999;
         double t14 = -8888;
         for (size_t i=0; i<t.fTdcHits.size(); i++) {
            if (t.fTdcHits[i].ch == 1 || t.fTdcHits[i].ch == 14) {
               t.fTdcHits[i].Print();
               if (t.fTdcHits[i].le) {
                  if (t.fTdcHits[i].ch == 1) {
                     t1 = t.fTdcHits[i].time_sec;
                     printf(" t1-t14 is %6.3f ns", sec_to_ns(t1-t14));
                  }
                  if (t.fTdcHits[i].ch == 14) {
                     t14 = t.fTdcHits[i].time_sec;
                     printf(" t1-t14 is %6.3f ns", sec_to_ns(t1-t14));
                  }
               }
               printf("\n");
            }
         }
      }

      ///////// pulser calibration histograms ///////////
         
      if (!fPulserMasterHave) {
         for (size_t itdc=0; itdc<=MAX_TDC_CHAN; itdc++) {
            //printf("tdc %d, up down %d %d\n", itdc, t.fHits[itdc].fUp, t.fHits[itdc].fDown);
            if (t.HaveCh(itdc)) {
               fPulserMaster = itdc;
               fPulserMasterHave = true;
               printf("Pulser master channel selected tdc %d\n", (int)itdc);
               break;
            }
         }
      }

      if (fPulserMasterHave) {
         for (size_t itdc=0; itdc<=MAX_TDC_CHAN; itdc++) {
            //printf("tdc %d, up down %d %d\n", itdc, t.fHits[itdc].fUp, t.fHits[itdc].fDown);
            if (t.HaveCh(itdc)) {
               if (itdc != fPulserMaster) {
                  fHpulserLe[itdc]->Fill(subtract_ns(t.GetCh(itdc).fLe, t.GetCh(fPulserMaster).fLe));
                  fHpulserTe[itdc]->Fill(subtract_ns(t.GetCh(itdc).fTe, t.GetCh(fPulserMaster).fTe));

                  fHpulserLeAll->Fill(subtract_ns(t.GetCh(itdc).fLe, t.GetCh(fPulserMaster).fLe));
                  fHpulserTeAll->Fill(subtract_ns(t.GetCh(itdc).fTe, t.GetCh(fPulserMaster).fTe));

                  fHpulserLePerTdcChan->Fill(itdc, subtract_ns(t.GetCh(itdc).fLe, t.GetCh(fPulserMaster).fLe));
                  fHpulserTePerTdcChan->Fill(itdc, subtract_ns(t.GetCh(itdc).fTe, t.GetCh(fPulserMaster).fTe));
               }

               fHpulserWiPerTdcChan->Fill(itdc, t.fHits[itdc].fWidthNs);
               fHpulserWiAll->Fill(t.fHits[itdc].fWidthNs);
               fHpulserWi[itdc]->Fill(t.fHits[itdc].fWidthNs);

               fHpulserPP[itdc]->Fill(t.fHits[itdc].fLe.phase, t.fHits[itdc].fTe.phase);
               fHpulserFF[itdc]->Fill(t.fHits[itdc].fLe.fine_ns, t.fHits[itdc].fTe.fine_ns);
            }
         }
      }

      ///////// create width calibration histograms ///////////

      for (size_t itdc=0; itdc<=MAX_TDC_CHAN; itdc++) {
         //printf("tdc %d, up down %d %d\n", itdc, t.fHits[itdc].fUp, t.fHits[itdc].fDown);
         if (!t.fHits[itdc].fUp && t.fHits[itdc].fDown) {
            fHwidth[itdc]->Fill(t.fHits[itdc].fWidthNs);
         }
      }

      ///////// plot unphysical time pairs /////////

      for (int i=0; i<=MAX_TDC_CHAN; i++) {
         if (t.HaveCh(i)) {
            fHunphysical_map->Fill(i);
            for (int j=i+1; j<=MAX_TDC_CHAN; j++) {
               if (t.HaveCh(j)) {
                  fHunphysical_map2->Fill(i,j);
                  fHunphysical_map2->Fill(j,i);
                  fHunphysical_ns[i*(MAX_TDC_CHAN+1) + j]->Fill(subtract_ns(t.GetCh(i).fLe, t.GetCh(j).fLe));
                  //fHunphysical_PP[i*(MAX_TDC_CHAN+1) + j]->Fill(t.GetCh(i).fLe.phase, t.GetCh(j).fLe.phase);
                  fHunphysical_FF[i*(MAX_TDC_CHAN+1) + j]->Fill(t.GetCh(i).fLe.fine_ns, t.GetCh(j).fLe.fine_ns);

                  //if (i==0 && j==1) {
                  //   printf("tdc00-01 %.0f\n", subtract_ns(t.GetCh(i).fLe, t.GetCh(j).fLe));
                  //   t.GetCh(i).Print();
                  //   t.GetCh(j).Print();
                  //}

                  if (i==1 && j==14) {
                     fHunphysical_ns_01_14_w01->Fill(subtract_ns(t.GetCh(i).fLe, t.GetCh(j).fLe), t.GetCh(i).fWidthNs);
                     fHunphysical_ns_01_14_w14->Fill(subtract_ns(t.GetCh(i).fLe, t.GetCh(j).fLe), t.GetCh(j).fWidthNs);

                     fHunphysical_ns_01_14_c01->Fill(subtract_ns(t.GetCh(i).fLe, t.GetCh(j).fLe), t.GetCh(i).fCount);
                     fHunphysical_ns_01_14_c14->Fill(subtract_ns(t.GetCh(i).fLe, t.GetCh(j).fLe), t.GetCh(j).fCount);

                     //printf("dt %6.3f ns\n", subtract_ns(t.GetCh(i).fLe, t.GetCh(j).fLe));
                     //t.GetCh(i).Print();
                     //t.GetCh(j).Print();
                  }
                  
                  //printf("count %d %d\n", t.GetCh(i).fCount, t.GetCh(j).fCount);
               }
            }
         }
      }
   }

   double fLastXtime = 0;

   bool fFirstTrig = true;

   std::vector<DlTdcHit> fPrevLe;

   std::list<DlTdcHit> fList;
   uint64_t fListEpoch = 0;
   uint32_t fLastCoarse = 0;
   bool     fEpochChanging = false;

   void PrintSorted() const
   {
      printf("Sorted list has %zu entries:\n", fList.size());
      for (const DlTdcHit& h: fList) {
         printf("coarse: 0x%016lx, ", h.coarse); h.Print(); printf("\n");
      }
      printf("Sorted list end!\n");
   }

   bool AddSorted(DlTdcHit& h)
   {
      bool debug_sort = false;

      if (!fEpochChanging) {
         if (h.coarse < h.kEpoch/4 && fLastCoarse > (h.kEpoch/4)*3) {
            printf("EPOCH!\n");
            fListEpoch += h.kEpoch;
            fEpochChanging = true;
         }
         
      }

      if (h.coarse > h.kEpoch/4) {
         fEpochChanging = false;
      }

      fLastCoarse = h.coarse;

      if (fEpochChanging) {
         if (h.coarse > (h.kEpoch/4)*3) { // previous epoch
            //printf("P ");
            h.coarse += fListEpoch - h.kEpoch;
         } else {
            //printf("N ");
            h.coarse += fListEpoch; // new epoch
         }
      } else {
         //printf("E ");
         h.coarse += fListEpoch;
      }

#if 0
      if (h.coarse == 0x0000000006ed96c2 && h.ch == 10 && h.le) {
         printf("HERE LE!!!\n");
      }

      if (h.coarse == 0x0000000006ed96c1 && h.ch == 10 && h.te) {
         printf("HERE TE!!!\n");
      }

      if (h.coarse == 0x0000000006ef4470 && h.ch == 10 && h.le) {
         printf("HERE LE!!!\n");
      }

      if (h.coarse == 0x0000000006ef446f && h.ch == 10 && h.te) {
         printf("HERE TE!!!\n");
      }

      if (h.ch == 10)
         debug_sort = true;
#endif

      if (debug_sort) {
         printf("coarse: 0x%016lx, ", h.coarse); h.Print(); printf(" ---> insert\n");
      }
         
      if (fList.empty()) {
         if (debug_sort)
            printf("FIRST!\n");
         fList.push_back(h);
      } else if (h.coarse > fList.back().coarse) {
         if (debug_sort)
            printf("BACK!\n");
         fList.push_back(h);
      } else {
         bool pushed = false;
         for (auto r=fList.rbegin(); r!=fList.rend(); r++) {
            if (debug_sort) {
               printf("coarse: 0x%016lx, ", r->coarse); r->Print(); printf("\n");
            }
            if ((h.coarse+1 == r->coarse) && h.te && r->le) {
               if (debug_sort)
                  printf("INSERT TE!\n");
               fList.insert(r.base(), h);
               pushed = true;
               break;
            } else 
            if (h.coarse >= r->coarse) {
               if ((h.coarse == r->coarse) && h.le && r->te) {
                  if (debug_sort)
                     printf("SKIP TE!\n");
               } else if ((h.coarse == r->coarse+1) && h.le && r->te) {
                  if (debug_sort)
                     printf("SKIP TE1!\n");
               } else {
                  if (debug_sort)
                     printf("INSERT!\n");
                  fList.insert(r.base(), h);
                  pushed = true;
                  break;
               }
            }
         }
         if (!pushed) {
            printf("PUSH FRONT!!!\n");
            fList.push_front(h);
         }
      }

      return debug_sort;
   }

   void ProcessSorted(TARunInfo* runinfo, bool flush = false)
   {
      if (fFlags->fPrintSortedHits) {
         printf("ProcessSorted: have %zu hits\n", fList.size());
      }

      while (!fList.empty()) {
         if (!flush) {
            if (fList.size() < 100)
               break;
         }

         DlTdcHit h = fList.front();
         fList.pop_front();
         
         fU->ComputeTimes(&h);

         if (fFlags->fPrintSortedHits) {
            h.Print(); printf("\n");
         }
         
         //static double prev_time_sec = 0;
         //
         //if (h.le) {
         //   if (fFlags->fDebug) {h.Print(); printf(", dt %4.0f ns\n", sec_to_ns(h.time_sec - prev_time_sec));}
         //   prev_time_sec = h.time_sec;
         //}
         
         if (fFlags->fCalibFineTime) {
            fU->fCalib[h.ch].AddHit(h);
         }
         
         if (h.le) {
            double dt = subtract_ns(h, fPrevLe[h.ch]);
            //printf("ch %d: dt %.3f ns\n", h.ch, dt);
            fHdt1Le[h.ch]->Fill(dt);
            fHdt2Le[h.ch]->Fill(dt);
            fHdt3Le[h.ch]->Fill(dt);
            fHdt4Le[h.ch]->Fill(dt);
            
            //fHdt_PP[h.ch]->Fill(h.phase, fPrevLe[h.ch].phase);
            fHdt_FF[h.ch]->Fill(h.fine_ns, fPrevLe[h.ch].fine_ns);
            
            fPrevLe[h.ch] = h;
         }
         
         double hit_dt_ns = sec_to_ns(h.time_sec - fCt->max_time_sec);
         
         fHhitdt1ns->Fill(hit_dt_ns);
         fHhitdt2ns->Fill(hit_dt_ns);
         fHhitdt3ns->Fill(hit_dt_ns);
         fHhitdt4ns->Fill(hit_dt_ns);
         
         //printf("TTX %.9f %.9f sec, first %.9f, last %.9f sec\n", fCt->min_time_sec, fCt->max_time_sec, fCt->first_time_sec, fCt->last_time_sec);
         
         //if (h.le && hit_dt_ns > 80.0) {
         
         if (hit_dt_ns > 80.0) {
            double event_dt_ns = sec_to_ns(h.time_sec - fCt->min_time_sec);
            
            if (fFlags->fDebug) {
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
               
               //FinishEventT(fPrevEventTimeSec, *fCt);

               if (fFlags->fPrintEventHits) {
                  printf("Event hits:\n");
                  fCt->Print();
               }
               
               runinfo->AddToFlowQueue(new DlTdcEventFlow(NULL, fCt));
               
               fCt = new DlTdcEvent();
               fCt->Init(MAX_TDC_CHAN+1);
               
               fPrevEventTimeSec = fCt->min_time_sec;
               
               fCt->Clear();
            }
         }
         
#ifdef HAVE_ROOT
         if (h.ch >= 0 && h.ch <= MAX_TDC_CHAN) {
            if (h.le) {
               if (h.phase < 0) {
                  fHphaseLePerTdcChan->Fill(h.ch, 0 - h.phase);
                  fHphaseLe[h.ch]->Fill(0 - h.phase);
                  fHfineLePerTdcChan->Fill(h.ch, h.fine_ns);
                  fHfineLe[h.ch]->Fill(h.fine_ns);
                  fHfineLeN[h.ch]->Fill(h.fine_ns);
               } else {
                  fHphaseLePerTdcChan->Fill(h.ch, 50 + h.phase);
                  fHphaseLe[h.ch]->Fill(50 + h.phase);
                  fHfineLePerTdcChan->Fill(h.ch, h.fine_ns);
                  fHfineLe[h.ch]->Fill(h.fine_ns);
                  fHfineLeP[h.ch]->Fill(h.fine_ns);
                  if (h.ch == 0) {
                     fHfineLeP0->Fill(h.fine_ns, h.time_sec);
                     //printf("X000 time %.0f\n", h.time_sec);
                  }
                  if (h.ch == 1) {
                     fHfineLeP1->Fill(h.fine_ns, h.time_sec);
                  }
               }
            }
            if (h.te) {
               if (h.phase < 0) {
                  fHphaseTePerTdcChan->Fill(h.ch, 0 - h.phase);
                  fHphaseTe[h.ch]->Fill(0 - h.phase);
                  fHfineTePerTdcChan->Fill(h.ch, h.fine_ns);
                  fHfineTe[h.ch]->Fill(h.fine_ns);
                  fHfineTeN[h.ch]->Fill(h.fine_ns);
                  if (h.ch == 0) {
                     fHfineTeN0->Fill(h.fine_ns, h.time_sec);
                  }
                  if (h.ch == 1) {
                     fHfineTeN1->Fill(h.fine_ns, h.time_sec);
                  }
               } else {
                  fHphaseTePerTdcChan->Fill(h.ch, 50 + h.phase);
                  fHphaseTe[h.ch]->Fill(50 + h.phase);
                  fHfineTePerTdcChan->Fill(h.ch, h.fine_ns);
                  fHfineTe[h.ch]->Fill(h.fine_ns);
                  fHfineTeP[h.ch]->Fill(h.fine_ns);
               }
            }
         }
#endif
         fCt->AddHit8(h);
      }
   }

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
            fCt->Init(MAX_TDC_CHAN+1);
         }

         if (fFlags->fPrintRawHits) {
            printf("TDC bank with %d hits\n", tdc_nw64);
         }

         bool debug_sort = false;

         DlTdcHit h;
     
	 for (int i=0; i<tdc_nw64; i++) {
            uint32_t wlo = tdc_data[i*2+0];
            uint32_t whi = tdc_data[i*2+1];

            if (wlo == 0 && whi == 0) {
               if (i <= 8192) 
                  printf("bank data at offset %d is zero\n", i);
               continue;
            }

            //printf("%2d: 0x%08x 0x%08x\n", i, whi, wlo);

            h.Unpack(wlo, whi);
            //fU->Unpack(&h, wlo, whi);

            if (fFlags->fPrintRawHits) {
               h.Print(); printf("\n");
            }

            debug_sort |= AddSorted(h);
         } // loop over data

         if (debug_sort)
            PrintSorted();

         if (fList.size() > 200) {
            ProcessSorted(runinfo);
         }

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

   void BookFineBins(int i)
   {
      const DlTdcFineCalib* c = &fU->fCalib[i];
      if (c) {
         char name[256];
         char title[256];
         
         sprintf(name,  "tdc%02d_fine_size_le", i);
         sprintf(title, "tdc%02d_LE fine time bin size, ns", i);

         fHfineBinsLe.push_back(new TProfile(name, title, 2*50, -0.5, 2*50-0.5));
         
         sprintf(name,  "tdc%02d_fine_size_te", i);
         sprintf(title, "tdc%02d_TE fine time bin size, ns", i);
         
         fHfineBinsTe.push_back(new TProfile(name, title, 2*50, -0.5, 2*50-0.5));
      } else {
         fHfineBinsLe.push_back(NULL);
         fHfineBinsTe.push_back(NULL);
      }
   }

   void PlotFineBins()
   {
      printf("PlotFineBins!\n");

      for (int i=0; i<(int)fU->fCalib.size(); i++) {
         const DlTdcFineCalib* c = &fU->fCalib[i];
         if (c) {
            TProfile* h = fHfineBinsLe[i];
            if (h) {
               for (int j=0; j<(int)c->leneg.fBinWidthNs.size(); j++) {
                  if (c->leneg.fBinWidthNs[j] > 0) {
                     h->Fill(  0+j, c->leneg.fBinWidthNs[j]);
                  }
               }

               for (int j=0; j<(int)c->lepos.fBinWidthNs.size(); j++) {
                  if (c->lepos.fBinWidthNs[j] > 0) {
                     h->Fill( 50+j, c->lepos.fBinWidthNs[j]);
                  }
               }
            }
 
            h = fHfineBinsTe[i];
            if (h) {
               for (int j=0; j<(int)c->teneg.fBinWidthNs.size(); j++) {
                  if (c->teneg.fBinWidthNs[j] > 0) {
                     h->Fill(  0+j, c->teneg.fBinWidthNs[j]);
                  }
               }

               for (int j=0; j<(int)c->tepos.fBinWidthNs.size(); j++) {
                  if (c->tepos.fBinWidthNs[j] > 0) {
                     h->Fill( 50+j, c->tepos.fBinWidthNs[j]);
                  }
               }
            }
         }
      }
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
      printf("--dltdc-print-raw-hits -- print hits in the TDC bank\n");
      printf("--dltdc-print-sorted-hits -- print hits sorted by time\n");
      printf("--dltdc-print-event-hits -- print hits in each event\n");
      printf("--dltdc-add40 -- add 40 ns to B cable\n");
      printf("--dltdc-sub40 -- subtract 40 ns from B cable\n");
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
         if (args[i] == "--dltdc-print-raw-hits") {
            fFlags.fPrintRawHits = true;
         }
         if (args[i] == "--dltdc-print-sorted-hits") {
            fFlags.fPrintSortedHits = true;
         }
         if (args[i] == "--dltdc-print-event-hits") {
            fFlags.fPrintEventHits = true;
         }
         if (args[i] == "--dltdc-add40") {
            fFlags.fAdd40 = true;
         }
         if (args[i] == "--dltdc-sub40") {
            fFlags.fSub40 = true;
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
