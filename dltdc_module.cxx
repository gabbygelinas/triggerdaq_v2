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
   bool fPrint = false;
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
   TH1D* fHphaseLe[MAX_TDC_CHAN+1];
   TH1D* fHphaseTe[MAX_TDC_CHAN+1];

   std::vector<TProfile*> fHfineBinsLe;
   std::vector<TProfile*> fHfineBinsTe;

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
      runinfo->fRoot->fOutputFile->cd(); // select correct ROOT directory
      TDirectory* dir = gDirectory->mkdir("dltdc");
      dir->cd(); // select correct ROOT directory

      dir->mkdir("fine_time")->cd();

      for (int i=0; i<=MAX_TDC_CHAN; i++) {
         char name[256];
         char title[256];

         sprintf(name,  "tdc%02d_fine_bin_le", i);
         sprintf(title, "tdc%02d_LE fine time bin occupancy", i);
         fHphaseLe[i] = new TH1D(name, title, 101, -50, 50);

         sprintf(name,  "tdc%02d_fine_bin_te", i);
         sprintf(title, "tdc%02d_TE fine time bin occupancy", i);
         fHphaseTe[i] = new TH1D(name, title, 101, -50, 50);

         BookFineBins(i);

         sprintf(name,  "tdc%02d_fine_time_le", i);
         sprintf(title, "tdc%02d_LE fine time distribution, ns", i);
         fHfineLe[i] = new TH1D(name, title, 200, -5, 15);

         sprintf(name,  "tdc%02d_fine_time_te", i);
         sprintf(title, "tdc%02d_TE fine time distribution, ns", i);
         fHfineTe[i] = new TH1D(name, title, 200, -5, 15);

         sprintf(name,  "tdc%02d_width_ns", i);
         sprintf(title, "tdc%02d pulse width, ns", i);
         fHwidth[i] = new TH1D(name, title, 100, -5, 5);
      }

      dir->mkdir("pulser")->cd();

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

      ///////// pulser calibration histograms ///////////
         
      if (!fPulserMasterHave) {
         for (size_t ch=0; ch<=MAX_TDC_CHAN; ch++) {
            //printf("ch %d, up down %d %d\n", ch, t.fHits[ch].fUp, t.fHits[ch].fDown);
            if (t.HaveCh(ch)) {
               fPulserMaster = ch;
               fPulserMasterHave = true;
               printf("Pulser master channel selected ch %d\n", (int)ch);
               break;
            }
         }
      }

      if (fPulserMasterHave) {
         for (size_t ch=0; ch<=MAX_TDC_CHAN; ch++) {
            //printf("ch %d, up down %d %d\n", ch, t.fHits[ch].fUp, t.fHits[ch].fDown);
            if (t.HaveCh(ch)) {
               if (ch != fPulserMaster) {
                  fHpulserLeAll->Fill(subtract_ns(t.GetCh(ch).fLe, t.GetCh(fPulserMaster).fLe));
                  fHpulserTeAll->Fill(subtract_ns(t.GetCh(ch).fTe, t.GetCh(fPulserMaster).fTe));

                  fHpulserLe[ch]->Fill(subtract_ns(t.GetCh(ch).fLe, t.GetCh(fPulserMaster).fLe));
                  fHpulserTe[ch]->Fill(subtract_ns(t.GetCh(ch).fTe, t.GetCh(fPulserMaster).fTe));
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

      ///////// plot unphysical time pairs /////////

      for (int i=0; i<=MAX_TDC_CHAN; i++) {
         for (int j=i+1; j<=MAX_TDC_CHAN; j++) {
            if (t.HaveCh(i) && t.HaveCh(j)) {
               fHunphysical_pair_ns[i*(MAX_TDC_CHAN+1) + j]->Fill(subtract_ns(t.GetCh(i).fLe, t.GetCh(j).fLe));
            }
         }
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
            fCt->Init(MAX_TDC_CHAN+1);
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
         
         sprintf(name,  "tdc%02d_fine_bins_le", i);
         sprintf(title, "tdc%02d_LE fine time bin size, ns", i);

         fHfineBinsLe.push_back(new TProfile(name, title, c->leneg.fBinWidthNs.size()+c->lepos.fBinWidthNs.size()-1, 0.5-c->leneg.fBinWidthNs.size(), +c->lepos.fBinWidthNs.size()-0.5));
         
         sprintf(name,  "tdc%02d_fine_bins_te", i);
         sprintf(title, "tdc%02d_TE fine time bin size, ns", i);
         
         fHfineBinsTe.push_back(new TProfile(name, title, c->teneg.fBinWidthNs.size()+c->tepos.fBinWidthNs.size()-1, 0.5-c->teneg.fBinWidthNs.size(), +c->tepos.fBinWidthNs.size()-0.5));
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
               for (int j=1; j<(int)c->leneg.fBinWidthNs.size(); j++) {
                  h->Fill(-j, c->leneg.fBinWidthNs[j]);
               }
               
               for (int j=1; j<(int)c->lepos.fBinWidthNs.size(); j++) {
                  h->Fill(+j, c->lepos.fBinWidthNs[j]);
               }
            }
 
            h = fHfineBinsTe[i];
            if (h) {
               for (int j=1; j<(int)c->teneg.fBinWidthNs.size(); j++) {
                  h->Fill(-j, c->teneg.fBinWidthNs[j]);
               }
               
               for (int j=1; j<(int)c->tepos.fBinWidthNs.size(); j++) {
                  h->Fill(+j, c->tepos.fBinWidthNs[j]);
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
      printf("--dltdc-print -- print events\n");
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
