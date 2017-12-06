//
// final_module.cxx
//
// final analysis of TPC data
//
// K.Olchanski
//

#include <stdio.h>

#include "manalyzer.h"
#include "midasio.h"

#include <assert.h> // assert()

#include <vector>
#include <deque>

#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TMath.h"
#include "TGraphPolar.h"

#include "Unpack.h"
#include "AgFlow.h"

#define NUM_AW 512
#define MAX_AW_AMP 16000
#define MAX_TIME 7000

#define MAX_PAD_AMP 4100

// number of pad columns
#define NUM_PC (8*4)

#define DELETE(x) if (x) { delete (x); (x) = NULL; }

#define MEMZERO(p) memset((p), 0, sizeof(p))

class FinalModule: public TARunObject
{
public:
   TCanvas* fC = NULL;
   TCanvas* fPH = NULL;

   TH1D* h_num_aw_hits;
   TH1D* h_num_pad_hits;
   TH2D* h_num_aw_pad_hits;

   TH1D* h_aw_time;
   TH1D* h_aw_amp;
   TH2D* h_aw_amp_time;

   TH1D* h_aw_map;
   TH2D* h_aw_map_time;
   TH2D* h_aw_map_amp;

   TH2D* h_aw_aw_hits;
   TH2D* h_aw_aw_time;
   TH2D* h_aw_aw_amp;

   TH1D* h_aw_286;
   TH1D* h_aw_287;
   TH1D* h_aw_288;
   TH1D* h_aw_289;
   TH1D* h_aw_290;

   TH1D* h_aw_299;
   TH1D* h_aw_300;
   TH1D* h_aw_301;

   TH1D* h_aw_310;
   TH1D* h_aw_320;
   TH1D* h_aw_330;
   TH1D* h_aw_340;

   TH1D* h_aw_352;

   TH1D* h_pad_time;
   TH1D* h_pad_amp;
   TH2D* h_pad_amp_time;

   TH2D *h_pad_amp_pad;
   TH2D *h_pad_time_pad;

   TH2D* h_aw_pad_hits;

   TH2D* h_aw_pad_time;
   TH2D* h_aw_pad_time_drift;

   TH2D* h_aw_pad_amp_pc;

   // time across to next channel
   TH1D* h_cal_adc05_0_1 = NULL;

   // time across to another 100 MHz ADC on the same module
   TH1D* h_cal_adc05_0_4 = NULL;
   TH1D* h_cal_adc05_0_8 = NULL;
   TH1D* h_cal_adc05_0_12 = NULL;

   // time across to the 62.5 MHz 32ch ADC on the same module
   TH1D* h_cal_adc05_0_16 = NULL;

   // time across to next channel 62.5 MHz ADC on the same module
   TH1D* h_cal_adc05_16_17 = NULL;
   TH1D* h_cal_adc06_16_20 = NULL;

   // time across to another 62.5 MHz ADC on the same module
   TH1D* h_cal_adc05_16_24 = NULL;
   TH1D* h_cal_adc05_16_32 = NULL;
   TH1D* h_cal_adc05_16_40 = NULL;

   // time across to the next module, 100 MHz and 62.5 MHz ADC
   TH1D* h_cal_adc_05_06_chan0 = NULL;
   TH1D* h_cal_adc_05_06_chan16 = NULL;

   // time across to next channel
   TH1D* h_cal_time_pos00_seqsca_04_05 = NULL;
   TH1D* h_cal_time_pos01_seqsca_04_05 = NULL;
   TH1D* h_cal_time_pos02_seqsca_04_05 = NULL;

   // time across to next SCA
   TH1D* h_cal_time_pos01_seqsca_04_84 = NULL;
   TH1D* h_cal_time_pos02_seqsca_04_84 = NULL;

   // time across to next ADC
   TH1D* h_cal_time_pos01_seqsca_04_164 = NULL;
   TH1D* h_cal_time_pos02_seqsca_04_164 = NULL;

   // time across 2 PWB boards
   TH1D* h_cal_time_pos_01_02_seqsca04 = NULL;
   TH1D* h_cal_time_pos_01_02_seqsca05 = NULL;
   TH1D* h_cal_time_pos_01_03_seqsca04 = NULL;
   TH1D* h_cal_time_pos_01_03_seqsca05 = NULL;

   // time across 100MHz ADC to PWB
   TH1D* h_cal_time_pos01_seqsca_04_minus_adc5_0 = NULL;
   TH1D* h_cal_time_pos01_seqsca_04_minus_adc5_16 = NULL;
   TH1D* h_cal_time_pos02_seqsca_04_minus_adc5_0 = NULL;
   TH1D* h_cal_time_pos02_seqsca_04_minus_adc5_16 = NULL;

   FinalModule(TARunInfo* runinfo)
      : TARunObject(runinfo)
   {
      printf("FinalModule::ctor!\n");

      fC = new TCanvas("fC", "FinalModule event display", 800, 800);
      //fPH = new TCanvas("fPH","Pulseheights",600,1000);
      if (fPH) {
         fPH->Divide(1,2);
      }
   }

   ~FinalModule()
   {
      printf("FinalModule::dtor!\n");
      DELETE(fC);
      DELETE(fPH);
   }

   void BeginRun(TARunInfo* runinfo)
   {
      printf("FinalModule::BeginRun, run %d, file %s\n", runinfo->fRunNo, runinfo->fFileName.c_str());
      time_t run_start_time = runinfo->fOdb->odbReadUint32("/Runinfo/Start time binary", 0, 0);
      printf("ODB Run start time: %d: %s", (int)run_start_time, ctime(&run_start_time));
      runinfo->fRoot->fOutputFile->cd(); // select correct ROOT directory

      runinfo->fRoot->fOutputFile->cd();
      TDirectory* dir = gDirectory->mkdir("final");
      dir->cd(); // select correct ROOT directory

      dir->mkdir("summary")->cd();

      h_num_aw_hits = new TH1D("h_num_aw_hits", "number of anode wire hits", 100, 0, 100);
      h_aw_time = new TH1D("h_aw_time", "aw hit time", 100, 0, MAX_TIME);
      h_aw_amp = new TH1D("h_aw_amp", "aw hit pulse height", 100, 0, MAX_AW_AMP);
      h_aw_amp_time = new TH2D("h_aw_amp_time", "aw p.h. vs time", 100, 0, MAX_TIME, 50, 0, MAX_AW_AMP);

      h_aw_map = new TH1D("h_aw_map", "aw hit occupancy", NUM_AW, -0.5, NUM_AW-0.5);
      h_aw_map_time = new TH2D("h_aw_map_time", "aw hit time vs wire", NUM_AW, -0.5, NUM_AW-0.5, 50, 0, MAX_TIME);
      h_aw_map_amp  = new TH2D("h_aw_map_amp", "aw hit p.h. vs wire", NUM_AW, -0.5, NUM_AW-0.5, 50, 0, MAX_AW_AMP);

      h_aw_aw_hits = new TH2D("h_aw_aw_hits", "hits in aw vs aw", NUM_AW, -0.5, NUM_AW-0.5, NUM_AW, -0.5, NUM_AW-0.5);
      h_aw_aw_time = new TH2D("h_aw_aw_time", "time in aw vs aw", 50, 0, MAX_TIME, 50, 0, MAX_TIME);
      h_aw_aw_amp  = new TH2D("h_aw_aw_amp",  "p.h. in aw vs aw", 50, 0, MAX_AW_AMP, 50, 0, MAX_AW_AMP);

      h_aw_286 = new TH1D("h_aw_286", "h_aw_286", NUM_AW, -0.5, NUM_AW-0.5);
      h_aw_287 = new TH1D("h_aw_287", "h_aw_287", NUM_AW, -0.5, NUM_AW-0.5);
      h_aw_288 = new TH1D("h_aw_288", "h_aw_288", NUM_AW, -0.5, NUM_AW-0.5);
      h_aw_289 = new TH1D("h_aw_289", "h_aw_289", NUM_AW, -0.5, NUM_AW-0.5);
      h_aw_290 = new TH1D("h_aw_290", "h_aw_290", NUM_AW, -0.5, NUM_AW-0.5);

      h_aw_299 = new TH1D("h_aw_299", "h_aw_299", NUM_AW, -0.5, NUM_AW-0.5);
      h_aw_300 = new TH1D("h_aw_300", "h_aw_300", NUM_AW, -0.5, NUM_AW-0.5);
      h_aw_301 = new TH1D("h_aw_301", "h_aw_301", NUM_AW, -0.5, NUM_AW-0.5);

      h_aw_310 = new TH1D("h_aw_310", "h_aw_310", NUM_AW, -0.5, NUM_AW-0.5);
      h_aw_320 = new TH1D("h_aw_320", "h_aw_320", NUM_AW, -0.5, NUM_AW-0.5);
      h_aw_330 = new TH1D("h_aw_330", "h_aw_330", NUM_AW, -0.5, NUM_AW-0.5);
      h_aw_340 = new TH1D("h_aw_340", "h_aw_340", NUM_AW, -0.5, NUM_AW-0.5);

      h_aw_352 = new TH1D("h_aw_352", "h_aw_352", NUM_AW, -0.5, NUM_AW-0.5);

      h_num_pad_hits = new TH1D("h_num_pad_hits", "number of cathode pad hits", 100, 0, 100);
      h_pad_time = new TH1D("h_pad_time", "pad hit time", 100, 0, MAX_TIME);
      h_pad_amp = new TH1D("h_pad_amp", "pad hit pulse height", 100, 0, MAX_PAD_AMP);
      h_pad_amp_time = new TH2D("h_pad_amp_time", "pad p.h vs time", 50, 0, MAX_TIME, 50, 0, MAX_PAD_AMP);
      int npads = MAX_FEAM*MAX_FEAM_PAD_COL*MAX_FEAM_PAD_ROWS;
      h_pad_amp_pad = new TH2D("h_pad_amp_pad", "pad p.h vs pad number", npads, -0.5, npads-0.5, 100, 0, MAX_PAD_AMP);
      h_pad_time_pad = new TH2D("h_pad_time_pad", "pad time vs pad number", npads, -0.5, npads-0.5, 100, 0, MAX_TIME);
      if (fPH) {
         fPH->cd(2);
         h_pad_amp_pad->Draw();
      }

      h_num_aw_pad_hits = new TH2D("h_num_aw_pad_hits", "number of aw vs pad hits", 50, 0, 100, 100, 0, 50);

      h_aw_pad_hits = new TH2D("h_aw_pad_hits", "hits in aw vs hits in pads", 4*8, -0.5, 4*8-0.5, 256, -0.5, 256-0.5);

      h_aw_pad_time = new TH2D("h_aw_pad_time", "time of hits in aw vs pads", 50, 0, MAX_TIME, 50, 0, MAX_TIME);

      if (fPH) {
         fPH->cd(1);
      }

      h_aw_pad_time_drift = new TH2D("h_aw_pad_time_drift", "time of hits in aw vs pads, drift region", 50, 0, 500, 50, 0, MAX_TIME);

      h_aw_pad_amp_pc = new TH2D("h_aw_pad_amp_pc", "p.h. of hits in aw vs pads, pc region", 50, 0, MAX_PAD_AMP, 50, 0, MAX_AW_AMP);

      dir->mkdir("pulser")->cd();

      // ADC timing

      h_cal_adc05_0_1 = new TH1D("h_cal_time_adc05_0_1", "calibration pulse time, adc05 chan 1-0; time, ns", 101, -50, 50);
      h_cal_adc05_0_4 = new TH1D("h_cal_time_adc05_0_4", "calibration pulse time, adc05 chan 4-0; time, ns", 101, -50, 50);
      h_cal_adc05_0_8 = new TH1D("h_cal_time_adc05_0_8", "calibration pulse time, adc05 chan 8-0; time, ns", 101, -50, 50);
      h_cal_adc05_0_12 = new TH1D("h_cal_time_adc05_0_12", "calibration pulse time, adc05 chan 12-0; time, ns", 101, -50, 50);
      h_cal_adc05_0_16 = new TH1D("h_cal_time_adc05_0_16", "calibration pulse time, adc05 chan 16-0; time, ns", 101, -50, 50);

      h_cal_adc_05_06_chan0 = new TH1D("h_cal_time_adc_05_06_chan0", "calibration pulse time, adc06 - adc05 chan 0; time, ns", 101, -50, 50);
      h_cal_adc_05_06_chan16 = new TH1D("h_cal_time_adc_05_06_chan16", "calibration pulse time, adc06 - adc05 chan 16; time, ns", 101, -50, 50);

      h_cal_adc05_16_17 = new TH1D("h_cal_time_adc05_16_17", "calibration pulse time, adc05 chan 17-16; time, ns", 101, -50, 50);
      h_cal_adc06_16_20 = new TH1D("h_cal_time_adc06_16_20", "calibration pulse time, adc06 chan 20-16; time, ns", 101, -50, 50);

      h_cal_adc05_16_24 = new TH1D("h_cal_time_adc05_16_24", "calibration pulse time, adc05 chan 24-16; time, ns", 101, -50, 50);
      h_cal_adc05_16_32 = new TH1D("h_cal_time_adc05_16_32", "calibration pulse time, adc05 chan 32-16; time, ns", 101, -50, 50);
      h_cal_adc05_16_40 = new TH1D("h_cal_time_adc05_16_40", "calibration pulse time, adc05 chan 40-16; time, ns", 101, -50, 50);


      // PWB timing

      h_cal_time_pos00_seqsca_04_05 = new TH1D("h_cal_time_pos00_seqsca_04_05", "calibration pulse time, pos00 seqsca 05-04", 101, -50, 50);

      h_cal_time_pos01_seqsca_04_05 = new TH1D("h_cal_time_pos01_seqsca_04_05", "calibration pulse time, pos01 seqsca 05-04", 101, -50, 50);
      h_cal_time_pos01_seqsca_04_84 = new TH1D("h_cal_time_pos01_seqsca_04_84", "calibration pulse time, pos01 seqsca 84-04", 101, -50, 50);
      h_cal_time_pos01_seqsca_04_164 = new TH1D("h_cal_time_pos01_seqsca_04_164", "calibration pulse time, pos01 seqsca 164-04", 101, -50, 50);

      h_cal_time_pos02_seqsca_04_05 = new TH1D("h_cal_time_pos02_seqsca_04_05", "calibration pulse time, pos02 seqsca 05-04", 101, -50, 50);
      h_cal_time_pos02_seqsca_04_84 = new TH1D("h_cal_time_pos02_seqsca_04_84", "calibration pulse time, pos02 seqsca 84-04", 101, -50, 50);
      h_cal_time_pos02_seqsca_04_164 = new TH1D("h_cal_time_pos02_seqsca_04_164", "calibration pulse time, pos02 seqsca 164-04", 101, -50, 50);

      h_cal_time_pos_01_02_seqsca04 = new TH1D("h_cal_time_pos_01_02_seqsca04", "calibration pulse time, pos02-pos01 seqsca 04; time, ns", 101, -50, 50);
      h_cal_time_pos_01_02_seqsca05 = new TH1D("h_cal_time_pos_01_02_seqsca05", "calibration pulse time, pos02-pos01 seqsca 05; time, ns", 101, -50, 50);

      h_cal_time_pos_01_03_seqsca04 = new TH1D("h_cal_time_pos_01_03_seqsca04", "calibration pulse time, pos03-pos01 seqsca 04; time, ns", 101, -50, 50);
      h_cal_time_pos_01_03_seqsca05 = new TH1D("h_cal_time_pos_01_03_seqsca05", "calibration pulse time, pos03-pos01 seqsca 05; time, ns", 101, -50, 50);

      // timing across PWB and Alpha16

      h_cal_time_pos01_seqsca_04_minus_adc5_0 = new TH1D("h_cal_time_pos01_seqsca_04_minus_adc5_0", "calibration pulse time, PWB pos01 - adc05 chan 0; time, ns", 101, -50, 50);
      h_cal_time_pos01_seqsca_04_minus_adc5_16 = new TH1D("h_cal_time_pos01_seqsca_04_minus_adc5_16", "calibration pulse time, PWB pos01 - adc05 chan 16; time, ns", 101, -50, 50);
      h_cal_time_pos02_seqsca_04_minus_adc5_0 = new TH1D("h_cal_time_pos02_seqsca_04_minus_adc5_0", "calibration pulse time, PWB pos02 - adc05 chan 0; time, ns", 101, -50, 50);
      h_cal_time_pos02_seqsca_04_minus_adc5_16 = new TH1D("h_cal_time_pos02_seqsca_04_minus_adc5_16", "calibration pulse time, PWB pos02 - adc05 chan 16; time, ns", 101, -50, 50);
   }

   void EndRun(TARunInfo* runinfo)
   {
      printf("FinalModule::EndRun, run %d\n", runinfo->fRunNo);
      //time_t run_stop_time = runinfo->fOdb->odbReadUint32("/Runinfo/Stop time binary", 0, 0);
      //printf("ODB Run stop time: %d: %s", (int)run_stop_time, ctime(&run_stop_time));
   }

   void PauseRun(TARunInfo* runinfo)
   {
      printf("FinalModule::PauseRun, run %d\n", runinfo->fRunNo);
   }

   void ResumeRun(TARunInfo* runinfo)
   {
      printf("FinalModule::ResumeRun, run %d\n", runinfo->fRunNo);
   }

   TAFlowEvent* Analyze(TARunInfo* runinfo, TMEvent* event, TAFlags* flags, TAFlowEvent* flow)
   {
      //printf("FinalModule::Analyze, run %d, event serno %d, id 0x%04x, data size %d\n", runinfo->fRunNo, event->serial_number, (int)event->event_id, event->data_size);

      AgEventFlow *ef = flow->Find<AgEventFlow>();

      if (!ef || !ef->fEvent)
         return flow;

      AgAwHitsFlow* eawh = flow->Find<AgAwHitsFlow>();
      AgPadHitsFlow* eph = flow->Find<AgPadHitsFlow>();

      //int force_plot = false;

      AgEvent* age = ef->fEvent;

      if (age->a16 && age->feam) {
         static bool first = true;
         static double a16_time0 = 0;
         static double feam_time0 = 0;
         static uint32_t a16_ts0 = 0;
         static uint32_t feam_ts0 = 0;

         double a16_time  = age->a16->time;
         double feam_time = age->feam->time;

         uint32_t a16_ts  = 0;
         uint32_t feam_ts = 0;

         if (1)
            a16_ts = age->a16->udp[0]->eventTimestamp;

         if (age->feam->modules[0])
            feam_ts = age->feam->modules[0]->ts_start;

         if (first) {
            first = false;
            a16_time0 = a16_time;
            feam_time0 = feam_time;
            a16_ts0 = a16_ts;
            feam_ts0 = feam_ts;
         }

         double atr = a16_time-a16_time0;
         double ftr = feam_time-feam_time0;
         double dtr = atr-ftr;

         uint32_t a16_tsr = a16_ts - a16_ts0;
         uint32_t feam_tsr = feam_ts - feam_ts0;

         uint32_t a16_tsr_ns = a16_tsr * 10.0;
         uint32_t feam_tsr_ns = feam_tsr * 8.0;

         int dts_ns = (a16_tsr_ns%0x80000000) - (feam_tsr_ns%0x80000000);

         printf("Have AgEvent: %d %d, %f %f, diff %f, ts 0x%08x 0x%08x, ns: %12d %12d, diff %d\n", age->a16->counter, age->feam->counter, atr, ftr, dtr, a16_tsr, feam_tsr, a16_tsr_ns, feam_tsr_ns, dts_ns);
      }

      double adc5_0  = -1010;
      double adc5_1  = -1020;
      double adc5_4  = -1030;
      double adc5_8  = -1040;
      double adc5_12 = -1050;
      double adc5_16 = -1150;
      double adc5_17 = -1160;
      double adc5_24 = -1170;
      double adc5_32 = -1180;
      double adc5_40 = -1190;
      double adc6_0  = -2060;
      double adc6_16 = -2150;
      double adc6_20 = -2160;

      if (eawh) {
         if (1) {
            printf("AW event %d, time %f, anode wire hits: %d\n", ef->fEvent->counter, ef->fEvent->time, (int)eawh->fAwHits.size());
         }

         h_num_aw_hits->Fill(eawh->fAwHits.size());

         for (unsigned j=0; j<eawh->fAwHits.size(); j++) {
            int wire = eawh->fAwHits[j].wire;
            double time = eawh->fAwHits[j].time;
            double amp = eawh->fAwHits[j].amp;

            h_aw_time->Fill(time);
            h_aw_amp->Fill(amp);
            h_aw_amp_time->Fill(time, amp);

            h_aw_map->Fill(wire);
            h_aw_map_time->Fill(wire, time);
            h_aw_map_amp->Fill(wire, amp);

            if (eawh->fAwHits[j].module == 5) {
               if (eawh->fAwHits[j].channel == 0) adc5_0 = time;
               if (eawh->fAwHits[j].channel == 1) adc5_1 = time;
               if (eawh->fAwHits[j].channel == 4) adc5_4 = time;
               if (eawh->fAwHits[j].channel == 8) adc5_8 = time;
               if (eawh->fAwHits[j].channel == 12) adc5_12 = time;
               if (eawh->fAwHits[j].channel == 16) adc5_16 = time;
               if (eawh->fAwHits[j].channel == 17) adc5_17 = time;
               if (eawh->fAwHits[j].channel == 24) adc5_24 = time;
               if (eawh->fAwHits[j].channel == 32) adc5_32 = time;
               if (eawh->fAwHits[j].channel == 40) adc5_40 = time;
            }

            if (eawh->fAwHits[j].module == 6) {
               if (eawh->fAwHits[j].channel == 0) adc6_0 = time;
               if (eawh->fAwHits[j].channel == 16) adc6_16 = time;
               if (eawh->fAwHits[j].channel == 20) adc6_20 = time;
            }

            for (unsigned k=0; k<eawh->fAwHits.size(); k++) {
               if (k==j)
                  continue;
               h_aw_aw_hits->Fill(eawh->fAwHits[j].wire, eawh->fAwHits[k].wire);
               h_aw_aw_time->Fill(eawh->fAwHits[j].time, eawh->fAwHits[k].time);
               h_aw_aw_amp->Fill(eawh->fAwHits[j].amp, eawh->fAwHits[k].amp);
               
               if (eawh->fAwHits[j].wire == 286) h_aw_286->Fill(eawh->fAwHits[k].wire);
               if (eawh->fAwHits[j].wire == 287) h_aw_287->Fill(eawh->fAwHits[k].wire);
               if (eawh->fAwHits[j].wire == 288) h_aw_288->Fill(eawh->fAwHits[k].wire);
               if (eawh->fAwHits[j].wire == 289) h_aw_289->Fill(eawh->fAwHits[k].wire);
               if (eawh->fAwHits[j].wire == 290) h_aw_290->Fill(eawh->fAwHits[k].wire);

               if (eawh->fAwHits[j].wire == 299) h_aw_299->Fill(eawh->fAwHits[k].wire);
               if (eawh->fAwHits[j].wire == 300) h_aw_300->Fill(eawh->fAwHits[k].wire);
               if (eawh->fAwHits[j].wire == 301) h_aw_301->Fill(eawh->fAwHits[k].wire);

               if (eawh->fAwHits[j].wire == 310) h_aw_310->Fill(eawh->fAwHits[k].wire);
               if (eawh->fAwHits[j].wire == 320) h_aw_320->Fill(eawh->fAwHits[k].wire);
               if (eawh->fAwHits[j].wire == 330) h_aw_330->Fill(eawh->fAwHits[k].wire);
               if (eawh->fAwHits[j].wire == 340) h_aw_340->Fill(eawh->fAwHits[k].wire);

               if (eawh->fAwHits[j].wire == 352) h_aw_352->Fill(eawh->fAwHits[k].wire);
            }
         }
      }

      double pos00_seqsca04 = -10;
      double pos00_seqsca05 = -20;

      double pos01_seqsca04 = -30;
      double pos01_seqsca05 = -40;
      double pos01_seqsca84 = -130;
      double pos01_seqsca164 = -140;

      double pos02_seqsca04 = -50;
      double pos02_seqsca05 = -60;
      double pos02_seqsca84 = -150;
      double pos02_seqsca164 = -160;

      double pos03_seqsca04 = -70;
      double pos03_seqsca05 = -80;

      if (eph) {
         if (1) {
            printf("PA event %d, time %f, pad hits: %d\n", ef->fEvent->counter, ef->fEvent->time, (int)eph->fPadHits.size());
         }

         h_num_pad_hits->Fill(eph->fPadHits.size());

         for (unsigned i=0; i<eph->fPadHits.size(); i++) {
            int col = eph->fPadHits[i].col;
            int row = eph->fPadHits[i].row;

            if (col < 0 || row < 0)
               continue;

            int pos = eph->fPadHits[i].pos;
            int seqsca = eph->fPadHits[i].seqsca;
            int seqpad = col*MAX_FEAM_PAD_ROWS + row;

            double time = eph->fPadHits[i].time;
            double amp = eph->fPadHits[i].amp;

            //printf("pos %d, seqsca %d, time %f\n", pos, seqsca, time);

            if (pos==0) {
               if (seqsca == 4) pos00_seqsca04 = time;
               if (seqsca == 5) pos00_seqsca05 = time;
            }

            if (pos==1) {
               if (seqsca == 4) pos01_seqsca04 = time;
               if (seqsca == 5) pos01_seqsca05 = time;
               if (seqsca == 84) pos01_seqsca84 = time;
               if (seqsca == 164) pos01_seqsca164 = time;
            }

            if (pos==2) {
               if (seqsca == 4) pos02_seqsca04 = time;
               if (seqsca == 5) pos02_seqsca05 = time;
               if (seqsca == 84) pos02_seqsca84 = time;
               if (seqsca == 164) pos02_seqsca164 = time;
            }

            if (pos==3) {
               if (seqsca == 4) pos03_seqsca04 = time;
               if (seqsca == 5) pos03_seqsca05 = time;
            }

            h_pad_time->Fill(time);
            h_pad_amp->Fill(amp);
            h_pad_amp_time->Fill(time, amp);
            h_pad_amp_pad->Fill(seqpad, amp);
            h_pad_time_pad->Fill(seqpad, time);
         }
      }

      if (eawh && eph) {
         if (1) {
            printf("AA event %d, time %f, anode wire hits: %d, pad hits: %d\n", ef->fEvent->counter, ef->fEvent->time, (int)eawh->fAwHits.size(), (int)eph->fPadHits.size());
         }

         h_num_aw_pad_hits->Fill(eph->fPadHits.size(), eawh->fAwHits.size());

         for (unsigned i=0; i<eph->fPadHits.size(); i++) {
            for (unsigned j=0; j<eawh->fAwHits.size(); j++) {
               int xcol = (eph->fPadHits[i].pos%8)*4 + eph->fPadHits[i].col;
               h_aw_pad_hits->Fill(xcol, eawh->fAwHits[j].wire);
               h_aw_pad_time->Fill(eph->fPadHits[i].time, eawh->fAwHits[j].time);

               if ((eawh->fAwHits[j].time > 200) && eph->fPadHits[i].time > 200) {
                  h_aw_pad_time_drift->Fill(eph->fPadHits[i].time, eawh->fAwHits[j].time);
               }

               if ((eawh->fAwHits[j].time < 200) && eph->fPadHits[i].time < 200) {
                  h_aw_pad_amp_pc->Fill(eph->fPadHits[i].amp, eawh->fAwHits[j].amp);
               }
            }
         }
      }

      if (1 || ((pos00_seqsca04 > 0)
                && (pos00_seqsca05 > 0)
                && (pos01_seqsca04 > 0)
                && (pos01_seqsca05 > 0)
                && (pos02_seqsca04 > 0)
                && (pos02_seqsca05 > 0))) {

#if 0
         printf("ADC times: adc05: %f %f %f %f %f, adc06: %f\n",
                adc5_0,
                adc5_1,
                adc5_4,
                adc5_8,
                adc5_12,
                adc6_0
                );
#endif

         if (adc5_0 > 0) {
            h_cal_adc05_0_1->Fill(adc5_1-adc5_0);
            h_cal_adc05_0_4->Fill(adc5_4-adc5_0);
            h_cal_adc05_0_8->Fill(adc5_8-adc5_0);
            h_cal_adc05_0_12->Fill(adc5_12-adc5_0);
            h_cal_adc05_0_16->Fill(adc5_16-adc5_0);

            h_cal_adc05_16_17->Fill(adc5_17-adc5_16);
            h_cal_adc05_16_24->Fill(adc5_24-adc5_16);
            h_cal_adc05_16_32->Fill(adc5_32-adc5_16);
            h_cal_adc05_16_40->Fill(adc5_40-adc5_16);

            h_cal_adc06_16_20->Fill(adc6_20-adc6_16);

            h_cal_adc_05_06_chan0->Fill(adc6_0-adc5_0);
            h_cal_adc_05_06_chan16->Fill(adc6_16-adc5_16);
         }

         double pulse_width = 5350.0 + 30.0;
         double xpad = pos01_seqsca04;
         double t = xpad - adc5_0 - pulse_width;
         printf("aw %.1f pad %.1f %.1f, diff %.1f\n", adc5_0, pos01_seqsca04, xpad, t);

         h_cal_time_pos00_seqsca_04_05->Fill(pos00_seqsca05-pos00_seqsca04);

         h_cal_time_pos01_seqsca_04_05->Fill(pos01_seqsca05-pos01_seqsca04);
         h_cal_time_pos01_seqsca_04_84->Fill(pos01_seqsca84-pos01_seqsca04);
         h_cal_time_pos01_seqsca_04_164->Fill(pos01_seqsca164-pos01_seqsca04);

         h_cal_time_pos02_seqsca_04_05->Fill(pos02_seqsca05-pos02_seqsca04);
         h_cal_time_pos02_seqsca_04_05->Fill(pos02_seqsca05-pos02_seqsca04);
         h_cal_time_pos02_seqsca_04_84->Fill(pos02_seqsca84-pos02_seqsca04);
         h_cal_time_pos02_seqsca_04_164->Fill(pos02_seqsca164-pos02_seqsca04);

         h_cal_time_pos_01_02_seqsca04->Fill(pos02_seqsca04-pos01_seqsca04);
         h_cal_time_pos_01_02_seqsca05->Fill(pos02_seqsca05-pos01_seqsca05);

         h_cal_time_pos_01_03_seqsca04->Fill(pos03_seqsca04-pos01_seqsca04);
         h_cal_time_pos_01_03_seqsca05->Fill(pos03_seqsca05-pos01_seqsca05);

         h_cal_time_pos01_seqsca_04_minus_adc5_0->Fill(pos01_seqsca04 - adc5_0 - pulse_width);
         h_cal_time_pos01_seqsca_04_minus_adc5_16->Fill(pos01_seqsca04 - adc5_16 - pulse_width);

         h_cal_time_pos02_seqsca_04_minus_adc5_0->Fill(pos02_seqsca04 - adc5_0 - pulse_width);
         h_cal_time_pos02_seqsca_04_minus_adc5_16->Fill(pos02_seqsca04 - adc5_16 - pulse_width);
      }


      //hamp[ichan]->Fill(wamp);

      bool do_plot = (runinfo->fRoot->fgApp != NULL);

      if (do_plot) {
         if (fPH) {
            fPH->GetPad(1)->Modified();
            fPH->GetPad(2)->Modified();
            fPH->Update();
         }

         if (fC) {
            fC->Clear();
            fC->Divide(2,2);

            fC->cd(1);
            TH1D* hh = new TH1D("hh", "hh", NUM_AW, -0.5, NUM_AW-0.5);
            hh->SetMinimum(0);
            hh->SetMaximum(7000);
            hh->Draw();

            fC->cd(2);
            TH1D* ha = new TH1D("ha", "ha", NUM_AW, -0.5, NUM_AW-0.5);
            ha->SetMinimum(0);
            ha->SetMaximum(66000);
            ha->Draw();

            fC->cd(4);
            TH1D* hp = new TH1D("hp", "hp", NUM_PC, -0.5, NUM_PC-0.5);
            hp->SetMinimum(0);
            hp->SetMaximum(512);
            hp->Draw();

            std::vector<Double_t> theta;
            std::vector<Double_t> radius;
            std::vector<Double_t> etheta;
            std::vector<Double_t> eradius;

            double rmin = 0.6;
            double rmax = 1.0;

            if (0) {
               for (int i=0; i<8; i++) {
                  theta.push_back((i+1)*(TMath::Pi()/4.));
                  radius.push_back((i+1)*0.05);
                  etheta.push_back(TMath::Pi()/8.);
                  eradius.push_back(0.05);
               }
            }

            if (eawh) {
               //h_num_aw_hits->Fill(eawh->fAwHits.size());

               for (unsigned j=0; j<eawh->fAwHits.size(); j++) {
                  //h_aw_map->Fill(eawh->fAwHits[j].wire);
                  //h_aw_time->Fill(eawh->fAwHits[j].time);
                  //h_aw_amp->Fill(eawh->fAwHits[j].amp);
                  hh->SetBinContent(1+eawh->fAwHits[j].wire, eawh->fAwHits[j].time);
                  ha->SetBinContent(1+eawh->fAwHits[j].wire, eawh->fAwHits[j].amp);

                  int num_wires = 256;

                  int iwire = eawh->fAwHits[j].wire%num_wires;
                  int itb = eawh->fAwHits[j].wire/num_wires;

                  double dist = (eawh->fAwHits[j].time - 1000.0)/4000.0;
                  if (dist < 0)
                     dist = 0;
                  if (dist > 1)
                     dist = 1;

                  double t = ((iwire-8.0)/(1.0*num_wires))*(2.0*TMath::Pi());
                  double r = rmax-dist*(rmax-rmin);

                  printf("hit %d, wire %d, tb %d, iwire %d, t %f (%f), r %f\n", j, eawh->fAwHits[j].wire, itb, iwire, t, t/TMath::Pi(), r);

                  theta.push_back(t+0.5*TMath::Pi());
                  radius.push_back(r);
                  etheta.push_back(2.0*TMath::Pi()/num_wires);
                  eradius.push_back(0.05);
               }
            }

            if (eph) {
               //h_num_aw_hits->Fill(eawh->fAwHits.size());

               for (unsigned i=0; i<eph->fPadHits.size(); i++) {
                  double time = eph->fPadHits[i].time;
                  double amp = eph->fPadHits[i].amp;
                  int pos = eph->fPadHits[i].pos;
                  int col = eph->fPadHits[i].col;
                  int row = eph->fPadHits[i].row;

                  if (pos < 0)
                     continue;
                  if (col < 0)
                     continue;

                  int pc = pos*4 + col;

                  printf("pad hit %d: pos %d col %d pc %d, row %d, time %f, amp %f\n", i, pos, col, pc, row, time, amp);

                  //h_aw_map->Fill(eawh->fAwHits[j].wire);
                  //h_aw_time->Fill(eawh->fAwHits[j].time);
                  //h_aw_amp->Fill(eawh->fAwHits[j].amp);
                  hp->SetBinContent(1+pc, time);
                  //ha->SetBinContent(eawh->fAwHits[j].wire, eawh->fAwHits[j].amp);

                  //double dist = (eawh->fAwHits[j].time - 1000.0)/4000.0;
                  //if (dist < 0)
                  //   dist = 0;
                  //if (dist > 1)
                  //   dist = 1;

                  double dist = -0.1;

                  double t = ((pc+0.5)/(1.0*NUM_PC))*(2.0*TMath::Pi());
                  double r = rmax-dist*(rmax-rmin);

                  //printf("hit %d, wire %d, tb %d, iwire %d, t %f (%f), r %f\n", j, eawh->fAwHits[j].wire, itb, iwire, t, t/TMath::Pi(), r);

                  theta.push_back(t+0.5*TMath::Pi());
                  radius.push_back(r);
                  etheta.push_back(2.0*TMath::Pi()/NUM_PC/2.0);
                  eradius.push_back(0.05);
                  
                  if (1) {
                     double dist = (time - 150.0)/300.0;
                     if (dist < 0)
                        dist = 0;
                     if (dist > 1)
                        dist = 1;
                     double r = rmax-dist*(rmax-rmin);
                     
                     theta.push_back(t+0.5*TMath::Pi());
                     radius.push_back(r);
                     etheta.push_back(2.0*TMath::Pi()/NUM_PC/2.0);
                     eradius.push_back(0.05);
                  }
               }
            }

            if (1) {
               theta.push_back(0+0.5*TMath::Pi());
               radius.push_back(0);
               etheta.push_back(TMath::Pi()/8.);
               eradius.push_back(0.10);
            }

            if (1) {
               theta.push_back(0+0.5*TMath::Pi());
               radius.push_back(rmin);
               etheta.push_back(2.0*TMath::Pi());
               eradius.push_back(0.10);
            }

            if (1) {
               theta.push_back(2.5*TMath::Pi());
               radius.push_back(rmax);
               etheta.push_back(2.0*TMath::Pi());
               eradius.push_back(0.10);
            }

            fC->cd(3);
            TGraphPolar * grP1 = new TGraphPolar(theta.size(), theta.data(), radius.data(), etheta.data(), eradius.data());
            grP1->SetTitle("TGraphPolar Example");
            grP1->SetMarkerStyle(20);
            grP1->SetMarkerSize(0.75);
            grP1->SetMarkerColor(4);
            grP1->SetLineColor(2);
            grP1->SetLineWidth(3);
            grP1->SetMinPolar(0);
            grP1->SetMaxPolar(2.0*TMath::Pi());
            grP1->SetMinRadial(0);
            grP1->SetMaxRadial(1.0);
            //grP1->SetMinimum(0);
            //grP1->SetMaximum(1);
            //grP1->Draw("PRE");
            grP1->Draw("PE");
            // Update, otherwise GetPolargram returns 0
            gPad->Update();
            grP1->GetPolargram()->SetToRadian();

            //hh->Modified();
            //ha->Modified();

            fC->Modified();
            fC->Draw();
            fC->Update();
         }

#if 0
         // plot waveforms

         fC->Clear();
         fC->Divide(2,3);

         if (1) {
            fC->cd(1);
            TH1D* hh = new TH1D("hh", "hh", nbins, 0, nbins);
            for (int ibin=0; ibin<nbins; ibin++) {
               hh->SetBinContent(ibin+1, e->adcs[0]->adc[0][0][ibin]);
            }
            hh->Draw();
         }

         if (1) {
            fC->cd(2);
            TH1D* hhh = new TH1D("hhh", "hhh", nbins, 0, nbins);
            for (int ibin=0; ibin<nbins; ibin++) {
               hhh->SetBinContent(ibin+1, e->adcs[0]->adc[0][0][ibin]);
            }
            hhh->SetMinimum(-33000);
            hhh->SetMaximum(+33000);
            hhh->Draw();
         }

         if (1) {
            fC->cd(3);
            int nbins = ww[39]->nsamples;
            TH1D* hhh = new TH1D("hhhh", "hhhh", nbins, 0, nbins);
            for (int ibin=0; ibin<nbins; ibin++) {
               hhh->SetBinContent(ibin+1, ww[39]->samples[ibin]);
            }
            hhh->SetMinimum(-9000);
            hhh->SetMaximum(+9000);
            hhh->Draw();
         }

         if (1) {
            fC->cd(4);
            int nbins = ww[iplot]->nsamples;
            TH1D* hhh = new TH1D("hhhhh", "hhhhh", nbins, 0, nbins);
            for (int ibin=0; ibin<nbins; ibin++) {
               hhh->SetBinContent(ibin+1, ww[iplot]->samples[ibin]);
            }
            hhh->SetMinimum(-9000);
            hhh->SetMaximum(+9000);
            hhh->Draw();
         }

         if (1) {
            fC->cd(5);
            int nbins = ww[33]->nsamples;
            TH1D* h33 = new TH1D("h33", "h33", nbins, 0, nbins);
            for (int ibin=0; ibin<nbins; ibin++) {
               h33->SetBinContent(ibin+1, ww[33]->samples[ibin]);
            }
            h33->SetMinimum(-9000);
            h33->SetMaximum(+9000);
            h33->Draw();
         }

         if (1) {
            fC->cd(6);
            int nbins = ww[34]->nsamples;
            TH1D* h34 = new TH1D("h34", "h34", nbins, 0, nbins);
            for (int ibin=0; ibin<nbins; ibin++) {
               h34->SetBinContent(ibin+1, ww[34]->samples[ibin]);
            }
            h34->SetMinimum(-9000);
            h34->SetMaximum(+9000);
            h34->Draw();
         }

         fC->Modified();
         fC->Draw();
         fC->Update();
#endif
      }

#if 0
      if (fModule->fPlotPad >= 0) {
         if (!fModule->fPlotPadCanvas)
            fModule->fPlotPadCanvas = new TCanvas("FEAM PAD", "FEAM PAD", 900, 650);

         TCanvas*c = fModule->fPlotPadCanvas;

         c->cd();

         int nbins = ww[fModule->fPlotPad]->nsamples;
         TH1D* h = new TH1D("h", "h", nbins, 0, nbins);
         for (int ibin=0; ibin<nbins; ibin++) {
            h->SetBinContent(ibin+1, ww[fModule->fPlotPad]->samples[ibin]);
         }

         h->SetMinimum(-9000);
         h->SetMaximum(+9000);
         h->Draw();

         c->Modified();
         c->Draw();
         c->Update();
      }
#endif

#if 0
      time_t now = time(NULL);

      if (force_plot) {
         static time_t plot_next = 0;
         if (now > plot_next) {
            //fATX->PlotA16Canvas();
            plot_next = time(NULL) + 15;
         }
      }

      static time_t t = 0;

      if (now - t > 15) {
         t = now;
         //fATX->Plot();
      }

      *flags |= TAFlag_DISPLAY;
#endif

      return flow;
   }

   void AnalyzeSpecialEvent(TARunInfo* runinfo, TMEvent* event)
   {
      printf("FinalModule::AnalyzeSpecialEvent, run %d, event serno %d, id 0x%04x, data size %d\n", runinfo->fRunNo, event->serial_number, (int)event->event_id, event->data_size);
   }
};

class FinalModuleFactory: public TAFactory
{
public:
   void Init(const std::vector<std::string> &args)
   {
      printf("FinalModuleFactory::Init!\n");
      
      //fDoPads = true;
      //fPlotPad = -1;
      //fPlotPadCanvas = NULL;
      
      for (unsigned i=0; i<args.size(); i++) {
         //if (args[i] == "--nopads")
         //   fDoPads = false;
         //if (args[i] == "--plot1")
         //   fPlotPad = atoi(args[i+1].c_str());
      }
   }

   void Finish()
   {
      printf("FinalModuleFactory::Finish!\n");
   }

   TARunObject* NewRunObject(TARunInfo* runinfo)
   {
      printf("FinalModule::NewRunObject, run %d, file %s\n", runinfo->fRunNo, runinfo->fFileName.c_str());
      return new FinalModule(runinfo);
   }
};

static TARegister tar(new FinalModuleFactory);

/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */
