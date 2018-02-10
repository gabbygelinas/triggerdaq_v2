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
#include "ko_limits.h"

// histogram limit for number of hits in aw and pads
#define MAX_HITS 250

#define DELETE(x) if (x) { delete (x); (x) = NULL; }

#define MEMZERO(p) memset((p), 0, sizeof(p))

class FinalModule: public TARunObject
{
public:
   TCanvas* fC = NULL;

   TH1D* h_time_between_events = NULL;
   TH1D* h_time_between_events_zoom = NULL;

   TH1D* h_aw_num_hits;

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

   TH1D* h_adc16_bits;
   TH2D* h_adc16_bits_vs_aw;

   TH1D* h_pad_num_hits;

   TH1D* h_pad_time;
   TH1D* h_pad_amp;
   TH2D* h_pad_amp_time;

   TH2D *h_pad_amp_pad;
   TH2D *h_pad_time_pad;

   //TH2D* h_pad_pad_num_hits;

   TH2D* h_aw_pad_num_hits;
   TH2D* h_aw_pad_hits;
   TH2D* h_aw_pad_time;

   //TH2D* h_aw_pad_time_drift;
   //TH2D* h_aw_pad_amp_pc;

   // time of calibration pulse relative to ADC trigger

   TH1D* h_cal_adc05_0_full_range = NULL;
   TH1D* h_cal_adc05_16_full_range = NULL;

   TH1D* h_cal_adc05_0 = NULL;
   TH1D* h_cal_adc05_16 = NULL;
   TH1D* h_cal_adc05_17 = NULL;

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

   // time across all 32 channels pf 62.5 MHz ADC
   TProfile* h_cal_adc05_16_xx = NULL;

   // time across to the next module, 100 MHz and 62.5 MHz ADC
   TH1D* h_cal_adc_05_06_chan0 = NULL;
   TH1D* h_cal_adc_05_06_chan16 = NULL;

   // time of PWB calibration pulse
   TH1D* h_cal_time_pos00_seqsca_04_full_range = NULL;
   TH1D* h_cal_time_pos01_seqsca_04_full_range = NULL;
   TH1D* h_cal_time_pos02_seqsca_04_full_range = NULL;
   TH1D* h_cal_time_pos00_seqsca_04 = NULL;
   TH1D* h_cal_time_pos01_seqsca_04 = NULL;
   TH1D* h_cal_time_pos02_seqsca_04 = NULL;

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

   bool fTrace = false;

   FinalModule(TARunInfo* runinfo)
      : TARunObject(runinfo)
   {
      if (fTrace)
         printf("FinalModule::ctor!\n");

      fC = new TCanvas("fC", "event display", 800, 800);
   }

   ~FinalModule()
   {
      if (fTrace)
         printf("FinalModule::dtor!\n");
      DELETE(fC);
   }

   void BeginRun(TARunInfo* runinfo)
   {
      if (fTrace)
         printf("FinalModule::BeginRun, run %d, file %s\n", runinfo->fRunNo, runinfo->fFileName.c_str());
      //time_t run_start_time = runinfo->fOdb->odbReadUint32("/Runinfo/Start time binary", 0, 0);
      //printf("ODB Run start time: %d: %s", (int)run_start_time, ctime(&run_start_time));

      runinfo->fRoot->fOutputFile->cd(); // select correct ROOT directory

      TDirectory* dir = gDirectory->mkdir("final");
      dir->cd(); // select correct ROOT directory

      dir->mkdir("summary")->cd();

      h_time_between_events = new TH1D("h_time_between_events", "time between events; time, sec", 100, 0, 3.0);
      h_time_between_events_zoom = new TH1D("h_time_between_events_zoom", "time between events, early part; time, sec", 100, 0, 1.0);

      h_aw_num_hits = new TH1D("h_aw_num_hits", "number of anode wire hits", 100, 0, MAX_HITS);
      h_aw_time = new TH1D("h_aw_time", "aw hit time; time, ns", 100, 0, MAX_TIME);
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

      h_adc16_bits = new TH1D("h_adc16_bits", "FPGA adc16_coinc_dff bits; link bit 0..15", 16, -0.5, 16-0.5);
      h_adc16_bits_vs_aw = new TH2D("h_adc16_bits_vs_aw", "FPGA adc16_coinc_dff bits vs AW tpc wire number; tpc wire number; link bit 0..15", NUM_AW, -0.5, NUM_AW-0.5, 16, -0.5, 16-0.5);

      h_pad_num_hits = new TH1D("h_pad_num_hits", "number of pad hits; number of hits in pads", 100, 0, MAX_HITS);
      h_pad_time = new TH1D("h_pad_time", "pad hit time; time, ns", 100, 0, MAX_TIME);
      h_pad_amp = new TH1D("h_pad_amp", "pad hit pulse height; adc counts", 100, 0, MAX_PAD_AMP);
      h_pad_amp_time = new TH2D("h_pad_amp_time", "pad p.h vs time; time, ns; adc counts", 50, 0, MAX_TIME, 50, 0, MAX_PAD_AMP);
      int npads = MAX_FEAM*MAX_FEAM_PAD_COL*MAX_FEAM_PAD_ROWS;
      //h_pad_pad_num_hits = new TH2D("h_pad_pad_num_hits", "pad number vs pad number; pad number, col*N+row; pad number, col*N+row", npads, -0.5, npads-0.5, npads, -0.5, npads-0.5);
      h_pad_amp_pad = new TH2D("h_pad_amp_pad", "pad p.h vs pad number; pad number, col*N+row; adc counts", npads, -0.5, npads-0.5, 100, 0, MAX_PAD_AMP);
      h_pad_time_pad = new TH2D("h_pad_time_pad", "pad time vs pad number; pad number, col*N+row; time, ns", npads, -0.5, npads-0.5, 100, 0, MAX_TIME);

      h_aw_pad_num_hits = new TH2D("h_aw_pad_num_hits", "number of aw vs pad hits; number if hits in aw; number of hits in pads", 50, 0, MAX_HITS, 50, 0, MAX_HITS);
      h_aw_pad_hits = new TH2D("h_aw_pad_hits", "hits in aw vs hits in pads; tpc wire; pad column", NUM_AW, -0.5, NUM_AW-0.5, NUM_PC, -0.5, NUM_PC);
      h_aw_pad_time = new TH2D("h_aw_pad_time", "time of hits in aw vs pads; time in aw, ns; time in pads, ns", 50, 0, MAX_TIME, 50, 0, MAX_TIME);

      //h_aw_pad_time_drift = new TH2D("h_aw_pad_time_drift", "time of hits in aw vs pads, drift region", 50, 0, 500, 50, 0, MAX_TIME);

      //h_aw_pad_amp_pc = new TH2D("h_aw_pad_amp_pc", "p.h. of hits in aw vs pads, pc region", 50, 0, MAX_PAD_AMP, 50, 0, MAX_AW_AMP);

      dir->mkdir("pulser")->cd();

      // ADC timing

      h_cal_adc05_0_full_range = new TH1D("h_cal_time_adc05_0_full_range", "calibration pulse time, adc05 chan 0, full time range; time, ns", 200, 0, MAX_TIME);
      h_cal_adc05_16_full_range = new TH1D("h_cal_time_adc05_16_ful_range", "calibration pulse time, adc05 chan 16, full time range; time, ns", 200, 0, MAX_TIME);

      h_cal_adc05_0 = new TH1D("h_cal_time_adc05_0", "calibration pulse time, adc05 chan 0; time, ns", 200, AW_PULSER_TIME_100-50, AW_PULSER_TIME_100+50);
      h_cal_adc05_16 = new TH1D("h_cal_time_adc05_16", "calibration pulse time, adc05 chan 16; time, ns", 200, AW_PULSER_TIME_625-50, AW_PULSER_TIME_625+50);
      h_cal_adc05_17 = new TH1D("h_cal_time_adc05_17", "calibration pulse time, adc05 chan 17; time, ns", 200, AW_PULSER_TIME_625-50, AW_PULSER_TIME_625+50);

      h_cal_adc05_0_1 = new TH1D("h_cal_time_adc05_0_1", "calibration pulse time, adc05 chan 1-0; time, ns", 201, -50, 50);
      h_cal_adc05_0_4 = new TH1D("h_cal_time_adc05_0_4", "calibration pulse time, adc05 chan 4-0; time, ns", 201, -50, 50);
      h_cal_adc05_0_8 = new TH1D("h_cal_time_adc05_0_8", "calibration pulse time, adc05 chan 8-0; time, ns", 201, -50, 50);
      h_cal_adc05_0_12 = new TH1D("h_cal_time_adc05_0_12", "calibration pulse time, adc05 chan 12-0; time, ns", 201, -50, 50);
      h_cal_adc05_0_16 = new TH1D("h_cal_time_adc05_0_16", "calibration pulse time, adc05 chan 16-0; time, ns", 201, -50, 50);

      h_cal_adc_05_06_chan0 = new TH1D("h_cal_time_adc_05_06_chan0", "calibration pulse time, adc06 - adc05 chan 0; time, ns", 201, -50, 50);
      h_cal_adc_05_06_chan16 = new TH1D("h_cal_time_adc_05_06_chan16", "calibration pulse time, adc06 - adc05 chan 16; time, ns", 201, -50, 50);

      h_cal_adc05_16_17 = new TH1D("h_cal_time_adc05_16_17", "calibration pulse time, adc05 chan 17-16; time, ns", 201, -50, 50);
      h_cal_adc06_16_20 = new TH1D("h_cal_time_adc06_16_20", "calibration pulse time, adc06 chan 20-16; time, ns", 201, -50, 50);

      h_cal_adc05_16_24 = new TH1D("h_cal_time_adc05_16_24", "calibration pulse time, adc05 chan 24-16; time, ns", 201, -50, 50);
      h_cal_adc05_16_32 = new TH1D("h_cal_time_adc05_16_32", "calibration pulse time, adc05 chan 32-16; time, ns", 201, -50, 50);
      h_cal_adc05_16_40 = new TH1D("h_cal_time_adc05_16_40", "calibration pulse time, adc05 chan 40-16; time, ns", 201, -50, 50);

      h_cal_adc05_16_xx = new TProfile("h_cal_time_adc05_16_xx", "calibration pulse time, adc05 chan 16..47-16; adc chan; time, ns", 32, 16-0.5, 48-0.5);


      // PWB timing

      h_cal_time_pos00_seqsca_04_full_range = new TH1D("h_cal_time_pos00_seqsca_04_full_range", "calibration pulse time, pos00 seqsca 04, full range; time, ns", 200, 0, MAX_TIME);
      h_cal_time_pos01_seqsca_04_full_range = new TH1D("h_cal_time_pos01_seqsca_04_full_range", "calibration pulse time, pos01 seqsca 04, full range; time, ns", 200, 0, MAX_TIME);
      h_cal_time_pos02_seqsca_04_full_range = new TH1D("h_cal_time_pos02_seqsca_04_full_range", "calibration pulse time, pos02 seqsca 04, full range; time, ns", 200, 0, MAX_TIME);

      h_cal_time_pos00_seqsca_04 = new TH1D("h_cal_time_pos00_seqsca_04", "calibration pulse time, pos00 seqsca 04; time, ns", 200, PAD_PULSER_TIME_625-50.0, PAD_PULSER_TIME_625+50.0);
      h_cal_time_pos01_seqsca_04 = new TH1D("h_cal_time_pos01_seqsca_04", "calibration pulse time, pos01 seqsca 04; time, ns", 200, PAD_PULSER_TIME_625-50.0, PAD_PULSER_TIME_625+50.0);
      h_cal_time_pos02_seqsca_04 = new TH1D("h_cal_time_pos02_seqsca_04", "calibration pulse time, pos02 seqsca 04; time, ns", 200, PAD_PULSER_TIME_625-50.0, PAD_PULSER_TIME_625+50.0);

      h_cal_time_pos00_seqsca_04_05 = new TH1D("h_cal_time_pos00_seqsca_04_05", "calibration pulse time, pos00 seqsca 05-04; time, ns", 101, -50, 50);

      h_cal_time_pos01_seqsca_04_05 = new TH1D("h_cal_time_pos01_seqsca_04_05", "calibration pulse time, pos01 seqsca 05-04; time, ns", 101, -50, 50);
      h_cal_time_pos01_seqsca_04_84 = new TH1D("h_cal_time_pos01_seqsca_04_84", "calibration pulse time, pos01 seqsca 84-04; time, ns", 101, -50, 50);
      h_cal_time_pos01_seqsca_04_164 = new TH1D("h_cal_time_pos01_seqsca_04_164", "calibration pulse time, pos01 seqsca 164-04; time, ns", 101, -50, 50);

      h_cal_time_pos02_seqsca_04_05 = new TH1D("h_cal_time_pos02_seqsca_04_05", "calibration pulse time, pos02 seqsca 05-04; time, ns", 101, -50, 50);
      h_cal_time_pos02_seqsca_04_84 = new TH1D("h_cal_time_pos02_seqsca_04_84", "calibration pulse time, pos02 seqsca 84-04; time, ns", 101, -50, 50);
      h_cal_time_pos02_seqsca_04_164 = new TH1D("h_cal_time_pos02_seqsca_04_164", "calibration pulse time, pos02 seqsca 164-04; time, ns", 101, -50, 50);

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
      if (fTrace)
         printf("FinalModule::EndRun, run %d\n", runinfo->fRunNo);
      //time_t run_stop_time = runinfo->fOdb->odbReadUint32("/Runinfo/Stop time binary", 0, 0);
      //printf("ODB Run stop time: %d: %s", (int)run_stop_time, ctime(&run_stop_time));
   }

   void PauseRun(TARunInfo* runinfo)
   {
      if (fTrace)
         printf("FinalModule::PauseRun, run %d\n", runinfo->fRunNo);
   }

   void ResumeRun(TARunInfo* runinfo)
   {
      if (fTrace)
         printf("FinalModule::ResumeRun, run %d\n", runinfo->fRunNo);
   }

   TAFlowEvent* AnalyzeFlowEvent(TARunInfo* runinfo, TAFlags* flags, TAFlowEvent* flow)
   {
      //printf("FinalModule::Analyze, run %d, event serno %d, id 0x%04x, data size %d\n", runinfo->fRunNo, event->serial_number, (int)event->event_id, event->data_size);

      AgEventFlow *ef = flow->Find<AgEventFlow>();

      if (!ef || !ef->fEvent)
         return flow;

      AgAwHitsFlow* eawh = flow->Find<AgAwHitsFlow>();
      AgPadHitsFlow* eph = flow->Find<AgPadHitsFlow>();

      //int force_plot = false;

      AgEvent* age = ef->fEvent;

      h_time_between_events->Fill(age->timeIncr);
      h_time_between_events_zoom->Fill(age->timeIncr);

      uint32_t adc16_coinc_dff = 0;

      if (age->trig && age->trig->udpData.size() > 7) {
         adc16_coinc_dff = (age->trig->udpData[6]>>8)&0xFFFF;
      }

      int trig_counter = -1;
      int adc_counter = -1;
      int pwb_counter = -1;

      if (age->trig) {
         trig_counter = age->trig->counter;
      }

      if (age->a16) {
         adc_counter = age->a16->counter;
      }

      if (age->feam) {
         pwb_counter = age->feam->counter;
      }

      printf("Have AgEvent: %d %d %d\n", trig_counter, adc_counter, pwb_counter);

      if (adc16_coinc_dff) {
         //printf("adc16_coinc_dff: 0x%04x\n", adc16_coinc_dff);
         for (int i=0; i<16; i++) {
            if (adc16_coinc_dff & (1<<i)) {
               h_adc16_bits->Fill(i);
            }
         }
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

         h_aw_num_hits->Fill(eawh->fAwHits.size());

         for (unsigned j=0; j<eawh->fAwHits.size(); j++) {
            int adc_module = eawh->fAwHits[j].adc_module;
            int adc_chan = eawh->fAwHits[j].adc_chan;
            int wire = eawh->fAwHits[j].wire;
            double time = eawh->fAwHits[j].time;
            double amp = eawh->fAwHits[j].amp;

            h_aw_time->Fill(time);
            h_aw_amp->Fill(amp);
            h_aw_amp_time->Fill(time, amp);

            h_aw_map->Fill(wire);
            h_aw_map_time->Fill(wire, time);
            h_aw_map_amp->Fill(wire, amp);

            if (adc16_coinc_dff) {
               //printf("adc16_coinc_dff: 0x%04x: ", adc16_coinc_dff);
               for (int i=0; i<16; i++) {
                  if (adc16_coinc_dff & (1<<i)) {
                     //printf(" link%d", i);
                     //h_adc16_bits->Fill(i);
                     h_adc16_bits_vs_aw->Fill(wire, i);
                  }
               }
               //printf("\n");
            }

            if (adc_module == 5) {
               if (adc_chan ==  0) adc5_0  = time;
               if (adc_chan ==  1) adc5_1  = time;
               if (adc_chan ==  4) adc5_4  = time;
               if (adc_chan ==  8) adc5_8  = time;
               if (adc_chan == 12) adc5_12 = time;
               if (adc_chan == 16) adc5_16 = time;
               if (adc_chan == 17) adc5_17 = time;
               if (adc_chan == 24) adc5_24 = time;
               if (adc_chan == 32) adc5_32 = time;
               if (adc_chan == 40) adc5_40 = time;
            }

            if (adc_module == 6) {
               if (adc_chan ==  0) adc6_0  = time;
               if (adc_chan == 16) adc6_16 = time;
               if (adc_chan == 20) adc6_20 = time;
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

         h_pad_num_hits->Fill(eph->fPadHits.size());

         for (unsigned i=0; i<eph->fPadHits.size(); i++) {
            int col = eph->fPadHits[i].col;
            int row = eph->fPadHits[i].row;

            if (col < 0 || row < 0)
               continue;

            int pos = eph->fPadHits[i].pos;
            int seqsca = eph->fPadHits[i].seqsca;
            int xcol = (pos%8)*4 + col;
            int seqpad = xcol*MAX_FEAM_PAD_ROWS + row;

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

#if 0
            for (unsigned ii=0; ii<eph->fPadHits.size(); ii++) {
               int col = eph->fPadHits[ii].col;
               int row = eph->fPadHits[ii].row;
               
               if (col < 0 || row < 0)
                  continue;
               
               int iipos = eph->fPadHits[ii].pos;
               int iixcol = (iipos%8)*4 + col;
               //int iiseqsca = eph->fPadHits[ii].seqsca;
               int iiseqpad = iixcol*MAX_FEAM_PAD_ROWS + row;

               h_pad_pad_num_hits->Fill(seqpad, iiseqpad);
            }
#endif
         }
      }

      if (eawh && eph) {
         if (0) {
            printf("AA event %d, time %f, anode wire hits: %d, pad hits: %d\n", ef->fEvent->counter, ef->fEvent->time, (int)eawh->fAwHits.size(), (int)eph->fPadHits.size());
         }

         h_aw_pad_num_hits->Fill(eawh->fAwHits.size(), eph->fPadHits.size());

         for (unsigned i=0; i<eph->fPadHits.size(); i++) {
            for (unsigned j=0; j<eawh->fAwHits.size(); j++) {
               int xcol = (eph->fPadHits[i].pos%8)*4 + eph->fPadHits[i].col;
               h_aw_pad_hits->Fill(eawh->fAwHits[j].wire, xcol);
               h_aw_pad_time->Fill(eawh->fAwHits[j].time, eph->fPadHits[i].time);

               //if ((eawh->fAwHits[j].time > 200) && eph->fPadHits[i].time > 200) {
               //   h_aw_pad_time_drift->Fill(eph->fPadHits[i].time, eawh->fAwHits[j].time);
               //}

               //if ((eawh->fAwHits[j].time < 200) && eph->fPadHits[i].time < 200) {
               //   h_aw_pad_amp_pc->Fill(eph->fPadHits[i].amp, eawh->fAwHits[j].amp);
               //}
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
         printf("ADC times: adc05: 100MHz: %f %f %f %f %f, 62.5MHz: %f %f, adc06: %f\n",
                adc5_0,
                adc5_1,
                adc5_4,
                adc5_8,
                adc5_12,
                adc5_16,
                adc5_17,
                adc6_0
                );
#endif

         if (adc5_0 > 0) {
            h_cal_adc05_0_full_range->Fill(adc5_0);
            h_cal_adc05_16_full_range->Fill(adc5_16);

            h_cal_adc05_0->Fill(adc5_0);
            h_cal_adc05_16->Fill(adc5_16);
            h_cal_adc05_17->Fill(adc5_17);

            h_cal_adc05_0_1->Fill(adc5_1-adc5_0);
            h_cal_adc05_0_4->Fill(adc5_4-adc5_0);
            h_cal_adc05_0_8->Fill(adc5_8-adc5_0);
            h_cal_adc05_0_12->Fill(adc5_12-adc5_0);
            h_cal_adc05_0_16->Fill(adc5_16-adc5_0);

            h_cal_adc05_16_17->Fill(adc5_17-adc5_16);
            h_cal_adc05_16_24->Fill(adc5_24-adc5_16);
            h_cal_adc05_16_32->Fill(adc5_32-adc5_16);
            h_cal_adc05_16_40->Fill(adc5_40-adc5_16);

            for (unsigned j=0; j<eawh->fAwHits.size(); j++) {
               int adc_module = eawh->fAwHits[j].adc_module;
               int adc_chan = eawh->fAwHits[j].adc_chan;
               //int wire = eawh->fAwHits[j].wire;
               double time = eawh->fAwHits[j].time;
               //double amp = eawh->fAwHits[j].amp;

               if (adc_module == 5) {
                  if (adc_chan >= 16 && adc_chan < 48) {
                     //printf("AAA ADC chan %d, %f %f, %f\n", adc_chan, time, adc5_16, time-adc5_16);
                     double dt = time-adc5_16;
                     if (fabs(dt) < 50.0) {
                        h_cal_adc05_16_xx->Fill(adc_chan, dt);
                     }
                  }
               }
            }

            h_cal_adc06_16_20->Fill(adc6_20-adc6_16);

            h_cal_adc_05_06_chan0->Fill(adc6_0-adc5_0);
            h_cal_adc_05_06_chan16->Fill(adc6_16-adc5_16);
         }

         double pulse_width = 5350.0 + 30.0 + 40.0;
         //double xpad = pos01_seqsca04;
         // double t = xpad - adc5_0 - pulse_width;
         //printf("aw %.1f pad %.1f %.1f, diff %.1f\n", adc5_0, pos01_seqsca04, xpad, t);

#if 0
         printf("PAD times: %f %f %f\n", pos00_seqsca04, pos01_seqsca04, pos02_seqsca04);
#endif

         h_cal_time_pos00_seqsca_04_full_range->Fill(pos00_seqsca04);
         h_cal_time_pos01_seqsca_04_full_range->Fill(pos01_seqsca04);
         h_cal_time_pos02_seqsca_04_full_range->Fill(pos02_seqsca04);

         h_cal_time_pos00_seqsca_04->Fill(pos00_seqsca04);
         h_cal_time_pos01_seqsca_04->Fill(pos01_seqsca04);
         h_cal_time_pos02_seqsca_04->Fill(pos02_seqsca04);

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

      bool do_plot = (runinfo->fRoot->fgApp != NULL);

      if (do_plot) {
         if (fC) {
            fC->Clear();
            fC->Divide(2,2);

            TVirtualPad *p = fC->cd(1);

            p->Divide(1, 2);

            p->cd(1);
            TH1D* hh = new TH1D("hh_aw_wire_time", "AW wire hit drift time; TPC wire number; drift time, ns", NUM_AW, -0.5, NUM_AW-0.5);
            hh->SetMinimum(0);
            hh->SetMaximum(MAX_TIME);
            hh->Draw();

            p->cd(2);
            TH1D* ha = new TH1D("hh_aw_wire_amp", "AW wire hit amplitude; TPC wire number; AW ADC counts", NUM_AW, -0.5, NUM_AW-0.5);
            ha->SetMinimum(0);
            ha->SetMaximum(66000);
            ha->Draw();

            TVirtualPad *p_pad = fC->cd(2);

            p_pad->Divide(1, 2);

            std::vector<Double_t> zpad_col[NUM_PC];
            std::vector<Double_t> zpad_row[NUM_PC];
            std::vector<Double_t> zpad_time[NUM_PC];
            std::vector<Double_t> zpad_amp[NUM_PC];
               
            std::vector<Double_t> pad_col;
            std::vector<Double_t> pad_row;
            std::vector<Double_t> pad_time;
            std::vector<Double_t> pad_amp;
               
            if (1) {
               std::vector<Double_t> theta;
               std::vector<Double_t> radius;
               std::vector<Double_t> etheta;
               std::vector<Double_t> eradius;
               
               std::vector<Double_t> aw_theta;
               std::vector<Double_t> aw_radius;
               std::vector<Double_t> aw_etheta;
               std::vector<Double_t> aw_eradius;
               
               std::vector<Double_t> pads_theta;
               std::vector<Double_t> pads_radius;
               std::vector<Double_t> pads_etheta;
               std::vector<Double_t> pads_eradius;
               
               bool preamp_hits[16];
               for (unsigned j=0; j<16; j++) {
                  preamp_hits[j] = 0;
               }

               std::vector<Double_t> preamp_theta;
               std::vector<Double_t> preamp_radius;
               std::vector<Double_t> preamp_etheta;
               std::vector<Double_t> preamp_eradius;
               
               double rmin = 0.6;
               double rmax = 1.0;
               double r_aw = 1.1;
               double r_preamp = 1.15;
               double r_pads = 1.2;
               
               if (0) {
                  for (int i=0; i<8; i++) {
                     theta.push_back((i+1)*(TMath::Pi()/4.));
                     radius.push_back((i+1)*0.05);
                     etheta.push_back(TMath::Pi()/8.);
                     eradius.push_back(0.05);
                  }
               }
               
               if (eawh) {
                  for (unsigned j=0; j<eawh->fAwHits.size(); j++) {
                     hh->SetBinContent(1+eawh->fAwHits[j].wire, eawh->fAwHits[j].time);
                     ha->SetBinContent(1+eawh->fAwHits[j].wire, eawh->fAwHits[j].amp);
                     
                     int num_wires = 256;
                     
                     int iwire = eawh->fAwHits[j].wire%num_wires;
                     int itb = eawh->fAwHits[j].wire/num_wires;
                     
                     int ipreamp = iwire/16;
                     
                     double dist = (eawh->fAwHits[j].time - 1000.0)/4000.0;
                     if (dist < 0)
                        dist = 0;
                     if (dist > 1)
                        dist = 1;
                     
                     double t = ((iwire-8.0)/(1.0*num_wires))*(2.0*TMath::Pi());
                     double r = rmax-dist*(rmax-rmin);
                     
                     printf("aw hit %3d, wire %3d, tb %d, preamp %2d, iwire %3d, t %f (%f), r %f\n", j, eawh->fAwHits[j].wire, itb, ipreamp, iwire, t, t/TMath::Pi(), r);

                     preamp_hits[ipreamp] = true;
                     
                     //theta.push_back(t+0.5*TMath::Pi());
                     //radius.push_back(r);
                     //etheta.push_back(2.0*TMath::Pi()/num_wires);
                     //eradius.push_back(0.05);

                     aw_theta.push_back(t+0.5*TMath::Pi());
                     aw_radius.push_back(r_aw);
                     aw_etheta.push_back(2.0*TMath::Pi()/num_wires);
                     aw_eradius.push_back(0.05);

                     aw_theta.push_back(t+0.5*TMath::Pi());
                     aw_radius.push_back(r);
                     aw_etheta.push_back(2.0*TMath::Pi()/num_wires);
                     aw_eradius.push_back(0.05);
                  }
               }

               if (eph) {
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
                     
                     int pos_ring = pos/8;
                     int pos_col = pos%8;
                     
                     int pc = pos_col*4 + col;
                     int pr = pos_ring*72+row;
                     
                     printf("pad hit %d: pos %d col %d pc %d, row %d, time %f, amp %f\n", i, pos, col, pc, row, time, amp);
                     
                     pad_col.push_back(pc);
                     pad_row.push_back(pr);
                     pad_time.push_back(time);
                     pad_amp.push_back(amp);
                     
                     assert(pc >= 0 && pc < NUM_PC);
                     
                     zpad_col[pc].push_back(pc);
                     zpad_row[pc].push_back(pr);
                     zpad_time[pc].push_back(time);
                     zpad_amp[pc].push_back(amp);
                     
                     //double dist = -0.2;
                     
                     double t = ((pc+0.5)/(1.0*NUM_PC))*(2.0*TMath::Pi());
                     //double r = rmax-dist*(rmax-rmin);
                     
                     //printf("hit %d, wire %d, tb %d, iwire %d, t %f (%f), r %f\n", j, eawh->fAwHits[j].wire, itb, iwire, t, t/TMath::Pi(), r);
                     
                     //theta.push_back(t+0.5*TMath::Pi());
                     //radius.push_back(r);
                     //etheta.push_back(2.0*TMath::Pi()/NUM_PC/2.0);
                     //eradius.push_back(0.05);

                     pads_theta.push_back(t+0.5*TMath::Pi());
                     pads_radius.push_back(r_pads);
                     pads_etheta.push_back(2.0*TMath::Pi()/NUM_PC/2.0);
                     pads_eradius.push_back(0.05);

                     if (1) {
                        double dist = (time - 2300.0)/4000.0;
                        if (dist < 0)
                           dist = 0;
                        if (dist > 1)
                           dist = 1;
                        double r = rmax-dist*(rmax-rmin);
                        
                        //theta.push_back(t+0.5*TMath::Pi());
                        //radius.push_back(r);
                        //etheta.push_back(2.0*TMath::Pi()/NUM_PC/2.0);
                        //eradius.push_back(0.05);

                        pads_theta.push_back(t+0.5*TMath::Pi());
                        pads_radius.push_back(r);
                        pads_etheta.push_back(2.0*TMath::Pi()/NUM_PC/2.0);
                        pads_eradius.push_back(0.05);
                     }
                  }
               }
               
               if (1) {
                  printf("preamp hits:");
                  for (unsigned ipreamp=0; ipreamp<16; ipreamp++) {
                     if (preamp_hits[ipreamp]) {
                        printf(" %d", ipreamp);
                        double t = ((ipreamp)/(1.0*16))*(2.0*TMath::Pi());
                        preamp_theta.push_back(t+0.5*TMath::Pi());
                        preamp_radius.push_back(r_preamp);
                        preamp_etheta.push_back(16*2.0*TMath::Pi()/256.0/2.0);
                        preamp_eradius.push_back(0.0);
                     }
                  }
                  printf("\n");
               }

               if (1) {
                  printf("adc16_coinc_dff: 0x%04x: link hits: ", adc16_coinc_dff);
                  for (int i=0; i<16; i++) {
                     if (adc16_coinc_dff & (1<<i)) {
                        printf(" %d", i);
                     }
                  }
                  printf("\n");
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
               grP1->SetTitle("TPC end view from T side, wire 0 is at pi/2");
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

               if (pads_theta.size() > 0) {
                  TGraphPolar* gpads = new TGraphPolar(pads_theta.size(), pads_theta.data(), pads_radius.data(), pads_etheta.data(), pads_eradius.data());
                  gpads->SetMarkerStyle(20);
                  gpads->SetMarkerSize(0.75);
                  gpads->SetMarkerColor(7);
                  gpads->SetLineColor(3);
                  gpads->SetLineWidth(3);
                  gpads->Draw("PE");
               }

               if (aw_theta.size() > 0) {
                  TGraphPolar* gaw = new TGraphPolar(aw_theta.size(), aw_theta.data(), aw_radius.data(), aw_etheta.data(), aw_eradius.data());
                  gaw->SetMarkerStyle(20);
                  gaw->SetMarkerSize(0.75);
                  gaw->SetMarkerColor(4);
                  gaw->SetLineColor(2);
                  gaw->SetLineWidth(3);
                  gaw->Draw("PE");
               }

               if (preamp_theta.size() > 0) {
                  TGraphPolar* gpreamp = new TGraphPolar(preamp_theta.size(), preamp_theta.data(), preamp_radius.data(), preamp_etheta.data(), preamp_eradius.data());
                  gpreamp->SetMarkerStyle(20);
                  gpreamp->SetMarkerSize(0.75);
                  gpreamp->SetMarkerColor(4);
                  gpreamp->SetLineColor(4);
                  gpreamp->SetLineWidth(3);
                  gpreamp->Draw("PE");
               }
            }

            TVirtualPad *p_pad_row = fC->cd(4);
            p_pad_row->Divide(1, 4);

            int zpad_colour[NUM_PC];
            int zpad_side[NUM_PC];

            if (1) {
               int col = 1;
               for (int i=0; i<NUM_PC; i++) {
                  zpad_side[i] = -1;
                  if (zpad_col[i].size() == 0) {
                     zpad_colour[i] = 0;
                  } else {
                     zpad_colour[i] = col;
                     col++;
                     if (col == 5)
                        col++;
                     if (col > 8) {
                        col = 1;
                     }
                  }
               }
            }


            if (1) { // split into two groups based on a gap
               int ifirst = -1;
               for (int i=0; i<NUM_PC; i++) {
                  if (zpad_col[i].size() == 0)
                     continue;
                  ifirst = i;
                  break;
               }
               
               zpad_side[ifirst] = 1;

               int igap = 0;
               for (int j=0; j<NUM_PC; j++) {
                  int i = (ifirst+j)%NUM_PC;
                  if (zpad_col[i].size() == 0) {
                     igap++;
                     printf("ifirst %d, j %d, i %d, igap %d\n", ifirst, j, i, igap);
                     if (igap>1)
                        break;
                     continue;
                  }
                  igap=0;
                  printf("ifirst %d, j %d, i %d, igap %d\n", ifirst, j, i, igap);
                  zpad_side[i] = 1;
               }

               igap = 0;
               for (int j=0; j<NUM_PC; j++) {
                  int i = (ifirst-j)%NUM_PC;
                  if (i<0)
                     i+=NUM_PC;
                  if (zpad_col[i].size() == 0) {
                     igap++;
                     printf("ifirst %d, j %d, i %d, igap %d\n", ifirst, j, i, igap);
                     if (igap>1)
                        break;
                     continue;
                  }
                  igap=0;
                  printf("ifirst %d, j %d, i %d, igap %d\n", ifirst, j, i, igap);
                  zpad_side[i] = 1;
               }
            }

            bool split_by_max_drift = false;

            if (1) {
               int count_plus = 0;
               int count_minus = 0;

               for (int i=0; i<NUM_PC; i++) {
                  if (zpad_col[i].size() == 0)
                     continue;
                  if (zpad_side[i]>0) count_plus++;
                  if (zpad_side[i]<0) count_minus++;
               }
               printf("counts: plus %d, minus %d\n", count_plus, count_minus);
               if (count_plus>0 && count_minus>0)
                  split_by_max_drift = false;
               else
                  split_by_max_drift = true;
            }

            if (split_by_max_drift) { // split into two groups separated by biggest drift time
               double max_time = -1;
               int max_pc = -1;
               for (int i=0; i<NUM_PC; i++) {
                  if (zpad_col[i].size() == 0)
                     continue;
                  for (unsigned j=0; j<zpad_col[i].size(); j++) {
                     if (zpad_time[i][j] > max_time) {
                        printf("max_pc %d->%d, time %f->%f\n", max_pc, i, max_time, zpad_time[i][j]);
                        max_time = zpad_time[i][j];
                        max_pc = i;
                     }
                  }
               }
               if (max_pc >= 0) {
                  for (int i=0; i<NUM_PC; i++) {
                     if (i<max_pc) {
                        zpad_side[i] = 1;
                     } else {
                        zpad_side[i] = -1;
                     }
                  }
               }
            }

            if (pad_col.size() > 0) {
               p_pad->cd(1);
               TGraph* hpct = new TGraph(pad_col.size(), pad_col.data(), pad_time.data());
               hpct->SetTitle("Pads hit time per pad column; pad column; hit time, ns");
               hpct->SetMarkerStyle(20);
               hpct->SetMarkerSize(0.75);
               hpct->SetMarkerColor(4);
               hpct->GetXaxis()->SetLimits(-0.5, NUM_PC-0.5);
               hpct->SetMinimum(0);
               hpct->SetMaximum(MAX_TIME);
               //hpct->Draw("AC*");
               hpct->Draw("A*");

               for (int i=0; i<NUM_PC; i++) {
                  if (zpad_col[i].size() == 0)
                     continue;

                  TGraph* g = new TGraph(zpad_col[i].size(), zpad_col[i].data(), zpad_time[i].data());
                  g->SetMarkerStyle(20);
                  g->SetMarkerSize(0.75);
                  g->SetMarkerColor(zpad_colour[i]);
                  g->GetXaxis()->SetLimits(-0.5, NUM_PC-0.5);
                  g->SetMinimum(0);
                  g->SetMaximum(MAX_TIME);
                  //hpct->Draw("AC*");
                  g->Draw("*");
               }
            }

            if (pad_col.size() > 0) {
               p_pad->cd(2);
               TGraph* hpca = new TGraph(pad_col.size(), pad_col.data(), pad_amp.data());
               hpca->SetTitle("Pads hit amplitude per pad column; pad column; hit amplitude, ADC counts");
               hpca->SetMarkerStyle(20);
               hpca->SetMarkerSize(0.75);
               hpca->SetMarkerColor(4);
               hpca->GetXaxis()->SetLimits(-0.5, NUM_PC-0.5);
               hpca->SetMinimum(0);
               hpca->SetMaximum(MAX_PAD_AMP);
               //hpct->Draw("AC*");
               hpca->Draw("A*");

               for (int i=0; i<NUM_PC; i++) {
                  if (zpad_col[i].size() == 0)
                     continue;

                  TGraph* g = new TGraph(zpad_col[i].size(), zpad_col[i].data(), zpad_amp[i].data());
                  g->SetTitle("Pads hit amplitude per pad row; pad row; hit amplitude, ADC counts");
                  g->SetMarkerStyle(20);
                  g->SetMarkerSize(0.75);
                  g->SetMarkerColor(zpad_colour[i]);
                  g->GetXaxis()->SetLimits(-0.5, NUM_PC-0.5);
                  g->SetMinimum(0);
                  g->SetMaximum(MAX_PAD_AMP);
                  //hpct->Draw("AC*");
                  g->Draw("*");
               }
            }

            if (1) {
               p_pad_row->cd(1);
               
               //TGraph* hprt = new TGraph();
               //hprt->SetMarkerStyle(20);
               //hprt->SetMarkerSize(0.75);
               //hprt->SetMarkerColor(4);
               //hprt->GetXaxis()->SetLimits(-0.5, NUM_PR-0.5);
               //hprt->SetMinimum(0);
               //hprt->SetMaximum(MAX_PAD_AMP);
               //hprt->Draw("AC*");
               //hprt->Draw("A*");

               bool first = true;

               for (int i=0; i<NUM_PC; i++) {
                  if (zpad_col[i].size() == 0)
                     continue;
                  if (zpad_side[i] != 1)
                     continue;

                  //printf("111 col %d, size %d\n", i, zpad_row[i].size());
                  TGraph* g = new TGraph(zpad_row[i].size(), zpad_row[i].data(), zpad_amp[i].data());
                  g->SetTitle("Pads hit time per pad row; pad row; hit time, ns");
                  g->SetMarkerStyle(20);
                  g->SetMarkerSize(0.75);
                  g->SetMarkerColor(zpad_colour[i]);
                  g->GetXaxis()->SetLimits(-0.5, NUM_PR-0.5);
                  g->SetMinimum(0);
                  g->SetMaximum(MAX_PAD_AMP);
                  if (first) {
                     g->Draw("A*");
                     first = false;
                  } else {
                     g->Draw("*");
                  }
               }
            }

            if (1) {
               p_pad_row->cd(2);

               //TGraph* hprt = new TGraph(pad_row.size(), pad_row.data(), pad_time.data());
               //hprt->SetTitle("XXX; aaa; bbb");
               //hprt->SetMarkerStyle(20);
               //hprt->SetMarkerSize(0.75);
               //hprt->SetMarkerColor(4);
               //hprt->GetXaxis()->SetLimits(-0.5, NUM_PR-0.5);
               //hprt->SetMinimum(-MAX_TIME);
               //hprt->SetMaximum(0);
               //hprt->Draw("AC*");
               //hprt->Draw("A*");

               bool first = true;

               for (int i=0; i<NUM_PC; i++) {
                  if (zpad_row[i].size() == 0)
                     continue;
                  if (zpad_side[i] != 1)
                     continue;

                  for (unsigned j=0; j<zpad_row[i].size(); j++) {
                     zpad_time[i][j] = -zpad_time[i][j];
                  }

                  //printf("222 col %d, size %d\n", i, zpad_row[i].size());
                  TGraph* g = new TGraph(zpad_row[i].size(), zpad_row[i].data(), zpad_time[i].data());
                  g->SetTitle("Pads hit time per pad row; pad row; hit time, ns");
                  g->SetMarkerStyle(20);
                  g->SetMarkerSize(0.75);
                  g->SetMarkerColor(zpad_colour[i]);
                  g->GetXaxis()->SetLimits(-0.5, NUM_PR-0.5);
                  g->SetMinimum(-MAX_TIME);
                  g->SetMaximum(0);

                  if (first) {
                     g->Draw("A*");
                     first = false;
                  } else {
                     g->Draw("*");
                  }
               }
            }

            if (1) {
               p_pad_row->cd(3);
               
               //TGraph* hprt = new TGraph(pad_row.size(), pad_row.data(), pad_time.data());
               //hprt->SetTitle("XXX; aaa; bbb");
               //hprt->SetMarkerStyle(20);
               //hprt->SetMarkerSize(0.75);
               //hprt->SetMarkerColor(4);
               //hprt->GetXaxis()->SetLimits(-0.5, NUM_PR-0.5);
               //hprt->SetMinimum(0);
               //hprt->SetMaximum(MAX_TIME);
               //hprt->Draw("AC*");
               //hprt->Draw("A*");

               bool first = true;

               for (int i=0; i<NUM_PC; i++) {
                  if (zpad_row[i].size() == 0)
                     continue;
                  if (zpad_side[i] != -1)
                     continue;

                  //printf("333 col %d, size %d\n", i, zpad_row[i].size());
                  TGraph* g = new TGraph(zpad_row[i].size(), zpad_row[i].data(), zpad_time[i].data());
                  g->SetTitle("Pads hit time per pad row; pad row; hit time, ns");
                  g->SetMarkerStyle(20);
                  g->SetMarkerSize(0.75);
                  g->SetMarkerColor(zpad_colour[i]);
                  g->GetXaxis()->SetLimits(-0.5, NUM_PR-0.5);
                  g->SetMinimum(0);
                  g->SetMaximum(MAX_TIME);
                  if (first) {
                     g->Draw("A*");
                     first = false;
                  } else {
                     g->Draw("*");
                  }
               }
            }

            if (1) {
               p_pad_row->cd(4);

               //TGraph* hprt = new TGraph(pad_row.size(), pad_row.data(), pad_amp.data());
               //hprt->SetMarkerStyle(20);
               //hprt->SetMarkerSize(0.75);
               //hprt->SetMarkerColor(4);
               //hprt->GetXaxis()->SetLimits(-0.5, NUM_PR-0.5);
               //hprt->SetMinimum(0);
               //hprt->SetMaximum(MAX_PAD_AMP);
               //hprt->Draw("AC*");
               //hprt->Draw("A*");

               bool first = true;

               for (int i=0; i<NUM_PC; i++) {
                  if (zpad_row[i].size() == 0)
                     continue;
                  if (zpad_side[i] != -1)
                     continue;

                  //printf("444 col %d, size %d\n", i, zpad_row[i].size());
                  TGraph* g = new TGraph(zpad_row[i].size(), zpad_row[i].data(), zpad_amp[i].data());
                  g->SetTitle("Pads hit amplitude per pad row; pad row; hit amplitude, ADC counts");
                  g->SetMarkerStyle(20);
                  g->SetMarkerSize(0.75);
                  g->SetMarkerColor(zpad_colour[i]);
                  g->GetXaxis()->SetLimits(-0.5, NUM_PR-0.5);
                  g->SetMinimum(0);
                  g->SetMaximum(MAX_PAD_AMP);
                  if (first) {
                     g->Draw("A*");
                     first = false;
                  } else {
                     g->Draw("*");
                  }
               }
            }

            fC->Modified();
            fC->Draw();
            fC->Update();
         }

         if (age->trig && age->trig->udpData.size() > 0) {
            for (unsigned i=0; i<age->trig->udpData.size(); i++) {
               printf("ATAT[%d]: 0x%08x (%d)\n", i, age->trig->udpData[i], age->trig->udpData[i]);
            }
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

      //*flags |= TAFlag_DISPLAY;
#endif

      return flow;
   }

   void AnalyzeSpecialEvent(TARunInfo* runinfo, TMEvent* event)
   {
      if (fTrace)
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
