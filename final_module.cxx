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

#define PLOT_MIN_TIME (900.0)
#define PLOT_MAX_TIME (5100.0)

#define NUM_PREAMP 32

#define DELETE(x) if (x) { delete (x); (x) = NULL; }

#define MEMZERO(p) memset((p), 0, sizeof(p))

class FinalModule: public TARunObject
{
public:
   TCanvas* fC = NULL;

   TH1D* h_time_between_events = NULL;
   TH1D* h_time_between_events_zoom_1sec = NULL;
   TH1D* h_time_between_events_zoom_01sec = NULL;

   TH1D* h_bsc_adc_num_hits;

   TH1D* h_bsc_adc_time;
   TH1D* h_bsc_adc_amp;
   TH2D* h_bsc_adc_amp_time;

   TH1D* h_bsc_adc_map;
   TH2D* h_bsc_adc_map_time;
   TH2D* h_bsc_adc_map_amp;

   TH2D* h_bsc_bsc_adc_hits;
   TH2D* h_bsc_bsc_adc_time;
   TH2D* h_bsc_bsc_adc_amp;

   TH1D* h_aw_num_hits;

   TH1D* h_aw_time;
   TH1D* h_aw_amp;
   TH2D* h_aw_amp_time;

   TH1D* h_aw_time_preamp_even_pc;
   TH1D* h_aw_amp_preamp_even_pc;

   TH1D* h_aw_time_preamp_odd_pc;
   TH1D* h_aw_amp_preamp_odd_pc;

   TH1D* h_aw_map;
   TH2D* h_aw_map_time;
   TH2D* h_aw_map_amp;

   TH1D* h_preamp_map;
   TH1D* h_preamp_map_pc;
   TH1D* h_preamp_map_dc;
   TH2D* h_preamp_map_time;
   TH2D* h_preamp_map_amp;
   TH2D* h_preamp_map_amp_pc;
   TProfile* h_preamp_map_amp_prof;
   TProfile* h_preamp_map_amp_prof_pc;

   TH1D* h_aw_map_early;
   TH1D* h_aw_map_pc;
   TH1D* h_aw_map_pc_8000;
   TH1D* h_aw_map_dc;
   TH1D* h_aw_map_late;

   TH2D* h_aw_aw_hits;
   TH2D* h_aw_aw_time;
   TH2D* h_aw_aw_amp;

#if 0
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
#endif

   //TH1D* h_adc16_bits;
   //TH2D* h_adc16_bits_vs_aw;

   TH1D* h_aw16_prompt_bits;
   TH2D* h_aw16_prompt_bits_vs_aw;

   TH2D* h_aw_bsc_adc_hits;

   TH1D* h_pad_num_hits;

   TH1D* h_pad_time;
   TH1D* h_pad_amp;
   TH2D* h_pad_amp_time;

   TH2D *h_pad_amp_pad;
   TH2D *h_pad_time_pad;

   //TH2D* h_pad_pad_num_hits;

   TH1D* h_pad_hits_per_column = NULL;
   TH1D* h_pad_hits_per_row = NULL;
   TH2D* h_pad_hits_per_row_column = NULL;

   TH2D* h_pad_time_per_column = NULL;
   TH2D* h_pad_time_per_row = NULL;

   TH1D* h_pad_time_pos[64];

   TH2D* h_aw_pad_num_hits;
   TH2D* h_aw_pad_hits;
   TH2D* h_aw_pad_time;

   TH2D* h_pad_bsc_adc_hits = NULL;

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
   bool fDoEventDisplay = false;

   FinalModule(TARunInfo* runinfo, bool do_event_display)
      : TARunObject(runinfo)
   {
      if (fTrace)
         printf("FinalModule::ctor!\n");

      fDoEventDisplay = do_event_display;

      if (fDoEventDisplay) {
         fC = new TCanvas("fC", "event display", 800, 800);
      }
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
      h_time_between_events_zoom_1sec = new TH1D("h_time_between_events_zoom_1sec", "time between events, zoom 1 sec; time, sec", 100, 0, 1.0);
      h_time_between_events_zoom_01sec = new TH1D("h_time_between_events_zoom_01sec", "time between events, zoom 0.1 sec; time, sec", 100, 0, 0.1);

      h_bsc_adc_num_hits = new TH1D("h_bsc_adc_num_hits", "BSC ADC number of hits", 20, 0-0.5, 20-0.5);
      h_bsc_adc_time = new TH1D("h_bsc_adc_time", "BSC ADC hit time; time, ns", 100, 0, MAX_TIME);
      h_bsc_adc_amp = new TH1D("h_bsc_adc_amp", "BSC ADC hit pulse height", 100, 0, MAX_AW_AMP);
      h_bsc_adc_amp_time = new TH2D("h_bsc_adc_amp_time", "BSC ADC p.h. vs time", 100, 0, MAX_TIME, 50, 0, MAX_AW_AMP);

      h_bsc_adc_map = new TH1D("h_bsc_adc_map", "BSC ADC bar occupancy; bar number", NUM_BSC, -0.5, NUM_BSC-0.5);
      h_bsc_adc_map_time = new TH2D("h_bsc_adc_map_time", "BSC ADC hit time vs bar; bar number; hit time, ns", NUM_BSC, -0.5, NUM_BSC-0.5, 50, 0, MAX_TIME);
      h_bsc_adc_map_amp  = new TH2D("h_bsc_adc_map_amp", "BSC ADC hit p.h. vs bar; bar number; hit p.h. adc units", NUM_BSC, -0.5, NUM_BSC-0.5, 50, 0, MAX_AW_AMP);

      h_bsc_bsc_adc_hits = new TH2D("h_bsc_bsc_adc_hits", "hits in bsc adc vs bsc adc", NUM_BSC, -0.5, NUM_BSC-0.5, NUM_BSC, -0.5, NUM_BSC-0.5);
      h_bsc_bsc_adc_time = new TH2D("h_bsc_bsc_adc_time", "time in bsc adc vs bsc adc", 50, 0, MAX_TIME, 50, 0, MAX_TIME);
      h_bsc_bsc_adc_amp  = new TH2D("h_bsc_bsc_adc_amp",  "p.h. in bsc adc vs bsc adc", 50, 0, MAX_AW_AMP, 50, 0, MAX_AW_AMP);

      h_aw_num_hits = new TH1D("h_aw_num_hits", "number of anode wire hits", 100, 0, MAX_HITS);
      h_aw_time = new TH1D("h_aw_time", "aw hit time; time, ns", 100, 0, MAX_TIME);
      h_aw_amp = new TH1D("h_aw_amp", "aw hit pulse height", 100, 0, MAX_AW_AMP);
      h_aw_amp_time = new TH2D("h_aw_amp_time", "aw p.h. vs time", 100, 0, MAX_TIME, 50, 0, MAX_AW_AMP);

      h_aw_time_preamp_even_pc = new TH1D("h_aw_time_preamp_even_pc", "aw hit time, even preamp, PC hits; time, ns", 100, 0, MAX_TIME);
      h_aw_amp_preamp_even_pc = new TH1D("h_aw_amp_preamp_even_pc", "aw hit pulse height, even preamp, PC hits", 100, 0, MAX_AW_AMP);

      h_aw_time_preamp_odd_pc = new TH1D("h_aw_time_preamp_odd_pc", "aw hit time, odd preamp, PC hits; time, ns", 100, 0, MAX_TIME);
      h_aw_amp_preamp_odd_pc = new TH1D("h_aw_amp_preamp_odd_pc", "aw hit pulse height, odd preamp, PC hits", 100, 0, MAX_AW_AMP);

      h_aw_map = new TH1D("h_aw_map", "aw hit occupancy", NUM_AW, -0.5, NUM_AW-0.5);
      h_aw_map_early = new TH1D("h_aw_map_early", "aw hit occupancy, early hits", NUM_AW, -0.5, NUM_AW-0.5);
      h_aw_map_pc = new TH1D("h_aw_map_pc", "aw hit occupancy, PC hits", NUM_AW, -0.5, NUM_AW-0.5);
      h_aw_map_pc_8000 = new TH1D("h_aw_map_pc_8000", "aw hit occupancy, PC hits with p.h. > 8000", NUM_AW, -0.5, NUM_AW-0.5);
      h_aw_map_dc = new TH1D("h_aw_map_dc", "aw hit occupancy, DC hits", NUM_AW, -0.5, NUM_AW-0.5);
      h_aw_map_late = new TH1D("h_aw_map_late", "aw hit occupancy, late hits", NUM_AW, -0.5, NUM_AW-0.5);

      h_aw_map_time = new TH2D("h_aw_map_time", "aw hit time vs wire", NUM_AW, -0.5, NUM_AW-0.5, 50, 0, MAX_TIME);
      h_aw_map_amp  = new TH2D("h_aw_map_amp", "aw hit p.h. vs wire", NUM_AW, -0.5, NUM_AW-0.5, 50, 0, MAX_AW_AMP);

      h_preamp_map = new TH1D("h_preamp_map", "preamp hit occupancy", NUM_PREAMP, -0.5, NUM_PREAMP-0.5);
      h_preamp_map_pc = new TH1D("h_preamp_map_pc", "preamp hit occupancy, PC hits", NUM_PREAMP, -0.5, NUM_PREAMP-0.5);
      h_preamp_map_dc = new TH1D("h_preamp_map_dc", "preamp hit occupancy, DC hits", NUM_PREAMP, -0.5, NUM_PREAMP-0.5);
      h_preamp_map_time = new TH2D("h_preamp_map_time", "aw hit time vs preamp", NUM_PREAMP, -0.5, NUM_PREAMP-0.5, 50, 0, MAX_TIME);
      h_preamp_map_amp  = new TH2D("h_preamp_map_amp", "aw hit p.h. vs preamp", NUM_PREAMP, -0.5, NUM_PREAMP-0.5, 50, 0, MAX_AW_AMP);
      h_preamp_map_amp_pc  = new TH2D("h_preamp_map_amp_pc", "aw hit p.h. vs preamp, PC hits", NUM_PREAMP, -0.5, NUM_PREAMP-0.5, 50, 0, MAX_AW_AMP);
      h_preamp_map_amp_prof = new TProfile("h_preamp_map_amp_prof", "aw hit p.h. vs preamp profile", NUM_PREAMP, -0.5, NUM_PREAMP-0.5);
      h_preamp_map_amp_prof_pc = new TProfile("h_preamp_map_amp_prof_pc", "aw hit p.h. vs preamp profile, PC hits", NUM_PREAMP, -0.5, NUM_PREAMP-0.5);

      h_aw_aw_hits = new TH2D("h_aw_aw_hits", "hits in aw vs aw", NUM_AW, -0.5, NUM_AW-0.5, NUM_AW, -0.5, NUM_AW-0.5);
      h_aw_aw_time = new TH2D("h_aw_aw_time", "time in aw vs aw", 50, 0, MAX_TIME, 50, 0, MAX_TIME);
      h_aw_aw_amp  = new TH2D("h_aw_aw_amp",  "p.h. in aw vs aw", 50, 0, MAX_AW_AMP, 50, 0, MAX_AW_AMP);

#if 0
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
#endif

      //h_adc16_bits = new TH1D("h_adc16_bits", "FPGA adc16_coinc_dff bits; link bit 0..15", 16+1, -0.5, 16-0.5+1);
      //h_adc16_bits_vs_aw = new TH2D("h_adc16_bits_vs_aw", "FPGA adc16_coinc_dff bits vs AW tpc wire number; tpc wire number; link bit 0..15", NUM_AW, -0.5, NUM_AW-0.5, 16, -0.5, 16-0.5);

      h_aw16_prompt_bits = new TH1D("h_aw16_prompt_bits", "FPGA aw16_prompt bits; link bit 0..15", 16+1, -0.5, 16-0.5+1);
      h_aw16_prompt_bits_vs_aw = new TH2D("h_aw16_prompt_bits_vs_aw", "FPGA aw16_prompt bits vs AW tpc wire number; tpc wire number; link bit 0..15", NUM_AW, -0.5, NUM_AW-0.5, 16, -0.5, 16-0.5);

      h_aw_bsc_adc_hits = new TH2D("h_aw_bsc_adc_hits", "hits in aw vs bsc adc", NUM_AW, -0.5, NUM_AW-0.5, NUM_BSC, -0.5, NUM_BSC-0.5);

      h_pad_num_hits = new TH1D("h_pad_num_hits", "number of pad hits; number of hits in pads", 100, 0, MAX_HITS);
      h_pad_time = new TH1D("h_pad_time", "pad hit time; time, ns", 100, 0, MAX_TIME);
      h_pad_amp = new TH1D("h_pad_amp", "pad hit pulse height; adc counts", 100, 0, MAX_PAD_AMP);
      h_pad_amp_time = new TH2D("h_pad_amp_time", "pad p.h vs time; time, ns; adc counts", 50, 0, MAX_TIME, 50, 0, MAX_PAD_AMP);

      h_pad_hits_per_column = new TH1D("h_pad_hits_per_column", "pad hits per TPC column; tpc column", 32, -0.5, 32-0.5);
      h_pad_hits_per_row = new TH1D("h_pad_hits_per_row", "pad hits per TPC row; tpc row", 8*4*18, -0.5, 8*4*18-0.5);

      h_pad_hits_per_row_column = new TH2D("h_pad_hits_per_row_column", "pad hits per TPC row and column; tpc row; tpc column", 8*4*18, -0.5, 8*4*18-0.5, 32, -0.5, 32-0.5);

      h_pad_time_per_column = new TH2D("h_pad_time_per_column", "pad hit time per column; tpc column; time, ns", 32, -0.5, 32-0.5, 100, 0, MAX_TIME);
      h_pad_time_per_row = new TH2D("h_pad_time_per_row", "pad hit time per row; tpc row; time, ns", 8*4*18, -0.5, 8*4*18-0.5, 100, 0, MAX_TIME);

      for (int i=0; i<64; i++) {
         char name[256];
         char title[256];
         sprintf(name, "h_pad_time_pwb_seqpos_%02d", i);
         sprintf(title, "pad hit time pwb seqpos %d; time, ns", i);
         h_pad_time_pos[i] = new TH1D(name, title, 100, 0, MAX_TIME);
      }

      //int npads = MAX_FEAM*MAX_FEAM_PAD_COL*MAX_FEAM_PAD_ROWS;
      //h_pad_pad_num_hits = new TH2D("h_pad_pad_num_hits", "pad number vs pad number; pad number, col*N+row; pad number, col*N+row", npads, -0.5, npads-0.5, npads, -0.5, npads-0.5);
      //h_pad_amp_pad = new TH2D("h_pad_amp_pad", "pad p.h vs pad number; pad number, col*N+row; adc counts", npads, -0.5, npads-0.5, 100, 0, MAX_PAD_AMP);
      //h_pad_time_pad = new TH2D("h_pad_time_pad", "pad time vs pad number; pad number, col*N+row; time, ns", npads, -0.5, npads-0.5, 100, 0, MAX_TIME);

      h_aw_pad_num_hits = new TH2D("h_aw_pad_num_hits", "number of aw vs pad hits; number if hits in aw; number of hits in pads", 50, 0, MAX_HITS, 50, 0, MAX_HITS);
      h_aw_pad_hits = new TH2D("h_aw_pad_hits", "hits in aw vs hits in pads; tpc wire; pad column", NUM_AW, -0.5, NUM_AW-0.5, NUM_PC, -0.5, NUM_PC);
      h_aw_pad_time = new TH2D("h_aw_pad_time", "time of hits in aw vs pads; time in aw, ns; time in pads, ns", 50, 0, MAX_TIME, 50, 0, MAX_TIME);

      h_pad_bsc_adc_hits = new TH2D("h_pad_bsc_adc_hits", "hits in pads vs hits in bsc adc; pad column; bsc bar", NUM_PC, -0.5, NUM_PC, NUM_BSC, -0.5, NUM_BSC);

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
      AgBscAdcHitsFlow* eba = flow->Find<AgBscAdcHitsFlow>();

      //int force_plot = false;

      AgEvent* age = ef->fEvent;

      h_time_between_events->Fill(age->timeIncr);
      h_time_between_events_zoom_1sec->Fill(age->timeIncr);
      h_time_between_events_zoom_01sec->Fill(age->timeIncr);

      //uint32_t adc16_coinc_dff = 0;
      uint32_t aw16_prompt = 0;
      uint32_t trig_bitmap  = 0;

      //if (age->trig && age->trig->udpData.size() > 7) {
      //   adc16_coinc_dff = (age->trig->udpData[6]>>8)&0xFFFF;
      //}

      if (age->trig && age->trig->udpData.size() > 9) {
         aw16_prompt = (age->trig->udpData[9])&0xFFFF;
      }

      if (age->trig && age->trig->udpData.size() > 7) {
         trig_bitmap = age->trig->udpData[6];
      }

      int trig_counter = -1;
      int adc_counter = -1;
      int pwb_counter = -1;
      int tdc_counter = -1;

      if (age->trig) {
         trig_counter = age->trig->counter;
      }

      if (age->a16) {
         adc_counter = age->a16->counter;
      }

      if (age->feam) {
         pwb_counter = age->feam->counter;
      }

      if (age->tdc) {
         tdc_counter = age->tdc->counter;
      }

      printf("Have AgEvent: %d %d %d %d\n", trig_counter, adc_counter, pwb_counter, tdc_counter);

      //if (adc16_coinc_dff) {
      //   //printf("adc16_coinc_dff: 0x%04x\n", adc16_coinc_dff);
      //   for (int i=0; i<16; i++) {
      //      if (adc16_coinc_dff & (1<<i)) {
      //         h_adc16_bits->Fill(i);
      //      }
      //   }
      //}

      if (aw16_prompt) {
         //printf("aw16_prompt: 0x%04x\n", aw16_prompt);
         for (int i=0; i<16; i++) {
            if (aw16_prompt & (1<<i)) {
               h_aw16_prompt_bits->Fill(i);
            }
         }
      }

      if (eba) {
         if (1) {
            printf("BA event %d, time %f, bsc adc hits: %d\n", ef->fEvent->counter, ef->fEvent->time, (int)eba->fBscAdcHits.size());
         }

         h_bsc_adc_num_hits->Fill(eba->fBscAdcHits.size());

         for (unsigned j=0; j<eba->fBscAdcHits.size(); j++) {
            //int adc_module = eba->fBscAdcHits[j].adc_module;
            //int adc_chan = eba->fBscAdcHits[j].adc_chan;
            //int preamp = eawh->fAwHits[j].preamp_pos;
            int bar = eba->fBscAdcHits[j].bar;
            double time = eba->fBscAdcHits[j].time;
            double amp = eba->fBscAdcHits[j].amp;

            h_bsc_adc_time->Fill(time);
            h_bsc_adc_amp->Fill(amp);
            h_bsc_adc_amp_time->Fill(time, amp);

            h_bsc_adc_map->Fill(bar);
            h_bsc_adc_map_time->Fill(bar, time);
            h_bsc_adc_map_amp->Fill(bar, amp);

            for (unsigned k=0; k<eba->fBscAdcHits.size(); k++) {
               if (k==j)
                  continue;
               h_bsc_bsc_adc_hits->Fill(eba->fBscAdcHits[j].bar, eba->fBscAdcHits[k].bar);
               h_bsc_bsc_adc_time->Fill(eba->fBscAdcHits[j].time, eba->fBscAdcHits[k].time);
               h_bsc_bsc_adc_amp->Fill(eba->fBscAdcHits[j].amp, eba->fBscAdcHits[k].amp);
            }
         }
      }

      if (eba && eawh) {
         for (unsigned j=0; j<eba->fBscAdcHits.size(); j++) {
            for (unsigned k=0; k<eawh->fAwHits.size(); k++) {
               h_aw_bsc_adc_hits->Fill(eawh->fAwHits[k].wire, eba->fBscAdcHits[j].bar);
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
            //int preamp = eawh->fAwHits[j].preamp_pos;
            int wire = eawh->fAwHits[j].wire;
            int preamp = wire/16;
            double time = eawh->fAwHits[j].time;
            double amp = eawh->fAwHits[j].amp;

            h_aw_time->Fill(time);
            h_aw_amp->Fill(amp);
            h_aw_amp_time->Fill(time, amp);

            h_aw_map->Fill(wire);
            h_aw_map_time->Fill(wire, time);
            h_aw_map_amp->Fill(wire, amp);

            h_preamp_map->Fill(preamp);
            h_preamp_map_time->Fill(preamp, time);
            h_preamp_map_amp->Fill(preamp, amp);
            if (amp < 10000) {
               h_preamp_map_amp_prof->Fill(preamp, amp);
            }

            bool aw_early = false;
            bool aw_pc = false;
            bool aw_dc = false;
            bool aw_late = false;

            if (time < 800) {
               aw_early = true;
               h_aw_map_early->Fill(wire);
            } else if (time < 1200) {
               aw_pc = true;
               h_aw_map_pc->Fill(wire);
               if (amp > 8000) {
                  h_aw_map_pc_8000->Fill(wire);
               }
               h_preamp_map_pc->Fill(preamp);
               h_preamp_map_amp_pc->Fill(preamp, amp);
               if (amp < 10000) {
                  h_preamp_map_amp_prof_pc->Fill(preamp, amp);
               }
               if (preamp == 16+2 /*preamp%2 == 0*/) {
                  h_aw_time_preamp_even_pc->Fill(time);
                  h_aw_amp_preamp_even_pc->Fill(amp);
               } else if (preamp == 16+3) {
                  h_aw_time_preamp_odd_pc->Fill(time);
                  h_aw_amp_preamp_odd_pc->Fill(amp);
               }
            } else if (time < 5000) {
               aw_dc = true;
               h_aw_map_dc->Fill(wire);
               h_preamp_map_dc->Fill(preamp);
            } else {
               aw_late = true;
               h_aw_map_late->Fill(wire);
            }

            if (aw16_prompt) {
               for (int i=0; i<16; i++) {
                  if (aw16_prompt & (1<<i)) {
                     if (aw_pc) {
                        h_aw16_prompt_bits_vs_aw->Fill(wire, i);
                     }
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

#if 0               
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
#endif
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
            int ipwb = eph->fPadHits[i].imodule;
            int seqsca = eph->fPadHits[i].seqsca;
            //int xcol = (pos%8)*4 + col;
            //int seqpad = xcol*MAX_FEAM_PAD_ROWS + row;
            int col = eph->fPadHits[i].tpc_col;
            int row = eph->fPadHits[i].tpc_row;

            double time = eph->fPadHits[i].time_ns;
            double amp = eph->fPadHits[i].amp;

            //printf("pos %d, seqsca %d, time %f\n", pos, seqsca, time);

            if (ipwb==0) {
               if (seqsca == 4) pos00_seqsca04 = time;
               if (seqsca == 5) pos00_seqsca05 = time;
            }

            if (ipwb==1) {
               if (seqsca == 4) pos01_seqsca04 = time;
               if (seqsca == 5) pos01_seqsca05 = time;
               if (seqsca == 84) pos01_seqsca84 = time;
               if (seqsca == 164) pos01_seqsca164 = time;
            }

            if (ipwb==2) {
               if (seqsca == 4) pos02_seqsca04 = time;
               if (seqsca == 5) pos02_seqsca05 = time;
               if (seqsca == 84) pos02_seqsca84 = time;
               if (seqsca == 164) pos02_seqsca164 = time;
            }

            if (ipwb==3) {
               if (seqsca == 4) pos03_seqsca04 = time;
               if (seqsca == 5) pos03_seqsca05 = time;
            }

            h_pad_time->Fill(time);
            h_pad_amp->Fill(amp);
            h_pad_amp_time->Fill(time, amp);
            //h_pad_amp_pad->Fill(seqpad, amp);
            //h_pad_time_pad->Fill(seqpad, time);

            h_pad_hits_per_column->Fill(col);
            h_pad_hits_per_row->Fill(row);

            h_pad_hits_per_row_column->Fill(row, col);

            h_pad_time_per_column->Fill(col, time);
            h_pad_time_per_row->Fill(row, time);

            int pwb_col = col/4;
            int pwb_ring = row/(4*18);
            int pwb_seqpos = pwb_col*8 + pwb_ring;
            assert(pwb_seqpos >= 0 && pwb_seqpos < 64);

            h_pad_time_pos[pwb_seqpos]->Fill(time);

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
               int col = eph->fPadHits[i].tpc_col;
               h_aw_pad_hits->Fill(eawh->fAwHits[j].wire, col);
               h_aw_pad_time->Fill(eawh->fAwHits[j].time, eph->fPadHits[i].time_ns);

               //if ((eawh->fAwHits[j].time > 200) && eph->fPadHits[i].time > 200) {
               //   h_aw_pad_time_drift->Fill(eph->fPadHits[i].time, eawh->fAwHits[j].time);
               //}

               //if ((eawh->fAwHits[j].time < 200) && eph->fPadHits[i].time < 200) {
               //   h_aw_pad_amp_pc->Fill(eph->fPadHits[i].amp, eawh->fAwHits[j].amp);
               //}
            }
         }
      }

      if (eba && eph) {
         for (unsigned j=0; j<eba->fBscAdcHits.size(); j++) {
            for (unsigned k=0; k<eph->fPadHits.size(); k++) {
               h_pad_bsc_adc_hits->Fill(eph->fPadHits[k].tpc_col, eba->fBscAdcHits[j].bar);
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

      if (fDoEventDisplay && do_plot) {
         if (fC) {
            TVirtualPad* save_gpad = gPad;

            
            fC->Clear();
            fC->Divide(2,2);

            TVirtualPad *p = fC->cd(1);

            p->Divide(1, 2);

            TDirectory* dir = runinfo->fRoot->fgDir;
            dir->cd();

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
                     int iwire = eawh->fAwHits[j].wire;
                     double time = eawh->fAwHits[j].time;
                     double amp = eawh->fAwHits[j].amp;

                     hh->SetBinContent(1+iwire, time);
                     ha->SetBinContent(1+iwire, amp);
                     
                     int num_wires = 256;
                     
                     int itbwire = iwire%num_wires;
                     int itb = iwire/num_wires;
                     
                     int ipreamp = itbwire/16;
                     
                     double dist = (time - 1000.0)/4000.0;
                     if (dist < 0)
                        dist = 0;
                     if (dist > 1)
                        dist = 1;
                     
                     double t = (itbwire-8.0)/(1.0*num_wires)*(2.0*TMath::Pi());
                     double r = rmax-dist*(rmax-rmin);
                     
                     printf("aw hit %3d, wire %3d, tb %d, preamp %2d, iwire %3d, time %f, amp %f, theta %f (%f), radius %f\n", j, iwire, itb, ipreamp, itbwire, time, amp, t, t/TMath::Pi(), r);

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
                     int imodule = eph->fPadHits[i].imodule;
                     double time = eph->fPadHits[i].time_ns;
                     double amp = eph->fPadHits[i].amp;
                     int col = eph->fPadHits[i].tpc_col;
                     int row = eph->fPadHits[i].tpc_row;

                     if (eph->fPadHits.size() < 5000) {
                        printf("pad hit %d: pwb%02d, col %d, row %d, time %f, amp %f\n", i, imodule, col, row, time, amp);
                     }

                     if (col < 0 || row < 0) {
                        continue;
                     }
                     
                     pad_col.push_back(col);
                     pad_row.push_back(row);
                     pad_time.push_back(time);
                     pad_amp.push_back(amp);
                     
                     assert(col >= 0 && col < NUM_PC);
                     
                     zpad_col[col].push_back(col);
                     zpad_row[col].push_back(row);
                     zpad_time[col].push_back(time);
                     zpad_amp[col].push_back(amp);
                     
                     //double dist = -0.2;
                     
                     double t = ((col+0.5)/(1.0*NUM_PC))*(2.0*TMath::Pi());
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
                        double dist = (time - 1000.0)/4000.0;
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
                  printf("aw16_prompt_bits: 0x%04x: link hits: ", aw16_prompt);
                  for (int i=0; i<16; i++) {
                     if (aw16_prompt & (1<<i)) {
                        printf(" %d", i);
                     }
                  }
                  printf("\n");
               }

               if (1) {
                  printf("(out of date!) trig_bitmap: 0x%08x: bits: ", trig_bitmap);
                  if (trig_bitmap & (1<<0)) printf("adc16_grand_or ");
                  if (trig_bitmap & (1<<1)) printf("adc32_grand_or ");
                  if (trig_bitmap & (1<<2)) printf("adc_grand_or ");
                  if (trig_bitmap & (1<<3)) printf("esata_nim_grand_or ");
                  if (trig_bitmap & (1<<7)) printf("MLU ");
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
                     //printf("ifirst %d, j %d, i %d, igap %d\n", ifirst, j, i, igap);
                     if (igap>1)
                        break;
                     continue;
                  }
                  igap=0;
                  //printf("ifirst %d, j %d, i %d, igap %d\n", ifirst, j, i, igap);
                  zpad_side[i] = 1;
               }

               igap = 0;
               for (int j=0; j<NUM_PC; j++) {
                  int i = (ifirst-j)%NUM_PC;
                  if (i<0)
                     i+=NUM_PC;
                  if (zpad_col[i].size() == 0) {
                     igap++;
                     //printf("ifirst %d, j %d, i %d, igap %d\n", ifirst, j, i, igap);
                     if (igap>1)
                        break;
                     continue;
                  }
                  igap=0;
                  //printf("ifirst %d, j %d, i %d, igap %d\n", ifirst, j, i, igap);
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
                        //printf("max_pc %d->%d, time %f->%f\n", max_pc, i, max_time, zpad_time[i][j]);
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
                  g->SetMinimum(-PLOT_MAX_TIME);
                  g->SetMaximum(-PLOT_MIN_TIME);

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
                  g->SetMinimum(PLOT_MIN_TIME);
                  g->SetMaximum(PLOT_MAX_TIME);
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
            
            save_gpad->cd();
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
   bool fDoEventDisplay = true;

public:
   void Usage()
   {
      printf("FinalModuleFactory::Usage:\n");
      printf("--nodisplay ## disable the event display\n");
   }

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

         if (args[i] == "--nodisplay")
            fDoEventDisplay = false;
      }
   }

   void Finish()
   {
      printf("FinalModuleFactory::Finish!\n");
   }

   TARunObject* NewRunObject(TARunInfo* runinfo)
   {
      printf("FinalModule::NewRunObject, run %d, file %s\n", runinfo->fRunNo, runinfo->fFileName.c_str());
      return new FinalModule(runinfo, fDoEventDisplay);
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
