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

#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
//#include "TMath.h"

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
   bool fPrint = false;

public:
   TH1D* h_time_between_events = NULL;
   TH1D* h_time_between_events_zoom_1sec = NULL;
   TH1D* h_time_between_events_zoom_01sec = NULL;
   TH1D* h_time_between_events_zoom_001sec = NULL;
   TH1D* h_time_between_events_zoom_0001sec = NULL;
   TH1D* h_time_between_events_zoom_00001sec = NULL;

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

   TH1D* h_bsc64_bits;
   TH2D* h_bsc64_vs_bsc_adc_hits;

   TH1D* h_aw_num_hits;

   TH1D* h_aw_time;
   TH1D* h_aw_amp;
   TH2D* h_aw_amp_time;

#if 0
   TH1D* h_aw_time_preamp_even_pc;
   TH1D* h_aw_amp_preamp_even_pc;

   TH1D* h_aw_time_preamp_odd_pc;
   TH1D* h_aw_amp_preamp_odd_pc;
#endif

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

   TH2D* h_aw_vs_bsc_adc_hits;
   TH2D* h_aw_pc_vs_bsc_adc_hits;

   TH2D* h_aw16_prompt_vs_bsc64 = NULL;
   TH2D* h_aw16_prompt_vs_bsc_adc_hits = NULL;

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

#if 0
   TH1D* h_pad_time_pos[64];
#endif

   TH2D* h_aw_pad_num_hits;
   TH2D* h_aw_pad_hits;
   TH2D* h_aw_pad_time;

   TH2D* h_pad_bsc_adc_hits = NULL;

   //TH2D* h_aw_pad_time_drift;
   //TH2D* h_aw_pad_amp_pc;

   bool fTrace = false;

   FinalModule(TARunInfo* runinfo, bool do_print)
      : TARunObject(runinfo)
   {
      if (fTrace)
         printf("FinalModule::ctor!\n");

      fModuleName = "final_module";

      fPrint = do_print;
   }

   ~FinalModule()
   {
      if (fTrace)
         printf("FinalModule::dtor!\n");
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
      h_time_between_events_zoom_001sec = new TH1D("h_time_between_events_zoom_001sec", "time between events, zoom 10 msec; time, sec", 100, 0, 0.01);
      h_time_between_events_zoom_0001sec = new TH1D("h_time_between_events_zoom_0001sec", "time between events, zoom 1 msec; time, sec", 100, 0, 0.001);
      h_time_between_events_zoom_00001sec = new TH1D("h_time_between_events_zoom_00001sec", "time between events, zoom 100 usec; time, sec", 100, 0, 0.0001);

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

      h_bsc64_bits = new TH1D("h_bsc64_bits", "FPGA bsc64 bits; bsc64 bit 0..63", 64+1, -0.5, 64-0.5+1);
      h_bsc64_vs_bsc_adc_hits = new TH2D("h_bsc64_vs_bsc_adc_hits", "FPGA bsc64 bits vs BSC ADC hit bar number; bsc64 bit 0..63; bar number", 64, -0.5, 64-0.5, NUM_BSC, -0.5, NUM_BSC-0.5);

      h_aw_num_hits = new TH1D("h_aw_num_hits", "number of anode wire hits", 100, 0, MAX_HITS);
      h_aw_time = new TH1D("h_aw_time", "aw hit time; time, ns", 100, 0, MAX_TIME);
      h_aw_amp = new TH1D("h_aw_amp", "aw hit pulse height", 100, 0, MAX_AW_AMP);
      h_aw_amp_time = new TH2D("h_aw_amp_time", "aw p.h. vs time", 100, 0, MAX_TIME, 50, 0, MAX_AW_AMP);

#if 0
      h_aw_time_preamp_even_pc = new TH1D("h_aw_time_preamp_even_pc", "aw hit time, even preamp, PC hits; time, ns", 100, 0, MAX_TIME);
      h_aw_amp_preamp_even_pc = new TH1D("h_aw_amp_preamp_even_pc", "aw hit pulse height, even preamp, PC hits", 100, 0, MAX_AW_AMP);

      h_aw_time_preamp_odd_pc = new TH1D("h_aw_time_preamp_odd_pc", "aw hit time, odd preamp, PC hits; time, ns", 100, 0, MAX_TIME);
      h_aw_amp_preamp_odd_pc = new TH1D("h_aw_amp_preamp_odd_pc", "aw hit pulse height, odd preamp, PC hits", 100, 0, MAX_AW_AMP);
#endif

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

      h_aw_vs_bsc_adc_hits = new TH2D("h_aw_vs_bsc_adc_hits", "hits in aw vs bsc adc; wire number; bar number", NUM_AW, -0.5, NUM_AW-0.5, NUM_BSC, -0.5, NUM_BSC-0.5);
      h_aw_pc_vs_bsc_adc_hits = new TH2D("h_aw_pc_vs_bsc_adc_hits", "hits in aw (PC region) vs bsc adc; wire number; bar number", NUM_AW, -0.5, NUM_AW-0.5, NUM_BSC, -0.5, NUM_BSC-0.5);

      h_aw16_prompt_vs_bsc64 = new TH2D("h_aw16_prompt_vs_bsc64", "FPGA aw16_prompt hits vs bsc64 hits; bsc64 bit; aw16_prompt bit", 64, -0.5, 64-0.5, 16, -0.5, 16-0.5);

      h_aw16_prompt_vs_bsc_adc_hits = new TH2D("h_aw16_prompt_vs_bsc_adc_hits", "FPGA aw16_prompt hits vs bsc adc hits; aw16_prompt bit; bsc bar number", 16, -0.5, 16-0.5, NUM_BSC, -0.5, NUM_BSC-0.5);

      h_pad_num_hits = new TH1D("h_pad_num_hits", "number of pad hits; number of hits in pads", 100, 0, MAX_HITS);
      h_pad_time = new TH1D("h_pad_time", "pad hit time; time, ns", 100, 0, MAX_TIME);
      h_pad_amp = new TH1D("h_pad_amp", "pad hit pulse height; adc counts", 100, 0, MAX_PAD_AMP);
      h_pad_amp_time = new TH2D("h_pad_amp_time", "pad p.h vs time; time, ns; adc counts", 50, 0, MAX_TIME, 50, 0, MAX_PAD_AMP);

      h_pad_hits_per_column = new TH1D("h_pad_hits_per_column", "pad hits per TPC column; tpc column", 32, -0.5, 32-0.5);
      h_pad_hits_per_row = new TH1D("h_pad_hits_per_row", "pad hits per TPC row; tpc row", 8*4*18, -0.5, 8*4*18-0.5);

      h_pad_hits_per_row_column = new TH2D("h_pad_hits_per_row_column", "pad hits per TPC row and column; tpc row; tpc column", 8*4*18, -0.5, 8*4*18-0.5, 32, -0.5, 32-0.5);

      h_pad_time_per_column = new TH2D("h_pad_time_per_column", "pad hit time per column; tpc column; time, ns", 32, -0.5, 32-0.5, 100, 0, MAX_TIME);
      h_pad_time_per_row = new TH2D("h_pad_time_per_row", "pad hit time per row; tpc row; time, ns", 8*4*18, -0.5, 8*4*18-0.5, 100, 0, MAX_TIME);

#if 0
      for (int icol=0; icol<8; icol++) {
         for (int iring=0; iring<8; iring++) {
            char name[256];
            char title[256];
            sprintf(name, "h_pad_time_pwb_c%dr%d", icol, iring);
            sprintf(title, "pad hit time pwb col %d ring %d; time, ns", icol, iring);
            h_pad_time_pos[icol*8+iring] = new TH1D(name, title, 100, 0, MAX_TIME);
         }
      }
#endif

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
      h_time_between_events_zoom_001sec->Fill(age->timeIncr);
      h_time_between_events_zoom_0001sec->Fill(age->timeIncr);
      h_time_between_events_zoom_00001sec->Fill(age->timeIncr);

      //uint32_t adc16_coinc_dff = 0;
      uint32_t aw16_prompt = 0;
      int aw16_mlu_out = 0;
      uint32_t trig_bitmap  = 0;
      uint32_t aw16_bus = 0;
      uint64_t bsc64_bus = 0;
      int bsc64_mult = 0;
      uint32_t coinc_latch = 0;

      //if (age->trig && age->trig->udpData.size() > 7) {
      //   adc16_coinc_dff = (age->trig->udpData[6]>>8)&0xFFFF;
      //}

      if (age->trig && age->trig->udpData.size() > 9) {
         aw16_prompt = (age->trig->udpData[9])&0xFFFF;
         aw16_mlu_out = ((age->trig->udpData[9])&0x80000000) != 0;
      }

      if (age->trig && age->trig->udpData.size() > 7) {
         trig_bitmap = age->trig->udpData[6];
      }

      if (age->trig && age->trig->udpData.size() > 13) {
         aw16_bus = age->trig->udpData[13] & 0xFFFF;
      }

      if (age->trig && age->trig->udpData.size() > 17) {
         uint64_t bsc64_bus_lo = age->trig->udpData[14];
         uint64_t bsc64_bus_hi = age->trig->udpData[15];
         bsc64_bus = (bsc64_bus_hi << 32) | bsc64_bus_lo;
         bsc64_mult = age->trig->udpData[16] & 0xFF;
      }

      if (age->trig && age->trig->udpData.size() > 17) {
         coinc_latch = age->trig->udpData[17] & 0xFF;
      }

      //printf("trig_bitmap 0x%08x, aw16_prompt 0x%04x, aw16_mlu_out %d, aw16_bus 0x%04x, bsc64_bus 0x%016lx, bsc64_mult %2d, coinc_latch 0x%02x\n", trig_bitmap, aw16_prompt, aw16_mlu_out, aw16_bus, bsc64_bus, bsc64_mult, coinc_latch);

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

      if (fPrint) {
         printf("Have AgEvent: %d %d %d %d\n", trig_counter, adc_counter, pwb_counter, tdc_counter);
      }

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

      if (bsc64_bus) {
         for (uint64_t i=0; i<64; i++) {
            if (bsc64_bus & (((uint64_t)1)<<i)) {
               h_bsc64_bits->Fill(i);
            }
         }
      }

      if (aw16_prompt && bsc64_bus) {
         for (int i=0; i<16; i++) {
            if (aw16_prompt & (1<<i)) {
               for (uint64_t j=0; j<64; j++) {
                  if (bsc64_bus & (((uint64_t)1)<<j)) {
                     h_aw16_prompt_vs_bsc64->Fill(j, i);
                  }
               }
            }
         }
      }

      if (eba) {
         if (fPrint) {
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

            if (bsc64_bus) {
               for (uint64_t i=0; i<64; i++) {
                  if (bsc64_bus & (((uint64_t)1)<<i)) {
                     h_bsc64_vs_bsc_adc_hits->Fill(i, bar);
                  }
               }
            }

            if (aw16_prompt) {
               //printf("aw16_prompt: 0x%04x\n", aw16_prompt);
               for (int i=0; i<16; i++) {
                  if (aw16_prompt & (1<<i)) {
                     h_aw16_prompt_vs_bsc_adc_hits->Fill(i, bar);
                  }
               }
            }

            //if (bar == 64+63) {
            //   printf("bar %3d, aw16_prompt 0x%04x, bsc64_bus 0x%016lx bit 0x%016lx %d\n", bar, aw16_prompt, bsc64_bus, bsc64_bus&(((uint64_t)1)<<bar), (bsc64_bus&(((uint64_t)1)<<bar))!=0);
            //}

            //int xbit = 32;
            //if (bsc64_bus&(((uint64_t)1)<<xbit)) {
            //   if (amp < 8000)
            //      continue;
            //   printf("bsc64 bit %2d, bar %3d, time %6.0f, amp %6.0f, aw16_prompt 0x%04x, bsc64_bus 0x%016lx bit 0x%016lx %d\n", xbit, bar, time, amp, aw16_prompt, bsc64_bus, bsc64_bus&(((uint64_t)1)<<bar), (bsc64_bus&(((uint64_t)1)<<bar))!=0);
            //
            //}
         }
      }

      if (eba && eawh) {
         for (unsigned j=0; j<eba->fBscAdcHits.size(); j++) {
            for (unsigned k=0; k<eawh->fAwHits.size(); k++) {
               double time = eawh->fAwHits[k].time;
               h_aw_vs_bsc_adc_hits->Fill(eawh->fAwHits[k].wire, eba->fBscAdcHits[j].bar);
               if (time >= 800 && time < 1200) {
                  h_aw_pc_vs_bsc_adc_hits->Fill(eawh->fAwHits[k].wire, eba->fBscAdcHits[j].bar);
               }
            }
         }
      }

      if (eawh) {
         if (fPrint) {
            printf("AW event %d, time %f, anode wire hits: %d\n", ef->fEvent->counter, ef->fEvent->time, (int)eawh->fAwHits.size());
         }

         h_aw_num_hits->Fill(eawh->fAwHits.size());

         for (unsigned j=0; j<eawh->fAwHits.size(); j++) {
            //int adc_module = eawh->fAwHits[j].adc_module;
            //int adc_chan = eawh->fAwHits[j].adc_chan;
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
#if 0
               if (preamp == 16+2 /*preamp%2 == 0*/) {
                  h_aw_time_preamp_even_pc->Fill(time);
                  h_aw_amp_preamp_even_pc->Fill(amp);
               } else if (preamp == 16+3) {
                  h_aw_time_preamp_odd_pc->Fill(time);
                  h_aw_amp_preamp_odd_pc->Fill(amp);
               }
#endif
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

      if (eph) {
         if (fPrint) {
            printf("PA event %d, time %f, pad hits: %d\n", ef->fEvent->counter, ef->fEvent->time, (int)eph->fPadHits.size());
         }

         h_pad_num_hits->Fill(eph->fPadHits.size());

         for (unsigned i=0; i<eph->fPadHits.size(); i++) {
            //int ipwb = eph->fPadHits[i].imodule;
            //int seqsca = eph->fPadHits[i].seqsca;
            //int xcol = (pos%8)*4 + col;
            //int seqpad = xcol*MAX_FEAM_PAD_ROWS + row;
            int col = eph->fPadHits[i].tpc_col;
            int row = eph->fPadHits[i].tpc_row;

            double time = eph->fPadHits[i].time_ns;
            double amp = eph->fPadHits[i].amp;

            //printf("pos %d, seqsca %d, time %f\n", pos, seqsca, time);

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

#if 0
            h_pad_time_pos[pwb_seqpos]->Fill(time);
#endif

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
   bool fDoPrint = false;

public:
   void Usage()
   {
      printf("FinalModuleFactory::Usage:\n");
      printf("--print ## be verbose\n");
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

         if (args[i] == "--print")
            fDoPrint = true;
      }
   }

   void Finish()
   {
      printf("FinalModuleFactory::Finish!\n");
   }

   TARunObject* NewRunObject(TARunInfo* runinfo)
   {
      printf("FinalModule::NewRunObject, run %d, file %s\n", runinfo->fRunNo, runinfo->fFileName.c_str());
      return new FinalModule(runinfo, fDoPrint);
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
