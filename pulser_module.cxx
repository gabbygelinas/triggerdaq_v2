//
// pulser_module.cxx
//
// pulser analysis of TPC data
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

class PulserModule: public TARunObject
{
public:
   bool fPrint = false;

public:
   // time of calibration pulse relative to ADC trigger

   TDirectory* hdir_pulser = NULL;

   bool  f_h_cal_adcxx_00_full_range = false;
   TH1D* h_cal_adcxx_00_full_range[ADC_MODULE_LAST+1];
   TH1D* h_cal_adcxx_01_full_range[ADC_MODULE_LAST+1];
   TH1D* h_cal_adcxx_04_full_range[ADC_MODULE_LAST+1];
   TH1D* h_cal_adcxx_08_full_range[ADC_MODULE_LAST+1];
   TH1D* h_cal_adcxx_12_full_range[ADC_MODULE_LAST+1];
   TH1D* h_cal_adcxx_xx_00[ADC_MODULE_LAST+1];

   TH1D* h_cal_adcxx_16_full_range[ADC_MODULE_LAST+1];
   TH1D* h_cal_adcxx_17_full_range[ADC_MODULE_LAST+1];
   TH1D* h_cal_adcxx_32_full_range[ADC_MODULE_LAST+1];
   TH1D* h_cal_adcxx_xx_16[ADC_MODULE_LAST+1];

   TH1D* h_cal_adcxx_00_16[ADC_MODULE_LAST+1];

   TH1D* h_cal_adcxx_17_16_zoom[ADC_MODULE_LAST+1];
   TH1D* h_cal_adcxx_32_16_zoom[ADC_MODULE_LAST+1];

   TProfile* h_cal_adcxx_profile_00[ADC_MODULE_LAST+1];
   TProfile* h_cal_adcxx_profile_16[ADC_MODULE_LAST+1];

   int h_first_adc_0 = 0;
   int h_first_adc_16 = 0;

   TH1D* h_cal_adcnn_00_all = NULL;
   TH1D* h_cal_adcnn_00[ADC_MODULE_LAST+1];
   TProfile* h_cal_adcnn_profile_00 = NULL;

   TH1D* h_cal_adcnn_16_all = NULL;
   TH1D* h_cal_adcnn_16[ADC_MODULE_LAST+1];
   TProfile* h_cal_adcnn_profile_16 = NULL;

   // time of TDC

   TH1D* h_cal_time_tdc_1 = NULL;
   TH1D* h_cal_time_tdc_2 = NULL;
   TH1D* h_cal_time_tdc_3 = NULL;
   TH1D* h_cal_time_tdc_2_1 = NULL;
   TH1D* h_cal_time_tdc_3_1 = NULL;
   TH1D* h_cal_time_tdc_2_1_zoom = NULL;
   TH1D* h_cal_time_tdc_3_1_zoom = NULL;

   TH1D* h_cal_time_tdc_fine = NULL;
   TH1D* h_cal_time_tdc_fine_0 = NULL;
   TH1D* h_cal_time_tdc_fine_1 = NULL;
   TH1D* h_cal_time_tdc_fine_2 = NULL;
   TH1D* h_cal_time_tdc_fine_3 = NULL;

   // time of PWB calibration pulse
   TH1D* h_cal_time_pwbaa_seqsca_04_full_range = NULL;
   TH1D* h_cal_time_pwbbb_seqsca_04_full_range = NULL;
   TH1D* h_cal_time_pwbcc_seqsca_04_full_range = NULL;
   TH1D* h_cal_time_pwbdd_seqsca_04_full_range = NULL;

   //TH1D* h_cal_time_pwbxx_seqsca_04_full_range[PWB_MODULE_LAST+1];

   TH1D* h_cal_time_pwbaa_seqsca_04 = NULL;
   TH1D* h_cal_time_pwbbb_seqsca_04 = NULL;
   TH1D* h_cal_time_pwbcc_seqsca_04 = NULL;
   TH1D* h_cal_time_pwbdd_seqsca_04 = NULL;

   TH1D* h_cal_time_pwbxx[PWB_MODULE_LAST+1];

   //TProfile* h_cal_time_pwbxx_seqsca_04_profile = NULL;

   // time across to next channel
   TH1D* h_cal_time_pwbaa_seqsca_04_05 = NULL;
   TH1D* h_cal_time_pwbbb_seqsca_04_05 = NULL;
   TH1D* h_cal_time_pwbcc_seqsca_04_05 = NULL;
   TH1D* h_cal_time_pwbdd_seqsca_04_05 = NULL;

   // time across to next SCA
   TH1D* h_cal_time_pwbaa_seqsca_04_84 = NULL;
   TH1D* h_cal_time_pwbbb_seqsca_04_84 = NULL;
   TH1D* h_cal_time_pwbcc_seqsca_04_84 = NULL;
   TH1D* h_cal_time_pwbdd_seqsca_04_84 = NULL;

   // time across to next ADC
   TH1D* h_cal_time_pwbaa_seqsca_04_164 = NULL;
   TH1D* h_cal_time_pwbbb_seqsca_04_164 = NULL;
   TH1D* h_cal_time_pwbcc_seqsca_04_164 = NULL;
   TH1D* h_cal_time_pwbdd_seqsca_04_164 = NULL;

   // time across 2 PWB boards
   TH1D* h_cal_time_pos_00_01_seqsca04 = NULL;
   TH1D* h_cal_time_pos_00_02_seqsca04 = NULL;
   TH1D* h_cal_time_pos_00_03_seqsca04 = NULL;
   TH1D* h_cal_time_pos_01_02_seqsca04 = NULL;
   TH1D* h_cal_time_pos_01_03_seqsca04 = NULL;
   TH1D* h_cal_time_pos_02_03_seqsca04 = NULL;

   // time across 100MHz ADC to PWB
   TH1D* h_cal_time_pwbaa_seqsca_04_minus_adc_0 = NULL;
   TH1D* h_cal_time_pwbbb_seqsca_04_minus_adc_0 = NULL;
   TH1D* h_cal_time_pwbcc_seqsca_04_minus_adc_0 = NULL;
   TH1D* h_cal_time_pwbdd_seqsca_04_minus_adc_0 = NULL;

   // time across 62.5MHz ADC to PWB
   TH1D* h_cal_time_pwbaa_seqsca_04_minus_adc_16 = NULL;
   TH1D* h_cal_time_pwbbb_seqsca_04_minus_adc_16 = NULL;
   TH1D* h_cal_time_pwbcc_seqsca_04_minus_adc_16 = NULL;
   TH1D* h_cal_time_pwbdd_seqsca_04_minus_adc_16 = NULL;

   bool fTrace = false;
   bool fPulser = false;
   
   int fPwbAA = 12; // pwb12 col0 row0 "M"
   int fPwbBB = 13; // pwb13 col0 row1 "M"
   int fPwbCC = 20; // pwb20 col1 row0 "STC"
   int fPwbDD = 21; // pwb21 col1 row1 "STC"

   PulserModule(TARunInfo* runinfo, bool do_print, bool pulser, int pwbaa, int pwbbb, int pwbcc, int pwbdd)
      : TARunObject(runinfo)
   {
      if (fTrace)
         printf("PulserModule::ctor!\n");

      fPrint = do_print;
      fPulser = pulser;

      if (pwbaa >= 0)
         fPwbAA = pwbaa;

      if (pwbbb >= 0)
         fPwbBB = pwbbb;

      if (pwbcc >= 0)
         fPwbCC = pwbcc;

      if (pwbdd >= 0)
         fPwbDD = pwbdd;
   }

   ~PulserModule()
   {
      if (fTrace)
         printf("PulserModule::dtor!\n");
   }

   void BeginRun(TARunInfo* runinfo)
   {
      if (fTrace)
         printf("PulserModule::BeginRun, run %d, file %s\n", runinfo->fRunNo, runinfo->fFileName.c_str());
      //time_t run_start_time = runinfo->fOdb->odbReadUint32("/Runinfo/Start time binary", 0, 0);
      //printf("ODB Run start time: %d: %s", (int)run_start_time, ctime(&run_start_time));

      bool pulser = fPulser;
      runinfo->fOdb->RB("Equipment/Ctrl/Settings/FwPulserEnable", &pulser);
      if (pulser)
         fPulser = true;
      printf("pulser mode: enabled: %d\n", fPulser);

      runinfo->fRoot->fOutputFile->cd(); // select correct ROOT directory

      hdir_pulser = gDirectory->mkdir("pulser");
      hdir_pulser->cd(); // select correct ROOT directory

      // ADC timing

      // PWB timing

      // TDC timing

      h_cal_time_tdc_1 = new TH1D("h_cal_time_tdc_1", "TDC time chan 1 - chan 0; time, ns", 1000, 100, 200);
      h_cal_time_tdc_2 = new TH1D("h_cal_time_tdc_2", "TDC time chan 2 - chan 0; time, ns", 1000, 100, 200);
      h_cal_time_tdc_3 = new TH1D("h_cal_time_tdc_3", "TDC time chan 3 - chan 0; time, ns", 1000, 100, 200);
      h_cal_time_tdc_2_1 = new TH1D("h_cal_time_tdc_2_1", "TDC time chan 2 vs chan 1; time, ns", 1000, -50, 50);
      h_cal_time_tdc_3_1 = new TH1D("h_cal_time_tdc_3_1", "TDC time chan 3 vs chan 1; time, ns", 1000, -50, 50);
      h_cal_time_tdc_2_1_zoom = new TH1D("h_cal_time_tdc_2_1_zoom", "TDC time chan 2 vs chan 1 zoom in; time, ns", 1000, -1, 1);
      h_cal_time_tdc_3_1_zoom = new TH1D("h_cal_time_tdc_3_1_zoom", "TDC time chan 3 vs chan 1 zoom in; time, ns", 1000, -1, 1);

      h_cal_time_tdc_fine = new TH1D("h_cal_time_tdc_fine", "TDC time fine time; fine time units", 500, 0, 500);
      h_cal_time_tdc_fine_0 = new TH1D("h_cal_time_tdc_fine_0", "TDC time fine time chan 0; fine time units", 50, 400, 450);
      h_cal_time_tdc_fine_1 = new TH1D("h_cal_time_tdc_fine_1", "TDC time fine time chan 1; fine time units", 500, 0, 500);
      h_cal_time_tdc_fine_2 = new TH1D("h_cal_time_tdc_fine_2", "TDC time fine time chan 2; fine time units", 500, 0, 500);
      h_cal_time_tdc_fine_3 = new TH1D("h_cal_time_tdc_fine_3", "TDC time fine time chan 3; fine time units", 500, 0, 500);
   }

   void EndRun(TARunInfo* runinfo)
   {
      if (fTrace)
         printf("PulserModule::EndRun, run %d\n", runinfo->fRunNo);
      //time_t run_stop_time = runinfo->fOdb->odbReadUint32("/Runinfo/Stop time binary", 0, 0);
      //printf("ODB Run stop time: %d: %s", (int)run_stop_time, ctime(&run_stop_time));
   }

   void PauseRun(TARunInfo* runinfo)
   {
      if (fTrace)
         printf("PulserModule::PauseRun, run %d\n", runinfo->fRunNo);
   }

   void ResumeRun(TARunInfo* runinfo)
   {
      if (fTrace)
         printf("PulserModule::ResumeRun, run %d\n", runinfo->fRunNo);
   }

   void CreatePwbHist(int pwbaa, int pwbbb, int pwbcc, int pwbdd)
   {
      hdir_pulser->cd(); // select correct ROOT directory

      h_cal_time_pwbaa_seqsca_04_full_range = new TH1D("h_cal_time_pwbaa_seqsca_04_full_range", "calibration pulse time, pwbaa seqsca 04, full range; time, ns", 200, 0, MAX_TIME);
      h_cal_time_pwbbb_seqsca_04_full_range = new TH1D("h_cal_time_pwbbb_seqsca_04_full_range", "calibration pulse time, pwbbb seqsca 04, full range; time, ns", 200, 0, MAX_TIME);
      h_cal_time_pwbcc_seqsca_04_full_range = new TH1D("h_cal_time_pwbcc_seqsca_04_full_range", "calibration pulse time, pwbcc seqsca 04, full range; time, ns", 200, 0, MAX_TIME);
      h_cal_time_pwbdd_seqsca_04_full_range = new TH1D("h_cal_time_pwbdd_seqsca_04_full_range", "calibration pulse time, pwbdd seqsca 04, full range; time, ns", 200, 0, MAX_TIME);

      h_cal_time_pwbaa_seqsca_04 = new TH1D("h_cal_time_pwbaa_seqsca_04", "calibration pulse time, pwbaa seqsca 04; time, ns", 200, PAD_PULSER_TIME_625-50.0, PAD_PULSER_TIME_625+50.0);
      h_cal_time_pwbbb_seqsca_04 = new TH1D("h_cal_time_pwbbb_seqsca_04", "calibration pulse time, pwbbb seqsca 04; time, ns", 200, PAD_PULSER_TIME_625-50.0, PAD_PULSER_TIME_625+50.0);
      h_cal_time_pwbcc_seqsca_04 = new TH1D("h_cal_time_pwbcc_seqsca_04", "calibration pulse time, pwbcc seqsca 04; time, ns", 200, PAD_PULSER_TIME_625-50.0, PAD_PULSER_TIME_625+50.0);
      h_cal_time_pwbdd_seqsca_04 = new TH1D("h_cal_time_pwbdd_seqsca_04", "calibration pulse time, pwbdd seqsca 04; time, ns", 200, PAD_PULSER_TIME_625-50.0, PAD_PULSER_TIME_625+50.0);

      for (int ipwb=0; ipwb<=PWB_MODULE_LAST; ipwb++) {
         char name[256];
         char title[256];
         //h_cal_time_pwbxx_seqsca_04_full_range[ipwb];
         sprintf(name, "h_cal_time_pwb%02d_seqsca_04", ipwb);
         sprintf(title, "calibration pulse time, pwb%02d; time, ns", ipwb);
         h_cal_time_pwbxx[ipwb] = new TH1D(name, title, 400, PAD_PULSER_TIME_625-100.0, PAD_PULSER_TIME_625+100.0);
      }
      
      //h_cal_time_pwbxx_seqsca_04_profile = new TProfile("h_cal_time_pwbxx_seqsca_04_profile", "calibration pulse time, pwbxx seqsca 04; pwbNN; time, ns", PWB_MODULE_LAST+1, -0.5, PWB_MODULE_LAST+0.5, PAD_PULSER_TIME_625-50.0, PAD_PULSER_TIME_625+50.0);

      h_cal_time_pwbaa_seqsca_04_05 = new TH1D("h_cal_time_pwbaa_seqsca_04_05", "calibration pulse time, pwbaa seqsca 05-04; time, ns", 101, -50, 50);
      h_cal_time_pwbaa_seqsca_04_84 = new TH1D("h_cal_time_pwbaa_seqsca_04_84", "calibration pulse time, pwbaa seqsca 84-04; time, ns", 101, -50, 50);
      h_cal_time_pwbaa_seqsca_04_164 = new TH1D("h_cal_time_pwbaa_seqsca_04_164", "calibration pulse time, pwbaa seqsca 164-04; time, ns", 101, -50, 50);

      h_cal_time_pwbbb_seqsca_04_05 = new TH1D("h_cal_time_pwbbb_seqsca_04_05", "calibration pulse time, pwbbb seqsca 05-04; time, ns", 101, -50, 50);
      h_cal_time_pwbbb_seqsca_04_84 = new TH1D("h_cal_time_pwbbb_seqsca_04_84", "calibration pulse time, pwbbb seqsca 84-04; time, ns", 101, -50, 50);
      h_cal_time_pwbbb_seqsca_04_164 = new TH1D("h_cal_time_pwbbb_seqsca_04_164", "calibration pulse time, pwbbb seqsca 164-04; time, ns", 101, -50, 50);

      h_cal_time_pwbcc_seqsca_04_05 = new TH1D("h_cal_time_pwbcc_seqsca_04_05", "calibration pulse time, pwbcc seqsca 05-04; time, ns", 101, -50, 50);
      h_cal_time_pwbcc_seqsca_04_84 = new TH1D("h_cal_time_pwbcc_seqsca_04_84", "calibration pulse time, pwbcc seqsca 84-04; time, ns", 101, -50, 50);
      h_cal_time_pwbcc_seqsca_04_164 = new TH1D("h_cal_time_pwbcc_seqsca_04_164", "calibration pulse time, pwbcc seqsca 164-04; time, ns", 101, -50, 50);

      h_cal_time_pwbdd_seqsca_04_05 = new TH1D("h_cal_time_pwbdd_seqsca_04_05", "calibration pulse time, pwbdd seqsca 05-04; time, ns", 101, -50, 50);
      h_cal_time_pwbdd_seqsca_04_84 = new TH1D("h_cal_time_pwbdd_seqsca_04_84", "calibration pulse time, pwbdd seqsca 84-04; time, ns", 101, -50, 50);
      h_cal_time_pwbdd_seqsca_04_164 = new TH1D("h_cal_time_pwbdd_seqsca_04_164", "calibration pulse time, pwbdd seqsca 164-04; time, ns", 101, -50, 50);

      h_cal_time_pos_00_01_seqsca04 = new TH1D("h_cal_time_pos_00_01_seqsca04", "calibration pulse time, pwbbb-pwbaa seqsca 04; time, ns", 101, -50, 50);
      h_cal_time_pos_00_02_seqsca04 = new TH1D("h_cal_time_pos_00_02_seqsca04", "calibration pulse time, pwbcc-pwbaa seqsca 04; time, ns", 101, -50, 50);
      h_cal_time_pos_00_03_seqsca04 = new TH1D("h_cal_time_pos_00_03_seqsca04", "calibration pulse time, pwbdd-pwbaa seqsca 04; time, ns", 101, -50, 50);
      h_cal_time_pos_01_02_seqsca04 = new TH1D("h_cal_time_pos_01_02_seqsca04", "calibration pulse time, pwbcc-pwbbb seqsca 04; time, ns", 101, -50, 50);
      h_cal_time_pos_01_03_seqsca04 = new TH1D("h_cal_time_pos_01_03_seqsca04", "calibration pulse time, pwbdd-pwbbb seqsca 04; time, ns", 101, -50, 50);
      h_cal_time_pos_02_03_seqsca04 = new TH1D("h_cal_time_pos_02_03_seqsca04", "calibration pulse time, pwbdd-pwbcc seqsca 04; time, ns", 101, -50, 50);

      // timing across PWB and Alpha16

      h_cal_time_pwbaa_seqsca_04_minus_adc_0 = new TH1D("h_cal_time_pwbaa_seqsca_04_minus_adc_0", "calibration pulse time, PWB pwbaa - adc chan 0; time, ns", 101, -50, 50);
      h_cal_time_pwbbb_seqsca_04_minus_adc_0 = new TH1D("h_cal_time_pwbbb_seqsca_04_minus_adc_0", "calibration pulse time, PWB pwbbb - adc chan 0; time, ns", 101, -50, 50);
      h_cal_time_pwbcc_seqsca_04_minus_adc_0 = new TH1D("h_cal_time_pwbcc_seqsca_04_minus_adc_0", "calibration pulse time, PWB pwbcc - adc chan 0; time, ns", 101, -50, 50);
      h_cal_time_pwbdd_seqsca_04_minus_adc_0 = new TH1D("h_cal_time_pwbdd_seqsca_04_minus_adc_0", "calibration pulse time, PWB pwbdd - adc chan 0; time, ns", 101, -50, 50);

      h_cal_time_pwbaa_seqsca_04_minus_adc_16 = new TH1D("h_cal_time_pwbaa_seqsca_04_minus_adc_16", "calibration pulse time, PWB pwbaa - adc chan 16; time, ns", 101, -50, 50);
      h_cal_time_pwbbb_seqsca_04_minus_adc_16 = new TH1D("h_cal_time_pwbbb_seqsca_04_minus_adc_16", "calibration pulse time, PWB pwbbb - adc chan 16; time, ns", 101, -50, 50);
      h_cal_time_pwbcc_seqsca_04_minus_adc_16 = new TH1D("h_cal_time_pwbcc_seqsca_04_minus_adc_16", "calibration pulse time, PWB pwbcc - adc chan 16; time, ns", 101, -50, 50);
      h_cal_time_pwbdd_seqsca_04_minus_adc_16 = new TH1D("h_cal_time_pwbdd_seqsca_04_minus_adc_16", "calibration pulse time, PWB pwbdd - adc chan 16; time, ns", 101, -50, 50);
   }

   TAFlowEvent* AnalyzeFlowEvent(TARunInfo* runinfo, TAFlags* flags, TAFlowEvent* flow)
   {
      //printf("PulserModule::Analyze, run %d, event serno %d, id 0x%04x, data size %d\n", runinfo->fRunNo, event->serial_number, (int)event->event_id, event->data_size);

      if (!fPulser)
         return flow;

      AgEventFlow *ef = flow->Find<AgEventFlow>();

      if (!ef || !ef->fEvent)
         return flow;

      AgAwHitsFlow* eawh = flow->Find<AgAwHitsFlow>();
      AgPadHitsFlow* eph = flow->Find<AgPadHitsFlow>();
      //AgBscAdcHitsFlow* eba = flow->Find<AgBscAdcHitsFlow>();

      AgEvent* age = ef->fEvent;

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

      double pwbaa_seqsca04 = 0;
      double pwbaa_seqsca05 = 0;
      double pwbaa_seqsca84 = 0;
      double pwbaa_seqsca164 = 0;

      double pwbbb_seqsca04 = 0;
      double pwbbb_seqsca05 = 0;
      double pwbbb_seqsca84 = 0;
      double pwbbb_seqsca164 = 0;

      double pwbcc_seqsca04 = 0;
      double pwbcc_seqsca05 = 0;
      double pwbcc_seqsca84 = 0;
      double pwbcc_seqsca164 = 0;

      double pwbdd_seqsca04 = 0;
      double pwbdd_seqsca05 = 0;
      double pwbdd_seqsca84 = 0;
      double pwbdd_seqsca164 = 0;

      if (eph) {
         if (fPrint) {
            printf("PA event %d, time %f, pad hits: %d\n", ef->fEvent->counter, ef->fEvent->time, (int)eph->fPadHits.size());
         }

         if (eph->fPadHits.size() > 0) {
            if (h_cal_time_pwbaa_seqsca_04_full_range == NULL) {
               CreatePwbHist(fPwbAA, fPwbBB, fPwbCC, fPwbDD);
            }
         }

         for (unsigned i=0; i<eph->fPadHits.size(); i++) {
            int ipwb = eph->fPadHits[i].imodule;
            int seqsca = eph->fPadHits[i].seqsca;
            //int xcol = (pos%8)*4 + col;
            //int seqpad = xcol*MAX_FEAM_PAD_ROWS + row;
            //int col = eph->fPadHits[i].tpc_col;
            //int row = eph->fPadHits[i].tpc_row;

            double time = eph->fPadHits[i].time_ns;
            //double amp = eph->fPadHits[i].amp;

            //printf("pos %d, seqsca %d, time %f\n", pos, seqsca, time);
            //printf("PWB ipwb %d\n", ipwb);

            //h_cal_time_pwbxx_seqsca_04_full_range[ipwb]->Fill(time);
            h_cal_time_pwbxx[ipwb]->Fill(time);
            //h_cal_time_pwbxx_seqsca_04_profile->Fill(ipwb, time);

            if (ipwb==fPwbAA) {
               if (seqsca == 4) pwbaa_seqsca04 = time;
               if (seqsca == 5) pwbaa_seqsca05 = time;
               if (seqsca == 84) pwbaa_seqsca84 = time;
               if (seqsca == 164) pwbaa_seqsca164 = time;
            }

            if (ipwb==fPwbBB) {
               if (seqsca == 4) pwbbb_seqsca04 = time;
               if (seqsca == 5) pwbbb_seqsca05 = time;
               if (seqsca == 84) pwbbb_seqsca84 = time;
               if (seqsca == 164) pwbbb_seqsca164 = time;
            }

            if (ipwb==fPwbCC) {
               if (seqsca == 4) pwbcc_seqsca04 = time;
               if (seqsca == 5) pwbcc_seqsca05 = time;
               if (seqsca == 84) pwbcc_seqsca84 = time;
               if (seqsca == 164) pwbcc_seqsca164 = time;
            }

            if (ipwb==fPwbDD) {
               if (seqsca == 4) pwbdd_seqsca04 = time;
               if (seqsca == 5) pwbdd_seqsca05 = time;
               if (seqsca == 84) pwbdd_seqsca84 = time;
               if (seqsca == 164) pwbdd_seqsca164 = time;
            }
         }
      }

      double adc_time[ADC_MODULE_LAST+1][48];

      for (int iadc=ADC_MODULE_FIRST; iadc<=ADC_MODULE_LAST; iadc++) {
         for (int ichan=0; ichan<48; ichan++) {
            adc_time[iadc][ichan] = 0;
         }
      }

      double tdc_time[17];

      for (int ichan=0; ichan<17; ichan++) {
         tdc_time[ichan] = 0;
      }

      if (eawh) {
         if (fPrint) {
            printf("AW event %d, time %f, anode wire hits: %d\n", ef->fEvent->counter, ef->fEvent->time, (int)eawh->fAwHits.size());
         }

         for (unsigned j=0; j<eawh->fAwHits.size(); j++) {
            int adc_module = eawh->fAwHits[j].adc_module;
            int adc_chan = eawh->fAwHits[j].adc_chan;
            //int preamp = eawh->fAwHits[j].preamp_pos;
            //int wire = eawh->fAwHits[j].wire;
            //int preamp = wire/16;
            double time = eawh->fAwHits[j].time;
            //double amp = eawh->fAwHits[j].amp;

            //printf("RRR adc_module %d, adc_chan %d, time %f\n", adc_module, adc_chan, time);

            adc_time[adc_module][adc_chan] = time;
         }
      }

      if (age->tdc) {
         //printf("TdcEvent: "); age->tdc->Print(1); printf("\n");
         
         double tdc_time_0 = 0;

         for (unsigned i=0; i<age->tdc->hits.size(); i++) {
            int ifpga  = age->tdc->hits[i]->fpga;
            int ichan  = age->tdc->hits[i]->chan;
            int re  = age->tdc->hits[i]->rising_edge;
            int coarse_time = age->tdc->hits[i]->coarse_time;
            int fine_time = age->tdc->hits[i]->fine_time;
            double time_ns = coarse_time/200e6*1e+9;

            double fine_time_ns = 0;
            if (ichan==0) {
               fine_time_ns = (fine_time-409.0)/(435.0-409.0) * 0.0;
            } else {
               fine_time_ns =  - (fine_time-17.0)/(450.0-17.0) * 5.0;
            }

            h_cal_time_tdc_fine->Fill(fine_time);

            if (re && ifpga == 1 && ichan >= 0 && ichan <= 16) {

               if (tdc_time_0 == 0) {
                  if (ichan == 0) {
                     h_cal_time_tdc_fine_0->Fill(fine_time);
                     tdc_time_0 = time_ns;
                  }
               }

               if (tdc_time[ichan] == 0) {
                  if (ichan==1)
                     h_cal_time_tdc_fine_1->Fill(fine_time);
                  if (ichan==2)
                     h_cal_time_tdc_fine_2->Fill(fine_time);
                  if (ichan==3)
                     h_cal_time_tdc_fine_3->Fill(fine_time);

                  tdc_time[ichan] = time_ns + fine_time_ns - tdc_time_0;

                  //printf("TDC chan %d, coarse %d, time %f ns, fine %d/%f, relative to %f is %f ns\n", ichan, coarse_time, time_ns, fine_time, fine_time_ns, tdc_time_0, tdc_time[ichan]);
               }
            }
         }

         h_cal_time_tdc_1->Fill(tdc_time[1]);
         h_cal_time_tdc_2->Fill(tdc_time[2]);
         h_cal_time_tdc_3->Fill(tdc_time[3]);
         h_cal_time_tdc_2_1->Fill(tdc_time[2] - tdc_time[1]);
         h_cal_time_tdc_3_1->Fill(tdc_time[3] - tdc_time[1]);
         h_cal_time_tdc_2_1_zoom->Fill(tdc_time[2] - tdc_time[1]);
         h_cal_time_tdc_3_1_zoom->Fill(tdc_time[3] - tdc_time[1]);
      }

      double first_adc_time_0 = 0;
      double first_adc_time_16 = 0;

      if (!f_h_cal_adcxx_00_full_range) {
         f_h_cal_adcxx_00_full_range = true;

         for (int iadc = ADC_MODULE_FIRST; iadc <= ADC_MODULE_LAST; iadc++) {
            if (1 || (adc_time[iadc][0] > 0) || (adc_time[iadc][16] > 0)) {
               char name[256];
               char title[256];

               hdir_pulser->cd();

               sprintf(name, "h_cal_time_adc%02d_00_full_range", iadc);
               sprintf(title, "calibration pulse time, adc%02d chan 00, full time range; time, ns", iadc);
               h_cal_adcxx_00_full_range[iadc] = new TH1D(name, title, 200, 0, MAX_TIME);

               sprintf(name, "h_cal_time_adc%02d_01_full_range", iadc);
               sprintf(title, "calibration pulse time, adc%02d chan 01, full time range; time, ns", iadc);
               h_cal_adcxx_01_full_range[iadc] = new TH1D(name, title, 200, 0, MAX_TIME);

               sprintf(name, "h_cal_time_adc%02d_16_full_range", iadc);
               sprintf(title, "calibration pulse time, adc%02d chan 16, full time range; time, ns", iadc);
               h_cal_adcxx_16_full_range[iadc] = new TH1D(name, title, 200, 0, MAX_TIME);

               sprintf(name, "h_cal_time_adc%02d_17_full_range", iadc);
               sprintf(title, "calibration pulse time, adc%02d chan 17, full time range; time, ns", iadc);
               h_cal_adcxx_17_full_range[iadc] = new TH1D(name, title, 200, 0, MAX_TIME);

               sprintf(name, "h_cal_time_adc%02d_32_full_range", iadc);
               sprintf(title, "calibration pulse time, adc%02d chan 32, full time range; time, ns", iadc);
               h_cal_adcxx_32_full_range[iadc] = new TH1D(name, title, 200, 0, MAX_TIME);

               sprintf(name, "h_cal_time_adc%02d_xx_00", iadc);
               sprintf(title, "calibration pulse time, adc%02d channels vs chan 0; time, ns", iadc);
               h_cal_adcxx_xx_00[iadc] = new TH1D(name, title, 200, -50, +50);

               sprintf(name, "h_cal_time_adc%02d_xx_16", iadc);
               sprintf(title, "calibration pulse time, adc%02d channels vs chan 16; time, ns", iadc);
               h_cal_adcxx_xx_16[iadc] = new TH1D(name, title, 200, -50, +50);

               sprintf(name, "h_cal_time_adc%02d_00_16", iadc);
               sprintf(title, "calibration pulse time, adc%02d chan 00 vs chan 16; time, ns", iadc);
               h_cal_adcxx_00_16[iadc] = new TH1D(name, title, 200, -20, +20);

               sprintf(name, "h_cal_time_adc%02d_17_16_zoom", iadc);
               sprintf(title, "calibration pulse time, adc%02d chan 17 vs chan 16; time, ns", iadc);
               h_cal_adcxx_17_16_zoom[iadc] = new TH1D(name, title, 200, -10, +10);

               sprintf(name, "h_cal_time_adc%02d_32_16_zoom", iadc);
               sprintf(title, "calibration pulse time, adc%02d chan 32 vs chan 16; time, ns", iadc);
               h_cal_adcxx_32_16_zoom[iadc] = new TH1D(name, title, 200, -10, +10);

               sprintf(name, "h_cal_time_adc%02d_profile_00", iadc);
               sprintf(title, "calibration pulse time, adc%02d chan 0..47 vs chan 0; ADC channel, 0..15 are 100MHz, 16..47 are 62.5MHz", iadc);
               h_cal_adcxx_profile_00[iadc] = new TProfile(name, title, 48, -0.5, 47.5, -10, 10);
               h_cal_adcxx_profile_00[iadc]->SetMaximum(+10.0);
               h_cal_adcxx_profile_00[iadc]->SetMinimum(-10.0);

               sprintf(name, "h_cal_time_adc%02d_profile_16", iadc);
               sprintf(title, "calibration pulse time, adc%02d chan 0..47 vs chan 16; ADC channel, 0..15 are 100MHz, 16..47 are 62.5MHz", iadc);
               h_cal_adcxx_profile_16[iadc] = new TProfile(name, title, 48, -0.5, 47.5, -10, 10);
               h_cal_adcxx_profile_16[iadc]->SetMaximum(+10.0);
               h_cal_adcxx_profile_16[iadc]->SetMinimum(-10.0);
            }
         }
      }

      for (int iadc = ADC_MODULE_FIRST; iadc <= ADC_MODULE_LAST; iadc++) {
         if ((adc_time[iadc][0] > 0) || (adc_time[iadc][16] > 0)) {

            if (h_first_adc_0 == 0)
               if (adc_time[iadc][0] > 0)
                  h_first_adc_0 = iadc;

            if (h_first_adc_16 == 0)
               if (adc_time[iadc][16] > 0)
                  h_first_adc_16 = iadc;

            //printf("adc%02d chan 16 %f, chan 17 %f, 17_16 %f\n", iadc, adc_time[iadc][16], adc_time[iadc][17], adc_time[iadc][17]-adc_time[iadc][16]);

            if (iadc == h_first_adc_0) {
               first_adc_time_0 = adc_time[iadc][0];
            }

            if (iadc == h_first_adc_16) {
               first_adc_time_16 = adc_time[iadc][16];
            }

            h_cal_adcxx_00_full_range[iadc]->Fill(adc_time[iadc][0]);
            h_cal_adcxx_01_full_range[iadc]->Fill(adc_time[iadc][1]);

            h_cal_adcxx_16_full_range[iadc]->Fill(adc_time[iadc][16]);
            h_cal_adcxx_17_full_range[iadc]->Fill(adc_time[iadc][17]);
            h_cal_adcxx_17_16_zoom[iadc]->Fill(adc_time[iadc][17] - adc_time[iadc][16]);

            h_cal_adcxx_32_full_range[iadc]->Fill(adc_time[iadc][32]);
            h_cal_adcxx_32_16_zoom[iadc]->Fill(adc_time[iadc][32] - adc_time[iadc][16]);

            for (int ichan=0; ichan<48; ichan++) {
               if (adc_time[iadc][ichan] > 0) {
                  h_cal_adcxx_xx_00[iadc]->Fill(adc_time[iadc][ichan] - adc_time[iadc][0]);
                  h_cal_adcxx_xx_16[iadc]->Fill(adc_time[iadc][ichan] - adc_time[iadc][16]);

                  h_cal_adcxx_profile_00[iadc]->Fill(ichan, adc_time[iadc][ichan] - adc_time[iadc][0]);
                  h_cal_adcxx_profile_16[iadc]->Fill(ichan, adc_time[iadc][ichan] - adc_time[iadc][16]);
               }
            }

            h_cal_adcxx_00_16[iadc]->Fill(adc_time[iadc][0] - adc_time[iadc][16] + 2200 + 50);
         }
      }

      if (h_cal_adcnn_00_all == NULL) {
         char name[256];
         char title[256];
         
         hdir_pulser->cd();
         
         sprintf(name, "h_cal_time_adcNN_00_all");
         sprintf(title, "calibration pulse time, adcNN chan 0 vs chan 0, all ADCs; time, ns");
         h_cal_adcnn_00_all = new TH1D(name, title, 200, -50, +50);
         
         for (int iadc = ADC_MODULE_FIRST; iadc <= ADC_MODULE_LAST; iadc++) {
            if (1 || adc_time[iadc][16] > 0) {
               sprintf(name, "h_cal_time_adc%02d_00_NN", iadc);
               sprintf(title, "calibration pulse time, adc%02d chan 0 vs chan 0; time, ns", iadc);
               h_cal_adcnn_00[iadc] = new TH1D(name, title, 200, -50, +50);
            }
         }
 
         sprintf(name, "h_cal_time_adcNN_profile_00");
         sprintf(title, "calibration pulse time, adcNN chan 0 vs chan 0; ADC module, adcNN");
         h_cal_adcnn_profile_00 = new TProfile(name, title, ADC_MODULE_LAST-ADC_MODULE_FIRST+1, ADC_MODULE_FIRST-0.5, ADC_MODULE_LAST+0.5, -50, 50);
         h_cal_adcnn_profile_00->SetMaximum(+50.0);
         h_cal_adcnn_profile_00->SetMinimum(-50.0);
      }

      if (h_cal_adcnn_16_all == NULL) {
         char name[256];
         char title[256];
         
         hdir_pulser->cd();
         
         sprintf(name, "h_cal_time_adcNN_16_all");
         sprintf(title, "calibration pulse time, adcNN chan 16 vs chan 16, all ADCs; time, ns");
         h_cal_adcnn_16_all = new TH1D(name, title, 200, -50, +50);
         
         for (int iadc = ADC_MODULE_FIRST; iadc <= ADC_MODULE_LAST; iadc++) {
            if (1 || adc_time[iadc][16] > 0) {
               sprintf(name, "h_cal_time_adc%02d_16_NN", iadc);
               sprintf(title, "calibration pulse time, adc%02d chan 16 vs chan 16; time, ns", iadc);
               h_cal_adcnn_16[iadc] = new TH1D(name, title, 200, -50, +50);
            }
         }
 
         sprintf(name, "h_cal_time_adcNN_profile_16");
         sprintf(title, "calibration pulse time, adcNN chan 16 vs chan 16; ADC module, adcNN");
         h_cal_adcnn_profile_16 = new TProfile(name, title, ADC_MODULE_LAST-ADC_MODULE_FIRST+1, ADC_MODULE_FIRST-0.5, ADC_MODULE_LAST+0.5, -50, 50);
         h_cal_adcnn_profile_16->SetMaximum(+50.0);
         h_cal_adcnn_profile_16->SetMinimum(-50.0);
      }

      for (int iadc = ADC_MODULE_FIRST; iadc <= ADC_MODULE_LAST; iadc++) {
         if (h_cal_adcnn_00[iadc] && adc_time[iadc][0] > 0) {
            if (first_adc_time_0 > 0) {
               //printf("adc%02d chan 00 %f, first %f, diff %f\n", iadc, adc_time[iadc][0], first_adc_time_0, adc_time[iadc][0]-first_adc_time_0);
               h_cal_adcnn_00_all->Fill(adc_time[iadc][0] - first_adc_time_0);
               h_cal_adcnn_00[iadc]->Fill(adc_time[iadc][0] - first_adc_time_0);
               h_cal_adcnn_profile_00->Fill(iadc, adc_time[iadc][0] - first_adc_time_0);
            }
         }
         if (h_cal_adcnn_16[iadc] && adc_time[iadc][16] > 0) {
            if (first_adc_time_16 > 0) {
               //printf("adc%02d chan 16 %f, first %f, diff %f\n", iadc, adc_time[iadc][16], first_adc_time_16, adc_time[iadc][16]-first_adc_time_16);
               h_cal_adcnn_16_all->Fill(adc_time[iadc][16] - first_adc_time_16);
               h_cal_adcnn_16[iadc]->Fill(adc_time[iadc][16] - first_adc_time_16);
               h_cal_adcnn_profile_16->Fill(iadc, adc_time[iadc][16] - first_adc_time_16);
            }
         }
      }

      if (1 || ((pwbaa_seqsca04 > 0)
                && (pwbaa_seqsca05 > 0)
                && (pwbbb_seqsca04 > 0)
                && (pwbbb_seqsca05 > 0)
                && (pwbcc_seqsca04 > 0)
                && (pwbcc_seqsca05 > 0))) {

         double offset_0 = 203;
         double offset_16 = -2040.0 +45.0;

#if 0
         printf("PAD times: %f %f %f\n", pwbaa_seqsca04, pwbbb_seqsca04, pwbcc_seqsca04);
#endif

         if (h_cal_time_pwbaa_seqsca_04_full_range) {

            h_cal_time_pwbaa_seqsca_04_full_range->Fill(pwbaa_seqsca04);
            h_cal_time_pwbbb_seqsca_04_full_range->Fill(pwbbb_seqsca04);
            h_cal_time_pwbcc_seqsca_04_full_range->Fill(pwbcc_seqsca04);
            h_cal_time_pwbdd_seqsca_04_full_range->Fill(pwbdd_seqsca04);
            
            h_cal_time_pwbaa_seqsca_04->Fill(pwbaa_seqsca04);
            h_cal_time_pwbbb_seqsca_04->Fill(pwbbb_seqsca04);
            h_cal_time_pwbcc_seqsca_04->Fill(pwbcc_seqsca04);
            h_cal_time_pwbdd_seqsca_04->Fill(pwbdd_seqsca04);
            
            h_cal_time_pwbaa_seqsca_04_05->Fill(pwbaa_seqsca05-pwbaa_seqsca04);
            h_cal_time_pwbaa_seqsca_04_84->Fill(pwbaa_seqsca84-pwbaa_seqsca04);
            h_cal_time_pwbaa_seqsca_04_164->Fill(pwbaa_seqsca164-pwbaa_seqsca04);
            
            h_cal_time_pwbbb_seqsca_04_05->Fill(pwbbb_seqsca05-pwbbb_seqsca04);
            h_cal_time_pwbbb_seqsca_04_84->Fill(pwbbb_seqsca84-pwbbb_seqsca04);
            h_cal_time_pwbbb_seqsca_04_164->Fill(pwbbb_seqsca164-pwbbb_seqsca04);
            
            h_cal_time_pwbcc_seqsca_04_05->Fill(pwbcc_seqsca05-pwbcc_seqsca04);
            h_cal_time_pwbcc_seqsca_04_84->Fill(pwbcc_seqsca84-pwbcc_seqsca04);
            h_cal_time_pwbcc_seqsca_04_164->Fill(pwbcc_seqsca164-pwbcc_seqsca04);
            
            h_cal_time_pwbdd_seqsca_04_05->Fill(pwbdd_seqsca05-pwbdd_seqsca04);
            h_cal_time_pwbdd_seqsca_04_84->Fill(pwbdd_seqsca84-pwbdd_seqsca04);
            h_cal_time_pwbdd_seqsca_04_164->Fill(pwbdd_seqsca164-pwbdd_seqsca04);
            
            h_cal_time_pos_00_01_seqsca04->Fill(pwbbb_seqsca04-pwbaa_seqsca04);
            h_cal_time_pos_00_02_seqsca04->Fill(pwbcc_seqsca04-pwbaa_seqsca04);
            h_cal_time_pos_00_03_seqsca04->Fill(pwbdd_seqsca04-pwbaa_seqsca04);
            h_cal_time_pos_01_02_seqsca04->Fill(pwbcc_seqsca04-pwbbb_seqsca04);
            h_cal_time_pos_01_03_seqsca04->Fill(pwbdd_seqsca04-pwbbb_seqsca04);
            h_cal_time_pos_02_03_seqsca04->Fill(pwbdd_seqsca04-pwbcc_seqsca04);
            
            if (first_adc_time_0 > 0) {
               //printf("PWB %f, ADC_00 %f, diff %f, plot %f\n", pwbaa_seqsca04, first_adc_time_0, pwbaa_seqsca04 - first_adc_time_0, pwbaa_seqsca04 - first_adc_time_0 - offset_0);
               h_cal_time_pwbaa_seqsca_04_minus_adc_0->Fill(pwbaa_seqsca04 - first_adc_time_0 - offset_0);
               h_cal_time_pwbbb_seqsca_04_minus_adc_0->Fill(pwbbb_seqsca04 - first_adc_time_0 - offset_0);
               h_cal_time_pwbcc_seqsca_04_minus_adc_0->Fill(pwbcc_seqsca04 - first_adc_time_0 - offset_0);
               h_cal_time_pwbdd_seqsca_04_minus_adc_0->Fill(pwbdd_seqsca04 - first_adc_time_0 - offset_0);
            }

            if (first_adc_time_16 > 0) {
               //printf("PWB %f, ADC_16 %f, diff %f, plot %f\n", pwbaa_seqsca04, first_adc_time_16, pwbaa_seqsca04 - first_adc_time_16, pwbaa_seqsca04 - first_adc_time_16 - offset_16);
               h_cal_time_pwbaa_seqsca_04_minus_adc_16->Fill(pwbaa_seqsca04 - first_adc_time_16 - offset_16);
               h_cal_time_pwbbb_seqsca_04_minus_adc_16->Fill(pwbbb_seqsca04 - first_adc_time_16 - offset_16);
               h_cal_time_pwbcc_seqsca_04_minus_adc_16->Fill(pwbcc_seqsca04 - first_adc_time_16 - offset_16);
               h_cal_time_pwbdd_seqsca_04_minus_adc_16->Fill(pwbdd_seqsca04 - first_adc_time_16 - offset_16);
            }
         }
      }

      return flow;
   }

   void AnalyzeSpecialEvent(TARunInfo* runinfo, TMEvent* event)
   {
      if (fTrace)
         printf("PulserModule::AnalyzeSpecialEvent, run %d, event serno %d, id 0x%04x, data size %d\n", runinfo->fRunNo, event->serial_number, (int)event->event_id, event->data_size);
   }
};

class PulserModuleFactory: public TAFactory
{
public:
   bool fDoPrint = false;
   bool fPulser  = false;

public:
   void Usage()
   {
      printf("PulserModuleFactory::Usage:\n");
      printf("--print ## be verbose\n");
      printf("--pulser # enable field-wire pulser analysis\n");
   }

   void Init(const std::vector<std::string> &args)
   {
      printf("PulserModuleFactory::Init!\n");
      
      for (unsigned i=0; i<args.size(); i++) {
         if (args[i] == "--print")
            fDoPrint = true;
         if (args[i] == "--pulser")
            fPulser = true;
      }
   }

   void Finish()
   {
      printf("PulserModuleFactory::Finish!\n");
   }

   TARunObject* NewRunObject(TARunInfo* runinfo)
   {
      printf("PulserModule::NewRunObject, run %d, file %s\n", runinfo->fRunNo, runinfo->fFileName.c_str());
      return new PulserModule(runinfo, fDoPrint, fPulser, -1, -1, -1, -1);
   }
};

static TARegister tar(new PulserModuleFactory);

/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */
