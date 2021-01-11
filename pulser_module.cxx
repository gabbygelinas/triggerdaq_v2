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

   TH1D* h_cal_adcxx_16_full_range[ADC_MODULE_LAST+1];
   TH1D* h_cal_adcxx_17_full_range[ADC_MODULE_LAST+1];
   TH1D* h_cal_adcxx_32_full_range[ADC_MODULE_LAST+1];
   TH1D* h_cal_adcxx_xx_16[ADC_MODULE_LAST+1];
   TH1D* h_cal_adcxx_17_16_zoom[ADC_MODULE_LAST+1];
   TH1D* h_cal_adcxx_32_16_zoom[ADC_MODULE_LAST+1];
   TProfile* h_cal_adcxx_profile_16[ADC_MODULE_LAST+1];

   int h_first_adc = 0;
   TH1D* h_cal_adcnn_16_all = NULL;
   TH1D* h_cal_adcnn_16[ADC_MODULE_LAST+1];
   TProfile* h_cal_adcnn_profile_16 = NULL;

   // time of PWB calibration pulse
   TH1D* h_cal_time_pos00_seqsca_04_full_range = NULL;
   TH1D* h_cal_time_pos01_seqsca_04_full_range = NULL;
   TH1D* h_cal_time_pos02_seqsca_04_full_range = NULL;
   TH1D* h_cal_time_pos03_seqsca_04_full_range = NULL;

   TH1D* h_cal_time_pos00_seqsca_04 = NULL;
   TH1D* h_cal_time_pos01_seqsca_04 = NULL;
   TH1D* h_cal_time_pos02_seqsca_04 = NULL;
   TH1D* h_cal_time_pos03_seqsca_04 = NULL;

   // time across to next channel
   TH1D* h_cal_time_pos00_seqsca_04_05 = NULL;
   TH1D* h_cal_time_pos01_seqsca_04_05 = NULL;
   TH1D* h_cal_time_pos02_seqsca_04_05 = NULL;
   TH1D* h_cal_time_pos03_seqsca_04_05 = NULL;

   // time across to next SCA
   TH1D* h_cal_time_pos00_seqsca_04_84 = NULL;
   TH1D* h_cal_time_pos01_seqsca_04_84 = NULL;
   TH1D* h_cal_time_pos02_seqsca_04_84 = NULL;
   TH1D* h_cal_time_pos03_seqsca_04_84 = NULL;

   // time across to next ADC
   TH1D* h_cal_time_pos00_seqsca_04_164 = NULL;
   TH1D* h_cal_time_pos01_seqsca_04_164 = NULL;
   TH1D* h_cal_time_pos02_seqsca_04_164 = NULL;
   TH1D* h_cal_time_pos03_seqsca_04_164 = NULL;

   // time across 2 PWB boards
   TH1D* h_cal_time_pos_01_02_seqsca04 = NULL;
   TH1D* h_cal_time_pos_01_02_seqsca05 = NULL;
   TH1D* h_cal_time_pos_01_03_seqsca04 = NULL;
   TH1D* h_cal_time_pos_01_03_seqsca05 = NULL;

   // time across 100MHz ADC to PWB
   TH1D* h_cal_time_pos00_seqsca_04_minus_adc_0 = NULL;
   TH1D* h_cal_time_pos01_seqsca_04_minus_adc_0 = NULL;
   TH1D* h_cal_time_pos02_seqsca_04_minus_adc_0 = NULL;
   TH1D* h_cal_time_pos03_seqsca_04_minus_adc_0 = NULL;

   // time across 62.5MHz ADC to PWB
   TH1D* h_cal_time_pos00_seqsca_04_minus_adc_16 = NULL;
   TH1D* h_cal_time_pos01_seqsca_04_minus_adc_16 = NULL;
   TH1D* h_cal_time_pos02_seqsca_04_minus_adc_16 = NULL;
   TH1D* h_cal_time_pos03_seqsca_04_minus_adc_16 = NULL;

   bool fTrace = false;

   PulserModule(TARunInfo* runinfo, bool do_print)
      : TARunObject(runinfo)
   {
      if (fTrace)
         printf("PulserModule::ctor!\n");

      fPrint = do_print;
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

      runinfo->fRoot->fOutputFile->cd(); // select correct ROOT directory

      hdir_pulser = gDirectory->mkdir("pulser");
      hdir_pulser->cd(); // select correct ROOT directory

      // ADC timing

      // PWB timing

      h_cal_time_pos00_seqsca_04_full_range = new TH1D("h_cal_time_pos00_seqsca_04_full_range", "calibration pulse time, pos00 seqsca 04, full range; time, ns", 200, 0, MAX_TIME);
      h_cal_time_pos01_seqsca_04_full_range = new TH1D("h_cal_time_pos01_seqsca_04_full_range", "calibration pulse time, pos01 seqsca 04, full range; time, ns", 200, 0, MAX_TIME);
      h_cal_time_pos02_seqsca_04_full_range = new TH1D("h_cal_time_pos02_seqsca_04_full_range", "calibration pulse time, pos02 seqsca 04, full range; time, ns", 200, 0, MAX_TIME);
      h_cal_time_pos03_seqsca_04_full_range = new TH1D("h_cal_time_pos03_seqsca_04_full_range", "calibration pulse time, pos03 seqsca 04, full range; time, ns", 200, 0, MAX_TIME);

      h_cal_time_pos00_seqsca_04 = new TH1D("h_cal_time_pos00_seqsca_04", "calibration pulse time, pos00 seqsca 04; time, ns", 200, PAD_PULSER_TIME_625-50.0, PAD_PULSER_TIME_625+50.0);
      h_cal_time_pos01_seqsca_04 = new TH1D("h_cal_time_pos01_seqsca_04", "calibration pulse time, pos01 seqsca 04; time, ns", 200, PAD_PULSER_TIME_625-50.0, PAD_PULSER_TIME_625+50.0);
      h_cal_time_pos02_seqsca_04 = new TH1D("h_cal_time_pos02_seqsca_04", "calibration pulse time, pos02 seqsca 04; time, ns", 200, PAD_PULSER_TIME_625-50.0, PAD_PULSER_TIME_625+50.0);
      h_cal_time_pos03_seqsca_04 = new TH1D("h_cal_time_pos03_seqsca_04", "calibration pulse time, pos03 seqsca 04; time, ns", 200, PAD_PULSER_TIME_625-50.0, PAD_PULSER_TIME_625+50.0);

      h_cal_time_pos00_seqsca_04_05 = new TH1D("h_cal_time_pos00_seqsca_04_05", "calibration pulse time, pos00 seqsca 05-04; time, ns", 101, -50, 50);
      h_cal_time_pos00_seqsca_04_84 = new TH1D("h_cal_time_pos00_seqsca_04_84", "calibration pulse time, pos00 seqsca 84-04; time, ns", 101, -50, 50);
      h_cal_time_pos00_seqsca_04_164 = new TH1D("h_cal_time_pos00_seqsca_04_164", "calibration pulse time, pos00 seqsca 164-04; time, ns", 101, -50, 50);

      h_cal_time_pos01_seqsca_04_05 = new TH1D("h_cal_time_pos01_seqsca_04_05", "calibration pulse time, pos01 seqsca 05-04; time, ns", 101, -50, 50);
      h_cal_time_pos01_seqsca_04_84 = new TH1D("h_cal_time_pos01_seqsca_04_84", "calibration pulse time, pos01 seqsca 84-04; time, ns", 101, -50, 50);
      h_cal_time_pos01_seqsca_04_164 = new TH1D("h_cal_time_pos01_seqsca_04_164", "calibration pulse time, pos01 seqsca 164-04; time, ns", 101, -50, 50);

      h_cal_time_pos02_seqsca_04_05 = new TH1D("h_cal_time_pos02_seqsca_04_05", "calibration pulse time, pos02 seqsca 05-04; time, ns", 101, -50, 50);
      h_cal_time_pos02_seqsca_04_84 = new TH1D("h_cal_time_pos02_seqsca_04_84", "calibration pulse time, pos02 seqsca 84-04; time, ns", 101, -50, 50);
      h_cal_time_pos02_seqsca_04_164 = new TH1D("h_cal_time_pos02_seqsca_04_164", "calibration pulse time, pos02 seqsca 164-04; time, ns", 101, -50, 50);

      h_cal_time_pos03_seqsca_04_05 = new TH1D("h_cal_time_pos03_seqsca_04_05", "calibration pulse time, pos03 seqsca 05-04; time, ns", 101, -50, 50);
      h_cal_time_pos03_seqsca_04_84 = new TH1D("h_cal_time_pos03_seqsca_04_84", "calibration pulse time, pos03 seqsca 84-04; time, ns", 101, -50, 50);
      h_cal_time_pos03_seqsca_04_164 = new TH1D("h_cal_time_pos03_seqsca_04_164", "calibration pulse time, pos03 seqsca 164-04; time, ns", 101, -50, 50);

      h_cal_time_pos_01_02_seqsca04 = new TH1D("h_cal_time_pos_01_02_seqsca04", "calibration pulse time, pos02-pos01 seqsca 04; time, ns", 101, -50, 50);
      h_cal_time_pos_01_02_seqsca05 = new TH1D("h_cal_time_pos_01_02_seqsca05", "calibration pulse time, pos02-pos01 seqsca 05; time, ns", 101, -50, 50);

      h_cal_time_pos_01_03_seqsca04 = new TH1D("h_cal_time_pos_01_03_seqsca04", "calibration pulse time, pos03-pos01 seqsca 04; time, ns", 101, -50, 50);
      h_cal_time_pos_01_03_seqsca05 = new TH1D("h_cal_time_pos_01_03_seqsca05", "calibration pulse time, pos03-pos01 seqsca 05; time, ns", 101, -50, 50);

      // timing across PWB and Alpha16

      h_cal_time_pos00_seqsca_04_minus_adc_0 = new TH1D("h_cal_time_pos00_seqsca_04_minus_adc_0", "calibration pulse time, PWB pos00 - adc chan 0; time, ns", 101, -50, 50);
      h_cal_time_pos01_seqsca_04_minus_adc_0 = new TH1D("h_cal_time_pos01_seqsca_04_minus_adc_0", "calibration pulse time, PWB pos01 - adc chan 0; time, ns", 101, -50, 50);
      h_cal_time_pos02_seqsca_04_minus_adc_0 = new TH1D("h_cal_time_pos02_seqsca_04_minus_adc_0", "calibration pulse time, PWB pos02 - adc chan 0; time, ns", 101, -50, 50);
      h_cal_time_pos03_seqsca_04_minus_adc_0 = new TH1D("h_cal_time_pos03_seqsca_04_minus_adc_0", "calibration pulse time, PWB pos03 - adc chan 0; time, ns", 101, -50, 50);

      h_cal_time_pos00_seqsca_04_minus_adc_16 = new TH1D("h_cal_time_pos00_seqsca_04_minus_adc_16", "calibration pulse time, PWB pos00 - adc chan 16; time, ns", 101, -50, 50);
      h_cal_time_pos01_seqsca_04_minus_adc_16 = new TH1D("h_cal_time_pos01_seqsca_04_minus_adc_16", "calibration pulse time, PWB pos01 - adc chan 16; time, ns", 101, -50, 50);
      h_cal_time_pos02_seqsca_04_minus_adc_16 = new TH1D("h_cal_time_pos02_seqsca_04_minus_adc_16", "calibration pulse time, PWB pos02 - adc chan 16; time, ns", 101, -50, 50);
      h_cal_time_pos03_seqsca_04_minus_adc_16 = new TH1D("h_cal_time_pos03_seqsca_04_minus_adc_16", "calibration pulse time, PWB pos03 - adc chan 16; time, ns", 101, -50, 50);
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

   TAFlowEvent* AnalyzeFlowEvent(TARunInfo* runinfo, TAFlags* flags, TAFlowEvent* flow)
   {
      //printf("PulserModule::Analyze, run %d, event serno %d, id 0x%04x, data size %d\n", runinfo->fRunNo, event->serial_number, (int)event->event_id, event->data_size);

      AgEventFlow *ef = flow->Find<AgEventFlow>();

      if (!ef || !ef->fEvent)
         return flow;

      AgAwHitsFlow* eawh = flow->Find<AgAwHitsFlow>();
      AgPadHitsFlow* eph = flow->Find<AgPadHitsFlow>();
      AgBscAdcHitsFlow* eba = flow->Find<AgBscAdcHitsFlow>();

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

      double pos00_seqsca04 = -10;
      double pos00_seqsca05 = -20;
      double pos00_seqsca84 = -130;
      double pos00_seqsca164 = -140;

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
      double pos03_seqsca84 = -150;
      double pos03_seqsca164 = -160;

      if (eph) {
         if (fPrint) {
            printf("PA event %d, time %f, pad hits: %d\n", ef->fEvent->counter, ef->fEvent->time, (int)eph->fPadHits.size());
         }

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
               if (seqsca == 84) pos00_seqsca84 = time;
               if (seqsca == 164) pos00_seqsca164 = time;
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
               if (seqsca == 84) pos03_seqsca84 = time;
               if (seqsca == 164) pos03_seqsca164 = time;
            }
         }
      }

      double adc_time[ADC_MODULE_LAST+1][48];

      for (int iadc=ADC_MODULE_FIRST; iadc<=ADC_MODULE_LAST; iadc++) {
         for (int ichan=0; ichan<48; ichan++) {
            adc_time[iadc][ichan] = 0;
         }
      }

      double first_adc_time_0 = 0;
      double first_adc_time_16 = 0;

      for (int iadc = ADC_MODULE_FIRST; iadc <= ADC_MODULE_LAST; iadc++) {
         if (adc_time[iadc][16] > 0) {
            if (h_cal_adcxx_16_full_range[iadc] == NULL) {
               char name[256];
               char title[256];

               hdir_pulser->cd();

               sprintf(name, "h_cal_time_adc%02d_16_full_range", iadc);
               sprintf(title, "calibration pulse time, adc%02d chan 16, full time range; time, ns", iadc);
               h_cal_adcxx_16_full_range[iadc] = new TH1D(name, title, 200, 0, MAX_TIME);

               sprintf(name, "h_cal_time_adc%02d_17_full_range", iadc);
               sprintf(title, "calibration pulse time, adc%02d chan 17, full time range; time, ns", iadc);
               h_cal_adcxx_17_full_range[iadc] = new TH1D(name, title, 200, 0, MAX_TIME);

               sprintf(name, "h_cal_time_adc%02d_32_full_range", iadc);
               sprintf(title, "calibration pulse time, adc%02d chan 32, full time range; time, ns", iadc);
               h_cal_adcxx_32_full_range[iadc] = new TH1D(name, title, 200, 0, MAX_TIME);

               sprintf(name, "h_cal_time_adc%02d_xx_16", iadc);
               sprintf(title, "calibration pulse time, adc%02d channels vs chan 16; time, ns", iadc);
               h_cal_adcxx_xx_16[iadc] = new TH1D(name, title, 200, -50, +50);

               sprintf(name, "h_cal_time_adc%02d_17_16_zoom", iadc);
               sprintf(title, "calibration pulse time, adc%02d chan 17 vs chan 16; time, ns", iadc);
               h_cal_adcxx_17_16_zoom[iadc] = new TH1D(name, title, 200, -5, +5);

               sprintf(name, "h_cal_time_adc%02d_32_16_zoom", iadc);
               sprintf(title, "calibration pulse time, adc%02d chan 32 vs chan 16; time, ns", iadc);
               h_cal_adcxx_32_16_zoom[iadc] = new TH1D(name, title, 200, -5, +5);

               sprintf(name, "h_cal_time_adc%02d_profile_16", iadc);
               sprintf(title, "calibration pulse time, adc%02d chan 0..47 vs chan 16; ADC channel, 0..15 are 100MHz, 16..47 are 62.5MHz", iadc);
               h_cal_adcxx_profile_16[iadc] = new TProfile(name, title, 48, -0.5, 47.5, -10, 10);
               h_cal_adcxx_profile_16[iadc]->SetMaximum(+5.0);
               h_cal_adcxx_profile_16[iadc]->SetMinimum(-5.0);

               if (h_first_adc == 0)
                  h_first_adc = iadc;
            }

            //printf("adc%02d chan 16 %f, chan 17 %f, 17_16 %f\n", iadc, adc_time[iadc][16], adc_time[iadc][17], adc_time[iadc][17]-adc_time[iadc][16]);

            if (iadc == h_first_adc) {
               first_adc_time_16 = adc_time[iadc][16];
               first_adc_time_0 = adc_time[iadc][0];
            }

            h_cal_adcxx_16_full_range[iadc]->Fill(adc_time[iadc][16]);
            h_cal_adcxx_17_full_range[iadc]->Fill(adc_time[iadc][17]);
            h_cal_adcxx_17_16_zoom[iadc]->Fill(adc_time[iadc][17] - adc_time[iadc][16]);

            h_cal_adcxx_32_full_range[iadc]->Fill(adc_time[iadc][32]);
            h_cal_adcxx_32_16_zoom[iadc]->Fill(adc_time[iadc][32] - adc_time[iadc][16]);

            for (int ichan=0; ichan<48; ichan++) {
               if (adc_time[iadc][ichan] > 0) {
                  h_cal_adcxx_xx_16[iadc]->Fill(adc_time[iadc][ichan] - adc_time[iadc][16]);
                  h_cal_adcxx_profile_16[iadc]->Fill(ichan, adc_time[iadc][ichan] - adc_time[iadc][16]);
               }
            }
         }
      }

      if (h_cal_adcnn_16_all == NULL) {
         char name[256];
         char title[256];
         
         hdir_pulser->cd();
         
         sprintf(name, "h_cal_time_adcNN_16_all");
         sprintf(title, "calibration pulse time, adcNN chan 16 vs chan 16, all ADCs; time, ns");
         h_cal_adcnn_16_all = new TH1D(name, title, 200, -50, +50);
         
         for (int iadc = ADC_MODULE_FIRST; iadc <= ADC_MODULE_LAST; iadc++) {
            if (adc_time[iadc][16] > 0) {
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
         if (adc_time[iadc][16] > 0) {
            if (first_adc_time_16 > 0) {
               printf("adc%02d chan 16 %f, first %f, diff %f\n", iadc, adc_time[iadc][16], first_adc_time_16, adc_time[iadc][16]-first_adc_time_16);
               h_cal_adcnn_16_all->Fill(adc_time[iadc][16] - first_adc_time_16);
               h_cal_adcnn_16[iadc]->Fill(adc_time[iadc][16] - first_adc_time_16);
               h_cal_adcnn_profile_16->Fill(iadc, adc_time[iadc][16] - first_adc_time_16);
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

#if 0
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
#endif

         double pulse_width = -2040.0;
         //double pulse_width = 5350.0 + 30.0 + 40.0;
         //double xpad = pos01_seqsca04;
         // double t = xpad - adc5_0 - pulse_width;
         //printf("aw %.1f pad %.1f %.1f, diff %.1f\n", adc5_0, pos01_seqsca04, xpad, t);

#if 0
         printf("PAD times: %f %f %f\n", pos00_seqsca04, pos01_seqsca04, pos02_seqsca04);
#endif

         h_cal_time_pos00_seqsca_04_full_range->Fill(pos00_seqsca04);
         h_cal_time_pos01_seqsca_04_full_range->Fill(pos01_seqsca04);
         h_cal_time_pos02_seqsca_04_full_range->Fill(pos02_seqsca04);
         h_cal_time_pos03_seqsca_04_full_range->Fill(pos03_seqsca04);

         h_cal_time_pos00_seqsca_04->Fill(pos00_seqsca04);
         h_cal_time_pos01_seqsca_04->Fill(pos01_seqsca04);
         h_cal_time_pos02_seqsca_04->Fill(pos02_seqsca04);
         h_cal_time_pos03_seqsca_04->Fill(pos03_seqsca04);

         h_cal_time_pos00_seqsca_04_05->Fill(pos00_seqsca05-pos00_seqsca04);
         h_cal_time_pos00_seqsca_04_84->Fill(pos00_seqsca84-pos00_seqsca04);
         h_cal_time_pos00_seqsca_04_164->Fill(pos00_seqsca164-pos00_seqsca04);

         h_cal_time_pos01_seqsca_04_05->Fill(pos01_seqsca05-pos01_seqsca04);
         h_cal_time_pos01_seqsca_04_84->Fill(pos01_seqsca84-pos01_seqsca04);
         h_cal_time_pos01_seqsca_04_164->Fill(pos01_seqsca164-pos01_seqsca04);

         h_cal_time_pos02_seqsca_04_05->Fill(pos02_seqsca05-pos02_seqsca04);
         h_cal_time_pos02_seqsca_04_84->Fill(pos02_seqsca84-pos02_seqsca04);
         h_cal_time_pos02_seqsca_04_164->Fill(pos02_seqsca164-pos02_seqsca04);

         h_cal_time_pos03_seqsca_04_05->Fill(pos03_seqsca05-pos03_seqsca04);
         h_cal_time_pos03_seqsca_04_84->Fill(pos03_seqsca84-pos03_seqsca04);
         h_cal_time_pos03_seqsca_04_164->Fill(pos03_seqsca164-pos03_seqsca04);

         h_cal_time_pos_01_02_seqsca04->Fill(pos02_seqsca04-pos01_seqsca04);
         h_cal_time_pos_01_02_seqsca05->Fill(pos02_seqsca05-pos01_seqsca05);

         h_cal_time_pos_01_03_seqsca04->Fill(pos03_seqsca04-pos01_seqsca04);
         h_cal_time_pos_01_03_seqsca05->Fill(pos03_seqsca05-pos01_seqsca05);

         if (first_adc_time_0 > 0) {
            h_cal_time_pos00_seqsca_04_minus_adc_0->Fill(pos00_seqsca04 - first_adc_time_0 - pulse_width);
            h_cal_time_pos01_seqsca_04_minus_adc_0->Fill(pos01_seqsca04 - first_adc_time_0 - pulse_width);
            h_cal_time_pos02_seqsca_04_minus_adc_0->Fill(pos02_seqsca04 - first_adc_time_0 - pulse_width);
            h_cal_time_pos03_seqsca_04_minus_adc_0->Fill(pos03_seqsca04 - first_adc_time_0 - pulse_width);
         }

         if (first_adc_time_16 > 0) {
            printf("PWB %f, ADC %f, diff %f, plot %f\n", pos00_seqsca04, first_adc_time_16, pos00_seqsca04 - first_adc_time_16, pos00_seqsca04 - first_adc_time_16 - pulse_width);
            h_cal_time_pos00_seqsca_04_minus_adc_16->Fill(pos00_seqsca04 - first_adc_time_16 - pulse_width);
            h_cal_time_pos01_seqsca_04_minus_adc_16->Fill(pos01_seqsca04 - first_adc_time_16 - pulse_width);
            h_cal_time_pos02_seqsca_04_minus_adc_16->Fill(pos02_seqsca04 - first_adc_time_16 - pulse_width);
            h_cal_time_pos03_seqsca_04_minus_adc_16->Fill(pos03_seqsca04 - first_adc_time_16 - pulse_width);
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

public:
   void Usage()
   {
      printf("PulserModuleFactory::Usage:\n");
      printf("--print ## be verbose\n");
   }

   void Init(const std::vector<std::string> &args)
   {
      printf("PulserModuleFactory::Init!\n");
      
      for (unsigned i=0; i<args.size(); i++) {
         if (args[i] == "--print")
            fDoPrint = true;
      }
   }

   void Finish()
   {
      printf("PulserModuleFactory::Finish!\n");
   }

   TARunObject* NewRunObject(TARunInfo* runinfo)
   {
      printf("PulserModule::NewRunObject, run %d, file %s\n", runinfo->fRunNo, runinfo->fFileName.c_str());
      return new PulserModule(runinfo, fDoPrint);
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
