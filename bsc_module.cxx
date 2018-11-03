//
// bsc_module.cxx - q.a. for the ADC BSC data
//
// K.Olchanski
//

#include "manalyzer.h"
#include "midasio.h"

#include <assert.h> // assert()

#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TVirtualFFT.h"

#include "Alpha16.h"
#include "Unpack.h"
#include "AgFlow.h"
#include "ko_limits.h"

#define MAX_BSC_BRMS 2000
#define MAX_AW_BRANGE 2000
// zoom on the ADC pedestal
#define MAX_BAR_PED 5000

static void compute_mean_rms(const int* aptr, int start, int end, double* xmean, double* xrms, double* xmin, double* xmax)
{
   double sum0 = 0;
   double sum1 = 0;
   double sum2 = 0;
   
   double bmin = aptr[start]; // baseline minimum
   double bmax = aptr[start]; // baseline maximum
   
   for (int i=start; i<end; i++) {
      double a = aptr[i];
      sum0 += 1;
      sum1 += a;
      sum2 += a*a;
      if (a < bmin)
         bmin = a;
      if (a > bmax)
         bmax = a;
   }
   
   double bmean = 0;
   double bvar = 0;
   double brms = 0;
   
   if (sum0 > 0) {
      bmean = sum1/sum0;
      bvar = sum2/sum0 - bmean*bmean;
      if (bvar>0)
         brms = sqrt(bvar);
   }

   if (xmean)
      *xmean = bmean;
   if (xrms)
      *xrms = brms;
   if (xmin)
      *xmin = bmin;
   if (xmax)
      *xmax = bmax;
}

struct PlotHistograms
{
   TH1D* fHbaselineMean;
   TH1D* fHbaselineRms;
   TH1D* fHbaselineRange;
   TH1D* fHbaselineBadRmsBscMap;
   TProfile* fHbaselineMeanBscMap;
   TProfile* fHbaselineRmsBscMap;
   TProfile* fHbaselineRangeBscMap;

   TH1D* fHrange;
   TH2D* fHrangeBscMap;

   TH1D* fBscMapPh500;
   TH1D* fBscMapPh1000;
   TH1D* fBscMapPh5000;
   TH1D* fBscMapPh10000;

   TH1D* fHph;
   TH1D* fHped;
   TH1D* fHphHit;
   TH1D* fHphHitBscMap;
   TH1D* fHle;
   TH1D* fHbscHitTime0;
   TH1D* fHbscHitTime;
   TH1D* fHbscHitAmp;
   TH1D* fHbscHitMap;

   PlotHistograms() // ctor
   {
      TDirectory *dir = gDirectory->mkdir("summary");
      dir->cd();

      fHbaselineMean = new TH1D("bsc_adc_baseline_mean", "waveform baseline mean, all wires; ADC counts", 100, -2000, 2000);
      fHbaselineRms = new TH1D("bsc_adc_baseline_rms", "waveform baseline rms, all wires; ADC counts", 100, 0, MAX_BSC_BRMS);
      fHbaselineRange = new TH1D("bsc_adc_baseline_range", "waveform baseline range, all wires; ADC counts", 100, 0, MAX_AW_BRANGE);

      fHbaselineBadRmsBscMap = new TH1D("bsc_adc_baseline_bad_rms_vs_bar", "BSC waveform with bad baseline vs wire number; TPC wire number", NUM_BSC, -0.5, NUM_BSC-0.5);

      fHbaselineMeanBscMap = new TProfile("bsc_adc_baseline_mean_vs_bar", "BSC waveform baseline mean vs wire number; TPC wire number; ADC counts", NUM_BSC, -0.5, NUM_BSC-0.5);
      fHbaselineRmsBscMap = new TProfile("bsc_adc_baseline_rms_vs_bar", "BSC waveform baseline rms vs wire number; TPC wire number; ADC counts", NUM_BSC, -0.5, NUM_BSC-0.5);
      fHbaselineRmsBscMap->SetMinimum(0);
      fHbaselineRangeBscMap = new TProfile("bsc_adc_baseline_range_vs_aw", "BSC waveform baseline range vs wire number; TPC wire number; ADC counts", NUM_BSC, -0.5, NUM_BSC-0.5);

      fHrange = new TH1D("bsc_adc_range", "BSC waveform range, max-min; ADC counts", 100, 0, MAX_AW_AMP);
      fHrangeBscMap = new TH2D("bsc_adc_range_vs_bar", "waveform range, max-min vs wire number; TPC wire number; ADC counts", NUM_BSC, -0.5, NUM_BSC-0.5,1000,0,10000);

      fBscMapPh500 = new TH1D("bsc_map_ph500_vs_bar", "BSC waveforms with ph > 500 vs wire number; TPC wire number", NUM_BSC, -0.5, NUM_BSC-0.5);
      fBscMapPh1000 = new TH1D("bsc_map_ph1000_vs_bar", "BSC waveforms with ph > 1000 vs wire number; TPC wire number", NUM_BSC, -0.5, NUM_BSC-0.5);
      fBscMapPh5000 = new TH1D("bsc_map_ph5000_vs_bar", "BSC waveforms with ph > 5000 vs wire number; TPC wire number", NUM_BSC, -0.5, NUM_BSC-0.5);
      fBscMapPh10000 = new TH1D("bsc_map_ph10000_vs_bar", "BSC waveforms with ph > 10000 vs wire number; TPC wire number", NUM_BSC, -0.5, NUM_BSC-0.5);

      fHph = new TH1D("bsc_adc_pulse_height", "BSC waveform pulse height; ADC counts", 100, 0, MAX_AW_AMP);

      fHped = new TH1D("bsc_adc_pedestal_pulse_height", "waveform pulse height, zoom on the pedestal; ADC counts", 100, 0, MAX_BAR_PED);
      
      fHphHit = new TH1D("bsc_adc_pulse_height_hit", "BSC pulse height, after ADC cut; ADC counts", 100, 0, MAX_AW_AMP);
      
      fHphHitBscMap = new TH1D("bsc_adc_pulse_height_hit_vs_bar", "BSC bar map, after ADC cut; BSC bar number", NUM_BSC, -0.5, NUM_BSC-0.5);

      fHle = new TH1D("bsc_adc_pulse_time", "BSC waveform pulse time, after ADC cut; ADC time bins", 200, 0, 1000);

      fHbscHitTime0 = new TH1D("bsc_adc_hit_time0", "BSC hit time, after ADC cut; time, ns", 200, 0, MAX_TIME);
      fHbscHitTime = new TH1D("bsc_adc_hit_time", "BSC hit time, after ADC and time cut; time, ns", 200, 0, MAX_TIME);
      fHbscHitAmp = new TH1D("bsc_adc_hit_amp", "BSC hit amplitude, after ADC and time cut; pulse height, ADC c.u.", 200, 0, MAX_AW_AMP);
      fHbscHitMap = new TH1D("bsc_adc_hit_map", "BSC hit map, after ADC and time cut; BSC bar", NUM_BSC, -0.5, NUM_BSC-0.5);
   }
};

class AnalyzeNoise
{
public:
   TDirectory* fDirTmp = NULL;
   TDirectory* fDirFile = NULL;
   TH1D* fWaveform = NULL;
   TH1* fFft = NULL;
   TH1D* fFftSum = NULL;
   int fCount = 0;

public:
   AnalyzeNoise(const char* prefix, TDirectory* file, TDirectory* tmp, int nbins = 701) // ctor
   {
      fDirFile = file;
      fDirFile->cd();

      fFftSum = new TH1D((std::string(prefix) + "_fftsum").c_str(), (std::string(prefix) + " FFT sum").c_str(), nbins, 0, nbins);

      fDirTmp = tmp;
      fDirTmp->cd();

      fWaveform = new TH1D((std::string(prefix) + "_waveform").c_str(), (std::string(prefix) + " waveform").c_str(), nbins, 0, nbins);
   }

   void AddWaveform(const std::vector<int>& w)
   {
      for(unsigned b = 0; b < w.size(); b++) {
         double a = w[b];
         //a += 1000*sin(b*2.0*M_PI/(700.0/20.0));
         //a += 100*sin(b*2.0*M_PI/(700.0/100.0));
         fWaveform->SetBinContent(b+1, a);
      }

      fDirTmp->cd();
      
      TVirtualFFT::SetTransform(0);
      fFft = fWaveform->FFT(fFft, "MAG");
      fFftSum->Add(fFft);
      fCount++;
   }

   void Finish()
   {
      if (fCount > 0) {
         fFftSum->Scale(1./(fCount*sqrt(701.)));
      }
   }
};

class PlotNoise
{
public:
   TCanvas* c = NULL;

   PlotNoise(const char* prefix) // ctor
   {
      std::string name = std::string(prefix) + "_cnoise";
      std::string title = std::string(prefix) + " noise analysis";
      c = new TCanvas(name.c_str(), title.c_str(), 700, 700);
   }

   ~PlotNoise() // dtor
   {
      if (c)
         delete c;
      c = NULL;
   }

   void Plot(AnalyzeNoise *n)
   {
      c->Clear();
      c->Divide(1,3);
      c->cd(1);
      n->fWaveform->Draw();

      c->cd(2);
      gPad->SetLogy();
      n->fFft->Draw();

      c->cd(3);
      gPad->SetLogy();
      n->fFftSum->Draw();
      
      c->Modified();
      c->Draw();
      c->Update();
   }
};

class A16ChanHistograms
{
public:
   std::string fNameBase;
   std::string fTitleBase;
   int fNbins = 0;

   TH1D* hbmean = NULL;
   TH1D* hbrms  = NULL;

public:
   A16ChanHistograms(const char* xname, const char* xtitle, int nbins) // ctor
   {
      //printf("Create name [%s] title [%s] with %d bins\n", xname, xtitle, nbins);

      fNameBase = xname;
      fTitleBase = xtitle;
      fNbins = nbins;

      //char name[256];
      //char title[256];
   }

   ~A16ChanHistograms() // dtor
   {

   }

};

class A16Flags
{
public:
   bool fPrint = false;
   bool fFft = false;
   bool fFilterWaveform = false;
   bool fInvertWaveform = false;
};

static double find_pulse_time(const int* adc, int nbins, double baseline, double gain, double threshold)
{
   for (int i=1; i<nbins; i++) {
      double v1 = (adc[i]-baseline)*gain;
      if (v1 > threshold) {
         double v0 = (adc[i-1]-baseline)*gain;
         if (!(v0 <= threshold))
            return 0;
         double ii = i-1+(v0-threshold)/(v0-v1);
         //printf("find_pulse_time: %f %f %f, bins %d %f %d\n", v0, threshold, v1, i-1, ii, i);
         return ii;
      }
   }

   return 0;
}

class BscModule: public TARunObject
{
public:
   A16Flags* fFlags = NULL;
   int fCounter = 0;
   bool fTrace = false;

   PlotHistograms* fH = NULL;
   PlotNoise* fPN16 = NULL;
   AnalyzeNoise* fAN16 = NULL;
   PlotNoise* fPN32 = NULL;
   AnalyzeNoise* fAN32 = NULL;

   std::vector<AnalyzeNoise> fAN16AWB;

public:
   BscModule(TARunInfo* runinfo, A16Flags* f)
      : TARunObject(runinfo)
   {
      if (fTrace)
         printf("BscModule::ctor!\n");
      fFlags = f;

      runinfo->fRoot->fOutputFile->cd();
      TDirectory* bsc = gDirectory->mkdir("bsc");
      bsc->cd(); // select correct ROOT directory

      fH = new PlotHistograms();

      TDirectory* fft_file = bsc->mkdir("noise_fft");
      TDirectory* fft_tmp = runinfo->fRoot->fgDir->mkdir("aw_noise_fft");

      if (fFlags->fFft) {
         fAN16 = new AnalyzeNoise("adc16", fft_file, fft_tmp, 701);
         fPN16 = new PlotNoise("adc16");

         for(int awb=0; awb<32; ++awb)
            {
               TString ANname = TString::Format("adc16AWB%02d",awb);
               fAN16AWB.emplace_back(ANname.Data(), fft_file, fft_tmp, 701);
            }

         fAN32 = new AnalyzeNoise("adc32", fft_file, fft_tmp, 511);
         fPN32 = new PlotNoise("adc32");
      }

      bsc->cd(); // select correct ROOT directory
   }

   ~BscModule()
   {
      if (fTrace)
         printf("BscModule::dtor!\n");
      if (fH) {
         delete fH;
         fH = NULL;
      }
   }

   void BeginRun(TARunInfo* runinfo)
   {
      if (fTrace)
         printf("BscModule::BeginRun, run %d, file %s\n", runinfo->fRunNo, runinfo->fFileName.c_str());
      //time_t run_start_time = runinfo->fOdb->odbReadUint32("/Runinfo/Start time binary", 0, 0);
      //printf("ODB Run start time: %d: %s", (int)run_start_time, ctime(&run_start_time));
      fCounter = 0;
      runinfo->fRoot->fOutputFile->cd(); // select correct ROOT directory
   }

   void EndRun(TARunInfo* runinfo)
   {
      if (fTrace)
         printf("BscModule::EndRun, run %d, events %d\n", runinfo->fRunNo, fCounter);
      if (fAN16)
         fAN16->Finish();
      if(fAN16AWB.size())
         {
            for(int awb=0; awb<32; ++awb)
            {
               fAN16AWB.at(awb).Finish();
            }
         }
      if (fAN32)
         fAN32->Finish();
      //time_t run_stop_time = runinfo->fOdb->odbReadUint32("/Runinfo/Stop time binary", 0, 0);
      //printf("ODB Run stop time: %d: %s", (int)run_stop_time, ctime(&run_stop_time));
   }

   void PauseRun(TARunInfo* runinfo)
   {
      if (fTrace)
         printf("BscModule::PauseRun, run %d\n", runinfo->fRunNo);
   }

   void ResumeRun(TARunInfo* runinfo)
   {
      if (fTrace)
         printf("BscModule::ResumeRun, run %d\n", runinfo->fRunNo);
   }

   void FilterWaveform(Alpha16Channel* hit)
   {
      std::vector<int> v = hit->adc_samples;
      int n=v.size();
      int k=20;
      for (int i=0; i<n; i++) {
         double sum0 = 0;
         double sum1 = 0;
         for (int j=i-k; j<i+k; j++) {
            if (j>=0 && j<n) {
               sum0 += 1;
               sum1 += v[j];
            }
         }
         hit->adc_samples[i] = sum1/sum0;
      }
   }

   void InvertWaveform(Alpha16Channel* hit)
   {
      int n = hit->adc_samples.size();
      for (int i=0; i<n; i++) {
         hit->adc_samples[i] = -hit->adc_samples[i];
      }
   }

   bool fft_first_adc16 = true;
   bool fft_first_adc32 = true;

   void AnalyzeHit(const TARunInfo* runinfo, const Alpha16Channel* hit, std::vector<AgBscAdcHit>* flow_hits)
   {
      //char xname[256];
      //char xtitle[256];
      //sprintf(xname, "m%02d_c%02d_w%03d", hit->adc_module, hit->adc_chan, hit->tpc_wire);
      //sprintf(xtitle, "AW Waveform ADC module %d, channel %d, tpc wire %d", hit->adc_module, hit->adc_chan, hit->tpc_wire);

      int ibar = hit->bsc_bar;

#if 0
      printf("hit: bank [%s], adc_module %d, adc_chan %d, preamp_pos %d, preamp_wire %d, tpc_wire %d, bsc_bar %d, nbins %d+%d, seqno %d, region %d\n",
             hit->bank.c_str(),
             hit->adc_module,
             hit->adc_chan,
             hit->preamp_pos,
             hit->preamp_wire,
             hit->tpc_wire,
             hit->bsc_bar,
             hit->first_bin,
             (int)hit->adc_samples.size(),
             ibar,
             r);
#endif

      if (ibar < 0)
         return;

      ////// Plot waveforms
      
      //while ((int)fHC.size() <= i) {
      //   fHC.push_back(NULL);
      //}

      //if (fHC[i] == NULL){
      //   fHC[i] = new A16ChanHistograms(xname, xtitle, hit->adc_samples.size());
      //}

      int max_fft_count = 30; // FFT is slow, limit number of events analyzed

      if (fft_first_adc16) {
         fft_first_adc16 = false;
         if (fAN16 && fAN16->fCount < max_fft_count) {
            fAN16->AddWaveform(hit->adc_samples);
            fAN16AWB.at(hit->preamp_pos).AddWaveform(hit->adc_samples);
            if (fPN16)
               fPN16->Plot(fAN16);
         }                  
      }
      
      int max_fft_awb_count = 30000;
      if (fAN16AWB.size() > 0 )
         {
            if ( fAN16AWB.at(hit->preamp_pos).fCount < max_fft_awb_count) 
               {
                  fAN16AWB.at(hit->preamp_pos).AddWaveform(hit->adc_samples);
               }                  
         }

      // analyze baseline
      
      int is_start = 0;
      int is_baseline = 100;

      double bmean, brms, bmin, bmax;

      compute_mean_rms(hit->adc_samples.data(), is_start, is_baseline, &bmean, &brms, &bmin, &bmax);

      double brange = bmax-bmin;
      
      double wmin = hit->adc_samples[0];
      double wmax = hit->adc_samples[0];

      for (unsigned s=0; s<hit->adc_samples.size(); s++) {
         double a = hit->adc_samples[s];
         if (a < wmin)
            wmin = a;
         if (a > wmax)
            wmax = a;
      }

      double wrange = wmax-wmin;
      
      bool have_baseline = false;
      bool have_pulse = false;
      bool have_hit = false;

      double ph = 0;
      double le = 0;
      double hit_time = 0;
      double hit_amp = 0;

      double cut_brms = 2000; // 300;

      if (brms > cut_brms) {
         if (fH) {
            fH->fHbaselineMean->Fill(0);

            if (brms < MAX_BSC_BRMS)
               fH->fHbaselineRms->Fill(brms);
            else
               fH->fHbaselineRms->Fill(MAX_BSC_BRMS-1);
            
            if (brange < MAX_AW_BRANGE)
               fH->fHbaselineRange->Fill(brange);
            else
               fH->fHbaselineRange->Fill(MAX_AW_BRANGE-1);

            fH->fHbaselineBadRmsBscMap->Fill(ibar);
         }
      } else {
         have_baseline = true;

         if (fH) {
            fH->fHbaselineMean->Fill(bmean);

            if (brms < MAX_BSC_BRMS)
               fH->fHbaselineRms->Fill(brms);
            else
               fH->fHbaselineRms->Fill(MAX_BSC_BRMS-1);

            if (brange < MAX_AW_BRANGE) {
               fH->fHbaselineRange->Fill(brange);
            } else {
               fH->fHbaselineRange->Fill(MAX_AW_BRANGE-1);
            }

            fH->fHbaselineMeanBscMap->Fill(ibar, bmean);
            fH->fHbaselineRmsBscMap->Fill(ibar, brms);
            fH->fHbaselineRangeBscMap->Fill(ibar, brange);

            if (wrange < MAX_AW_AMP) {
               fH->fHrange->Fill(wrange);
               fH->fHrangeBscMap->Fill(ibar, wrange);
            } else {
               fH->fHrange->Fill(MAX_AW_AMP-1);
               fH->fHrangeBscMap->Fill(ibar, MAX_AW_AMP-1);
            }
         }
      
         ph = wmax - bmean;

         if (ph > 500)
            fH->fBscMapPh500->Fill(ibar);
         if (ph > 1000)
            fH->fBscMapPh1000->Fill(ibar);
         if (ph > 5000)
            fH->fBscMapPh5000->Fill(ibar);
         if (ph > 10000)
            fH->fBscMapPh10000->Fill(ibar);

         double cfd_thr = 0.75*ph;

         if (wmax == 32764.0) {
            ph = MAX_AW_AMP-1;
            cfd_thr = 10000;
         }

         //if (ph > MAX_AW_AMP) {
         //   ph = MAX_AW_AMP-1;
         //}

         if (fH) {
            fH->fHph->Fill(ph);
            fH->fHped->Fill(ph);
         }

         double ph_hit_thr = 1000;

         if (brms > 300) {
            ph_hit_thr = 4000;
         }

#if 0
         if (0) {
         } else if (runinfo->fRunNo < 2181) {
            ph_hit_thr_adc16 =  2000000; // adc16 not connected
            ph_hit_thr_adc32 =  2500;
         } else if (runinfo->fRunNo < 2202) {
            ph_hit_thr_adc16 =  2000000; // adc16 not connected
            ph_hit_thr_adc32 =  750;
         } else if (runinfo->fRunNo < 902340) {
            ph_hit_thr_adc16 =  1000;
            ph_hit_thr_adc32 =  750;
         } else if (runinfo->fRunNo < 999999) {
            ph_hit_thr_adc16 =  1000;
            ph_hit_thr_adc32 =  99999; // not used
         }
#endif

         double time_bin = 1000.0/100.0; // 100 MHz ADC

         double pulse_time_middle = 125; // ADC time bins

#if 0
         if (runinfo->fRunNo < 458) {
            pulse_time_middle_adc16 = 160; // ADC time bins
            pulse_time_middle_adc32 = 160; // ADC time bins
         } else if (runinfo->fRunNo < 1244) {
            pulse_time_middle_adc16 = 165; // ADC time bins
            pulse_time_middle_adc32 = 165; // ADC time bins
         } else if (runinfo->fRunNo < 2164) {
            pulse_time_middle_adc16 = 147; // ADC time bins
            pulse_time_middle_adc32 = 150; // ADC time bins
         } else if (runinfo->fRunNo < 9999) {
            pulse_time_middle_adc16 = 147; // ADC time bins
            pulse_time_middle_adc32 = 135; // ADC time bins
         }
#endif

         double time_pc = 1000.0; // place PC drift times at 1000ns.

         double time_offset = time_pc - pulse_time_middle*time_bin;

         double adc_offset = 0.0;
         
         double adc_gain = 1.0;

         if (ph > ph_hit_thr) {
            have_pulse = true;

            if (fH) {
               fH->fHphHit->Fill(ph);
               fH->fHphHitBscMap->Fill(ibar);
            }

            //int le = led(w, b, -1.0, cfd_thr);
            le = find_pulse_time(hit->adc_samples.data(), hit->adc_samples.size(), bmean, 1.0, cfd_thr);

            if (fH) {
               fH->fHle->Fill(le);
            }
            
            hit_time = le * time_bin + time_offset;
            hit_amp = ph * adc_gain + adc_offset;

            fH->fHbscHitTime0->Fill(hit_time);
            
            if (1 || (hit_time > 700 && hit_time < 6000)) {
               have_hit = true;

               fH->fHbscHitTime->Fill(hit_time);
               fH->fHbscHitAmp->Fill(hit_amp);
               fH->fHbscHitMap->Fill(ibar);

               AgBscAdcHit h;
               h.adc_module = hit->adc_module;
               h.adc_chan = hit->adc_chan;
               h.bar  = ibar;
               h.time = hit_time;
               h.amp = hit_amp;

               flow_hits->push_back(h);
            }
         }
      }
      
      if (fFlags->fPrint) {
         printf("bar %3d: baseline mean %8.1f, rms %4.1f, range %8.1f %6.1f, pulse %6.1f, le %6.1f, time %5.0f, amp %6.1f, flags: baseline %d, pulse %d, hit %d\n", ibar, bmean, brms, wmin, wmax, ph, le, hit_time, hit_amp, have_baseline, have_pulse, have_hit);
      }
   }
   
   TAFlowEvent* AnalyzeFlowEvent(TARunInfo* runinfo, TAFlags* flags, TAFlowEvent* flow)
   {
      //if (fTrace)
      //printf("BscModule::Analyze, run %d, event serno %d, id 0x%04x, data size %d\n", runinfo->fRunNo, event->serial_number, (int)event->event_id, event->data_size);

      AgEventFlow *ef = flow->Find<AgEventFlow>();

      if (!ef || !ef->fEvent)
         return flow;

      Alpha16Event* e = ef->fEvent->a16;

      if (!e) {
         return flow;
      }

      if (1) {
         printf("Have ADC event:  ");
         e->Print();
         printf("\n");
      }

      AgBscAdcHitsFlow* flow_hits = new AgBscAdcHitsFlow(flow);
      flow = flow_hits;

      fft_first_adc16 = true;
      fft_first_adc32 = true;

      for (unsigned i=0; i<e->hits.size(); i++) {
         if (fFlags->fFilterWaveform) {
            FilterWaveform(e->hits[i]);
         }
         if (fFlags->fInvertWaveform) {
            InvertWaveform(e->hits[i]);
         }
         AnalyzeHit(runinfo, e->hits[i], &flow_hits->fBscAdcHits);
      }

      //*flags |= TAFlag_DISPLAY;

      fCounter++;

      return flow;
   }

   void AnalyzeSpecialEvent(TARunInfo* runinfo, TMEvent* event)
   {
      if (fTrace)
         printf("AnalyzeSpecialEvent, run %d, event serno %d, id 0x%04x, data size %d\n", runinfo->fRunNo, event->serial_number, (int)event->event_id, event->data_size);
   }
};

class BscModuleFactory: public TAFactory
{
public:
   A16Flags fFlags;
   
public:
   void Init(const std::vector<std::string> &args)
   {
      printf("BscModuleFactory::Init!\n");
      
      for (unsigned i=0; i<args.size(); i++) {
         if (args[i] == "--bscprint")
            fFlags.fPrint = true;
         if (args[i] == "--bscfft")
            fFlags.fFft = true;
         if (args[i] == "--bscfwf") {
            fFlags.fFilterWaveform = true;
            i++;
         }
         if (args[i] == "--bscinv") {
            fFlags.fInvertWaveform = true;
            i++;
         }
      }
      
      TARootHelper::fgDir->cd(); // select correct ROOT directory
   }
   
   void Finish()
   {
      printf("BscModuleFactory::Finish!\n");
   }
   
   TARunObject* NewRunObject(TARunInfo* runinfo)
   {
      printf("BscModuleFactory::NewRunObject, run %d, file %s\n", runinfo->fRunNo, runinfo->fFileName.c_str());
      return new BscModule(runinfo, &fFlags);
   }

public:
};

static TARegister tar(new BscModuleFactory);

/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */
