//
// MIDAS analyzer example 2: ROOT analyzer
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

#include "Alpha16.h"
#include "Unpack.h"
#include "AgFlow.h"
#include "ko_limits.h"

#define MAX_AW_BRMS 600
#define MAX_AW_BRANGE 2000
// zoom on the ADC pedestal
#define MAX_AW_PED 2000

#include <iostream>
#include "TVirtualFFT.h"

#define DELETE(x) if (x) { delete (x); (x) = NULL; }

#define MEMZERO(p) memset((p), 0, sizeof(p))

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
   TProfile* fHbaselineMeanAwMap;
   TProfile* fHbaselineRmsAwMap;

   TH1D* fHrange;

   TH1D* fHph;
   TH1D* fHph_adc16;
   TH1D* fHph_adc32;

   TH1D* fHped;
   TH1D* fHped_adc16;
   TH1D* fHped_adc32;

   TH1D* fHphHit;
   TH1D* fHphHit_adc16;
   TH1D* fHphHit_adc32;

   TH1D* fHphHitAwMap;

   TH1D* fHle;
   TH1D* fHle_adc16;
   TH1D* fHle_adc32;

   TH1D* fHawHitTime0;
   TH1D* fHawHitTime;
   TH1D* fHawHitAmp;
   TH1D* fHawHitMap;

   TH1D* fHphCal;
   TH1D* fHleCal;
   TH1D* fHoccCal;
   TProfile* fHphVsChanCal;
   TProfile* fHleVsChanCal;

   PlotHistograms() // ctor
   {
      TDirectory *dir = gDirectory->mkdir("summary");
      dir->cd();

      fHbaselineMean = new TH1D("adc_baseline_mean", "waveform baseline mean; ADC counts", 100, -2000, 2000);
      fHbaselineRms = new TH1D("adc_baseline_rms", "waveform baseline rms; ADC counts", 100, 0, MAX_AW_BRMS);
      fHbaselineRange = new TH1D("adc_baseline_range", "waveform baseline range; ADC counts", 100, 0, MAX_AW_BRANGE);

      fHbaselineMeanAwMap = new TProfile("adc_baseline_mean_vs_aw", "waveform baseline mean; TPC wire number; ADC counts", NUM_AW, -0.5, NUM_AW-0.5);
      fHbaselineRmsAwMap = new TProfile("adc_baseline_rms_vs_aw", "waveform baseline rms vs wire number; TPC wire number; ADC counts", NUM_AW, -0.5, NUM_AW-0.5);
      fHbaselineRmsAwMap->SetMinimum(0);

      fHrange = new TH1D("adc_range", "waveform range, max-min; ADC counts", 100, 0, MAX_AW_AMP);

      fHph = new TH1D("adc_pulse_height", "waveform pulse height; ADC counts", 100, 0, MAX_AW_AMP);

      fHph_adc16 = new TH1D("adc16_pulse_height", "adc16 pulse height (100MHz); ADC counts", 100, 0, MAX_AW_AMP);
      fHph_adc32 = new TH1D("adc32_pulse_height", "adc32 pulse height (62.5MHz); ADC counts", 100, 0, MAX_AW_AMP);
      
      fHped = new TH1D("adc_pedstal_pulse_height", "waveform pulse height, zoom on the pedestal; ADC counts", 100, 0, MAX_AW_PED);

      fHped_adc16 = new TH1D("adc16_pedestal_pulse_height", "adc16 pulse height (100MHz), zoom on the pedestal; ADC counts", 100, 0, MAX_AW_PED);
      fHped_adc32 = new TH1D("adc32_pedestal_pulse_height", "adc32 pulse height (62.5MHz), zoom on the pedestal; ADC counts", 100, 0, MAX_AW_PED);
      
      fHphHit = new TH1D("adc_pulse_height_hit", "pulse height, after ADC cut; ADC counts", 100, 0, MAX_AW_AMP);
      
      fHphHit_adc16 = new TH1D("adc16_pulse_height_hit", "adc16 pulse height (100MHz), after ADC cut; ADC counts", 100, 0, MAX_AW_AMP);
      fHphHit_adc32 = new TH1D("adc32_pulse_height_hit", "adc32 pulse height (62.5MHz), after ADC cut; ADC counts", 100, 0, MAX_AW_AMP);
      
      fHphHitAwMap = new TH1D("adc_pulse_height_hit_vs_aw", "TPC wire map, after ADC cut; TPC wire number", NUM_AW, -0.5, NUM_AW-0.5);

      fHle = new TH1D("adc_pulse_time", "waveform pulse time, after ADC cut; ADC time bins", 200, 0, 1000);

      fHle_adc16 = new TH1D("adc16_pulse_time", "adc16 pulse time (100MHz); ADC time bins", 200, 0, 1000);
      fHle_adc32 = new TH1D("adc32_pulse_time", "adc32 pulse time (62.5MHz); ADC time bins", 200, 0, 1000);

      fHawHitTime0 = new TH1D("aw_hit_time0", "hit time, after ADC cut; time, ns", 200, 0, MAX_TIME);
      fHawHitTime = new TH1D("aw_hit_time", "hit time, after ADC and time cut; time, ns", 200, 0, MAX_TIME);
      fHawHitAmp = new TH1D("aw_hit_amp", "hit amplitude, after ADC and time cut; pulse height, ADC c.u.", 200, 0, MAX_AW_AMP);
      fHawHitMap = new TH1D("aw_hit_map", "hit map, after ADC and time cut; TPC wire", NUM_AW, -0.5, NUM_AW-0.5);

      fHphCal = new TH1D("adccal_pulse_height", "pulse_height_cal; ADC counts", 100, 0, MAX_AW_AMP);
      fHleCal = new TH1D("adccal_pulse_time", "pulse_time_cal; ADC time bins", 100, 100, 200);
      fHoccCal = new TH1D("adccal_channel_occupancy", "channel_occupancy_cal; TPC wire number", NUM_AW, -0.5, NUM_AW-0.5);
      fHoccCal->SetMinimum(0);
      fHphVsChanCal = new TProfile("adccal_pulse_height_vs_wire", "pulse_height_vs_chan_cal; TPC wire number; ADC counts", NUM_AW, -0.5, NUM_AW-0.5);
      fHleVsChanCal = new TProfile("adccal_pulse_time_vs_wire", "pulse_time_vs_chan_cal; TPC wire number; ADC time bins", NUM_AW, -0.5, NUM_AW-0.5);
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
      printf("Create name [%s] title [%s] with %d bins\n", xname, xtitle, nbins);

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

struct PlotA16
{
   TCanvas* fCanvas = NULL;
   std::vector<TH1D*> fH;

   int fModule = 0;
   int fSection = 0;

   PlotA16(TCanvas* c, int module, int section) // ctor
   {
      if (!c) {
         char title[256];
         sprintf(title, "ADC mod%d section %d", module, section);
         c = new TCanvas(title, title, 900, 650);
         //c = new TCanvas("ALPHA16 ADC", "ALPHA16 ADC", 900, 650);
         if (!(c->GetShowEventStatus()))
            c->ToggleEventStatus();
         if (!(c->GetShowToolBar()))
            c->ToggleToolBar();
      }

      fCanvas = c;

      fCanvas->cd();
      fCanvas->Divide(4,4);

      fModule = module;
      fSection = section;

      for (int i=0; i<16; i++) {
         fH.push_back(NULL);
      }
   }

   ~PlotA16()
   {
      if (fCanvas) {
         delete fCanvas;
         fCanvas = NULL;
      }
   }

   void DrawChannel(int i, const Alpha16Channel* c)
   {
      assert(i>=0 && i<16);

      fCanvas->cd(i+1);

      int color = 1;
      
      if (!fH[i]) {
         char name[256];
         sprintf(name, "mod%dch%d", c->adc_module, c->adc_chan);
         
         fH[i] = new TH1D(name, name, c->adc_samples.size(), 0, c->adc_samples.size());
         
         fH[i]->Draw();
         fH[i]->SetMinimum(-(1<<15));
         fH[i]->SetMaximum(1<<15);
         //fH[i]->SetMinimum(-2000);
         //fH[i]->SetMaximum(2000);
         fH[i]->GetYaxis()->SetLabelSize(0.10);
         //fH[i]->SetLineColor(color);
      }

      fH[i]->SetLineColor(color);

      for (unsigned s=0; s<c->adc_samples.size(); s++) {
         int v = c->adc_samples[s];
         //printf("bin %d, v %d\n", i, v);
         fH[i]->SetBinContent(s+1, v);
      }
   }

   void ClearChannel(int i)
   {
      if (fH[i]) {
         fH[i]->SetLineColor(2);
         //fH[i]->Clear();
      }
   }
   
   void Draw(const Alpha16Event* e)
   {
      // colors:
      // 0 = white
      // 1 = black
      // 2 = red
      // 3 = green
      // 4 = blue

      for (int i=0; i<16; i++) {
         ClearChannel(i);
      }

      for (unsigned i=0; i<e->hits.size(); i++) {
         Alpha16Channel* c = e->hits[i];
         if (c->adc_module != fModule)
            continue;
         if (fSection == 0) { // 16 onboard ADCs
            if (c->adc_chan < 0 || c->adc_chan >= 16)
               continue;
            DrawChannel(c->adc_chan, c);
         } else if (fSection == 1) {
            if (c->adc_chan < 16 || c->adc_chan >= 32)
               continue;
            DrawChannel(c->adc_chan-16, c);
         } else if (fSection == 2) {
            if (c->adc_chan < 32 || c->adc_chan >= 48)
               continue;
            DrawChannel(c->adc_chan-32, c);
         }
      }

      fCanvas->Modified();
      fCanvas->Draw();
      fCanvas->Update();
   }
};

#if 0
class AlphaTpcX
{
public:
   int fRunNo;

   int fCountEarlyBad;
   int fCountGood;
   int fCountBad;

   TDirectory *dnoise = NULL;
   TCanvas *c;
   TH1D *hwf0, *hwf1, *hfftsum0, *hfftsum1, *hbase0, *hbase1, *hRMS0, *hRMS1;

   int entries = 0;

      dnoise = gDirectory->mkdir("noise");
      dnoise->cd();
      c = new TCanvas("cnoise","noise analysis",800,1200);
      c->Divide(1,6);
      hwf0 = new TH1D("hwf0","waveform 0",701,0,701);
      c->cd(1);
      hwf0->Draw();
      c->GetPad(2)->SetLogy();

      hwf1 = new TH1D("hwf1","waveform 1",701,0,701);
      c->cd(3);
      hwf1->Draw();
      c->GetPad(4)->SetLogy();

      hfftsum0 = new TH1D("hfftsum0","FFT sum",701,0,701);
      c->cd(5);
      hfftsum0->Draw();
      c->GetPad(5)->SetLogy();

      hfftsum1 = new TH1D("hfftsum1","FFT sum",701,0,701);
      c->cd(6);
      hfftsum1->Draw();
      c->GetPad(6)->SetLogy();
   }

   void EndRun()
   {
      printf("AlphaTpcX::EndRun: early bad events: %d, good events: %d, bad events: %d, total %d events\n", fCountEarlyBad, fCountGood, fCountBad, fCountEarlyBad + fCountGood + fCountBad);
      hfftsum0->Scale(1./double(hbase0->GetEntries()*sqrt(701.)));
      hfftsum1->Scale(1./double(hbase1->GetEntries()*sqrt(701.)));
   }

   void AnalyzeNoise(Waveform *w, short i){
      TH1D *h, *hRMS, *hbase, *hfftsum;
      short index;
      switch(i){
      case 96: h = hwf0; hbase = hbase0; hRMS = hRMS0; hfftsum = hfftsum0; index = 2; break;
      case 127: h = hwf1; hbase = hbase1; hRMS = hRMS1; hfftsum = hfftsum1; index = 4; break;
      default: std::cerr << "Something's wrong!" << std::endl;
      }
      for(int b = 0; b < w->nsamples; b++)
         h->SetBinContent(b+1, w->samples[b]);

      double b, brms;
      b = baseline(w, 0, w->nsamples, NULL, &brms);
      hbase->Fill(b);
      hRMS->Fill(brms);

      TH1 *hfft = 0;
      TVirtualFFT::SetTransform(0);
      hfft = h->FFT(hfft, "MAG");
      c->cd(index);
      hfft->Draw();
      hfftsum->Add(hfft);
      entries++;
      for(int p = 1; p <= 6; p++)
         c->GetPad(p)->Modified();
      c->Update();
   }
};
#endif

class A16Flags
{
public:
   bool fPrint = false;
   bool fPlotWF = false;
   bool fDoPlotAll = false;
   bool fExportWaveforms = false;
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

class A16Module: public TARunObject
{
public:
   A16Flags* fFlags = NULL;
   int fCounter = 0;
   bool fTrace = false;

   PlotHistograms* fH = NULL;
   std::vector<A16ChanHistograms*> fHC;
   std::vector<PlotA16*> fPlotA16;

   TDirectory *dnoise = NULL;

public:
   A16Module(TARunInfo* runinfo, A16Flags* f)
      : TARunObject(runinfo)
   {
      fTrace = false;
      if (fTrace)
         printf("A16Module::ctor!\n");
      fFlags = f;

      runinfo->fRoot->fOutputFile->cd();
      TDirectory* aw = gDirectory->mkdir("aw");
      aw->cd(); // select correct ROOT directory

      fH = new PlotHistograms();

      dnoise = gDirectory->mkdir("noise");
   }

   ~A16Module()
   {
      if (fTrace)
         printf("A16Module::dtor!\n");
      if (fH) {
         delete fH;
         fH = NULL;
      }
   }

   void BeginRun(TARunInfo* runinfo)
   {
      printf("A16Module::BeginRun, run %d, file %s\n", runinfo->fRunNo, runinfo->fFileName.c_str());
      time_t run_start_time = runinfo->fOdb->odbReadUint32("/Runinfo/Start time binary", 0, 0);
      printf("ODB Run start time: %d: %s", (int)run_start_time, ctime(&run_start_time));
      fCounter = 0;
      runinfo->fRoot->fOutputFile->cd(); // select correct ROOT directory

      if (0) {
         fPlotA16.push_back(new PlotA16(NULL, 1, 0));
         fPlotA16.push_back(new PlotA16(NULL, 2, 0));
         fPlotA16.push_back(new PlotA16(NULL, 3, 0));
         fPlotA16.push_back(new PlotA16(NULL, 4, 0));
         fPlotA16.push_back(new PlotA16(NULL, 5, 0));
         fPlotA16.push_back(new PlotA16(NULL, 5, 1));
         fPlotA16.push_back(new PlotA16(NULL, 5, 2));
         fPlotA16.push_back(new PlotA16(NULL, 6, 0));
         fPlotA16.push_back(new PlotA16(NULL, 6, 1));
         fPlotA16.push_back(new PlotA16(NULL, 6, 2));
         fPlotA16.push_back(new PlotA16(NULL, 7, 0));
         fPlotA16.push_back(new PlotA16(NULL, 8, 0));
         
         fPlotA16.push_back(new PlotA16(NULL, 9, 0));
         fPlotA16.push_back(new PlotA16(NULL, 10, 0));
         fPlotA16.push_back(new PlotA16(NULL, 11, 0));
         fPlotA16.push_back(new PlotA16(NULL, 12, 0));
         fPlotA16.push_back(new PlotA16(NULL, 13, 0));
         fPlotA16.push_back(new PlotA16(NULL, 14, 0));
         fPlotA16.push_back(new PlotA16(NULL, 15, 0));
         fPlotA16.push_back(new PlotA16(NULL, 16, 0));
         fPlotA16.push_back(new PlotA16(NULL, 17, 0));
         fPlotA16.push_back(new PlotA16(NULL, 18, 0));
      }
   }

   void EndRun(TARunInfo* runinfo)
   {
      printf("A16Module::EndRun, run %d, events %d\n", runinfo->fRunNo, fCounter);
      time_t run_stop_time = runinfo->fOdb->odbReadUint32("/Runinfo/Stop time binary", 0, 0);
      printf("ODB Run stop time: %d: %s", (int)run_stop_time, ctime(&run_stop_time));
   }

   void PauseRun(TARunInfo* runinfo)
   {
      printf("A16Module::PauseRun, run %d\n", runinfo->fRunNo);
   }

   void ResumeRun(TARunInfo* runinfo)
   {
      printf("A16Module::ResumeRun, run %d\n", runinfo->fRunNo);
   }

   void AnalyzeHit(const TARunInfo* runinfo, const Alpha16Channel* hit, std::vector<AgAwHit>* flow_hits)
   {
      char xname[256];
      char xtitle[256];
      sprintf(xname, "m%02d_c%02d_w%03d", hit->adc_module, hit->adc_chan, hit->tpc_wire);
      sprintf(xtitle, "AW Waveform ADC module %d, channel %d, tpc wire %d", hit->adc_module, hit->adc_chan, hit->tpc_wire);

      if (fFlags->fExportWaveforms) {
         TDirectory* dir = runinfo->fRoot->fgDir;
         const char* dirname = "AW waveforms";

         if (!dir->cd(dirname)) {
            TDirectory* awdir = dir->mkdir(dirname);
            awdir->cd();
         }

         dir = dir->CurrentDirectory();

         std::string wname = std::string(xname) + "_waveform";
         std::string wtitle = std::string(xtitle) + " current waveform";

         TH1D* hwf = (TH1D*)dir->FindObject(wname.c_str());
         if (!hwf) {
            hwf = new TH1D(wname.c_str(), wtitle.c_str(), hit->adc_samples.size(), 0, hit->adc_samples.size());
         }
         
         for (unsigned i=0; i< hit->adc_samples.size(); i++) {
            hwf->SetBinContent(i+1, hit->adc_samples[i]);
         }
      }

      int i = hit->tpc_wire;
      int r = 1;

      if (hit->adc_module == 5 && hit->adc_chan >= 16) {
         r = 2;
      }

      if (hit->adc_module == 6 && hit->adc_chan >= 16) {
         r = 3;
      }

#if 0      
      printf("hit: bank [%s], adc_module %d, adc_chan %d, preamp_pos %d, preamp_wire %d, tpc_wire %d. nbins %d+%d, seqno %d, region %d\n",
             hit->bank.c_str(),
             hit->adc_module,
             hit->adc_chan,
             hit->preamp_pos,
             hit->preamp_wire,
             hit->tpc_wire,
             hit->first_bin,
             (int)hit->adc_samples.size(),
             i,
             r);
#endif

      if (i < 0)
         return;

      int iwire = i;

      bool is_aw = true;
      bool is_adc16 = (hit->adc_chan < 16);
      bool is_adc32 = (hit->adc_chan >= 16);

      ////// Plot waveforms
      
      while ((int)fHC.size() <= i) {
         fHC.push_back(NULL);
      }

      if (fHC[i] == NULL){
         fHC[i] = new A16ChanHistograms(xname, xtitle, hit->adc_samples.size());
      }

#if 0
      if(i == 96 || i == 127)
         AnalyzeNoise(w, i);
#endif

      // analyze baseline
      
      int is_start = 0;
      int is_baseline = 100;
      int calStart = 160;
      int calEnd = 200;

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

      double cut_brms = 500;

      if (brms > cut_brms) {
         if (fH) {
            fH->fHbaselineMean->Fill(0);

            if (brms < MAX_AW_BRMS)
               fH->fHbaselineRms->Fill(brms);
            else
               fH->fHbaselineRms->Fill(MAX_AW_BRMS-1);
            
            if (brange < MAX_AW_BRANGE)
               fH->fHbaselineRange->Fill(brange);
            else
               fH->fHbaselineRange->Fill(MAX_AW_BRANGE-1);
         }
      } else {
         have_baseline = true;

         if (fH) {
            fH->fHbaselineMean->Fill(bmean);

            if (brms < MAX_AW_BRMS)
               fH->fHbaselineRms->Fill(brms);
            else
               fH->fHbaselineRms->Fill(MAX_AW_BRMS-1);

            if (brange < MAX_AW_BRANGE)
               fH->fHbaselineRange->Fill(brange);
            else
               fH->fHbaselineRange->Fill(MAX_AW_BRANGE-1);

            if (is_aw) {
               fH->fHbaselineMeanAwMap->Fill(iwire, bmean);
               fH->fHbaselineRmsAwMap->Fill(iwire, brms);
            }

            if (wrange < MAX_AW_AMP)
               fH->fHrange->Fill(wrange);
            else
               fH->fHrange->Fill(MAX_AW_AMP-1);
         }
      
         ph = bmean - wmin;

         double cfd_thr = ph/2.0;

         if (wmin == -32768.0) {
            ph = MAX_AW_AMP-1;
            cfd_thr = 10000;
         }

         if (ph > MAX_AW_AMP) {
            ph = MAX_AW_AMP-1;
         }

         if (fH) {
            fH->fHph->Fill(ph);
            fH->fHped->Fill(ph);
         
            if (is_adc16)
               fH->fHph_adc16->Fill(ph);
            else if (is_adc32)
               fH->fHph_adc32->Fill(ph);

            if (is_adc16)
               fH->fHped_adc16->Fill(ph);
            else if (is_adc32)
               fH->fHped_adc32->Fill(ph);
         }

         double ph_hit_thr = 250;

         if (is_adc16)
            ph_hit_thr = 1200;
         else if (is_adc32)
            ph_hit_thr = 600;

         double time_bin = 0;
         double time_offset = 0;

         if (is_adc16)
            time_bin = 1000.0/100.0; // 100 MHz ADC
         else if (is_adc32)
            time_bin = 1000.0/62.5; // 62.5 MHz ADC

         double pulse_time_middle_adc16 = 160; // ADC time bins

         double time_pc = 1000.0; // place PC drift times at 1000ns.

         if (is_adc16)
            time_offset = time_pc - pulse_time_middle_adc16*time_bin;
         else if (is_adc32)
            time_offset = time_pc - 2410;

         double adc_gain = 1.0;
         double adc_offset = 0.0;
         
         if (is_adc16)
            adc_gain = 1.0;
         else if (r == 2)
            adc_gain = 4.0;
         else if (r == 3)
            adc_gain = 2.0;

         if (ph > ph_hit_thr) {
            have_pulse = true;

            if (fH) {
               fH->fHphHit->Fill(ph);

               if (is_aw)
                  fH->fHphHitAwMap->Fill(iwire);
               
               if (is_adc16)
                  fH->fHphHit_adc16->Fill(ph);
               else if (is_adc32)
                  fH->fHphHit_adc32->Fill(ph);
            }

            //int le = led(w, b, -1.0, cfd_thr);
            le = find_pulse_time(hit->adc_samples.data(), hit->adc_samples.size(), bmean, -1.0, cfd_thr);

            if (fH) {
               fH->fHle->Fill(le);
            
               if (is_adc16)
                  fH->fHle_adc16->Fill(le);
               else if (is_adc32)
                  fH->fHle_adc32->Fill(le);
            }
            
            if (le > calStart && le < calEnd) {
               fH->fHleCal->Fill(le);
               fH->fHphCal->Fill(ph);
               if (is_aw) {
                  fH->fHoccCal->Fill(iwire);
                  fH->fHleVsChanCal->Fill(iwire, le);
                  fH->fHphVsChanCal->Fill(iwire, ph);
               }
            }
            
            //if (ph > 4000) {
            //nhits++;
            //printf("samples %d %d, ", e->waveform[i].size(), w->nsamples);
            //}

            hit_time = le * time_bin + time_offset;
            hit_amp = ph * adc_gain + adc_offset;

            fH->fHawHitTime0->Fill(hit_time);
            
            if (hit_time > 800 && hit_time < 6000) {
               have_hit = true;

               fH->fHawHitTime->Fill(hit_time);
               fH->fHawHitAmp->Fill(hit_amp);
               fH->fHawHitMap->Fill(iwire);

               AgAwHit h;
               h.adc_module = hit->adc_module;
               h.adc_chan = hit->adc_chan;
               h.wire = iwire;
               h.time = hit_time;
               h.amp = hit_amp;

               flow_hits->push_back(h);
            }
         }
      }
      
      if (fFlags->fPrint) {
         printf("wire %3d: baseline mean %8.1f, rms %4.1f, range %8.1f %6.1f, pulse %6.1f, le %6.1f, time %5.0f, amp %6.1f, flags: baseline %d, pulse %d, hit %d\n", i, bmean, brms, wmin, wmax, ph, le, hit_time, hit_amp, have_baseline, have_pulse, have_hit);
      }
   }
   
   TAFlowEvent* Analyze(TARunInfo* runinfo, TMEvent* event, TAFlags* flags, TAFlowEvent* flow)
   {
      if (fTrace)
         printf("A16Module::Analyze, run %d, event serno %d, id 0x%04x, data size %d\n", runinfo->fRunNo, event->serial_number, (int)event->event_id, event->data_size);

      AgEventFlow *ef = flow->Find<AgEventFlow>();

      if (!ef || !ef->fEvent)
         return flow;

      Alpha16Event* e = ef->fEvent->a16;

      if (!e) {
         return flow;
      }

      if (1) {
         printf("Have A16 event: ");
         e->Print();
         printf("\n");
      }

      for (unsigned i=0; i<fPlotA16.size(); i++) {
         fPlotA16[i]->Draw(e);
      }

      AgAwHitsFlow* flow_hits = new AgAwHitsFlow(flow);
      flow = flow_hits;

      for (unsigned i=0; i<e->hits.size(); i++) {
         AnalyzeHit(runinfo, e->hits[i], &flow_hits->fAwHits);
      }

      //*flags |= TAFlag_DISPLAY;

      fCounter++;

      return flow;
   }

   void AnalyzeSpecialEvent(TARunInfo* runinfo, TMEvent* event)
   {
      printf("AnalyzeSpecialEvent, run %d, event serno %d, id 0x%04x, data size %d\n", runinfo->fRunNo, event->serial_number, (int)event->event_id, event->data_size);
   }
};

class A16ModuleFactory: public TAFactory
{
public:
   A16Flags fFlags;
   
public:
   void Init(const std::vector<std::string> &args)
   {
      printf("A16ModuleFactory::Init!\n");
      
      fFlags.fPlotWF = false;
      
      for (unsigned i=0; i<args.size(); i++) {
         if (args[i] == "--adcprint")
            fFlags.fPrint = true;
         if (args[i] == "--wf")
            fFlags.fPlotWF = true;
         if (args[i] == "--wfall") {
            fFlags.fDoPlotAll = true;
            fFlags.fPlotWF = true;
         }
         if (args[i] == "--wfexport") {
            fFlags.fExportWaveforms = true;
         }
      }
      
      TARootHelper::fgDir->cd(); // select correct ROOT directory
   }
   
   void Finish()
   {
      printf("A16ModuleFactory::Finish!\n");
   }
   
   TARunObject* NewRunObject(TARunInfo* runinfo)
   {
      printf("A16ModuleFactory::NewRunObject, run %d, file %s\n", runinfo->fRunNo, runinfo->fFileName.c_str());
      return new A16Module(runinfo, &fFlags);
   }

public:
};

static TARegister tar(new A16ModuleFactory);

/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */
