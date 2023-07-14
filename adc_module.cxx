//
// adc_module.cxx - q.a. for the ADC data
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
#include "AgFlow.h"
#include "ko_limits.h"
#include "ncfm.h"

#include "wfsuppress_adc.h"

#define MAX_AW_BRMS 600
#define MAX_AW_BRANGE 2000
// zoom on the ADC pedestal
#define MAX_AW_PED 3000

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
   TH1D* fHbaselineBadRmsAwMap;
   TProfile* fHbaselineMeanAwMap;
   TProfile* fHbaselineRmsAwMap;
   TProfile* fHbaselineRangeAwMap;

   TH1D* fHrange;
   // TProfile* fHrangeAwMap;
   TH2D* fHrangeAwMap;

   TH1D* fAwMapPh500;
   TH1D* fAwMapPh1000;
   TH1D* fAwMapPh5000;
   TH1D* fAwMapPh10000;

   TH1D* fHph;
   TH1D* fHph_adc16;
   TH1D* fHph_adc32;

   TH1D* fHped;
   TH1D* fHped_adc16;
   TH1D* fHped_adc32;
   TH1D* fHped_adc32_p0;
   TH1D* fHped_adc32_pN;

   TH1D* fHphHit;
   TH1D* fHphHit_adc16;
   TH1D* fHphHit_adc32;

   TH1D* fHphHitAwMap;

   TH1D* fHle;
   TH1D* fHle_adc16;
   TH1D* fHle_adc32;

   TH1D* fHle_aw282 = NULL;

   TH1D* fHwi = NULL;
   TH1D* fHwi_adc16 = NULL;
   TH1D* fHwi_adc32 = NULL;

   TH1D* fHawHitTime0;
   TH1D* fHawHitTime;
   TH1D* fHawHitAmp;
   TH1D* fHawHitMap;

   PlotHistograms(TDirectory* dir) // ctor
   {
      dir->cd();

      fHbaselineMean = new TH1D("adc_baseline_mean", "waveform baseline mean, all wires; ADC counts", 100, -4000, 4000);
      fHbaselineRms = new TH1D("adc_baseline_rms", "waveform baseline rms, all wires; ADC counts", 100, 0, MAX_AW_BRMS);
      fHbaselineRange = new TH1D("adc_baseline_range", "waveform baseline range, all wires; ADC counts", 100, 0, MAX_AW_BRANGE);

      fHbaselineBadRmsAwMap = new TH1D("adc_baseline_bad_rms_vs_aw", "waveform with bad baseline vs wire number; TPC wire number", NUM_AW, -0.5, NUM_AW-0.5);

      fHbaselineMeanAwMap = new TProfile("adc_baseline_mean_vs_aw", "waveform baseline mean vs wire number; TPC wire number; ADC counts", NUM_AW, -0.5, NUM_AW-0.5);
      fHbaselineRmsAwMap = new TProfile("adc_baseline_rms_vs_aw", "waveform baseline rms vs wire number; TPC wire number; ADC counts", NUM_AW, -0.5, NUM_AW-0.5);
      fHbaselineRmsAwMap->SetMinimum(0);
      fHbaselineRangeAwMap = new TProfile("adc_baseline_range_vs_aw", "waveform baseline range vs wire number; TPC wire number; ADC counts", NUM_AW, -0.5, NUM_AW-0.5);

      fHrange = new TH1D("adc_range", "waveform range, max-min; ADC counts", 100, 0, MAX_AW_AMP);
      // fHrangeAwMap = new TProfile("adc_range_vs_aw", "waveform range, max-min vs wire number; TPC wire number; ADC counts", NUM_AW, -0.5, NUM_AW-0.5);
      fHrangeAwMap = new TH2D("adc_range_vs_aw", "waveform range, max-min vs wire number; TPC wire number; ADC counts", NUM_AW, -0.5, NUM_AW-0.5,1000,0,10000);

      fAwMapPh500 = new TH1D("aw_map_ph500_vs_aw", "waveforms with ph > 500 vs wire number; TPC wire number", NUM_AW, -0.5, NUM_AW-0.5);
      fAwMapPh1000 = new TH1D("aw_map_ph1000_vs_aw", "waveforms with ph > 1000 vs wire number; TPC wire number", NUM_AW, -0.5, NUM_AW-0.5);
      fAwMapPh5000 = new TH1D("aw_map_ph5000_vs_aw", "waveforms with ph > 5000 vs wire number; TPC wire number", NUM_AW, -0.5, NUM_AW-0.5);
      fAwMapPh10000 = new TH1D("aw_map_ph10000_vs_aw", "waveforms with ph > 10000 vs wire number; TPC wire number", NUM_AW, -0.5, NUM_AW-0.5);

      fHph = new TH1D("adc_pulse_height", "waveform pulse height; ADC counts", 100, 0, MAX_AW_AMP);

      fHph_adc16 = new TH1D("adc16_pulse_height", "adc16 pulse height (100MHz); ADC counts", 100, 0, MAX_AW_AMP);
      fHph_adc32 = new TH1D("adc32_pulse_height", "adc32 pulse height (62.5MHz); ADC counts", 100, 0, MAX_AW_AMP);
      
      fHped = new TH1D("adc_pedestal_pulse_height", "waveform pulse height, zoom on the pedestal; ADC counts", 100, 0, MAX_AW_PED);

      fHped_adc16 = new TH1D("adc16_pedestal_pulse_height", "adc16 pulse height (100MHz), zoom on the pedestal; ADC counts", 100, 0, MAX_AW_PED);
      fHped_adc32 = new TH1D("adc32_pedestal_pulse_height", "adc32 pulse height (62.5MHz), zoom on the pedestal; ADC counts", 100, 0, MAX_AW_PED);
      fHped_adc32_p0 = new TH1D("adc32_pedestal_pulse_height_preamp0", "adc32 pulse height (62.5MHz) in preamp 0, zoom on the pedestal; ADC counts", 100, 0, MAX_AW_PED);
      fHped_adc32_pN = new TH1D("adc32_pedestal_pulse_height_preampN", "adc32 pulse height (62.5MHz) in preamp 1..15, zoom on the pedestal; ADC counts", 100, 0, MAX_AW_PED);
      
      fHphHit = new TH1D("adc_pulse_height_hit", "pulse height, after ADC cut; ADC counts", 100, 0, MAX_AW_AMP);
      
      fHphHit_adc16 = new TH1D("adc16_pulse_height_hit", "adc16 pulse height (100MHz), after ADC cut; ADC counts", 100, 0, MAX_AW_AMP);
      fHphHit_adc32 = new TH1D("adc32_pulse_height_hit", "adc32 pulse height (62.5MHz), after ADC cut; ADC counts", 100, 0, MAX_AW_AMP);
      
      fHphHitAwMap = new TH1D("adc_pulse_height_hit_vs_aw", "TPC wire map, after ADC cut; TPC wire number", NUM_AW, -0.5, NUM_AW-0.5);

      fHle = new TH1D("adc_pulse_time", "waveform pulse time, after ADC cut; ADC time bins", 200, 0, 1000);

      fHle_adc16 = new TH1D("adc16_pulse_time", "adc16 pulse time (100MHz); ADC time bins", 200, 0, 1000);
      fHle_adc32 = new TH1D("adc32_pulse_time", "adc32 pulse time (62.5MHz); ADC time bins", 200, 0, 1000);

      fHle_aw282 = new TH1D("aw282_pulse_time", "aw282 pulse time (62.5MHz); ADC time bins", 100, 300-0.5, 400-0.5);

      fHwi = new TH1D("adc_pulse_width", "waveform pulse width, after ADC cut; ADC time bins", 100, 0, 100);

      fHwi_adc16 = new TH1D("adc16_pulse_width", "adc16 pulse width (100MHz); ADC time bins", 100, 0, 100);
      fHwi_adc32 = new TH1D("adc32_pulse_width", "adc32 pulse width (62.5MHz); ADC time bins", 100, 0, 100);

      fHawHitTime0 = new TH1D("aw_hit_time0", "hit time, after ADC cut; time, ns", 200, 0, MAX_TIME);
      fHawHitTime = new TH1D("aw_hit_time", "hit time, after ADC and time cut; time, ns", 200, 0, MAX_TIME);
      fHawHitAmp = new TH1D("aw_hit_amp", "hit amplitude, after ADC and time cut; pulse height, ADC c.u.", 200, 0, MAX_AW_AMP);
      fHawHitMap = new TH1D("aw_hit_map", "hit map, after ADC and time cut; TPC wire", NUM_AW, -0.5, NUM_AW-0.5);
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
   int fNbins;

public:
   AnalyzeNoise(const char* prefix, TDirectory* file, TDirectory* tmp, int nbins = 701):fNbins(nbins) // ctor
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
         fFftSum->Scale(1./(fCount*sqrt(double(fNbins))));
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
         sprintf(title, "adc%02d section %d", module, section);
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
         sprintf(name, "adc%02dch%d", c->adc_module, c->adc_chan);
         
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

struct PlotAwWaveforms
{
   TCanvas* fCanvas = NULL;
   std::vector<int> fAwWires;
   std::vector<TH1D*> fH;

   int fModule = 0;
   int fSection = 0;

   PlotAwWaveforms(TCanvas* c, const std::vector<int>& aw_wires) // ctor
   {
      if (!c) {
         char title[256];
         sprintf(title, "AW waveforms");
         c = new TCanvas(title, title, 900, 650);
         if (!(c->GetShowEventStatus()))
            c->ToggleEventStatus();
         if (!(c->GetShowToolBar()))
            c->ToggleToolBar();
      }

      fCanvas = c;

      fCanvas->cd();

      fAwWires = aw_wires;

      for (unsigned i=0; i<fAwWires.size(); i++) {
         fH.push_back(NULL);
      }
   }

   ~PlotAwWaveforms()
   {
      if (fCanvas) {
         fCanvas->Clear();
         delete fCanvas;
         fCanvas = NULL;
      }

      for (unsigned i=0; i<fAwWires.size(); i++) {
         if (!fH[i])
            continue;
         delete fH[i];
         fH[i] = NULL;
      }
   }

   void Draw(TARunInfo* runinfo, const std::vector<Alpha16Channel*>& hits)
   {
      fCanvas->cd();
      fCanvas->Clear();

      for (unsigned ihit=0; ihit<hits.size(); ihit++) {
         for (unsigned i=0; i<fAwWires.size(); i++) {
            if (hits[ihit]->tpc_wire == fAwWires[i]) {
               if (!fH[i]) {
                  TDirectory* dir = runinfo->fRoot->fgDir;
                  dir->cd();

                  char name[256];
                  sprintf(name, "aw%d_waveform", hits[ihit]->tpc_wire);
         
                  fH[i] = new TH1D(name, name, hits[ihit]->adc_samples.size(), 0, hits[ihit]->adc_samples.size());
               }

               fH[i]->Clear();

               for (unsigned s=0; s<hits[ihit]->adc_samples.size(); s++) {
                  int v = hits[ihit]->adc_samples[s];
                  //printf("bin %d, v %d\n", i, v);
                  fH[i]->SetBinContent(s+1, v);
               }
            }
         }
      }

      int color = 1;
      bool first = true;
      for (unsigned i=0; i<fAwWires.size(); i++) {
         if (!fH[i])
            continue;

         fH[i]->SetLineColor(color);
         color++;

         if (first) {
            first = false;
            fH[i]->SetMinimum(-(1<<15));
            fH[i]->SetMaximum(1<<15);
            fH[i]->GetYaxis()->SetLabelSize(0.10);
            fH[i]->Draw();
         } else {
            fH[i]->Draw("SAME");
         }
      }

      fCanvas->Modified();
      fCanvas->Draw();
      fCanvas->Update();
   }
};

class AdcFlags
{
public:
   bool fPrint = false;
   bool fPulser = false;
   bool fFft = false;
   bool fFilterWaveform = false;
   bool fInvertWaveform = false;
   std::vector<int> fPlotAdc16;
   std::vector<int> fPlotAdc32;
   int fAdcPlotScaledown = 1;
   std::vector<int> fAwWires;
   bool fWfSuppress = false;
   bool fWfVerbose  = false;
   int  fWfThresholdAdc16 = 0;
   int  fWfThresholdAdc32 = 0;
   int  fWfKeepMoreAdc16 = 0;
   int  fWfKeepMoreAdc32 = 0;
   bool fWfSaveBad = false;
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

static double find_pulse_width(int iwire, const int* adc, int nbins, double baseline, double gain, double threshold)
{
   double istart = -1;
   int iend = -1;
   for (int i=1; i<nbins; i++) {
      double v1 = (adc[i]-baseline)*gain;
      if (iwire==262) {}
      //   printf("wire %d, bin %d, adc %d, v1 %f, thr %f, pulse %d %d\n", iwire, i, adc[i], v1, threshold, istart, iend);
      if (v1 > threshold) {
         if (istart < 0) {
            double v0 = (adc[i-1]-baseline)*gain;
            istart = i-1+(v0-threshold)/(v0-v1);
         }
      } else if (v1 < threshold) {
         if (istart > 0) {
            if (iend < 0) {
               double v0 = (adc[i-1]-baseline)*gain;
               iend = i-1+(v0-threshold)/(v0-v1);
               break;
            }
         }
      }
   }

   //if (iwire==262) {
   //   printf("find_pulse_width: bins %d %d\n", istart, iend);
   //}

   return iend-istart;
}

class AdcModule: public TARunObject
{
public:
   AdcFlags* fFlags = NULL;
   int fCounter = 0;
   bool fTrace = false;
   Ncfm* fCfm = NULL;
   NcfmParser* fCfmCuts = NULL;

   int fPulserStart = 0;
   int fPulserEnd   = 0;

   TDirectory* fDirSummary = NULL;
   TDirectory* fDirPulser = NULL;
   TDirectory* fDirWfSuppress = NULL;

   PlotHistograms* fH = NULL;
   PlotNoise* fPN16 = NULL;
   AnalyzeNoise* fAN16 = NULL;
   PlotNoise* fPN32 = NULL;
   AnalyzeNoise* fAN32 = NULL;
   std::vector<PlotA16*> fPlotA16;

   std::vector<AnalyzeNoise> fAN16AWB;
   std::vector<AnalyzeNoise> fAN32AWB;

   PlotAwWaveforms* fPlotAwWaveforms = NULL;

   // pulser analysis plots

   TH1D* fHpulserPh = NULL;
   TH1D* fHpulserLe = NULL;
   TH1D* fHpulserWi = NULL;
   TH1D* fHpulserAwMap = NULL;
   TH1D* fHpulserAwPh  = NULL;
   TH1D* fHpulserAwLe  = NULL;
   TH1D* fHpulserAwWi  = NULL;
   TProfile* fHpulserAwPhMap = NULL;
   TProfile* fHpulserAwLeMap = NULL;
   TProfile* fHpulserAwWiMap = NULL;
   TH2D* fHpulserAwWiPh = NULL;

   // data suppression plots

   int fWfThresholdAdc16 = 0;
   int fWfThresholdAdc32 = 0;
   int fWfKeepMoreAdc16  = 0;
   int fWfKeepMoreAdc32  = 0;

   TH1D* fWfSuppressAdc16Amp = NULL;
   TH1D* fWfSuppressAdc16AmpPos = NULL;
   TH1D* fWfSuppressAdc16AmpNeg = NULL;
   TH1D* fWfSuppressAdc16AmpCumulKeepAll = NULL;
   TH1D* fWfSuppressAdc16AmpCumulDropAll = NULL;
   TH2D* fWfSuppressAdc16AmpCumulKeepMap = NULL;
   TH2D* fWfSuppressAdc16AmpCumulKeepXMap = NULL;
   TH2D* fWfSuppressAdc16AmpCumulDropMap = NULL;
   std::vector<TH1D*> fWfSuppressAdc16AmpCumulKeep;
   std::vector<TH1D*> fWfSuppressAdc16AmpCumulDrop;

   std::vector<TH1D*> fWfSuppressAdc16AmpCumulKeep16;

   TH1D* fWfSuppressAdc32Amp = NULL;
   TH1D* fWfSuppressAdc32AmpPos = NULL;
   TH1D* fWfSuppressAdc32AmpNeg = NULL;
   TH1D* fWfSuppressAdc32AmpCumulKeepAll = NULL;
   TH1D* fWfSuppressAdc32AmpCumulDropAll = NULL;
   TH2D* fWfSuppressAdc32AmpCumulKeepMap = NULL;
   TH2D* fWfSuppressAdc32AmpCumulKeepXMap = NULL;
   TH2D* fWfSuppressAdc32AmpCumulDropMap = NULL;
   std::vector<TH1D*> fWfSuppressAdc32AmpCumulKeep;
   std::vector<TH1D*> fWfSuppressAdc32AmpCumulDrop;

   std::vector<TH1D*> fWfSuppressAdc32AmpCumulKeep32;

public:
   AdcModule(TARunInfo* runinfo, AdcFlags* f)
      : TARunObject(runinfo)
   {
      if (fTrace)
         printf("AdcModule::ctor!\n");

      fModuleName = "adc_module";
      fFlags = f;

      std::string agcfmdb_path="agcfmdb";
      if (getenv("AGRELEASE"))
      {
         agcfmdb_path=getenv("AGRELEASE");
         agcfmdb_path+="/agana/agcfmdb";
      }
      fCfm     = new Ncfm(agcfmdb_path.c_str());

      runinfo->fRoot->fOutputFile->cd();
      TDirectory* aw = gDirectory->mkdir("adc");
      aw->cd(); // select correct ROOT directory

      fDirSummary = aw->mkdir("summary");
      fDirPulser  = aw->mkdir("pulser");
      fDirWfSuppress = aw->mkdir("wfsuppress");

      fH = new PlotHistograms(fDirSummary);

      if (fFlags->fFft) {
         TDirectory* fft_file = aw->mkdir("noise_fft");

         TDirectory* fft_tmp = runinfo->fRoot->fgDir->mkdir("bv_noise_fft");

         fAN16 = new AnalyzeNoise("adc16", fft_file, fft_tmp, 701);
         fPN16 = new PlotNoise("adc16");

         for(int awb=0; awb<32; ++awb)
            {
               TString ANname = TString::Format("adc16AWB%02d",awb);
               fAN16AWB.emplace_back(ANname.Data(), fft_file, fft_tmp, 701);
            }

         fft_tmp = runinfo->fRoot->fgDir->mkdir("aw_noise_fft");

         fAN32 = new AnalyzeNoise("adc32", fft_file, fft_tmp, 511);
         fPN32 = new PlotNoise("adc32");

         for(int awb=0; awb<32; ++awb)
            {
               TString ANname = TString::Format("adc32AWB%02d",awb);
               fAN32AWB.emplace_back(ANname.Data(), fft_file, fft_tmp, 511);
            }
      }

      aw->cd(); // select correct ROOT directory

      if (fFlags->fAwWires.size() > 0) {
         fPlotAwWaveforms = new PlotAwWaveforms(NULL, fFlags->fAwWires);
      }
   }

   ~AdcModule()
   {
      if (fTrace)
         printf("AdcModule::dtor!\n");
      if (fH) {
         delete fH;
         fH = NULL;
      }
      if (fCfmCuts) {
         delete fCfmCuts;
         fCfmCuts = NULL;
      }
      if (fCfm) {
         delete fCfm;
         fCfm = NULL;
      }
   }

   void BeginRun(TARunInfo* runinfo)
   {
      if (fTrace)
         printf("AdcModule::BeginRun, run %d, file %s\n", runinfo->fRunNo, runinfo->fFileName.c_str());
      //time_t run_start_time = runinfo->fOdb->odbReadUint32("/Runinfo/Start time binary", 0, 0);
      //printf("ODB Run start time: %d: %s", (int)run_start_time, ctime(&run_start_time));

      runinfo->fOdb->RB("Equipment/Ctrl/Settings/FwPulserEnable", &fFlags->fPulser);
      printf("pulser mode: enabled: %d\n", fFlags->fPulser);

      fCfmCuts = fCfm->ParseFile("adc", "cuts", runinfo->fRunNo);

      fPulserStart = fCfmCuts->GetInt("adc_bin_pulser_start", 330); // 160;
      fPulserEnd   = fCfmCuts->GetInt("adc_bin_pulser_end", 350); // 200;

      fCounter = 0;
      runinfo->fRoot->fOutputFile->cd(); // select correct ROOT directory

      for (unsigned i=0; i<fFlags->fPlotAdc16.size(); i++) {
         fPlotA16.push_back(new PlotA16(NULL, fFlags->fPlotAdc16[i], 0));
      }

      for (unsigned i=0; i<fFlags->fPlotAdc32.size(); i++) {
         fPlotA16.push_back(new PlotA16(NULL, fFlags->fPlotAdc32[i], 1));
         fPlotA16.push_back(new PlotA16(NULL, fFlags->fPlotAdc32[i], 2));
      }

      if (fFlags->fPulser) {
         fDirPulser->cd();

         double wi_start = 0;
         double wi_end =  30;
         double wi_bins = 100.0;

         fHpulserPh = new TH1D("adc_pulser_ph",    "ADC pulser pulse height; ADC counts", 100, 0, MAX_AW_AMP);
         fHpulserLe = new TH1D("adc_pulser_time",  "ADC pulser pulse time; ADC time bins", 100.0, fPulserStart, fPulserEnd);
         fHpulserWi = new TH1D("adc_pulser_width", "ADC pulser pulse width; ADC time bins", wi_bins, wi_start, wi_end);


         fHpulserAwMap = new TH1D("adc_pulser_aw_map", "AW pulser map; AW wire number", NUM_AW, -0.5, NUM_AW-0.5);
         fHpulserAwMap->SetMinimum(0);
         fHpulserAwPh = new TH1D("adc_pulser_aw_ph",    "AW pulser pulse height; ADC counts", 100, 0, MAX_AW_AMP);
         fHpulserAwLe = new TH1D("adc_pulser_aw_time",  "AW pulser pulse time; ADC time bins", 100.0, fPulserStart, fPulserEnd);
         fHpulserAwWi = new TH1D("adc_pulser_aw_width", "AW pulser pulse width; ADC time bins", wi_bins, wi_start, wi_end);

         fHpulserAwPhMap = new TProfile("adc_pulse_aw_ph_map", "AW pulser pulse height vs wire number; AW wire number; ADC counts", NUM_AW, -0.5, NUM_AW-0.5);
         fHpulserAwLeMap = new TProfile("adc_pulse_aw_time_map", "AW pulser pulse time vs wire number; AW wire number; ADC time bins", NUM_AW, -0.5, NUM_AW-0.5, fPulserStart, fPulserEnd);
         fHpulserAwLeMap->SetMinimum(fPulserStart);
         fHpulserAwLeMap->SetMaximum(fPulserEnd);
         fHpulserAwWiMap = new TProfile("adc_pulse_aw_width_map", "AW pulser pulse width vs wire number; AW wire number; ADC time bins", NUM_AW, -0.5, NUM_AW-0.5, wi_start, wi_end);
         fHpulserAwWiPh  = new TH2D("adc_pulse_aw_width_vs_ph", "AW pulser pulse width vs pulse height; pulse height, ADC counts; pulse width, ADC time bins", 100, 0, MAX_AW_AMP, wi_bins, wi_start, wi_end);
      }

      if (fFlags->fWfSuppress) {

         if (fFlags->fWfThresholdAdc16)
            fWfThresholdAdc16 = fFlags->fWfThresholdAdc16;
         else
            runinfo->fOdb->RI("Equipment/Ctrl/Settings/ADC/adc16_sthreshold", &fWfThresholdAdc16);

         if (fFlags->fWfThresholdAdc32)
            fWfThresholdAdc32 = fFlags->fWfThresholdAdc32;
         else
            runinfo->fOdb->RI("Equipment/Ctrl/Settings/ADC/adc32_sthreshold", &fWfThresholdAdc32);

         if (fFlags->fWfKeepMoreAdc16)
            fWfKeepMoreAdc16 = fFlags->fWfKeepMoreAdc16;
         else
            runinfo->fOdb->RI("Equipment/Ctrl/Settings/ADC/adc16_keep", &fWfKeepMoreAdc16);

         if (fFlags->fWfKeepMoreAdc32)
            fWfKeepMoreAdc32 = fFlags->fWfKeepMoreAdc32;
         else
            runinfo->fOdb->RI("Equipment/Ctrl/Settings/ADC/adc32_keep", &fWfKeepMoreAdc32);

         fprintf(stderr, "ADC data suppression: adc16 threshold %d, keep_more %d\n", fWfThresholdAdc16, fWfKeepMoreAdc16);
         fprintf(stderr, "ADC data suppression: adc32 threshold %d, keep_more %d\n", fWfThresholdAdc32, fWfKeepMoreAdc32);

         fDirWfSuppress->cd();

         //int adc_amp = 32000;
         double adc_amp = 4000;
         double adc_bin_size = 50;

         // ADC16 plots

         if (fWfSuppressAdc16Amp == NULL) {
            int min = 0;
            int max = adc_amp;
            fWfSuppressAdc16Amp = new TH1D("adc16_suppress_adc_amp", "WfSuppress ADC amp; adc counts", (max-min)/adc_bin_size, min, max);
         }
         
         if (fWfSuppressAdc16AmpPos == NULL) {
            int min = 0;
            int max = adc_amp;
            fWfSuppressAdc16AmpPos = new TH1D("adc16_suppress_adc_amp_pos", "WfSuppress ADC amp pos; adc counts", (max-min)/adc_bin_size, min, max);
         }
         
         if (fWfSuppressAdc16AmpNeg == NULL) {
            int min = -adc_amp;
            int max = 0;
            fWfSuppressAdc16AmpNeg = new TH1D("adc16_suppress_adc_amp_neg", "WfSuppress ADC amp neg; adc counts", (max-min)/adc_bin_size, min, max);
         }
         
         if (fWfSuppressAdc16AmpCumulKeepAll == NULL) {
            int min = -1;
            int max = adc_amp;
            fWfSuppressAdc16AmpCumulKeepAll = new TH1D("adc16_suppress_cumul_keep", "WfSuppress cumulative kept channels; ch_threshold, adc counts", (max-min+1)/adc_bin_size, min-0.5, max+0.5);
         }
         
         if (fWfSuppressAdc16AmpCumulDropAll == NULL) {
            int min = -1;
            int max = adc_amp;
            fWfSuppressAdc16AmpCumulDropAll = new TH1D("adc16_suppress_cumul_drop", "WfSuppress cumulative dropped channels; ch_threshold, adc_counts", (max-min+1)/adc_bin_size, min-0.5, max+0.5);
         }

         if (fWfSuppressAdc16AmpCumulKeepMap == NULL) {
            int min = -1;
            int max = adc_amp;
            fWfSuppressAdc16AmpCumulKeepMap = new TH2D("adc16_suppress_cumul_keep_map", "WfSuppress cumulative kept channels; adcNN; ch_threshold, adc counts", ADC_MODULE_LAST+1, 0-0.5, ADC_MODULE_LAST+1-0.5, (max-min+1)/adc_bin_size, min-0.5, max+0.5);
         }
         
         if (fWfSuppressAdc16AmpCumulKeepXMap == NULL) {
            int min = -1;
            int max = adc_amp;
            fWfSuppressAdc16AmpCumulKeepXMap = new TH2D("adc16_suppress_cumul_keep_xmap", "WfSuppress cumulative kept channels; xmodule; ch_threshold, adc counts", 2*ADC_MODULE_LAST+1, 0-0.5, 2*ADC_MODULE_LAST+1-0.5, (max-min+1)/adc_bin_size, min-0.5, max+0.5);
         }
         
         if (fWfSuppressAdc16AmpCumulDropMap == NULL) {
            int min = -1;
            int max = adc_amp;
            fWfSuppressAdc16AmpCumulDropMap = new TH2D("adc16_suppress_cumul_drop_map", "WfSuppress cumulative dropped channels; adcNN; ch_threshold, adc_counts", ADC_MODULE_LAST+1, 0-0.5, ADC_MODULE_LAST+1-0.5, (max-min+1)/adc_bin_size, min-0.5, max+0.5);
         }

         for (int i=0; i<=ADC_MODULE_LAST; i++) {
            char name[100];
            char title[100];

            sprintf(name,  "adc16_adc%02d_keep", i);
            sprintf(title, "adc%02d cumulative kept channels; ch_threshold, adc counts", i);

            int min = -1;
            int max = adc_amp;

            TH1D* h = new TH1D(name, title, max-min+1, min-0.5, max+0.5);

            fWfSuppressAdc16AmpCumulKeep.push_back(h);
         }
         
         for (int i=0; i<16; i++) {
            char name[100];
            char title[100];

            sprintf(name,  "adc16_adcxx_keep%02d", i);
            sprintf(title, "adcxx_keep%02d cumulative kept channels; ch_threshold, adc counts", i);

            int min = -1;
            int max = adc_amp;

            TH1D* h = new TH1D(name, title, max-min+1, min-0.5, max+0.5);

            fWfSuppressAdc16AmpCumulKeep16.push_back(h);
         }
         
         for (int i=0; i<=ADC_MODULE_LAST; i++) {
            char name[100];
            char title[100];

            sprintf(name,  "adc16_adc%02d_drop", i);
            sprintf(title, "adc%02d cumulative dropped channels; ch_threshold, adc counts", i);

            int min = -1;
            int max = adc_amp;

            TH1D* h = new TH1D(name, title, max-min+1, min-0.5, max+0.5);

            fWfSuppressAdc16AmpCumulDrop.push_back(h);
         }

         // ADC32 plots

         if (fWfSuppressAdc32Amp == NULL) {
            int min = 0;
            int max = adc_amp;
            fWfSuppressAdc32Amp = new TH1D("adc32_suppress_adc_amp", "WfSuppress ADC amp; adc counts", (max-min)/adc_bin_size, min, max);
         }
         
         if (fWfSuppressAdc32AmpPos == NULL) {
            int min = 0;
            int max = adc_amp;
            fWfSuppressAdc32AmpPos = new TH1D("adc32_suppress_adc_amp_pos", "WfSuppress ADC amp pos; adc counts", (max-min)/adc_bin_size, min, max);
         }
         
         if (fWfSuppressAdc32AmpNeg == NULL) {
            int min = -adc_amp;
            int max = 0;
            fWfSuppressAdc32AmpNeg = new TH1D("adc32_suppress_adc_amp_neg", "WfSuppress ADC amp neg; adc counts", (max-min)/adc_bin_size, min, max);
         }
         
         if (fWfSuppressAdc32AmpCumulKeepAll == NULL) {
            int min = -1;
            int max = adc_amp;
            fWfSuppressAdc32AmpCumulKeepAll = new TH1D("adc32_suppress_cumul_keep", "WfSuppress cumulative kept channels; ch_threshold, adc counts", (max-min+1)/adc_bin_size, min-0.5, max+0.5);
         }
         
         if (fWfSuppressAdc32AmpCumulDropAll == NULL) {
            int min = -1;
            int max = adc_amp;
            fWfSuppressAdc32AmpCumulDropAll = new TH1D("adc32_suppress_cumul_drop", "WfSuppress cumulative dropped channels; ch_threshold, adc_counts", (max-min+1)/adc_bin_size, min-0.5, max+0.5);
         }

         if (fWfSuppressAdc32AmpCumulKeepMap == NULL) {
            int min = -1;
            int max = adc_amp;
            fWfSuppressAdc32AmpCumulKeepMap = new TH2D("adc32_suppress_cumul_keep_map", "WfSuppress cumulative kept channels; adcNN; ch_threshold, adc counts", ADC_MODULE_LAST+1, 0-0.5, ADC_MODULE_LAST+1-0.5, (max-min+1)/adc_bin_size, min-0.5, max+0.5);
         }
         
         if (fWfSuppressAdc32AmpCumulKeepXMap == NULL) {
            int min = -1;
            int max = adc_amp;
            fWfSuppressAdc32AmpCumulKeepXMap = new TH2D("adc32_suppress_cumul_keep_xmap", "WfSuppress cumulative kept channels; xmodule; ch_threshold, adc counts", 2*ADC_MODULE_LAST+1, 0-0.5, 2*ADC_MODULE_LAST+1-0.5, (max-min+1)/adc_bin_size, min-0.5, max+0.5);
         }
         
         if (fWfSuppressAdc32AmpCumulDropMap == NULL) {
            int min = -1;
            int max = adc_amp;
            fWfSuppressAdc32AmpCumulDropMap = new TH2D("adc32_suppress_cumul_drop_map", "WfSuppress cumulative dropped channels; adcNN; ch_threshold, adc_counts", ADC_MODULE_LAST+1, 0-0.5, ADC_MODULE_LAST+1-0.5, (max-min+1)/adc_bin_size, min-0.5, max+0.5);
         }

         for (int i=0; i<=ADC_MODULE_LAST; i++) {
            char name[100];
            char title[100];

            sprintf(name,  "adc32_adc%02da_keep", i);
            sprintf(title, "adc%02da cumulative kept channels; ch_threshold, adc counts", i);

            int min = -1;
            int max = adc_amp;

            TH1D* h = new TH1D(name, title, max-min+1, min-0.5, max+0.5);

            fWfSuppressAdc32AmpCumulKeep.push_back(h);

            sprintf(name,  "adc32_adc%02db_keep", i);
            sprintf(title, "adc%02db cumulative kept channels; ch_threshold, adc counts", i);

            h = new TH1D(name, title, max-min+1, min-0.5, max+0.5);

            fWfSuppressAdc32AmpCumulKeep.push_back(h);
         }
         
         for (int i=0; i<32; i++) {
            char name[100];
            char title[100];

            sprintf(name,  "adc32_adcxx_keep%02d", i);
            sprintf(title, "adcxx_keep%02d cumulative kept channels; ch_threshold, adc counts", i);

            int min = -1;
            int max = adc_amp;

            TH1D* h = new TH1D(name, title, max-min+1, min-0.5, max+0.5);

            fWfSuppressAdc32AmpCumulKeep32.push_back(h);
         }
         
         for (int i=0; i<=ADC_MODULE_LAST; i++) {
            char name[100];
            char title[100];

            sprintf(name,  "adc32_adc%02da_drop", i);
            sprintf(title, "adc%02da cumulative dropped channels; ch_threshold, adc counts", i);

            int min = -1;
            int max = adc_amp;

            TH1D* h = new TH1D(name, title, max-min+1, min-0.5, max+0.5);

            fWfSuppressAdc32AmpCumulDrop.push_back(h);

            sprintf(name,  "adc32_adc%02db_drop", i);
            sprintf(title, "adc%02db cumulative dropped channels; ch_threshold, adc counts", i);

            h = new TH1D(name, title, max-min+1, min-0.5, max+0.5);

            fWfSuppressAdc32AmpCumulDrop.push_back(h);
         }

         fDirSummary->cd();
      }
   }

   void EndRun(TARunInfo* runinfo)
   {
      if (fTrace)
         printf("AdcModule::EndRun, run %d, events %d\n", runinfo->fRunNo, fCounter);
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
         printf("AdcModule::PauseRun, run %d\n", runinfo->fRunNo);
   }

   void ResumeRun(TARunInfo* runinfo)
   {
      if (fTrace)
         printf("AdcModule::ResumeRun, run %d\n", runinfo->fRunNo);
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

   bool AnalyzeHit(const TARunInfo* runinfo, const Alpha16Packet* p, const Alpha16Channel* hit, std::vector<AgAwHit>* flow_hits)
   {
      //char xname[256];
      //char xtitle[256];
      //sprintf(xname, "m%02d_c%02d_w%03d", hit->adc_module, hit->adc_chan, hit->tpc_wire);
      //sprintf(xtitle, "AW Waveform ADC module %d, channel %d, tpc wire %d", hit->adc_module, hit->adc_chan, hit->tpc_wire);

      int i = hit->tpc_wire;
      int r = 1;

      if (runinfo->fRunNo < 1694) {
         if (hit->adc_module == 5 && hit->adc_chan >= 16) {
            r = 2; // fmc-adc32-rev0 gain 1
         }
         
         if (hit->adc_module == 6 && hit->adc_chan >= 16) {
            r = 3; // fmc-adc32-rev0 gain 2
         }
      } else {
         r = 4; // fmc-adc32-rev1 gain 3
      }

#if 0
      // blank off bad FMC-ADC32 channels
      if (runinfo->fRunNo >= 2120 && runinfo->fRunNo < 99999) {
         if (hit->adc_chan == 16) {
            return;
         }
         if (hit->adc_chan == 25) {
            return;
         }
         if (hit->adc_chan == 34) {
            return;
         }
         if (hit->adc_chan == 41) {
            return;
         }
      }
#endif

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

      //if (i < 0)
      //   return false;

      int iwire = i;

      bool is_aw = (iwire >= 0);
      bool is_bsc = true;

      bool is_adc16 = (hit->adc_chan < 16);
      bool is_adc32 = (hit->adc_chan >= 16);

      int max_fft_count = 30; // FFT is slow, limit number of events analyzed

      if (is_adc16 && fft_first_adc16) {
         fft_first_adc16 = false;
         if (fAN16 && fAN16->fCount < max_fft_count) {
            fAN16->AddWaveform(hit->adc_samples);
            fAN16AWB.at(hit->preamp_pos).AddWaveform(hit->adc_samples);
            if (fPN16)
               fPN16->Plot(fAN16);
         }                  
      }
      
      int max_fft_awb_count = 30000;
      if (is_adc16 && fAN16AWB.size() > 0 )
         {
            if ( fAN16AWB.at(hit->preamp_pos).fCount < max_fft_awb_count) 
               {
                  fAN16AWB.at(hit->preamp_pos).AddWaveform(hit->adc_samples);
               }                  
         }

      if (is_adc32 && fft_first_adc32) {
         fft_first_adc32 = false;
         if (fAN32 && fAN32->fCount < max_fft_count) {
            fAN32->AddWaveform(hit->adc_samples);
            if (fPN32)
               fPN32->Plot(fAN32);
         }
      }


      if (!is_adc16 && fAN32AWB.size() > 0 )
         {
            if ( fAN32AWB.at(hit->preamp_pos).fCount < max_fft_awb_count) 
               {
                  fAN32AWB.at(hit->preamp_pos).AddWaveform(hit->adc_samples);
               }                  
         }

      // analyze data suppression

      if (fFlags->fWfSuppress && (is_adc16 || is_adc32)) {
         int imodule = hit->adc_module;
         int ichan = hit->adc_chan;
         int xmodule = -1;
         int keep16 = -1;
         int keep32 = -1;
         
         //printf("imodule %d, size %d\n", imodule, (int)fWfSuppress.size());

         WfSuppressAdc ch_supp;
         WfSuppressAdc *s = &ch_supp;

         // chan 0..15 is adc16
         // chan 16..31 is adc32a
         // chan 32..47 is adc32b

         if (is_adc16) {
            s->Config(fWfThresholdAdc16, fWfKeepMoreAdc16, p->nsamples);
            xmodule = imodule;

            if (imodule == 12) {
               keep16 = ichan;
               assert(keep16>=0 && keep16<16);
            }
         }
         
         if (is_adc32) {
            s->Config(fWfThresholdAdc32, fWfKeepMoreAdc32, p->nsamples);
            xmodule = imodule*2;
            // chan 32..47 is 2nd section of adc32
            if (ichan >= 32)
               xmodule += 1;

            //xmodule = imodule*4;
            // chan 32..47 is 2nd section of adc32
            //if (ichan >= 32)
            //   xmodule += 2;
            //if (ichan % 2)
            //   xmodule += 1;

            //if (ichan >= 8*3+16) // 40+8
            //   xmodule += 3;
            //else if (ichan >= 8*2+16) // 32+8
            //   xmodule += 2;
            //else if (ichan >= 8*1+16) // 24+8
            //   xmodule += 1;
            //else // 16+8
            //   xmodule += 0;

            if (imodule == 12) {
               keep32 = ichan - 16;
               assert(keep32>=0 && keep32<32);
            }
         }

         s->Reset();
         
         int nsamples = hit->adc_samples.size();
         
         // do not pass fake samples added by unpacking to the data supression code
         if (p->nsamples_supp)
            nsamples = p->nsamples_supp;

         //if (nsamples > p->nsamples - 6)
         //   nsamples = p->nsamples - 6;

         for (int i=0; i<nsamples; i++) {
            s->Add(hit->adc_samples[i]);
            //if (0) {
            //   printf("adc %02d, chan %2d: bin %3d, adc %d, base %d, amp %4d, keep %d %d, state: %s\n", imodule, ichan, i, hit->adc_samples[i], base, amp, k, keep, s->PrintToString().c_str());
            //}
         }

         double xamp = -s->fAdcMaxNeg;
         if (s->fAdcMaxPos > xamp)
            xamp = s->fAdcMaxPos;

         if (is_adc16) {
            fWfSuppressAdc16Amp->Fill(xamp);
            fWfSuppressAdc16AmpPos->Fill(s->fAdcMaxPos);
            fWfSuppressAdc16AmpNeg->Fill(s->fAdcMaxNeg);
            
            for (int i=0; i<xamp; i++) {
               fWfSuppressAdc16AmpCumulKeepAll->Fill(i);
               fWfSuppressAdc16AmpCumulKeep[xmodule]->Fill(i);
               if (keep16>=0)
                  fWfSuppressAdc16AmpCumulKeep16[keep16]->Fill(i);
               fWfSuppressAdc16AmpCumulKeepMap->Fill(imodule, i);
               fWfSuppressAdc16AmpCumulKeepXMap->Fill(xmodule, i);
            }
         
            for (int i=xamp; i<32000; i++) {
               fWfSuppressAdc16AmpCumulDropAll->Fill(i);
               fWfSuppressAdc16AmpCumulDrop[xmodule]->Fill(i);
               fWfSuppressAdc16AmpCumulDropMap->Fill(imodule, i);
            }
         }
         
         if (is_adc32) {
            fWfSuppressAdc32Amp->Fill(xamp);
            fWfSuppressAdc32AmpPos->Fill(s->fAdcMaxPos);
            fWfSuppressAdc32AmpNeg->Fill(s->fAdcMaxNeg);
            
            for (int i=0; i<xamp; i++) {
               fWfSuppressAdc32AmpCumulKeepAll->Fill(i);
               fWfSuppressAdc32AmpCumulKeep[xmodule]->Fill(i);
               if (keep32>=0)
                  fWfSuppressAdc32AmpCumulKeep32[keep32]->Fill(i);
               fWfSuppressAdc32AmpCumulKeepMap->Fill(imodule, i);
               fWfSuppressAdc32AmpCumulKeepXMap->Fill(xmodule, i);
            }
         
            for (int i=xamp; i<32000; i++) {
               fWfSuppressAdc32AmpCumulDropAll->Fill(i);
               fWfSuppressAdc32AmpCumulDrop[xmodule]->Fill(i);
               fWfSuppressAdc32AmpCumulDropMap->Fill(imodule, i);
            }
         }
         
         //int range = ampmax - ampmin;

         bool bad_baseline  = false;
         bool bad_keep_bit  = false;
         bool bad_keep_last = false;
         bool bad_kept = false;
         bool bad = false;

         if (s->fBaselineReady && p->baseline && (s->fBaseline != p->baseline)) {
            if ((s->fBaseline&0xFFFF) != (p->baseline&0xFFFF) &&
                ((s->fBaseline+1)&0xFFFF) != (p->baseline&0xFFFF) &&
                (s->fBaseline&0xFFFF) != ((p->baseline+1)&0xFFFF)
                ) {
               printf("maybe bad baseline 0x%04x 0x%04x\n", s->fBaseline&0xFFFF, p->baseline&0xFFFF);
               bad_baseline = true;
               bad = true;
            }
         }
         
         if (s->fBaselineReady && p->baseline && (s->fKeepBit != p->keep_bit)) {
            bad_keep_bit = true;
            bad = true;
         }
         
         if (s->fBaselineReady && p->baseline && (s->fKeepLast != p->keep_last)) {
            bad_keep_last = true;
            bad = true;
         }

         int kept = 2*s->fKeepLast + 2*s->fKeepMore + 4;
         if (kept > p->nsamples - 2)
            kept = p->nsamples - 2;

         if (p->packetVersion == 3) {
            if (p->supp_enabled && !p->keep_bit) {
               // keep_bit is not set, this is a short packet without any adc samples.
               kept = 0;
            }
         }
         
         if (s->fBaselineReady && p->baseline && p->supp_enabled && (kept != p->nsamples_supp)) {
            bad_kept = true;
            bad = true;
         }
         
         if (fFlags->fWfVerbose || bad) {
            //printf("adc %02d, chan %2d: wfsuppress: %s, keep: %d, xamp %4d, kept %3d, keep_bit %d, keep_last %d, baseline %d", imodule, ichan, s->PrintToString().c_str(), keep, (int)xamp, p->nsamples_supp, p->keep_bit, p->keep_last, p->baseline);
            printf("adc %02d, chan %2d: samples %3d/%3d, baseline %d %4d vs %4d, range %5d..%4d, keep_bit %d vs %d, keep_last %3d vs %3d, adc_last %5d, supp %d, kept %3d vs %3d", imodule, ichan, (int)hit->adc_samples.size(), p->nsamples_supp, s->fBaselineReady, s->fBaseline, p->baseline, s->fAdcMaxNeg, s->fAdcMaxPos, s->fKeepBit, p->keep_bit, s->fKeepLast, p->keep_last, s->fAdcLast, p->supp_enabled, p->nsamples_supp, kept);
            if (bad_baseline)
               printf(" AAA bad baseline %4d should be %4d", s->fBaseline, p->baseline);
            if (bad_keep_bit)
               printf(" AAA bad keep_bit %d should be %d", s->fKeepBit, p->keep_bit);
            if (bad_keep_last)
               printf(" AAA bad keep_last %3d should be %3d", s->fKeepLast, p->keep_last);
            if (bad_kept)
               printf(" AAA bad kept %3d should be %3d", kept, p->nsamples_supp);
            printf("\n");
         }

#if 0
         if (bad) {
            printf("AAA samples: ");
            for (unsigned i=0; i<hit->adc_samples.size(); i++) {
               printf(" %4d", hit->adc_samples[i] - s->fBaseline);
            }
            printf("\n");
            exit(123);
         }
#endif

#if 0
         if (!keep) {
            printf("TTTT: ");
            printf("pwb %02d, sca %d, chan %2d: wfsuppress: %s, , keep: %d, xamp %d, mismatch!\n", imodule, isca, ichan, s->PrintToString().c_str(), keep, (int)xamp);
            
            static int count = 0;
            
            if (count < 100) {
               char name[256];
               char title[256];
               
               fDirWfSuppress->cd();
               
               sprintf(name, "pwb_%02d_sca_%d_chan_%02d_range_%d_count_%d", imodule, isca, ichan, range, count);
               sprintf(title, "pwb %02d, sca %d, chan %2d, range %d, count %d", imodule, isca, ichan, range, count);
               WfToTH1D(name, title, c->adc_samples);
               
               //sprintf(name, "pwb_%02d_sca_%d_chan_%02d_range_%d_count_%d_base", imodule, isca, ichan, range, count);
               //sprintf(title, "pwb %02d, sca %d, chan %2d, range %d, count %d, baseline", imodule, isca, ichan, range, count);
               //WfToTH1D(name, title, samples_base);
               
               sprintf(name, "pwb_%02d_sca_%d_chan_%02d_range_%d_count_%d_amp", imodule, isca, ichan, range, count);
               sprintf(title, "pwb %02d, sca %d, chan %2d, range %d, count %d, amplitude", imodule, isca, ichan, range, count);
               WfToTH1D(name, title, samples_amp);
               
               count++;
            }
         }
#endif
      }

      // analyze baseline
      
      int is_start = fCfmCuts->GetInt("adc_bin_baseline_start", 0);
      int is_baseline = fCfmCuts->GetInt("adc_bin_baseline_end", 100);

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
      double wi = 0;
      double hit_time = 0;
      double hit_amp = 0;

      double cut_brms = fCfmCuts->GetDouble("cut_baseline_rms", 500);

      bool special = false;

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

            if (is_aw) {
               fH->fHbaselineBadRmsAwMap->Fill(iwire);
            }
         }
      } else {
         have_baseline = true;

         if (fH) {
            fH->fHbaselineMean->Fill(bmean);

            if (brms < MAX_AW_BRMS)
               fH->fHbaselineRms->Fill(brms);
            else
               fH->fHbaselineRms->Fill(MAX_AW_BRMS-1);

            if (brange < MAX_AW_BRANGE) {
               fH->fHbaselineRange->Fill(brange);
            } else {
               fH->fHbaselineRange->Fill(MAX_AW_BRANGE-1);
            }

            if (is_aw) {
               fH->fHbaselineMeanAwMap->Fill(iwire, bmean);
               fH->fHbaselineRmsAwMap->Fill(iwire, brms);
               fH->fHbaselineRangeAwMap->Fill(iwire, brange);
            }

            if (wrange < MAX_AW_AMP) {
               fH->fHrange->Fill(wrange);
               if (is_aw) {
                  fH->fHrangeAwMap->Fill(iwire, wrange);
               }
            } else {
               fH->fHrange->Fill(MAX_AW_AMP-1);
               if (is_aw) {
                  fH->fHrangeAwMap->Fill(iwire, MAX_AW_AMP-1);
               }
            }
         }
      
         ph = bmean - wmin;

         if (is_aw) {
            if (ph > 500)
               fH->fAwMapPh500->Fill(iwire);
            if (ph > 1000)
               fH->fAwMapPh1000->Fill(iwire);
            if (ph > 5000)
               fH->fAwMapPh5000->Fill(iwire);
            if (ph > 10000)
               fH->fAwMapPh10000->Fill(iwire);
         }

         double cfd_thr_pct = fCfmCuts->GetDouble("cfd_thr_pct", 0.75);
         double cfd_thr = cfd_thr_pct*ph;

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
            else if (is_adc32) {
               fH->fHped_adc32->Fill(ph);
               if (hit->preamp_pos == 0 || hit->preamp_pos == 16)
                  fH->fHped_adc32_p0->Fill(ph);
               else
                  fH->fHped_adc32_pN->Fill(ph);
            }
         }

         double ph_hit_thr_adc16 = 0;
         double ph_hit_thr_adc32 = 0;

         if (runinfo->fRunNo < 1244) {
            ph_hit_thr_adc16 = 1200;
            ph_hit_thr_adc32 =  600;
         } else if (runinfo->fRunNo < 1405) {
            ph_hit_thr_adc16 =  500;
            ph_hit_thr_adc32 =  500;
         } else if (runinfo->fRunNo < 1450) {
            ph_hit_thr_adc16 =  1000;
            ph_hit_thr_adc32 =  1000;
         } else if (runinfo->fRunNo < 1694) {
            ph_hit_thr_adc16 =  2500;
            ph_hit_thr_adc32 =  2500;
         } else if (runinfo->fRunNo < 2028) {
            ph_hit_thr_adc16 =  1000;
            ph_hit_thr_adc32 =   800;
         } else if (runinfo->fRunNo < 2164) {
            ph_hit_thr_adc16 =  2000;
            ph_hit_thr_adc32 =  1500;
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
            ph_hit_thr_adc16 =  fCfmCuts->GetDouble("ph_hit_thr_adc16", 1000);
            ph_hit_thr_adc32 =  fCfmCuts->GetDouble("ph_hit_thr_adc32", 750);
         }

         double ph_hit_thr = 0;

         if (is_adc16)
            ph_hit_thr = ph_hit_thr_adc16;
         else if (is_adc32)
            ph_hit_thr = ph_hit_thr_adc32;

         double time_bin = 0;
         double time_offset = 0;

         if (is_adc16)
            time_bin = 1000.0/100.0; // 100 MHz ADC
         else if (is_adc32)
            time_bin = 1000.0/62.5; // 62.5 MHz ADC

         double pulse_time_middle_adc16 = 0; // ADC time bins
         double pulse_time_middle_adc32 = 0; // ADC time bins

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
         } else {
            pulse_time_middle_adc16 = fCfmCuts->GetInt("adc_bin_pc_middle_adc16", 147); // ADC time bins
            pulse_time_middle_adc32 = fCfmCuts->GetInt("adc_bin_pc_middle_adc32", 110); // ADC time bins
         }

         double time_pc = 1000.0; // place PC drift times at 1000ns.

         if (is_adc16)
            time_offset = time_pc - pulse_time_middle_adc16*time_bin;
         else if (is_adc32)
            time_offset = time_pc - pulse_time_middle_adc32*time_bin;

         double adc_gain = 1.0;
         double adc_offset = 0.0;
         
         if (is_adc16)
            adc_gain = 1.0;
         else if (r == 2)
            adc_gain = 4.0;
         else if (r == 3)
            adc_gain = 2.0;
         else if (r == 4)
            adc_gain = 4.0/3.0; // fmc-adc32-rev1 with gain 3

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

            wi = find_pulse_width(iwire, hit->adc_samples.data(), hit->adc_samples.size(), bmean, -1.0, 2000.0);

            if (fH) {
               fH->fHle->Fill(le);
               fH->fHwi->Fill(wi);

               if (is_aw) {
                  if (iwire == 282) {
                     fH->fHle_aw282->Fill(le);
                  }
               }
            
               if (is_adc16) {
                  fH->fHle_adc16->Fill(le);
                  fH->fHwi_adc16->Fill(wi);
               } else if (is_adc32) {
                  fH->fHle_adc32->Fill(le);
                  fH->fHwi_adc32->Fill(wi);
               }
            }
            
            if (fFlags->fPulser) {
               if (le > fPulserStart && le < fPulserEnd) {
                  fHpulserPh->Fill(ph);
                  fHpulserLe->Fill(le);
                  fHpulserWi->Fill(wi);

                  if (is_aw) {
                     fHpulserAwMap->Fill(iwire);
                     fHpulserAwPh->Fill(ph);
                     fHpulserAwLe->Fill(le);
                     fHpulserAwWi->Fill(wi);
                     fHpulserAwPhMap->Fill(iwire, ph);
                     fHpulserAwLeMap->Fill(iwire, le);
                     fHpulserAwWiMap->Fill(iwire, wi);
                     fHpulserAwWiPh->Fill(ph, wi);
                  }

                  if (is_bsc) {

                  }
               }
            }
            
            //if (ph > 4000) {
            //nhits++;
            //printf("samples %d %d, ", e->waveform[i].size(), w->nsamples);
            //}

            hit_time = le * time_bin + time_offset;
            hit_amp = ph * adc_gain + adc_offset;

            fH->fHawHitTime0->Fill(hit_time);

            //if (le < 342) {
            //   special = true;
            //}

            if (1) {
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
      
      if (fFlags->fPrint || special) {
         printf("wire %3d: baseline mean %8.1f, rms %4.1f, range %8.1f %6.1f, pulse %6.1f, le %6.1f, wi %3.1f, time %5.0f, amp %6.1f, flags: baseline %d, pulse %d, hit %d\n", i, bmean, brms, wmin, wmax, ph, le, wi, hit_time, hit_amp, have_baseline, have_pulse, have_hit);
      }

      return special;
   }
   
   TAFlowEvent* AnalyzeFlowEvent(TARunInfo* runinfo, TAFlags* flags, TAFlowEvent* flow)
   {
      //if (fTrace)
      //printf("AdcModule::Analyze, run %d, event serno %d, id 0x%04x, data size %d\n", runinfo->fRunNo, event->serial_number, (int)event->event_id, event->data_size);

      AgEventFlow *ef = flow->Find<AgEventFlow>();

      if (!ef || !ef->fEvent)
         return flow;

      Alpha16Event* e = ef->fEvent->a16;

      if (!e) {
         return flow;
      }

      if (fFlags->fPrint) {
         printf("Have ADC event:  ");
         e->Print();
         printf("\n");
      }

      if (fPlotAwWaveforms) {
         fPlotAwWaveforms->Draw(runinfo, e->hits);
      }

      if ((fFlags->fAdcPlotScaledown == 1) || ((ef->fEvent->counter % fFlags->fAdcPlotScaledown) == 0)) {
         for (unsigned i=0; i<fPlotA16.size(); i++) {
            fPlotA16[i]->Draw(e);
         }
      }

      AgAwHitsFlow* flow_hits = new AgAwHitsFlow(flow);
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
         bool special = AnalyzeHit(runinfo, e->udp[i], e->hits[i], &flow_hits->fAwHits);
         if (special)
            *flags |= TAFlag_DISPLAY;
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

class AdcModuleFactory: public TAFactory
{
public:
   AdcFlags fFlags;
   
public:
   void Usage()
   {
      printf("AdcModuleFactory flags:\n");
      printf("--pulser # enable field-wire pulser analysis\n");
      printf("--adcprint\n");
      printf("--adcfft\n");
      printf("--adc16 <imodule> # plot adc16 waveforms\n");
      printf("--adc32 <imodule> # plot adc32 waveforms\n");
      printf("--adcsd # adc plot scale down\n");
      printf("--adcfwf # filter waveform\n");
      printf("--adcinv # invert waveform\n");
      printf("--aw NNN # plot waveform for anode wire NNN\n");
      printf("--adc-wf-verbose  -- enable waveform suppression code verbose output\n");
      printf("--adc-wf-suppress -- enable waveform suppression code\n");
      printf("--adc-wf-threshold -- set the waveform suppression threshold\n");
      printf("--adc-wf-save-bad  -- write bad waveforms to the root output file\n");
   }

   void Init(const std::vector<std::string> &args)
   {
      printf("AdcModuleFactory::Init!\n");
      
      for (unsigned i=0; i<args.size(); i++) {
         if (args[i] == "--adcprint")
            fFlags.fPrint = true;
         if (args[i] == "--pulser")
            fFlags.fPulser = true;
         if (args[i] == "--adcfft")
            fFlags.fFft = true;
         if (args[i] == "--adc16") {
            int imodule = atoi(args[i+1].c_str());
            i++;
            fFlags.fPlotAdc16.push_back(imodule);
         }
         if (args[i] == "--adc32") {
            int imodule = atoi(args[i+1].c_str());
            i++;
            fFlags.fPlotAdc32.push_back(imodule);
         }
         if (args[i] == "--adcsd") {
            fFlags.fAdcPlotScaledown = atoi(args[i+1].c_str());
            i++;
         }
         if (args[i] == "--adcfwf") {
            fFlags.fFilterWaveform = true;
            i++;
         }
         if (args[i] == "--adcinv") {
            fFlags.fInvertWaveform = true;
            i++;
         }
         if (args[i] == "--aw") {
            int iwire = atoi(args[i+1].c_str());
            i++;
            fFlags.fAwWires.push_back(iwire);
         }
         if (args[i] == "--adc-wf-suppress")
            fFlags.fWfSuppress = true;
         if (args[i] == "--adc-wf-verbose")
            fFlags.fWfVerbose = true;
         if (args[i] == "--adc16-wf-threshold")
            fFlags.fWfThresholdAdc16 = atoi(args[i+1].c_str());
         if (args[i] == "--adc32-wf-threshold")
            fFlags.fWfThresholdAdc32 = atoi(args[i+1].c_str());
         if (args[i] == "--adc16-wf-keep-more")
            fFlags.fWfKeepMoreAdc16 = atoi(args[i+1].c_str());
         if (args[i] == "--adc32-wf-keep-more")
            fFlags.fWfKeepMoreAdc32 = atoi(args[i+1].c_str());
         if (args[i] == "--adc-wf-save-bad")
            fFlags.fWfSaveBad = true;
      }
      
      TARootHelper::fgDir->cd(); // select correct ROOT directory
   }
   
   void Finish()
   {
      printf("AdcModuleFactory::Finish!\n");
   }
   
   TARunObject* NewRunObject(TARunInfo* runinfo)
   {
      printf("AdcModuleFactory::NewRunObject, run %d, file %s\n", runinfo->fRunNo, runinfo->fFileName.c_str());
      return new AdcModule(runinfo, &fFlags);
   }

public:
};

static TARegister tar(new AdcModuleFactory);

/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */
