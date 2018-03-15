//
// MIDAS analyzer example 2: ROOT analyzer
//
// K.Olchanski
//

#include <stdio.h>

#include "manalyzer.h"
#include "midasio.h"

#include <assert.h> // assert()

#include <vector>
#include <deque>
#include <iostream>

#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"

#include "FeamEVB.h"
#include "Unpack.h"
#include "AgFlow.h"

#define DELETE(x) if (x) { delete (x); (x) = NULL; }

#define MEMZERO(p) memset((p), 0, sizeof(p))

#define ADC_MIN -33000
#define ADC_MAX  33000

// adjusted for 12-bit range 0xFFF = 4095
#define ADC_BINS 410
#define ADC_RANGE 4100
#define ADC_RANGE_PED 200

#define ADC_MIN_ADC -2048
#define ADC_OVERFLOW 4099

#define ADC_BINS_PULSER 100
#define ADC_RANGE_PULSER 300

#define ADC_RANGE_RMS 50

#define ADC_RMS_FPN_MIN 0.100
#define ADC_RMS_FPN_MAX 6.000

#define ADC_PULSER_TIME 450

#define NUM_SEQSCA (3*80+79)

#define NUM_TIME_BINS 512
#define MAX_TIME_BINS 512
#define MAX_TIME_NS 8200

class FeamHistograms
{
public:
   TProfile* hbmean_prof  = NULL;
   TProfile* hbrms_prof   = NULL;
   TH1D*     hbrms_pads   = NULL;
   TH1D*     hbrms_fpn    = NULL;
   TProfile* hbrange_prof = NULL;
   TH1D*     h_fpn_shift[4] = { NULL, NULL, NULL, NULL };
   TH1D*     h_amp = NULL;
   TH1D* hnhitchan = NULL;
   TH1D* h_nhitchan_seqsca = NULL;
   TH1D* h_spike_seqsca = NULL;
   TH1D* h_nhits_seqsca = NULL;
   TH1D* h_nhits_seqpad = NULL;
   TH1D* hnhits_pad_nospike = NULL;
   TH1D* hnhits_pad_drift = NULL;
   TH2D* h_hit_time_seqsca = NULL;
   TH2D* h_hit_amp_seqsca = NULL;
   TH2D* h_hit_amp_seqpad = NULL;
   TProfile* h_amp_seqsca = NULL;
   TProfile* h_amp_seqpad = NULL;
   TH1D* h_pulser_hit_amp = NULL;
   TH1D* h_pulser_hit_time = NULL;
   TH1D* h_pulser_hit_time_zoom = NULL;
   TH1D* h_pulser_hit_time_seqsca4_zoom = NULL;
   TProfile* h_pulser_hit_amp_seqpad = NULL;
   TProfile* h_pulser_hit_time_seqpad = NULL;
   TProfile* h_pulser_hit_time_seqsca = NULL;

public:
   FeamHistograms()
   {
   };

   void CreateHistograms(int position, int nbins, bool pulser)
   {
      char name[256];
      char title[256];

      sprintf(name,  "pos%02d_baseline_mean_prof", position);
      sprintf(title, "feam pos %2d baseline mean vs (SCA*80 + readout index)", position);
      hbmean_prof = new TProfile(name, title, NUM_SEQSCA, 0.5, NUM_SEQSCA+0.5);

      sprintf(name,  "pos%02d_baseline_rms_prof", position);
      sprintf(title, "feam pos %2d baseline rms vs (SCA*80 +  readout index)", position);
      hbrms_prof  = new TProfile(name, title,  NUM_SEQSCA, 0.5, NUM_SEQSCA+0.5);

      sprintf(name,  "pos%02d_baseline_rms_pads", position);
      sprintf(title, "feam pos %2d baseline rms for pad channels", position);
      hbrms_pads  = new TH1D(name, title,  100, 0, ADC_RANGE_RMS);

      sprintf(name,  "pos%02d_baseline_rms_fpn", position);
      sprintf(title, "feam pos %2d baseline rms fpr fpn channels", position);
      hbrms_fpn   = new TH1D(name, title,  100, 0, ADC_RANGE_RMS);

      sprintf(name,  "pos%02d_baseline_range_prof", position);
      sprintf(title, "feam pos %2d baseline range (max-min) vs (SCA*80 +  readout index)", position);
      hbrange_prof  = new TProfile(name, title,  NUM_SEQSCA, 0.5, NUM_SEQSCA+0.5);

      sprintf(name,  "pos%02d_fpn_shift_sca0", position);
      sprintf(title, "feam pos %2d fpn shift, sca 0", position);
      h_fpn_shift[0]  = new TH1D(name, title,  41, -20, 20);

      sprintf(name,  "pos%02d_fpn_shift_sca1", position);
      sprintf(title, "feam pos %2d fpn shift, sca 1", position);
      h_fpn_shift[1]  = new TH1D(name, title,  41, -20, 20);

      sprintf(name,  "pos%02d_fpn_shift_sca2", position);
      sprintf(title, "feam pos %2d fpn shift, sca 2", position);
      h_fpn_shift[2]  = new TH1D(name, title,  41, -20, 20);

      sprintf(name,  "pos%02d_fpn_shift_sca3", position);
      sprintf(title, "feam pos %2d fpn shift, sca 3", position);
      h_fpn_shift[3]  = new TH1D(name, title,  41, -20, 20);

      sprintf(name,  "pos%02d_amp", position);
      sprintf(title, "feam pos %2d waveform amplitude from baseline to minimum", position);
      h_amp = new TH1D(name, title,  ADC_BINS, 0, ADC_RANGE);

      sprintf(name,  "pos%02d_nhitchan", position);
      sprintf(title, "feam pos %2d number of hit channels", position);
      hnhitchan = new TH1D(name, title, 100, 0, 100);

      sprintf(name,  "pos%02d_nhitchan_map", position);
      sprintf(title, "feam pos %2d hit channels vs (SCA*80 + readout index)", position);
      h_nhitchan_seqsca = new TH1D(name, title, NUM_SEQSCA, 0.5, NUM_SEQSCA+0.5);

      sprintf(name,  "pos%02d_hit_map", position);
      sprintf(title, "feam pos %2d hits vs (SCA*80 + readout index)", position);
      h_nhits_seqsca = new TH1D(name, title, NUM_SEQSCA, 0.5, NUM_SEQSCA+0.5);

      sprintf(name,  "pos%02d_spike_seqsca", position);
      sprintf(title, "feam pos %2d spikes vs (SCA*80 + readout index)", position);
      h_spike_seqsca = new TH1D(name, title, NUM_SEQSCA, 0.5, NUM_SEQSCA+0.5);

      sprintf(name,  "pos%02d_hit_map_pads", position);
      sprintf(title, "feam pos %2d hits vs TPC seq.pad (col*4*18+row)", position);
      h_nhits_seqpad  = new TH1D(name, title, MAX_FEAM_PAD_ROWS*MAX_FEAM_PAD_COL, -0.5, MAX_FEAM_PAD_ROWS*MAX_FEAM_PAD_COL-0.5);

      sprintf(name,  "pos%02d_hit_map_pads_nospike", position);
      sprintf(title, "feam pos %2d hits with spikes removed vs TPC seq.pad (col*4*18+row)", position);
      hnhits_pad_nospike = new TH1D(name, title, MAX_FEAM_PAD_ROWS*MAX_FEAM_PAD_COL, -0.5, MAX_FEAM_PAD_ROWS*MAX_FEAM_PAD_COL-0.5);

      sprintf(name,  "pos%02d_hit_map_pads_drift", position);
      sprintf(title, "feam pos %2d hits in drift region vs TPC seq.pad (col*4*18+row)", position);
      hnhits_pad_drift = new TH1D(name, title, MAX_FEAM_PAD_ROWS*MAX_FEAM_PAD_COL, -0.5, MAX_FEAM_PAD_ROWS*MAX_FEAM_PAD_COL-0.5);

      sprintf(name,  "pos%02d_hit_time_seqsca", position);
      sprintf(title, "feam pos %2d hit time vs (SCA*80 + readout index)", position);
      h_hit_time_seqsca = new TH2D(name, title, NUM_SEQSCA, 0.5, NUM_SEQSCA+0.5, 50, 0, 500);

      sprintf(name,  "pos%02d_hit_amp_seqsca", position);
      sprintf(title, "feam pos %2d hit p.h. vs (SCA*80 + readout index)", position);
      h_hit_amp_seqsca = new TH2D(name, title, NUM_SEQSCA, 0.5, NUM_SEQSCA+0.5, 50, 0, ADC_RANGE);

      sprintf(name,  "pos%02d_hit_amp_seqpad", position);
      sprintf(title, "feam pos %2d hit p.h. vs TPC seq.pad (col*4*18+row)", position);
      h_hit_amp_seqpad = new TH2D(name, title, MAX_FEAM_PAD_ROWS*MAX_FEAM_PAD_COL, -0.5, MAX_FEAM_PAD_ROWS*MAX_FEAM_PAD_COL-0.5, 50, 0, ADC_RANGE);

      sprintf(name,  "pos%02d_amp_seqsca", position);
      sprintf(title, "feam pos %2d hit p.h. profile, cut 10000..40000 vs (SCA*80 + readout index)", position);
      h_amp_seqsca = new TProfile(name, title, NUM_SEQSCA, -0.5, NUM_SEQSCA-0.5);

      sprintf(name,  "pos%02d_amp_seqpad", position);
      sprintf(title, "feam pos %2d hit p.h. profile, cut 10000..40000 vs TPC seq.pad (col*4*18+row)", position);
      h_amp_seqpad = new TProfile(name, title, MAX_FEAM_PAD_ROWS*MAX_FEAM_PAD_COL+1, -1.5, MAX_FEAM_PAD_ROWS*MAX_FEAM_PAD_COL-0.5);

      if (pulser) {
         sprintf(name,  "pos%02d_pulser_hit_amp", position);
         sprintf(title, "feam pos %2d pulser hit p.h.", position);
         h_pulser_hit_amp = new TH1D(name, title, ADC_BINS_PULSER, 0, ADC_RANGE_PULSER);

         sprintf(name,  "pos%02d_pulser_hit_time", position);
         sprintf(title, "feam pos %2d pulser hit time", position);
         h_pulser_hit_time = new TH1D(name, title, nbins, 0, nbins);

         sprintf(name,  "pos%02d_pulser_hit_time_zoom", position);
         sprintf(title, "feam pos %2d pulser hit time zoom", position);
         h_pulser_hit_time_zoom = new TH1D(name, title, 100, ADC_PULSER_TIME-10, ADC_PULSER_TIME+10);

         sprintf(name,  "pos%02d_pulser_hit_time_seqsca4_zoom", position);
         sprintf(title, "feam pos %2d pulser hit time seqsca 4 zoom", position);
         h_pulser_hit_time_seqsca4_zoom = new TH1D(name, title, 100, ADC_PULSER_TIME-10, ADC_PULSER_TIME+10);

         sprintf(name,  "pos%02d_pulser_hit_amp_seqpad", position);
         sprintf(title, "feam pos %2d pulser hit p.h. vs TPC seq.pad (col*4*18+row)", position);
         h_pulser_hit_amp_seqpad = new TProfile(name, title, MAX_FEAM_PAD_ROWS*MAX_FEAM_PAD_COL+1, -1.5, MAX_FEAM_PAD_ROWS*MAX_FEAM_PAD_COL-0.5);

         sprintf(name,  "pos%02d_pulser_hit_time_seqpad", position);
         sprintf(title, "feam pos %2d pulser hit time vs TPC seq.pad (col*4*18+row)", position);
         h_pulser_hit_time_seqpad = new TProfile(name, title, MAX_FEAM_PAD_ROWS*MAX_FEAM_PAD_COL, -0.5, MAX_FEAM_PAD_ROWS*MAX_FEAM_PAD_COL-0.5);

         sprintf(name,  "pos%02d_pulser_hit_time_seqsca", position);
         sprintf(title, "feam pos %2d pulser hit time vs (SCA*80 + readout index)", position);
         h_pulser_hit_time_seqsca = new TProfile(name, title, NUM_SEQSCA, -0.5, NUM_SEQSCA-0.5);
      }
   }
};

static int find_pulse(const int* adc, int nbins, double baseline, double gain, double threshold)
{
   for (int i=0; i<nbins; i++) {
      if ((adc[i]-baseline)*gain > threshold) {
         return i;
      }
   }

   return 0;
}

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

class ChanHistograms
{
public:
   std::string fNameBase;
   std::string fTitleBase;
   int fNbins = 0;
   TDirectory* fDirBad = NULL;

   TH1D* hbmean = NULL;
   TH1D* hbrms  = NULL;

   TH1D* hwaveform_first = NULL;
   TH1D* hwaveform_bad   = NULL;
   TH1D* hwaveform_max   = NULL;
   TH1D* hwaveform_max_drift = NULL;
   TH1D* hwaveform_avg   = NULL;
   TH1D* hwaveform_avg_drift = NULL;

   int nwf = 0;
   int nwf_drift = 0;

   double fMaxWamp = 0;
   double fMaxWampDrift = 0;


public:
   ChanHistograms(const char* xname, const char* xtitle, TDirectory* dir, int nbins) // ctor
   {
      TDirectory* dir_first = dir->GetDirectory("pchan_waveform_first");
      if(!dir_first) dir_first = dir->mkdir("pchan_waveform_first");
      TDirectory* dir_bad = dir->GetDirectory("pchan_waveform_bad");
      if(!dir_bad) dir_bad = dir->mkdir("pchan_waveform_bad");
      TDirectory* dir_max = dir->GetDirectory("pchan_waveform_max");
      if(!dir_max) dir_max = dir->mkdir("pchan_waveform_max");
      TDirectory* dir_max_drift = dir->GetDirectory("pchan_waveform_max_drift");
      if(!dir_max_drift) dir_max_drift = dir->mkdir("pchan_waveform_max_drift");
      TDirectory* dir_avg = dir->GetDirectory("pchan_waveform_avg");
      if(!dir_avg) dir_avg = dir->mkdir("pchan_waveform_avg");
      TDirectory* dir_avg_drift = dir->GetDirectory("pchan_waveform_avg_drift");
      if(!dir_avg_drift) dir_avg_drift = dir->mkdir("pchan_waveform_avg_drift");

      fNameBase = xname;
      fTitleBase = xtitle;
      fNbins = nbins;
      fDirBad = dir_bad;

      char name[256];
      char title[256];

      sprintf(name, "hpwf_first_%s", xname);
      sprintf(title, "%s first waveform", xtitle);

      dir_first->cd();
      hwaveform_first = new TH1D(name, title, nbins, -0.5, nbins-0.5);

      sprintf(name, "hpwf_max_%s", xname);
      sprintf(title, "%s biggest waveform", xtitle);
      dir_max->cd();
      hwaveform_max = new TH1D(name, title, nbins, -0.5, nbins-0.5);

      sprintf(name, "hpwf_max_drift_%s", xname);
      sprintf(title, "%s biggest waveform, drift region", xtitle);
      dir_max_drift->cd();
      hwaveform_max_drift = new TH1D(name, title, nbins, -0.5, nbins-0.5);

      sprintf(name, "hpwf_avg_%s", xname);
      sprintf(title, "%s average waveform", xtitle);
      dir_avg->cd();
      hwaveform_avg = new TH1D(name, title, nbins, -0.5, nbins-0.5);

      sprintf(name, "hpwf_avg_drift_%s", xname);
      sprintf(title, "%s average waveform, drift region", xtitle);
      dir_avg_drift->cd();
      hwaveform_avg_drift = new TH1D(name, title, nbins, -0.5, nbins-0.5);
   }

   ~ChanHistograms() // dtor
   {

   }

   bool SaveBad(int nbins, const int* adc)
   {
      if (hwaveform_bad == NULL) {
         char name[256];
         char title[256];
         sprintf(name, "hpwf_bad_%s", fNameBase.c_str());
         sprintf(title, "%s bad waveform", fTitleBase.c_str());
         fDirBad->cd();
         hwaveform_bad = new TH1D(name, title, nbins, -0.5, nbins-0.5);

         for (int i=0; i<nbins; i++)
            hwaveform_bad->SetBinContent(i+1, adc[i]);

         return true;
      }

      return false;
   }
};

class FeamFlags
{
public:
   bool fDoPads = true;
   int  fPlotPad = -1;
   TCanvas* fPlotPadCanvas = NULL;
   bool fExportWaveforms = false;

public:
   FeamFlags() // ctor
   {
   }

   ~FeamFlags() // dtor
   {
      DELETE(fPlotPadCanvas);
   }
};

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

static double compute_rms(const int* aptr, int start, int end)
{
   double mean, rms;
   compute_mean_rms(aptr, start, end, &mean, &rms, NULL, NULL);
   return rms;
}

bool fpn_rms_ok(int ichan, double brms)
{
   if (ichan < 4)
      return true;

   if (brms > ADC_RMS_FPN_MIN && brms < ADC_RMS_FPN_MAX)
      return true;

   return false;
}

int fpn_wrap(int ifpn)
{
   while (ifpn < 0)
      ifpn += 80;

   while (ifpn >= 80)
      ifpn -= 80;

   return ifpn;
}

class FeamModule: public TARunObject
{
private:
   const padMap padMapper;
public:
   FeamFlags* fFlags = NULL;
   //FILE *fin = NULL;
   //TCanvas* fC = NULL;

   TH1D* h_spike_diff = NULL;
   TH1D* h_spike_diff_max = NULL;
   TH1D* h_spike_num = NULL;

   TH1D* hbmean_all;
   TH1D* hbrms_all;
   TH1D* hbrms_all_pads;
   TH1D* hbrms_all_fpn;

   TH1D* h_adc_range_all = NULL;
   TH1D* h_adc_range_baseline = NULL;
   TH1D* h_adc_range_drift = NULL;

   TH1D* hamp_pad;
   TH1D* hamp_pad_pedestal;
   TH1D* hamp_pad_above_pedestal;
   TH1D* hled_pad_amp;
   TH1D* hled_pad_amp_ns;

   TH1D* hled_all_hits;
   TH1D* hamp_all_hits;

   TH2D* h2led2amp;

   TH1D* hdrift_amp_all;
   TH1D* hdrift_amp_all_pedestal;
   TH1D* hdrift_amp_all_above_pedestal;
   TH1D* hdrift_led_all;
   TH2D* hdrift_led2amp;

   TH1D* hnhits;
   TH1D* hled_hit;
   TH1D* hamp_hit;
   TH2D* h_amp_hit_col = NULL;

   bool  fPulser = true;
   TH1D* h_pulser_led_hit = NULL;

   TH1D* hnhitchan = NULL;

   TH2D* hpadmap;

   TDirectory* hdir_summary;
   TDirectory* hdir_feam;
   TDirectory* hdir_pads;
   std::vector<FeamHistograms*> fHF;
   std::vector<ChanHistograms*> fHC;

   int fCountTestScaEvents = 0;
   int fCountBadScaEvents = 0;
   int fCountBadSca = 0;
   int fCountGoodFpn = 0;
   int fCountBadFpn = 0;

   bool fTrace = false;

   FeamModule(TARunInfo* runinfo, FeamFlags* f)
      : TARunObject(runinfo)
   {
      if (fTrace)
         printf("FeamModule::ctor!\n");

      fFlags = f;

      hbmean_all = NULL;

      if (!fFlags->fDoPads)
         return;

      //fC = new TCanvas();
   }

   ~FeamModule()
   {
      if (fTrace)
         printf("FeamModule::dtor!\n");
      //DELETE(fC);
      for (unsigned i=0; i<fHF.size(); i++) {
         DELETE(fHF[i]);
      }
      for (unsigned i=0; i<fHC.size(); i++) {
         DELETE(fHC[i]);
      }
   }

   void BeginRun(TARunInfo* runinfo)
   {
      if (fTrace)
         printf("BeginRun, run %d, file %s\n", runinfo->fRunNo, runinfo->fFileName.c_str());
      //time_t run_start_time = runinfo->fOdb->odbReadUint32("/Runinfo/Start time binary", 0, 0);
      //printf("ODB Run start time: %d: %s", (int)run_start_time, ctime(&run_start_time));
      runinfo->fRoot->fOutputFile->cd(); // select correct ROOT directory

      if (!fFlags->fDoPads)
         return;
   }

   void CreateHistograms(TARunInfo* runinfo, int nfeam, int nchan_feam, int nchan, int nbins)
   {
      if (hbmean_all) // already created
         return;

      runinfo->fRoot->fOutputFile->cd();
      hdir_pads = gDirectory->mkdir("pads");
      hdir_pads->cd(); // select correct ROOT directory

      hdir_summary = hdir_pads->mkdir("summary");
      hdir_feam = hdir_pads->mkdir("feam");

      hdir_summary->cd();

      h_spike_diff = new TH1D("spike_diff", "channel spike finder", 100, 0, 1000);
      h_spike_diff_max = new TH1D("spike_diff_max", "channel spike finder, max", 100, 0, 1000);
      h_spike_num = new TH1D("spike_num", "channel spike finder, num", 100, 0-0.5, 100-0.5);

      hbmean_all = new TH1D("hbaseline_mean", "baseline mean", 100, ADC_MIN, ADC_MAX);
      hbrms_all  = new TH1D("hbaseline_rms",  "baseline rms",  100, 0, ADC_RANGE_RMS);
      hbrms_all_pads  = new TH1D("hbaseline_rms_pads",  "baseline rms, tpc pad channels",  100, 0, ADC_RANGE_RMS);
      hbrms_all_fpn   = new TH1D("hbaseline_rms_fpn",  "baseline rms, fpn channels",  100, 0, ADC_RANGE_RMS);

      h_adc_range_all      = new TH1D("adc_range_all",      "waveform range (max-min)",  100, 0, ADC_RANGE_PED);
      h_adc_range_baseline = new TH1D("adc_range_baseline", "waveform range (max-min), baseline region",  100, 0, ADC_RANGE_PED);
      h_adc_range_drift    = new TH1D("adc_range_drift",    "waveform range (max-min), drift region",  100, 0, ADC_RANGE_PED);

      hamp_pad         = new TH1D("hamp_pad",   "pad channels pulse height; adc counts", 100, 0, ADC_RANGE);
      hamp_pad_pedestal = new TH1D("hamp_pad_pedestal", "pad channels pulse height, zoom on pedestal area; adc counts", 100, 0, ADC_RANGE_PED);
      hamp_pad_above_pedestal = new TH1D("hamp_pad_above_pedestal", "pad channels pulse height, away from pedestal area; adc counts", 100, ADC_RANGE_PED, ADC_RANGE);

      hled_pad_amp = new TH1D("hled_pad_amp",   "pad channels above threshold, pulse leading edge; adc time bins", 100, 0, MAX_TIME_BINS);
      hled_pad_amp_ns = new TH1D("hled_pad_amp_ns",   "pad channels above threshold, pulse leading edge; time, ns", 100, 0, MAX_TIME_NS);

      h2led2amp  = new TH2D("h2led2amp", "pulse amp vs time, adc time bins", 100, 0, MAX_TIME_BINS, 100, 0, ADC_RANGE);

      hled_all_hits = new TH1D("hled_all_hits",   "pulse leading edge, adc time bins, with p.h. cut", 100, 0, nbins);
      hamp_all_hits = new TH1D("hamp_all_hits",   "pulse height, with time cut", 100, 0, ADC_RANGE);

      hdrift_amp_all = new TH1D("drift_amp", "drift region pulse height", 100, 0, ADC_RANGE);
      hdrift_amp_all_pedestal = new TH1D("drift_amp_pedestal", "drift region pulse height, zoom on pedestal area", 100, 0, ADC_RANGE_PED);
      hdrift_amp_all_above_pedestal = new TH1D("drift_amp_above_pedestal", "drift region pulse height, away from pedestal area", 100, ADC_RANGE_PED, ADC_RANGE);
      hdrift_led_all = new TH1D("drift_led", "drift region pulse leading edge, adc time bins, above pedestal", 100, 0, nbins);
      hdrift_led2amp = new TH2D("drift_led2amp", "drift region pulse amp vs time, adc time bins, above pedestal", 100, 0, nbins, 100, 0, ADC_RANGE);

      hnhits = new TH1D("hnhits", "hits per channel", nchan, -0.5, nchan-0.5);
      hled_hit = new TH1D("hled_hit", "hit time, adc time bins", 100, 0, nbins);
      hamp_hit = new TH1D("hamp_hit", "hit pulse height", 100, 0, ADC_RANGE);
      h_amp_hit_col = new TH2D("hamp_hit_col", "hit pulse height vs column (ifeam*4+col)", 8*4, -0.5, 8*4-0.5, 100, 0, ADC_RANGE);

      hnhitchan = new TH1D("hnhitchan", "number of hit channels per event", 100, 0, 1000);

      if (fPulser) {
         h_pulser_led_hit = new TH1D("pulser_led_hit", "pulser time, adc bins", 30, 180-0.5, 210-0.5);
      }

      hpadmap = new TH2D("hpadmap", "map from TPC pad number (col*4*18+row) to SCA readout channel (sca*80+chan)", 4*4*18, -0.5, 4*4*18-0.5, NUM_SEQSCA, 0.5, NUM_SEQSCA+0.5);

      hdir_feam->cd();

      for (int i=0; i<nfeam; i++) {
         FeamHistograms *fh = new FeamHistograms();
         fh->CreateHistograms(i, nbins, fPulser);
         fHF.push_back(fh);
      }
   }

   void EndRun(TARunInfo* runinfo)
   {
      if (fTrace)
         printf("FeamModule::EndRun, run %d\n", runinfo->fRunNo);
      //time_t run_stop_time = runinfo->fOdb->odbReadUint32("/Runinfo/Stop time binary", 0, 0);
      //printf("ODB Run stop time: %d: %s", (int)run_stop_time, ctime(&run_stop_time));
      for(auto *hc: fHC){
         if(hc){
            if(hc->nwf) hc->hwaveform_avg->Scale(1./double(hc->nwf));
            if(hc->nwf_drift) hc->hwaveform_avg_drift->Scale(1./double(hc->nwf_drift));
         }
      }

      printf("FeamModule::EndRun: test for bad SCA: total events %d, bad events %d, bad sca %d, bad fpn %d, good fpn %d\n", fCountTestScaEvents, fCountBadScaEvents, fCountBadSca, fCountBadFpn, fCountGoodFpn);
   }

   void PauseRun(TARunInfo* runinfo)
   {
      if (fTrace)
         printf("PauseRun, run %d\n", runinfo->fRunNo);
   }

   void ResumeRun(TARunInfo* runinfo)
   {
      if (fTrace)
         printf("ResumeRun, run %d\n", runinfo->fRunNo);
   }

   TAFlowEvent* AnalyzeFlowEvent(TARunInfo* runinfo, TAFlags* flags, TAFlowEvent* flow)
   {
      bool verbose = false;
      
      //printf("Analyze, run %d, event serno %d, id 0x%04x, data size %d\n", runinfo->fRunNo, event->serial_number, (int)event->event_id, event->data_size);

      if (!fFlags->fDoPads)
         return flow;

      AgEventFlow *ef = flow->Find<AgEventFlow>();

      if (!ef || !ef->fEvent)
         return flow;

      FeamEvent* e = ef->fEvent->feam;

      if (!e) {
         return flow;
      }

      //int force_plot = false;

#if 0
      if (event->event_id != 1 && event->event_id != 2)
         return flow;
#endif

#if 0
      int adc[80][5120];

      const int xbins = 829;
      const int xchan = 79;
#endif

      if (!e) {
         return flow;
      }

      if (e) {
         if (1) {
            printf("Have FEAM event: ");
            e->Print();
            printf("\n");
         }

         if (0) {
            for (unsigned i=0; i<e->modules.size(); i++) {
               printf("FeamEvent slot %d: ", i);
               if (!e->modules[i]) {
                  printf("null\n");
                  continue;
               }
               FeamModuleData* m = e->modules[i];
               printf("pwb%2d, cnt %4d, ts_trig: 0x%08x %14.3f usec, ts_incr %14.3f usec\n", m->fModule, m->cnt, m->ts_trig, m->fTime*1e6, m->fTimeIncr*1e6);
            }
         }

         if (e->error) {
            //delete e;
            return flow;
         }

         if (0) {
            //delete e;
            return flow;
         }

         if (0 && !e->complete) {
            return flow;
         }

         //assert(a->next_n == 256);
         //assert(a->fSize == 310688);

         //MEMZERO(adc);


         //for (int ibin = 511; ibin < xbins; ibin++) {
         //   for (int ichan = 0; ichan < 76; ichan++) {
         //      adc[ichan][ibin] = adc[ichan][ibin-511];
         //   }
         //}

         //delete e;

#if 0
      } else if (event->event_id == 2) {

         const char* banks[] = { "YP01", "YP02", NULL };
         int itpc = -1;
         unsigned short *samples = NULL;

         for (int i=0; banks[i]; i++) {
            TMBank* b = event->FindBank(banks[i]);
            if (b) {
               samples = (unsigned short*)event->GetBankData(b);
               if (samples) {
                  itpc = i;
                  break;
               }
            }
         }

         printf("itpc %d, samples 0x%p\n", itpc, samples);

         if (itpc < 0 || samples == NULL) {
            return flow;
         }

#if 0
         int count = 0;
         for (int ibin=0; ibin<xbins; ibin++) {
            for (int ichan=0; ichan<xchan; ichan++) {
               adc[ichan][ibin] = samples[count];
               count++;
            }
         }
         printf("got %d samples\n", count);
#endif

      } else if (fin) {
         // good stuff goes here

         char buf[4*1024*1024];

         char *s = fgets(buf, sizeof(buf), fin);
         if (s == NULL) {
            *flags |= TAFlag_QUIT;
            return flow;
         }
         printf("read %d\n", (int)strlen(s));

         int event_no = strtoul(s, &s, 0);
         int t0 = strtoul(s, &s, 0);
         int t1 = strtoul(s, &s, 0);
         int t2 = strtoul(s, &s, 0);

         printf("event %d, t %d %d %d\n", event_no, t0, t1, t2);

#if 0
         int count = 0;
         for (int ibin=0; ibin<xbins; ibin++) {
            for (int ichan=0; ichan<xchan; ichan++) {
               count++;
               adc[ichan][ibin] = strtoul(s, &s, 0);
            }
         }

         printf("got %d samples\n", count);

         for (int i=0; ; i++) {
            if (!*s)
               break;
            int v = strtoul(s, &s, 0);
            if (v == 0)
               break;
            count++;
         }

         printf("total %d samples before zeros\n", count);

         for (int i=0; ; i++) {
            if (s[0]==0)
               break;
            if (s[0]=='\n')
               break;
            if (s[0]=='\r')
               break;
            int v = strtoul(s, &s, 0);
            if (v != 0)
               break;
            count++;
            if (*s == '+')
               s++;
         }

         printf("total %d samples with zeros\n", count);
#endif

         s[100] = 0;
         printf("pads data: [%s]\n", s);

         char buf1[1024];

         char *s1 = fgets(buf1, sizeof(buf1), fin);
         if (s1 == NULL) {
            *flags |= TAFlag_QUIT;
            return flow;
         }
         printf("read %d [%s]\n", (int)strlen(s1), s1);

         //int event_no = strtoul(s, &s, 0);
         //int t0 = strtoul(s, &s, 0);
         //int t1 = strtoul(s, &s, 0);
         //int t2 = strtoul(s, &s, 0);
         //printf("event %d, t %d %d %d\n", event_no, t0, t1, t2);
#endif
      }
      
      //

      bool doPrint = false;

      // got all the data here

      int nfeam = e->adcs.size();
      int nchan_feam = 0;
      int nbins = 0;

      for (unsigned i=0; i<e->adcs.size(); i++) {
         if (e->adcs[i]) {
            nchan_feam = e->adcs[i]->nsca * e->adcs[i]->nchan;
            nbins = e->adcs[i]->nbins;
            break;
         }
      }

      int nchan = nfeam * nchan_feam; // MAX_FEAM_SCA * ...;

      if (nbins == 0 || nchan == 0)
         return flow;

      // create histograms

      CreateHistograms(runinfo, nfeam, nchan_feam, nchan, nbins);

      // create pad hits flow event

      AgPadHitsFlow* hits = new AgPadHitsFlow(flow);
      flow = hits;

      // check for bad sca data

      bool bad_sca = false;

      fCountTestScaEvents++;

      for (unsigned ifeam=0; ifeam<e->adcs.size(); ifeam++) {
         FeamAdcData* aaa = e->adcs[ifeam];
         if (!aaa)
            continue;

         for (int isca=0; isca<aaa->nsca; isca++) {
            int first_chan = 4;
            int v = aaa->adc[isca][first_chan][0];

            bool bad = (v == 0x7FFC);

            if (bad) {
               printf("Error: FEAM pos %02d sca %d has gibberish data: ", ifeam, isca);
               e->Print();
               printf("\n");
            }

            if (bad) {
               fCountBadSca++;
               bad_sca = true;
            }
         }
      }

      if (bad_sca) {
         fCountBadScaEvents++;
         e->error = true;
         return flow;
      }

      // loop over all waveforms

      //int iplot = 0;
      double zmax = 0;
      bool first_zero_range = true;

      int ibaseline_start = 10;
      int ibaseline_end = 100;

      int iwire_start = 130;
      int iwire_end = 160;

      int idrift_start = iwire_end;
      int idrift_cut = iwire_end;
      int idrift_end = 410;

      int ipulser_start = 400;
      int ipulser_end   = 500;

      double hit_amp_threshold = 100;

      int nhitchan = 0;

      // check FPN channels

      for (unsigned ifeam=0; ifeam<e->adcs.size(); ifeam++) {
         FeamAdcData* aaa = e->adcs[ifeam];
         if (!aaa)
            continue;

         bool trace = true;
         bool trace_shift = false;

         for (int isca=0; isca<aaa->nsca; isca++) {

            int fpn_shift = 0;

            double rms_fpn1 = compute_rms(aaa->adc[isca][16], ibaseline_start, ibaseline_end);
            double rms_fpn2 = compute_rms(aaa->adc[isca][29], ibaseline_start, ibaseline_end);
            double rms_fpn3 = compute_rms(aaa->adc[isca][54], ibaseline_start, ibaseline_end);
            double rms_fpn4 = compute_rms(aaa->adc[isca][67], ibaseline_start, ibaseline_end);

            if (fpn_rms_ok(16, rms_fpn1)
                && fpn_rms_ok(29, rms_fpn2) 
                && fpn_rms_ok(54, rms_fpn3) 
                && fpn_rms_ok(67, rms_fpn4)) {

               if (trace) {
                  printf("XXX good fpn pos %2d, sca %d, fpn rms: %5.1f %5.1f %5.1f %5.1f\n", ifeam, isca, rms_fpn1, rms_fpn2, rms_fpn3, rms_fpn4);
               }

               fpn_shift = 0;
            } else {
               if (trace) {
                  printf("XXX bad  fpn pos %2d, sca %d, fpn rms: %5.1f %5.1f %5.1f %5.1f\n", ifeam, isca, rms_fpn1, rms_fpn2, rms_fpn3, rms_fpn4);
               }

               for (int i=0; i>-30; i--) {
                  int ifpn1 = fpn_wrap(i+16);
                  int ifpn2 = fpn_wrap(i+29);
                  int ifpn3 = fpn_wrap(i+54);
                  int ifpn4 = fpn_wrap(i+67);

                  double rms_fpn1 = compute_rms(aaa->adc[isca][ifpn1], ibaseline_start, ibaseline_end);
                  double rms_fpn2 = compute_rms(aaa->adc[isca][ifpn2], ibaseline_start, ibaseline_end);
                  double rms_fpn3 = compute_rms(aaa->adc[isca][ifpn3], ibaseline_start, ibaseline_end);
                  double rms_fpn4 = compute_rms(aaa->adc[isca][ifpn4], ibaseline_start, ibaseline_end);

                  if (trace_shift) {
                     printf("XXX shift %3d fpn pos %d, sca %d, fpn rms: %5.1f %5.1f %5.1f %5.1f, fpn chan %2d %2d %2d %2d", i, ifeam, isca, rms_fpn1, rms_fpn2, rms_fpn3, rms_fpn4, ifpn1, ifpn2, ifpn3, ifpn4);
                  }

                  if (fpn_rms_ok(ifpn1, rms_fpn1)
                      && fpn_rms_ok(ifpn2, rms_fpn2) 
                      && fpn_rms_ok(ifpn3, rms_fpn3) 
                      && fpn_rms_ok(ifpn4, rms_fpn4)) {
                     if (trace_shift) {
                        printf(", fpn ok!!!!\n");
                     }
                     fpn_shift = i;
                     break;
                  } else {
                     if (trace_shift) {
                        printf(", fpn bad\n");
                     }
                  }
               }
            }
            
            if (trace_shift) {
               if (fpn_shift != 0) {
                  printf("XXX pos %2d, sca %d, fpn_shift %d\n", ifeam, isca, fpn_shift);
               }
            }

            fHF[ifeam]->h_fpn_shift[isca]->Fill(fpn_shift);

            if (fpn_shift < 0) {
               int buf[MAX_FEAM_READOUT][MAX_FEAM_BINS];

               char* asrc = (char*)&aaa->adc[isca][0];
               char* adst = (char*)&buf[0][0];

               const int s = sizeof(int)*MAX_FEAM_BINS;

               for (int i=0; i<MAX_FEAM_READOUT; i++) {
                  int j = (i+fpn_shift+MAX_FEAM_READOUT)%MAX_FEAM_READOUT;

                  //printf("fpn_shift %d, copy %d from %d, size %d, ptr %p from %p\n", fpn_shift, i, j, s, adst+i*s, asrc+j*s);
                  
                  memcpy(adst+i*s, asrc+j*s, s);
               }

               memcpy(asrc, adst, sizeof(buf));

               if (1) {
                  bool trace = true;

                  double rms_fpn1 = compute_rms(aaa->adc[isca][16], ibaseline_start, ibaseline_end);
                  double rms_fpn2 = compute_rms(aaa->adc[isca][29], ibaseline_start, ibaseline_end);
                  double rms_fpn3 = compute_rms(aaa->adc[isca][54], ibaseline_start, ibaseline_end);
                  double rms_fpn4 = compute_rms(aaa->adc[isca][67], ibaseline_start, ibaseline_end);
                  
                  if (fpn_rms_ok(16, rms_fpn1)
                      && fpn_rms_ok(29, rms_fpn2) 
                      && fpn_rms_ok(54, rms_fpn3) 
                      && fpn_rms_ok(67, rms_fpn4)) {
                     
                     if (trace) {
                        printf("XXX good fpn pos %d, sca %d, fpn rms: %5.1f %5.1f %5.1f %5.1f, fpn_shift: %3d\n", ifeam, isca, rms_fpn1, rms_fpn2, rms_fpn3, rms_fpn4, fpn_shift);
                     }
                  } else {
                     if (trace) {
                        printf("XXX bad  fpn pos %d, sca %d, fpn rms: %5.1f %5.1f %5.1f %5.1f, fpn_shift: %3d\n", ifeam, isca, rms_fpn1, rms_fpn2, rms_fpn3, rms_fpn4, fpn_shift);
                     }
                  }
               }
            }
         }
      }

      for (unsigned ifeam=0; ifeam<e->adcs.size(); ifeam++) {
         FeamAdcData* aaa = e->adcs[ifeam];
         if (!aaa)
            continue;

         int nhitchan_feam = 0;

         bool fpn_is_ok = true;

         for (int isca=0; isca<aaa->nsca; isca++) {
            for (int ichan=1; ichan<=aaa->nchan; ichan++) {
               const int* aptr = &aaa->adc[isca][ichan][0];

               unsigned seqchan = ifeam*(aaa->nsca*aaa->nchan) + isca*aaa->nchan + ichan;
               int seqsca = isca*80 + ichan;

               // consult the pad map

               static bool once = true;
               if (once) {
                  once = false;
                  printf("Pad map:\n");
                  printf("  sca chan: ");
                  for (int i=0; i<=79; i++)
                     printf("%d ", padMapper.channel[i]);
                  printf("\n");
                  for (int sca=0; sca<4; sca++) {
                     printf("sca %d:\n", sca);
                     printf("  tpc col: ");
                     for (int i=0; i<=72; i++)
                        printf("%d ", padMapper.padcol[sca][i]);
                     printf("\n");
                     printf("  tpc row: ");
                     for (int i=0; i<=72; i++)
                        printf("%d ", padMapper.padrow[sca][i]);
                     printf("\n");
                  }

                  int test[4*4*18];
                  for (int i=0; i<4*4*18; i++)
                     test[i] = 0;

                  bool map_ok = true;

                  for (int sca=0; sca<4; sca++) {
                     for (int i=0; i<=79; i++) {
                        int seqsca = sca*80+i;
                        int chan = padMapper.channel[i];
                        if (chan > 0) {
                           int col = padMapper.padcol[sca][chan];
                           int row = padMapper.padrow[sca][chan];
                           int seqpad = col*4*18+row;
                           hpadmap->Fill(seqpad, seqsca);
                           if (test[seqpad] != 0) {
                              printf("pad map error: col %d, row %d, seqpad %d: duplicate mapping seqsca %d and %d\n", col, row, seqpad, seqsca, test[seqpad]);
                              map_ok = false;
                           } else {
                              test[seqpad] = seqsca;
                           }
                        }
                     }
                  }

                  for (int i=0; i<4*4*18; i++) {
                     if (test[i] == 0) {
                        printf("pad map error: seqpad %d is not mapped to sca channel!\n", i);
                        map_ok = false;
                     }
                  }

                  if (map_ok) {
                     printf("pad map is ok.\n");
                  } else {
                     printf("pad map has errors!\n");
                  }
               }

               int scachan = padMapper.channel[ichan];
               int col = -1; // TPC pad column
               int row = -1; // TPC pad row
               int seqpad = -1; // TPC sequential pad number col*4*72+row

               bool scachan_is_pad = (scachan > 0);
               bool scachan_is_fpn = (scachan >= -4) && (scachan <= -1);

               if (scachan_is_pad) {
                  col = padMapper.padcol[isca][scachan];
                  row = padMapper.padrow[isca][scachan];
                  //printf("isca %d, ichan %d, scachan %d, col %d, row %d\n", isca, ichan, scachan, col, row);
                  assert(col>=0 && col<4);
                  assert(row>=0 && row<4*72);
                  seqpad = col*MAX_FEAM_PAD_ROWS + row;
               } else {
                  row = scachan; // special channel
               }

               char xname[256];
               char xtitle[256];

               if (scachan_is_pad) {
                  sprintf(xname, "pos%02d_%03d_sca%d_chan%02d_scachan%02d_col%02d_row%02d", ifeam, seqsca, isca, ichan, scachan, col, row);
                  sprintf(xtitle, "FEAM pos %d, sca %d, readout chan %d, sca chan %d, col %d, row %d", ifeam, isca, ichan, scachan, col, row);
               } else if (scachan_is_fpn) {
                  sprintf(xname, "pos%02d_%03d_sca%d_chan%02d_fpn%d", ifeam, seqsca, isca, ichan, -scachan);
                  sprintf(xtitle, "FEAM pos %d, sca %d, readout chan %d, fpn %d", ifeam, isca, ichan, -scachan);
               } else {
                  sprintf(xname, "pos%02d_%03d_sca%d_chan%02d", ifeam, seqsca, isca, ichan);
                  sprintf(xtitle, "FEAM pos %d, sca %d, readout chan %d", ifeam, isca, ichan);
               }

               if (fFlags->fExportWaveforms) {
                  TDirectory* dir = runinfo->fRoot->fgDir;
                  const char* dirname = "Pad waveforms";
                  
                  if (!dir->cd(dirname)) {
                     TDirectory* awdir = dir->mkdir(dirname);
                     awdir->cd();
                  }
                  
                  dir = dir->CurrentDirectory();
                  
                  std::string wname = std::string(xname) + "_waveform";
                  std::string wtitle = std::string(xtitle) + " current waveform";
                  
                  TH1D* hwf = (TH1D*)dir->FindObject(wname.c_str());
                  if (!hwf) {
                     hwf = new TH1D(wname.c_str(), wtitle.c_str(), nbins, 0, nbins);
                  }
                  
                  for (int i=0; i<nbins; i++) {
                     hwf->SetBinContent(i+1, aptr[i]);
                  }
               }

               // create per-channel data

               if (seqchan >= fHC.size()) {
                  for (unsigned i=fHC.size(); i<=seqchan; i++)
                     fHC.push_back(NULL);
               }

               if (fHC[seqchan] == NULL) {
                  fHC[seqchan] = new ChanHistograms(xname, xtitle, hdir_pads, nbins);
               }

               // check for spikes

               bool spike = false;

               double spike_max = 0;
               int spike_num = 0;
               for (int i=1; i<nbins-1; i++) {
                  double a0 = aptr[i-1];
                  double a1 = aptr[i];
                  double a2 = aptr[i+1];
                  if (a0 <= a1 && a1 <= a2)
                     continue;
                  if (a0 >= a1 && a1 >= a2)
                     continue;
                  double aa = (a0+a2)/2.0;
                  double da = fabs(a1 - aa);
                  if (da > spike_max)
                     spike_max = da;
                  if (da > 300)
                     spike_num++;
                  h_spike_diff->Fill(da);
               }
               h_spike_diff_max->Fill(spike_max);
               h_spike_num->Fill(spike_num);

               if (spike_max > 500 && spike_num > 10) {
                  spike = true;
               }

               if (spike) {
                  spike = true;
                  if (fHC[seqchan]->SaveBad(nbins, aptr)) {
                     if (verbose)
                        printf("BBB feam pos %d, seqsca %d, spike %f %d\n", ifeam, seqsca, spike_max, spike_num);

                     for (int i=1; i<nbins-1; i++) {
                        double a0 = aptr[i-1];
                        double a1 = aptr[i];
                        double a2 = aptr[i+1];
                        if (a0 <= a1 && a1 <= a2)
                           continue;
                        if (a0 >= a1 && a1 >= a2)
                           continue;
                        double aa = (a0+a2)/2.0;
                        double da = fabs(a1 - aa);
                        if (da > spike_max)
                           spike_max = da;
                        if (da > 300) {
                           if (verbose)
                              printf("bin %d, %.0f %.0f %.0f, aa %.0f, da %.0f\n", i, a0, a1, a2, aa, da);
                        }
                     }
                  }
               }

               // compute baseline

               double sum0 = 0;
               double sum1 = 0;
               double sum2 = 0;

               double bmin = aptr[ibaseline_start]; // baseline minimum
               double bmax = aptr[ibaseline_start]; // baseline maximum

               for (int i=ibaseline_start; i<ibaseline_end; i++) {
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

               // scan the whole waveform

               double wmin = aptr[0]; // waveform minimum
               double wmax = aptr[0]; // waveform maximum

               for (int i=0; i<nbins; i++) {
                  double a = aptr[i];
                  if (a < wmin)
                     wmin = a;
                  if (a > wmax)
                     wmax = a;
               }

               // scan the drift time region of the waveform

               double dmin = aptr[idrift_start]; // waveform minimum
               double dmax = aptr[idrift_start]; // waveform maximum

               for (int i=idrift_start; i<idrift_end; i++) {
                  double a = aptr[i];
                  if (a < dmin)
                     dmin = a;
                  if (a > dmax)
                     dmax = a;
               }

               // diagnostics

               if (scachan_is_fpn) {
                  if (brms > ADC_RMS_FPN_MIN && brms < ADC_RMS_FPN_MAX) {
                  } else {
                     printf("XXX bad fpn, feam %d, sca %d, readout %d, scachan %d, col %d, row %d, bmin %f, bmax %f, in hex 0x%04x, brms %f\n", ifeam, isca, ichan, scachan, col, row, bmin, bmax, (uint16_t)bmin, brms);
                     fpn_is_ok = false;
                  }
               }

               // diagnostics

#if 0
               if (scachan_is_fpn) {
                  printf("XXX fpn, feam %d, sca %d, readout %d, scachan %d, col %d, row %d, bmin %f, bmax %f, in hex 0x%04x, brms %f\n", ifeam, isca, ichan, scachan, col, row, bmin, bmax, (uint16_t)bmin, brms);
               }
#endif

#if 1
               if (scachan_is_pad || scachan_is_fpn) {
                  if (bmax-bmin == 0) {
                     if (first_zero_range) {
                        first_zero_range = false;
                        printf("XXX zero baseline range, feam %d, sca %d, readout %d, scachan %d, col %d, row %d, bmin %f, bmax %f, in hex 0x%04x\n", ifeam, isca, ichan, scachan, col, row, bmin, bmax, (uint16_t)bmin);
                     }
#if 0
                     //FeamAdcData* aaa = e->adcs[ifeam];
                     e->Print();
                     printf("\n");
                     e->modules[ifeam]->Print(1);
                     printf("\n");
                     abort();
#endif
                  }
               }
#endif

               // find pulses

               double wamp = bmean - wmin;

               if (wmin == ADC_MIN_ADC)
                  wamp = ADC_OVERFLOW;
               
               if (0) {
                  static double xwmin = 0;
                  if (wmin < xwmin) {
                     xwmin = wmin;
                     printf("MMM --- mean %f, wmin %f, wamp %f\n", bmean, wmin, wamp);
                  }
               }

               fHF[ifeam]->h_amp->Fill(wamp);

               //int wpos = find_pulse(aptr, nbins, bmean, -1.0, wamp/2.0);
               double wpos = find_pulse_time(aptr, nbins, bmean, -1.0, wamp/2.0);

               double wpos_ns = wpos*16.0;

               double damp = bmean - dmin;
               int dpos = idrift_start + find_pulse(aptr+idrift_start, idrift_end-idrift_start, bmean, -1.0, damp/2.0);

               // decide if we have a hit

               bool hit_time = false;
               bool hit_amp = false;
               bool hit = false;

               if ((wpos > iwire_start) && (wpos < idrift_end)) {
                  hit_time = true;
               }

               if (fPulser) {
                  if ((wpos > ipulser_start) && (wpos < ipulser_end)) {
                     hit_time = true;
                  }
               }

               if (wamp > hit_amp_threshold) {
                  hit_amp = true;
               }

               hit = hit_time && hit_amp;

               if (1 || hit_amp) {

                  FeamChannel* c = new FeamChannel;
                  //c->bank = e->modules[ifeam]->fBank;
                  c->imodule    = e->modules[ifeam]->fModule;
                  c->pwb_column = e->modules[ifeam]->fColumn;
                  c->pwb_ring   = e->modules[ifeam]->fRing;
                  c->sca = isca;
                  c->sca_readout = ichan;
                  c->sca_chan = scachan;
                  c->pad_col = col;
                  c->pad_row = row;
                  c->first_bin = 0;
                  c->adc_samples.reserve(nbins);
                  for (int i=0; i<nbins; i++)
                     c->adc_samples.push_back(aptr[i]);

                  //printf("adc samples: size %d, capacity %d\n", c->adc_samples.size(), c->adc_samples.capacity());

                  e->hits.push_back(c);

                  nhitchan++;
                  nhitchan_feam++;

                  fHF[ifeam]->h_nhitchan_seqsca->Fill(seqsca);

                  if (hit) {
                     AgPadHit h;
                     h.pos = ifeam;
                     h.seqsca = seqsca;
                     h.col = col;
                     h.row = row;
                     h.time = wpos_ns;
                     h.amp  = wamp;
                     hits->fPadHits.push_back(h);
                  }
               }

               if (doPrint) {
                  printf("chan %3d: baseline %8.1f, rms %8.1f, min %8.1f, max %8.1f, amp %8.1f, wpos %5.1f, hit %d\n", ichan, bmean, brms, wmin, wmax, wamp, wpos, hit);
               }

               if (1 || (wpos > 0 && wpos < 4000 && wamp > 1000)) {
                  if (wamp > zmax) {
                     if (doPrint)
                        printf("plot this one.\n");
                     //iplot = ichan;
                     zmax = wamp;
                  }
               }

               // save first waveform

               if (fHC[seqchan]->hwaveform_first->GetEntries() == 0) {
                  if (doPrint)
                     printf("saving first waveform %d\n", seqchan);
                  for (int i=0; i<nbins; i++)
                     fHC[seqchan]->hwaveform_first->SetBinContent(i+1, aptr[i]);
               }

               // save biggest waveform

               if (wamp > fHC[seqchan]->fMaxWamp) {
                  fHC[seqchan]->fMaxWamp = wamp;
                  if (doPrint)
                     printf("saving biggest waveform %d\n", seqchan);
                  for (int i=0; i<nbins; i++)
                     fHC[seqchan]->hwaveform_max->SetBinContent(i+1, aptr[i]);
               }

               // add to average waveform

               for (int j=0; j< nbins; j++)
                  fHC[seqchan]->hwaveform_avg->AddBinContent(j+1, aptr[j]);
               fHC[seqchan]->nwf++;

               // save biggest drift region waveform

               if (dpos > idrift_cut){
                  if(damp > fHC[seqchan]->fMaxWampDrift) {
                     fHC[seqchan]->fMaxWampDrift = damp;
                     if (doPrint)
                        printf("saving biggest drift waveform %d\n", seqchan);
                     for (int i=0; i<nbins; i++)
                        fHC[seqchan]->hwaveform_max_drift->SetBinContent(i+1, aptr[i]);
                  }

                  // add to average waveform

                  for (int j=0; j< nbins; j++)
                     fHC[seqchan]->hwaveform_avg_drift->AddBinContent(j+1, aptr[j]);
                  fHC[seqchan]->nwf_drift++;

               }

               if (scachan_is_pad || scachan_is_fpn) {
                  hbmean_all->Fill(bmean);
                  hbrms_all->Fill(brms);
               }

               if (scachan_is_pad) {
                  hamp_pad->Fill(wamp);
                  hamp_pad_pedestal->Fill(wamp);
                  hamp_pad_above_pedestal->Fill(wamp);
                  if (hit_amp) {
                     hled_pad_amp->Fill(wpos);
                     hled_pad_amp_ns->Fill(wpos_ns);
                  }
               }

               if (scachan_is_pad) {
                  hbrms_all_pads->Fill(brms);
               } else if (scachan_is_fpn) {
                  hbrms_all_fpn->Fill(brms);
               }

               if (scachan_is_pad || scachan_is_fpn) {
                  h_adc_range_all->Fill(wmax-wmin);
                  h_adc_range_baseline->Fill(bmax-bmin);
                  h_adc_range_drift->Fill(dmax-dmin);

                  fHF[ifeam]->hbmean_prof->Fill(seqsca, bmean);
                  fHF[ifeam]->hbrms_prof->Fill(seqsca, brms);
                  fHF[ifeam]->hbrange_prof->Fill(seqsca, bmax-bmin);

                  if (scachan_is_pad) {
                     fHF[ifeam]->hbrms_pads->Fill(brms);
                  }

                  if (scachan_is_fpn) {
                     fHF[ifeam]->hbrms_fpn->Fill(brms);
                  }

                  h2led2amp->Fill(wpos, wamp);
               }

               // plots for hits

               if (hit_amp) {
                  hled_all_hits->Fill(wpos);
               }

               if (hit_time) {
                  hamp_all_hits->Fill(wamp);
               }

               // plots for the drift region

               if (dpos > idrift_cut) {
                  hdrift_amp_all->Fill(damp);
                  hdrift_amp_all_pedestal->Fill(damp);
                  hdrift_amp_all_above_pedestal->Fill(damp);

                  if (damp > hit_amp_threshold) {
                     hdrift_led_all->Fill(dpos);
                     hdrift_led2amp->Fill(dpos, damp);
                     if (seqpad >= 0) {
                        fHF[ifeam]->hnhits_pad_drift->Fill(seqpad);
                     }
                  }
               }

               if (hit) {
                  hnhits->Fill(seqchan);
                  hled_hit->Fill(wpos);
                  hamp_hit->Fill(wamp);

                  if (fPulser) {
                     h_pulser_led_hit->Fill(wpos);
                     fHF[ifeam]->h_pulser_hit_amp_seqpad->Fill(-1, 0); // force plot to start from 0
                     fHF[ifeam]->h_pulser_hit_amp_seqpad->Fill(seqpad, wamp);
                     fHF[ifeam]->h_pulser_hit_time_seqpad->Fill(seqpad, wpos);
                     fHF[ifeam]->h_pulser_hit_time_seqsca->Fill(seqsca, wpos);
                     fHF[ifeam]->h_pulser_hit_amp->Fill(wamp);
                     fHF[ifeam]->h_pulser_hit_time->Fill(wpos);
                     fHF[ifeam]->h_pulser_hit_time_zoom->Fill(wpos);
                     if (seqsca == 4)
                        fHF[ifeam]->h_pulser_hit_time_seqsca4_zoom->Fill(wpos);
                  }

                  fHF[ifeam]->h_nhits_seqsca->Fill(seqsca);
                  fHF[ifeam]->h_hit_time_seqsca->Fill(seqsca, wpos);
                  fHF[ifeam]->h_hit_amp_seqsca->Fill(seqsca, wamp);
                  fHF[ifeam]->h_hit_amp_seqpad->Fill(seqpad, wamp);

                  if (wamp >= 10000 && wamp <= 40000) {
                     fHF[ifeam]->h_amp_seqsca->Fill(seqsca, wamp);
                     fHF[ifeam]->h_amp_seqpad->Fill(seqpad, wamp);
                  }

                  //h_amp_hit_col->Fill((ifeam*4 + col)%(MAX_FEAM_PAD_COL*MAX_FEAM), wamp);

                  if (seqpad >= 0) {
                     fHF[ifeam]->h_nhits_seqpad->Fill(seqpad);
                     if (!spike) {
                        fHF[ifeam]->hnhits_pad_nospike->Fill(seqpad);
                     }
                  }
               }

               if (spike) {
                  fHF[ifeam]->h_spike_seqsca->Fill(seqsca);
               }
            }
         }

         if (fpn_is_ok) {
            fCountGoodFpn ++;
            printf("XXX good fpn count %d\n", fCountGoodFpn);
         } else {
            fCountBadFpn ++;
            printf("XXX bad fpn count %d\n", fCountBadFpn);
         }

         fHF[ifeam]->h_spike_seqsca->Fill(1); // event counter marker
         fHF[ifeam]->hnhitchan->Fill(nhitchan_feam);
      }

      hnhitchan->Fill(nhitchan);

#if 0
      bool do_plot = (runinfo->fRoot->fgApp != NULL);

      if (do_plot) {
         // plot waveforms

         fC->Clear();
         fC->Divide(2,3);

         // if (1) {
         //    fC->cd(1);
         //    TH1D* hh = new TH1D("hh", "hh", nbins, 0, nbins);
         //    for (int ibin=0; ibin<nbins; ibin++) {
         //       hh->SetBinContent(ibin+1, e->adcs[0]->adc[0][0][ibin]);
         //    }
         //    hh->Draw();
         // }

         // if (1) {
         //    fC->cd(2);
         //    TH1D* hhh = new TH1D("hhh", "hhh", nbins, 0, nbins);
         //    for (int ibin=0; ibin<nbins; ibin++) {
         //       hhh->SetBinContent(ibin+1, e->adcs[0]->adc[0][0][ibin]);
         //    }
         //    hhh->SetMinimum(-33000);
         //    hhh->SetMaximum(+33000);
         //    hhh->Draw();
         // }

#if 0
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
#endif

         fC->Modified();
         fC->Draw();
         fC->Update();
      }

      if (fFlags->fPlotPad >= 0) {
         if (!fFlags->fPlotPadCanvas)
            fFlags->fPlotPadCanvas = new TCanvas("FEAM PAD", "FEAM PAD", 900, 650);

         TCanvas*c = fFlags->fPlotPadCanvas;

         c->cd();

#if 0
         int nbins = ww[fFlags->fPlotPad]->nsamples;
         TH1D* h = new TH1D("h", "h", nbins, 0, nbins);
         for (int ibin=0; ibin<nbins; ibin++) {
            h->SetBinContent(ibin+1, ww[fFlags->fPlotPad]->samples[ibin]);
         }

         h->SetMinimum(-9000);
         h->SetMaximum(+9000);
         h->Draw();
#endif

         c->Modified();
         c->Draw();
         c->Update();
      }

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
         printf("AnalyzeSpecialEvent, run %d, event serno %d, id 0x%04x, data size %d\n", runinfo->fRunNo, event->serial_number, (int)event->event_id, event->data_size);
   }
};

class FeamModuleFactory: public TAFactory
{
public:
   FeamFlags fFlags;
   
public:
   void Init(const std::vector<std::string> &args)
   {
      printf("FeamModuleFactory::Init!\n");

      for (unsigned i=0; i<args.size(); i++) {
         if (args[i] == "--nopads")
            fFlags.fDoPads = false;
         if (args[i] == "--plot1")
            fFlags.fPlotPad = atoi(args[i+1].c_str());
         if (args[i] == "--wfexport")
            fFlags.fExportWaveforms = true;
      }
   }

   void Finish()
   {
      printf("FeamModuleFactory::Finish!\n");
   }

   TARunObject* NewRunObject(TARunInfo* runinfo)
   {
      printf("FeamModuleFactory::NewRunObject, run %d, file %s\n", runinfo->fRunNo, runinfo->fFileName.c_str());
      return new FeamModule(runinfo, &fFlags);
   }
};

static TARegister tar(new FeamModuleFactory);

/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */
