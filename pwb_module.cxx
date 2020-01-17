//
// pwb_module.cxx - create PWB diagnostic histograms
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

#include "AgFlow.h"

#include "wfsuppress.h"

#define DELETE(x) if (x) { delete (x); (x) = NULL; }

#define MEMZERO(p) memset((p), 0, sizeof(p))

#define ADC_MIN -2100
#define ADC_MAX  2100

#define NUM_PWB (8*8)

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
#define ADC_RMS_FPN_MAX 4.000

#define ADC_RMS_PAD_MIN 2.500
#define ADC_RMS_PAD_MAX 25.0

#define ADC_PULSER_TIME 450

#define NUM_SEQSCA (3*80+79)

#define NUM_TIME_BINS 512
#define MAX_TIME_BINS 512
#define MAX_TIME_NS 8200

class PwbHistograms
{
public:
   //TProfile* hbmean_prof  = NULL;
   //TProfile* hbrms_prof   = NULL;
   TProfile* hbmean_bis_prof  = NULL;
   TProfile* hbrms_bis_prof   = NULL;
   TH1D*     hbrms_pads   = NULL;
   TH1D*     hbrms_fpn    = NULL;
   //TProfile* hbrange_prof = NULL;
   TProfile* hbrange_bis_prof = NULL;
   TH1D*     h_fpn_shift[4] = { NULL, NULL, NULL, NULL };
   TH1D*     h_amp = NULL;
   TH1D* hnhitchan = NULL;
   TH1D* h_nhitchan_seqsca = NULL;
   TH1D* h_spike_seqsca = NULL;
   TH1D* h_nhits_seqsca = NULL;
   TH1D* h_nhits_seqpad = NULL;
   TH1D* hnhits_pad_nospike = NULL;
   //TH1D* hnhits_pad_drift = NULL;
   TH2D* h_hit_time_seqsca = NULL;
   TH2D* h_hit_amp_seqsca = NULL;
   TH2D* h_hit_amp_seqpad = NULL;
   //TProfile* h_amp_seqsca_prof = NULL;
   //TProfile* h_amp_seqpad_prof = NULL;
   TH1D* h_pulser_hit_amp = NULL;
   TH1D* h_pulser_hit_time = NULL;
   TH1D* h_pulser_hit_time_zoom = NULL;
   //TH1D* h_pulser_hit_time_seqsca4_zoom = NULL;
   TProfile* h_pulser_hit_amp_seqpad_prof = NULL;
   TProfile* h_pulser_hit_time_seqpad_prof = NULL;
   TProfile* h_pulser_hit_time_seqsca_prof = NULL;

public:
   PwbHistograms()
   {
   };

   void CreateHistograms(TDirectory* hdir, int imodule, int icolumn, int iring, bool pulser)
   {
      char xname[256];
      char xtitle[256];

      sprintf(xname,  "pwb%02d_c%dr%d", imodule, icolumn, iring);
      sprintf(xtitle, "pwb %02d, col %d, ring %d", imodule, icolumn, iring);

      hdir->cd();

      TDirectory* xdir = hdir->mkdir(xname);
      xdir->cd();

      char name[256];
      char title[256];

      //sprintf(name,  "%s_baseline_mean_prof", xname);
      //sprintf(title, "%s baseline mean vs (SCA*80 + readout index)", xtitle);
      //hbmean_prof = new TProfile(name, title, NUM_SEQSCA, 0.5, NUM_SEQSCA+0.5);

      //sprintf(name,  "%s_baseline_rms_prof", xname);
      //sprintf(title, "%s baseline rms vs (SCA*80 +  readout index)", xtitle);
      //hbrms_prof  = new TProfile(name, title,  NUM_SEQSCA, 0.5, NUM_SEQSCA+0.5);

      sprintf(name,  "%s_baseline_mean_bis_prof", xname);
      sprintf(title, "%s baseline mean vs seqsca; SCA*80 + readout index; mean, adc counts", xtitle);
      hbmean_bis_prof = new TProfile(name, title, NUM_SEQSCA, 0.5, NUM_SEQSCA+0.5);

      sprintf(name,  "%s_baseline_rms_bis_prof", xname);
      sprintf(title, "%s baseline rms vs seqsca; SCA*80 + readout index; rms, adc counts", xtitle);
      hbrms_bis_prof  = new TProfile(name, title,  NUM_SEQSCA, 0.5, NUM_SEQSCA+0.5);

      sprintf(name,  "%s_baseline_rms_pads", xname);
      sprintf(title, "%s baseline rms for pad channels; rms, adc counts", xtitle);
      hbrms_pads  = new TH1D(name, title,  100, 0, ADC_RANGE_RMS);

      sprintf(name,  "%s_baseline_rms_fpn", xname);
      sprintf(title, "%s baseline rms fpr fpn channels; rms, adc counts", xtitle);
      hbrms_fpn   = new TH1D(name, title,  100, 0, ADC_RANGE_RMS);

      //sprintf(name,  "%s_baseline_range_prof", xname);
      //sprintf(title, "%s baseline range (max-min) vs (SCA*80 +  readout index)", xtitle);
      //hbrange_prof  = new TProfile(name, title,  NUM_SEQSCA, 0.5, NUM_SEQSCA+0.5);

      sprintf(name,  "%s_baseline_range_bis_prof", xname);
      sprintf(title, "%s baseline range (max-min) vs seqsca; SCA*80 + readout index; max-min, adc counts", xtitle);
      hbrange_bis_prof  = new TProfile(name, title,  NUM_SEQSCA, 0.5, NUM_SEQSCA+0.5);

      sprintf(name,  "%s_fpn_shift_sca0", xname);
      sprintf(title, "%s fpn shift, sca 0; sca time bins", xtitle);
      h_fpn_shift[0]  = new TH1D(name, title,  41, -20, 20);

      sprintf(name,  "%s_fpn_shift_sca1", xname);
      sprintf(title, "%s fpn shift, sca 1; sca time bins", xtitle);
      h_fpn_shift[1]  = new TH1D(name, title,  41, -20, 20);

      sprintf(name,  "%s_fpn_shift_sca2", xname);
      sprintf(title, "%s fpn shift, sca 2; sca time bins", xtitle);
      h_fpn_shift[2]  = new TH1D(name, title,  41, -20, 20);

      sprintf(name,  "%s_fpn_shift_sca3", xname);
      sprintf(title, "%s fpn shift, sca 3; sca time bins", xtitle);
      h_fpn_shift[3]  = new TH1D(name, title,  41, -20, 20);

      sprintf(name,  "%s_amp", xname);
      sprintf(title, "%s waveform amplitude from baseline to minimum; adc counts", xtitle);
      h_amp = new TH1D(name, title,  ADC_BINS, 0, ADC_RANGE);

      sprintf(name,  "%s_nhitchan", xname);
      sprintf(title, "%s number of hit channels; number of hits", xtitle);
      hnhitchan = new TH1D(name, title, 100, 0, 100);

      sprintf(name,  "%s_nhitchan_map", xname);
      sprintf(title, "%s hit channels vs seqsca; SCA*80 + readout index; number of hits", xtitle);
      h_nhitchan_seqsca = new TH1D(name, title, NUM_SEQSCA, 0.5, NUM_SEQSCA+0.5);

      sprintf(name,  "%s_hit_map", xname);
      sprintf(title, "%s hits vs seqsca; SCA*80 + readout index", xtitle);
      h_nhits_seqsca = new TH1D(name, title, NUM_SEQSCA, 0.5, NUM_SEQSCA+0.5);

      sprintf(name,  "%s_spike_seqsca", xname);
      sprintf(title, "%s spikes vs seqsca; SCA*80 + readout index", xtitle);
      h_spike_seqsca = new TH1D(name, title, NUM_SEQSCA, 0.5, NUM_SEQSCA+0.5);

      sprintf(name,  "%s_hit_map_pads", xname);
      sprintf(title, "%s hits vs seqpad; col*4*18+row", xtitle);
      h_nhits_seqpad  = new TH1D(name, title, MAX_FEAM_PAD_ROWS*MAX_FEAM_PAD_COL, -0.5, MAX_FEAM_PAD_ROWS*MAX_FEAM_PAD_COL-0.5);

      sprintf(name,  "%s_hit_map_pads_nospike", xname);
      sprintf(title, "%s hits with spikes removed vs seqpad; col*4*18+row", xtitle);
      hnhits_pad_nospike = new TH1D(name, title, MAX_FEAM_PAD_ROWS*MAX_FEAM_PAD_COL, -0.5, MAX_FEAM_PAD_ROWS*MAX_FEAM_PAD_COL-0.5);

      //sprintf(name,  "%s_hit_map_pads_drift", xname);
      //sprintf(title, "%s hits in drift region vs TPC seq.pad (col*4*18+row)", xtitle);
      //hnhits_pad_drift = new TH1D(name, title, MAX_FEAM_PAD_ROWS*MAX_FEAM_PAD_COL, -0.5, MAX_FEAM_PAD_ROWS*MAX_FEAM_PAD_COL-0.5);

      sprintf(name,  "%s_hit_time_seqsca", xname);
      sprintf(title, "%s hit time vs seqsca; SCA*80 + readout index; hit time, sca time bins", xtitle);
      h_hit_time_seqsca = new TH2D(name, title, NUM_SEQSCA, 0.5, NUM_SEQSCA+0.5, 50, 0, 500);

      sprintf(name,  "%s_hit_amp_seqsca", xname);
      sprintf(title, "%s hit p.h. vs seqsca; SCA*80 + readout index; hit amp, adc counts", xtitle);
      h_hit_amp_seqsca = new TH2D(name, title, NUM_SEQSCA, 0.5, NUM_SEQSCA+0.5, 50, 0, ADC_RANGE);

      sprintf(name,  "%s_hit_amp_seqpad", xname);
      sprintf(title, "%s hit p.h. vs seqpad; col*4*18+row; hit amp, adc counts", xtitle);
      h_hit_amp_seqpad = new TH2D(name, title, MAX_FEAM_PAD_ROWS*MAX_FEAM_PAD_COL, -0.5, MAX_FEAM_PAD_ROWS*MAX_FEAM_PAD_COL-0.5, 50, 0, ADC_RANGE);

      //sprintf(name,  "%s_amp_seqsca_prof", xname);
      //sprintf(title, "%s hit p.h. profile, cut 10000..40000 vs (SCA*80 + readout index)", xtitle);
      //h_amp_seqsca_prof = new TProfile(name, title, NUM_SEQSCA, -0.5, NUM_SEQSCA-0.5);

      //sprintf(name,  "%s_amp_seqpad_prof", xname);
      //sprintf(title, "%s hit p.h. profile, cut 10000..40000 vs TPC seq.pad (col*4*18+row)", xtitle);
      //h_amp_seqpad_prof = new TProfile(name, title, MAX_FEAM_PAD_ROWS*MAX_FEAM_PAD_COL+1, -1.5, MAX_FEAM_PAD_ROWS*MAX_FEAM_PAD_COL-0.5);

      if (pulser) {
         sprintf(name,  "%s_pulser_hit_amp", xname);
         sprintf(title, "%s pulser hit p.h.", xtitle);
         h_pulser_hit_amp = new TH1D(name, title, ADC_BINS_PULSER, 0, ADC_RANGE_PULSER);

         sprintf(name,  "%s_pulser_hit_time", xname);
         sprintf(title, "%s pulser hit time", xtitle);
         h_pulser_hit_time = new TH1D(name, title, NUM_TIME_BINS, 0, NUM_TIME_BINS);

         sprintf(name,  "%s_pulser_hit_time_zoom", xname);
         sprintf(title, "%s pulser hit time zoom", xtitle);
         h_pulser_hit_time_zoom = new TH1D(name, title, 100, ADC_PULSER_TIME-10, ADC_PULSER_TIME+10);

         //sprintf(name,  "%s_pulser_hit_time_seqsca4_zoom", xname);
         //sprintf(title, "%s pulser hit time seqsca 4 zoom", xtitle);
         //h_pulser_hit_time_seqsca4_zoom = new TH1D(name, title, 100, ADC_PULSER_TIME-10, ADC_PULSER_TIME+10);

         sprintf(name,  "%s_pulser_hit_amp_seqpad_prof", xname);
         sprintf(title, "%s pulser hit p.h. vs seqpad (col*4*18+row)", xtitle);
         h_pulser_hit_amp_seqpad_prof = new TProfile(name, title, MAX_FEAM_PAD_ROWS*MAX_FEAM_PAD_COL+1, -1.5, MAX_FEAM_PAD_ROWS*MAX_FEAM_PAD_COL-0.5);

         sprintf(name,  "%s_pulser_hit_time_seqpad_prof", xname);
         sprintf(title, "%s pulser hit time vs seqpad; col*4*18+row; hit time, sca time bins", xtitle);
         h_pulser_hit_time_seqpad_prof = new TProfile(name, title, MAX_FEAM_PAD_ROWS*MAX_FEAM_PAD_COL, -0.5, MAX_FEAM_PAD_ROWS*MAX_FEAM_PAD_COL-0.5);

         sprintf(name,  "%s_pulser_hit_time_seqsca_prof", xname);
         sprintf(title, "%s pulser hit time vs seqsca; SCA*80 + readout index; hit time, sca time bins", xtitle);
         h_pulser_hit_time_seqsca_prof = new TProfile(name, title, NUM_SEQSCA, -0.5, NUM_SEQSCA-0.5);
      }
   }
};

//static int find_pulse(const std::vector<int>& adc, int istart, int iend, double baseline, double gain, double threshold)
//{
//   for (int i=istart; i<iend; i++) {
//      if ((adc[i]-baseline)*gain > threshold) {
//         return i;
//      }
//   }
//
//   return 0;
//}

static double find_pulse_time(const std::vector<int>& adc, int nbins, double baseline, double gain, double threshold)
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

   bool SaveBad(int nbins, const std::vector<int> &adc)
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

class PwbFlags
{
public:
   int  fPlotPad = -1;
   TCanvas* fPlotPadCanvas = NULL;
   bool fWfSuppress = false;

public:
   PwbFlags() // ctor
   {
   }

   ~PwbFlags() // dtor
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

static bool fpn_rms_ok(int ichan, double brms)
{
   if (ichan < 4)
      return true;

   if (brms > ADC_RMS_FPN_MIN && brms < ADC_RMS_FPN_MAX)
      return true;

   return false;
}

static int fpn_wrap(int ifpn)
{
   while (ifpn < 0)
      ifpn += 80;

   while (ifpn >= 80)
      ifpn -= 80;

   return ifpn;
}

class PwbModule: public TARunObject
{
public:
   PwbFlags* fFlags = NULL;

   std::vector<std::vector<std::vector<WfSuppress*>>> fWfSuppress;
   std::vector<TH1D*> fWfSuppressAdcMin;
   TH1D* fWfSuppressAdcAmp = NULL;
   TH1D* fWfSuppressAdcAmpPos = NULL;
   TH1D* fWfSuppressAdcAmpNeg = NULL;
   TH1D* fWfSuppressAdcAmpCumulA = NULL;
   TH1D* fWfSuppressAdcAmpCumulB = NULL;

   TH1D* h_all_fpn_count = NULL;

   //TProfile* h_all_fpn_rms_prof = NULL;
   TProfile* h_all_fpn_mean_bis_prof = NULL;
   TProfile* h_all_fpn_rms_bis_prof = NULL;

   std::vector<TProfile*> h_all_fpn_mean_per_col_prof;
   std::vector<TProfile*> h_all_fpn_rms_per_col_prof;

   TH1D* h_all_pad_baseline_count = NULL;
   TH1D* h_all_pad_baseline_good_count = NULL;
   TProfile* h_all_pad_baseline_mean_prof = NULL;
   TProfile* h_all_pad_baseline_rms_prof = NULL;

   //TH1D* h_spike_diff = NULL;
   //TH1D* h_spike_diff_max = NULL;
   //TH1D* h_spike_num = NULL;

   TH1D* hbmean_all = NULL;
   TH1D* hbrms_all = NULL;
   TH1D* hbrms_all_pads = NULL;
   TH1D* hbrms_all_fpn = NULL;

   TProfile* hbmean_pwb_prof = NULL;

   //TH1D* h_adc_range_all = NULL;
   //TH1D* h_adc_range_baseline = NULL;
   //TH1D* h_adc_range_drift = NULL;

   TH1D* hamp_pad;
   TH1D* hamp_pad_pedestal;
   TH1D* hamp_pad_above_pedestal;
   TH1D* hled_pad_amp;
   TH1D* hled_pad_amp_ns;

   TH1D* hled_all_hits;
   TH1D* hamp_all_hits;

   TH2D* h2led2amp;

   //TH1D* hdrift_amp_all;
   //TH1D* hdrift_amp_all_pedestal;
   //TH1D* hdrift_amp_all_above_pedestal;
   //TH1D* hdrift_led_all;
   //TH2D* hdrift_led2amp;

   //TH1D* hnhits;
   TH1D* hled_hit;
   TH1D* hamp_hit;
   //TH2D* h_amp_hit_col = NULL;

   bool  fPulser = true;
   TH1D* h_pulser_led_hit = NULL;

   TH1D* hnhitchan = NULL;

   //TH2D* hpadmap;

   TDirectory* hdir_summary = NULL;
   TDirectory* hdir_wfsuppress = NULL;
   TDirectory* hdir_pwb  = NULL;
   TDirectory* hdir_pads = NULL;
   std::vector<PwbHistograms*> fHF;
   std::vector<ChanHistograms*> fHC;

   int fCountTestScaEvents = 0;
   int fCountBadScaEvents = 0;
   int fCountBadSca = 0;
   int fCountGoodFpn = 0;
   int fCountBadFpn = 0;

   int fCountBadPad = 0;

   bool fTrace = false;

   PwbModule(TARunInfo* runinfo, PwbFlags* f)
      : TARunObject(runinfo)
   {
      if (fTrace)
         printf("PwbModule::ctor!\n");

      fFlags = f;
   }

   ~PwbModule()
   {
      if (fTrace)
         printf("PwbModule::dtor!\n");
      for (unsigned i=0; i<fHF.size(); i++) {
         DELETE(fHF[i]);
      }
      for (unsigned i=0; i<fHC.size(); i++) {
         DELETE(fHC[i]);
      }

      //for (unsigned i=0; i<h_all_fpn_mean_per_col.size(); i++) {
      //   DELETE(h_all_fpn_mean_per_col[i]);
      //}

      //for (unsigned i=0; i<h_all_fpn_rms_per_col.size(); i++) {
      //   DELETE(h_all_fpn_rms_per_col[i]);
      //}
   }

   void BeginRun(TARunInfo* runinfo)
   {
      if (fTrace)
         printf("BeginRun, run %d, file %s\n", runinfo->fRunNo, runinfo->fFileName.c_str());
      //time_t run_start_time = runinfo->fOdb->odbReadUint32("/Runinfo/Start time binary", 0, 0);
      //printf("ODB Run start time: %d: %s", (int)run_start_time, ctime(&run_start_time));
      runinfo->fRoot->fOutputFile->cd(); // select correct ROOT directory
   }

   void CreateHistograms(TARunInfo* runinfo)
   {
      if (hbmean_all) // already created
         return;

      runinfo->fRoot->fOutputFile->cd();
      hdir_pads = gDirectory->mkdir("xpads");
      hdir_pads->cd(); // select correct ROOT directory

      hdir_summary = hdir_pads->mkdir("summary");
      hdir_wfsuppress = hdir_pads->mkdir("wfsuppress");
      hdir_pwb = hdir_pads->mkdir("pwb");

      hdir_summary->cd();

      //int max_fpn = NUM_PWB*MAX_FEAM_SCA*4;
      int max_pad = NUM_PWB*MAX_FEAM_SCA;

      h_all_fpn_count = new TH1D("all_fpn_count", "count of all fpn channels; sca+4*(ring+8*column); fpn count", max_pad, -0.5, max_pad-0.5);

      //h_all_fpn_rms_prof  = new TProfile("all_fpn_rms_prof", "rms of all fpn channels; fpn+4*(sca+4*(ring+8*column)); rms, adc counts", max_fpn, -0.5, max_fpn-0.5);

      h_all_fpn_mean_bis_prof  = new TProfile("all_fpn_mean_bis_prof", "mean of all fpn channels; sca+4*(ring+8*column)", max_pad, -0.5, max_pad-0.5);
      h_all_fpn_rms_bis_prof  = new TProfile("all_fpn_rms_bis_prof", "rms of all fpn channels; sca+4*(ring+8*column)", max_pad, -0.5, max_pad-0.5);

      for (unsigned i=0; i<8; i++) {
         char name[256];
         char title[256];

         int num = 8*MAX_FEAM_SCA;

         sprintf(name, "all_fpn_mean_col%d_prof", i);
         sprintf(title, "mean of fpn channels column %d; 4*ring+sca", i);
         h_all_fpn_mean_per_col_prof.push_back(new TProfile(name, title, num, -0.5, num-0.5));

         sprintf(name, "all_fpn_rms_col%d_prof", i);
         sprintf(title, "rms of fpn channels column %d; 4*ring+sca", i);
         h_all_fpn_rms_per_col_prof.push_back(new TProfile(name, title, num, -0.5, num-0.5));
      }

      h_all_pad_baseline_count = new TH1D("all_pad_baseline_count", "count of all baselines; sca+4*(ring+8*column); pad count", max_pad, -0.5, max_pad-0.5);

      h_all_pad_baseline_good_count = new TH1D("all_pad_baseline_good_count", "count of good baselines; sca+4*(ring+8*column); pad count", max_pad, -0.5, max_pad-0.5);

      h_all_pad_baseline_mean_prof = new TProfile("all_pad_baseline_mean_prof", "baseline mean of all pad channels; sca+4*(ring+8*column); mean, adc counts", max_pad, -0.5, max_pad-0.5);
      h_all_pad_baseline_rms_prof  = new TProfile("all_pad_baseline_rms_prof", "baseline rms of all pad channels; sca+4*(ring+8*column); rms, adc counts", max_pad, -0.5, max_pad-0.5);

      //h_spike_diff = new TH1D("spike_diff", "channel spike finder", 100, 0, 1000);
      //h_spike_diff_max = new TH1D("spike_diff_max", "channel spike finder, max", 100, 0, 1000);
      //h_spike_num = new TH1D("spike_num", "channel spike finder, num", 100, 0-0.5, 100-0.5);

      hbmean_all = new TH1D("all_baseline_mean", "baseline mean of all channels; mean, adc counts", 100, ADC_MIN, ADC_MAX);
      hbrms_all  = new TH1D("all_baseline_rms",  "baseline rms of all channels; rms, adc counts",  100, 0, ADC_RANGE_RMS);
      hbrms_all_pads  = new TH1D("all_baseline_rms_pads",  "baseline rms of pad channels; rms, adc counts",  100, 0, ADC_RANGE_RMS);
      hbrms_all_fpn   = new TH1D("all_baseline_rms_fpn",  "baseline rms of fpn channels; rms, adc counts",  100, 0, ADC_RANGE_RMS);

      hbmean_pwb_prof = new TProfile("baseline_mean_pwb_prof", "baseline mean of all PWBs; seqpwb (column*8+ring); mean, adc counts", NUM_PWB, -0.5, NUM_PWB-0.5, ADC_MIN, ADC_MAX);

      //h_adc_range_all      = new TH1D("adc_range_all",      "waveform range (max-min)",  100, 0, ADC_RANGE_PED);
      //h_adc_range_baseline = new TH1D("adc_range_baseline", "waveform range (max-min), baseline region",  100, 0, ADC_RANGE_PED);
      //h_adc_range_drift    = new TH1D("adc_range_drift",    "waveform range (max-min), drift region",  100, 0, ADC_RANGE_PED);

      hamp_pad         = new TH1D("hamp_pad",   "pad channels pulse height; adc counts", 100, 0, ADC_RANGE);
      hamp_pad_pedestal = new TH1D("hamp_pad_pedestal", "pad channels pulse height, zoom on pedestal area; adc counts", 100, 0, ADC_RANGE_PED);
      hamp_pad_above_pedestal = new TH1D("hamp_pad_above_pedestal", "pad channels pulse height, away from pedestal area; adc counts", 100, ADC_RANGE_PED, ADC_RANGE);

      hled_pad_amp = new TH1D("hled_pad_amp",   "pad channels above threshold, pulse leading edge; adc time bins", 100, 0, MAX_TIME_BINS);
      hled_pad_amp_ns = new TH1D("hled_pad_amp_ns",   "pad channels above threshold, pulse leading edge; time, ns", 100, 0, MAX_TIME_NS);

      h2led2amp  = new TH2D("h2led2amp", "pulse amp vs time, adc time bins", 100, 0, MAX_TIME_BINS, 100, 0, ADC_RANGE);

      hled_all_hits = new TH1D("hled_all_hits",   "pulse leading edge, adc time bins, with p.h. cut", 100, 0, NUM_TIME_BINS);
      hamp_all_hits = new TH1D("hamp_all_hits",   "pulse height, with time cut", 100, 0, ADC_RANGE);

      //hdrift_amp_all = new TH1D("drift_amp", "drift region pulse height", 100, 0, ADC_RANGE);
      //hdrift_amp_all_pedestal = new TH1D("drift_amp_pedestal", "drift region pulse height, zoom on pedestal area", 100, 0, ADC_RANGE_PED);
      //hdrift_amp_all_above_pedestal = new TH1D("drift_amp_above_pedestal", "drift region pulse height, away from pedestal area", 100, ADC_RANGE_PED, ADC_RANGE);
      //hdrift_led_all = new TH1D("drift_led", "drift region pulse leading edge, adc time bins, above pedestal", 100, 0, NUM_TIME_BINS);
      //hdrift_led2amp = new TH2D("drift_led2amp", "drift region pulse amp vs time, adc time bins, above pedestal", 100, 0, NUM_TIME_BINS, 100, 0, ADC_RANGE);

      //hnhits = new TH1D("hnhits", "hits per channel", nchan, -0.5, nchan-0.5);
      hled_hit = new TH1D("hled_hit", "hit time, adc time bins", 100, 0, NUM_TIME_BINS);
      hamp_hit = new TH1D("hamp_hit", "hit pulse height", 100, 0, ADC_RANGE);
      //h_amp_hit_col = new TH2D("hamp_hit_col", "hit pulse height vs column (ifeam*4+col)", 8*4, -0.5, 8*4-0.5, 100, 0, ADC_RANGE);

      hnhitchan = new TH1D("hnhitchan", "number of hit channels per event", 100, 0, 1000);

      if (fPulser) {
         h_pulser_led_hit = new TH1D("pulser_led_hit", "pulser time, adc bins", 30, 180-0.5, 210-0.5);
      }

      //hpadmap = new TH2D("hpadmap", "map from TPC pad number (col*4*18+row) to SCA readout channel (sca*80+chan)", 4*4*18, -0.5, 4*4*18-0.5, NUM_SEQSCA, 0.5, NUM_SEQSCA+0.5);
   }

   void EndRun(TARunInfo* runinfo)
   {
      if (fTrace)
         printf("PwbModule::EndRun, run %d\n", runinfo->fRunNo);
      //time_t run_stop_time = runinfo->fOdb->odbReadUint32("/Runinfo/Stop time binary", 0, 0);
      //printf("ODB Run stop time: %d: %s", (int)run_stop_time, ctime(&run_stop_time));
      for(auto *hc: fHC){
         if(hc){
            if(hc->nwf) hc->hwaveform_avg->Scale(1./double(hc->nwf));
            if(hc->nwf_drift) hc->hwaveform_avg_drift->Scale(1./double(hc->nwf_drift));
         }
      }

      printf("PwbModule::EndRun: test for bad SCA: total events %d, bad events %d, bad sca %d, bad fpn %d, good fpn %d, bad pad %d\n", fCountTestScaEvents, fCountBadScaEvents, fCountBadSca, fCountBadFpn, fCountGoodFpn, fCountBadPad);
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

   // check for gibberish ADC data

#if 0
   bool TestBadSca(const FeamEvent* e)
   {
      bool bad_sca = false;

      fCountTestScaEvents++;

      for (unsigned ifeam=0; ifeam<e->adcs.size(); ifeam++) {
         FeamAdcData* aaa = e->adcs[ifeam];
         if (!aaa)
            continue;

         int imodule    = e->modules[ifeam]->fModule;

         for (int isca=0; isca<aaa->nsca; isca++) {
            int first_chan = 4;
            int v = aaa->adc[isca][first_chan][0];

            bool bad = (v == 0x7FFC);

            if (bad) {
               printf("Error: pwb%02d sca %d has gibberish data: ", imodule, isca);
               e->Print();
               printf("\n");
            }

            if (bad) {
               fCountBadSca++;
               bad_sca = true;
            }
         }
      }

      return bad_sca;
   }
#endif

#if 0
   const FeamChannel* FindChannel(const FeamEvent* e, int ipwb, int isca, int iri)
   {
      for (unsigned ii=0; ii<e->hits.size(); ii++) {
         FeamChannel* c = e->hits[ii];
         if (!c)
            continue;

         if (c->imodule == ipwb && c->sca == isca && c->sca_readout == iri) {
            return c;
         }
      }

      return NULL;
   }
#endif

   // check FPN channels and shifted channels

   void CheckAndShiftFpn(FeamEvent* e)
   {
      int ibaseline_start = 10;
      int ibaseline_end = 100;

      struct PwbPtr {
         int counter = 0;
         FeamChannel* ptr[4][80];
         PwbPtr() { for (int i=0; i<80; i++) { ptr[0][i]=NULL;ptr[1][i]=NULL;ptr[2][i]=NULL;ptr[3][i]=NULL; } };
      };

      std::vector<PwbPtr> pwb_ptr;

      for (unsigned ii=0; ii<e->hits.size(); ii++) {
         FeamChannel* c = e->hits[ii];
         if (!c)
            continue;

         unsigned imodule    = c->imodule;
         unsigned isca       = c->sca;
         unsigned iri        = c->sca_readout;

         assert(isca>=0 && isca<4);
         assert(iri>=0 && iri<80);

         while (imodule >= pwb_ptr.size()) {
            pwb_ptr.push_back(PwbPtr());
         }

         //printf("module %d sca %d ri %d, counter %d\n", imodule, isca, iri, pwb_ptr[imodule].counter);
         pwb_ptr[imodule].counter++;
         pwb_ptr[imodule].ptr[isca][iri] = c;
      }

      for (unsigned imodule = 0; imodule < pwb_ptr.size(); imodule++) {
         //printf("module %d, counter %d\n", imodule, pwb_ptr[imodule].counter);
         if (pwb_ptr[imodule].counter < 1) {
            continue;
         }

         PwbHistograms* hf = fHF[imodule];

         for (int isca = 0; isca < 4; isca++) {
            bool trace = false;
            bool trace_shift = false;

            int fpn_shift = 0;

            const FeamChannel *c16 = pwb_ptr[imodule].ptr[isca][16];
            const FeamChannel *c29 = pwb_ptr[imodule].ptr[isca][29];
            const FeamChannel *c54 = pwb_ptr[imodule].ptr[isca][54];
            const FeamChannel *c67 = pwb_ptr[imodule].ptr[isca][67];

            //printf("module %d sca %d, ptr %p %p %p %p\n", imodule, isca, c16, c29, c54, c67);

            if (!c16)
               continue;
            if (!c29)
               continue;
            if (!c54)
               continue;
            if (!c67)
               continue;

            //int pwb_column = c16->pwb_column;
            //int pwb_ring   = c16->pwb_ring;
            //int seqpwb = 0;
            //if (pwb_column >= 0) {
            //   seqpwb = pwb_column*8 + pwb_ring;
            //}
            //int seqpwbsca = seqpwb*4 + isca;
            //int seqpwbscafpn = seqpwbsca*4;

            double rms_fpn1 = compute_rms(c16->adc_samples.data(), ibaseline_start, ibaseline_end);
            double rms_fpn2 = compute_rms(c29->adc_samples.data(), ibaseline_start, ibaseline_end);
            double rms_fpn3 = compute_rms(c54->adc_samples.data(), ibaseline_start, ibaseline_end);
            double rms_fpn4 = compute_rms(c67->adc_samples.data(), ibaseline_start, ibaseline_end);

            //h_all_fpn_count->Fill(seqpwbscafpn + 0, 1);
            //h_all_fpn_count->Fill(seqpwbscafpn + 1, 1);
            //h_all_fpn_count->Fill(seqpwbscafpn + 2, 1);
            //h_all_fpn_count->Fill(seqpwbscafpn + 3, 1);

            //h_all_fpn_rms_prof->Fill(seqpwbscafpn + 0, rms_fpn1);
            //h_all_fpn_rms_prof->Fill(seqpwbscafpn + 1, rms_fpn2);
            //h_all_fpn_rms_prof->Fill(seqpwbscafpn + 2, rms_fpn3);
            //h_all_fpn_rms_prof->Fill(seqpwbscafpn + 3, rms_fpn4);

            if (fpn_rms_ok(16, rms_fpn1)
                && fpn_rms_ok(29, rms_fpn2) 
                && fpn_rms_ok(54, rms_fpn3) 
                && fpn_rms_ok(67, rms_fpn4)) {

               if (trace) {
                  printf("CheckAndShiftFpn: good fpn pwb%02d, sca %d, fpn rms: %5.1f %5.1f %5.1f %5.1f\n", imodule, isca, rms_fpn1, rms_fpn2, rms_fpn3, rms_fpn4);
               }

               fpn_shift = 0;
            } else {
               if (trace) {
                  printf("CheckAndShiftFpn: bad  fpn pwb%02d, sca %d, fpn rms: %5.1f %5.1f %5.1f %5.1f\n", imodule, isca, rms_fpn1, rms_fpn2, rms_fpn3, rms_fpn4);
               }

               for (int i=0; i>-30; i--) {
                  int ifpn1 = fpn_wrap(i+16);
                  int ifpn2 = fpn_wrap(i+29);
                  int ifpn3 = fpn_wrap(i+54);
                  int ifpn4 = fpn_wrap(i+67);

                  const FeamChannel *c16 = pwb_ptr[imodule].ptr[isca][ifpn1];
                  const FeamChannel *c29 = pwb_ptr[imodule].ptr[isca][ifpn2];
                  const FeamChannel *c54 = pwb_ptr[imodule].ptr[isca][ifpn3];
                  const FeamChannel *c67 = pwb_ptr[imodule].ptr[isca][ifpn4];
                  
                  if (!c16)
                     continue;
                  if (!c29)
                     continue;
                  if (!c54)
                     continue;
                  if (!c67)
                     continue;

                  double rms_fpn1 = compute_rms(c16->adc_samples.data(), ibaseline_start, ibaseline_end);
                  double rms_fpn2 = compute_rms(c29->adc_samples.data(), ibaseline_start, ibaseline_end);
                  double rms_fpn3 = compute_rms(c54->adc_samples.data(), ibaseline_start, ibaseline_end);
                  double rms_fpn4 = compute_rms(c67->adc_samples.data(), ibaseline_start, ibaseline_end);

                  if (trace_shift) {
                     printf("CheckAndShiftFpn: shift %3d fpn pwb%02d, sca %d, fpn rms: %5.1f %5.1f %5.1f %5.1f, fpn chan %2d %2d %2d %2d", i, imodule, isca, rms_fpn1, rms_fpn2, rms_fpn3, rms_fpn4, ifpn1, ifpn2, ifpn3, ifpn4);
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
            
            if (1||trace_shift) {
               if (fpn_shift != 0) {
                  printf("CheckAndShiftFpn: pwb%02d, sca %d, fpn_shift %d\n", imodule, isca, fpn_shift);
               }
            }

            hf->h_fpn_shift[isca]->Fill(fpn_shift);

            if (fpn_shift < 0) {
               int sca_readout[80];
               int sca_chan[80];
               int pad_col[80];
               int pad_row[80];

               for (int i=0; i<80; i++) {
                  sca_readout[i] = -1;
                  sca_chan[i] = -1;
                  pad_col[i] = -1;
                  pad_row[i] = -1;
               }

               for (int rr=0; rr<80; rr++) {
                  FeamChannel*c = pwb_ptr[imodule].ptr[isca][rr];
                  if (c) {
                     int ri = c->sca_readout;
                     assert(ri>=0 && ri<80);
                     sca_readout[ri] = c->sca_readout;
                     sca_chan[ri] = c->sca_chan;
                     pad_col[ri] = c->pad_col;
                     pad_row[ri] = c->pad_row;
                  }
               }

               for (int rr=0; rr<80; rr++) {
                  FeamChannel*c = pwb_ptr[imodule].ptr[isca][rr];
                  if (c) {
                     int ri = c->sca_readout - fpn_shift;
                     if (ri < 1)
                        ri += 79;
                     if (ri > 79)
                        ri -= 79;
                     c->sca_readout = ri;
                     if (sca_readout[ri] == -1) {
                        c->sca_chan = -1;
                        c->pad_col = -1;
                        c->pad_row = -1;
                     } else {
                        //if (ri != sca_readout[ri]) {
                        //printf("ri %d, %d\n", ri, sca_readout[ri]);
                        //}
                        assert(ri == sca_readout[ri]);
                        c->sca_chan = sca_chan[ri];
                        c->pad_col = pad_col[ri];
                        c->pad_row = pad_row[ri];
                     }
                  }
               }
            }
         }
      }
   }

   // Create per PWB histograms

   void CreatePwbHistograms(const FeamEvent* e)
   {
      for (unsigned i=0; i<e->hits.size(); i++) {
         if (!e->hits[i])
            continue;

         unsigned imodule = e->hits[i]->imodule;

         while (imodule >= fHF.size()) {
            fHF.push_back(NULL);
         }

         if (!fHF[imodule]) {
            fHF[imodule] = new PwbHistograms();
            fHF[imodule]->CreateHistograms(hdir_pwb, imodule, e->hits[i]->pwb_column, e->hits[i]->pwb_ring, fPulser);
         }
      }
   }

   TAFlowEvent* AnalyzeFlowEvent(TARunInfo* runinfo, TAFlags* flags, TAFlowEvent* flow)
   {
      //bool verbose = false;
      
      //printf("Analyze, run %d, event serno %d, id 0x%04x, data size %d\n", runinfo->fRunNo, event->serial_number, (int)event->event_id, event->data_size);

      AgEventFlow *ef = flow->Find<AgEventFlow>();

      if (!ef || !ef->fEvent)
         return flow;

      FeamEvent* e = ef->fEvent->feam;

      if (!e) {
         return flow;
      }

      if (1) {
         printf("Have PWB  event: ");
         e->Print();
         printf("\n");
      }

      if (e->error) {
         //delete e;
         return flow;
      }
      
      //

      bool doPrint = false;

      // create histograms

      CreateHistograms(runinfo);

      // create pad hits flow event

      AgPadHitsFlow* hits = new AgPadHitsFlow(flow);
      flow = hits;

      // check for bad sca data

#if 0
      bool bad_sca = TestBadSca(e);

      if (bad_sca) {
         fCountBadScaEvents++;
         e->error = true;
         return flow;
      }
#endif

      // loop over all waveforms

      //int iplot = 0;
      bool first_zero_range = true;

      int ibaseline_start = 10;
      int ibaseline_end = 100;

      //int iwire_start = 130;
      int iwire_end = 160;

      int idrift_start = iwire_end;
      //int idrift_cut = iwire_end;
      int idrift_end = 410;

      int ipulser_start = 400;
      int ipulser_end   = 500;

      double hit_amp_threshold = 100;

      int nhitchan = 0;

      CreatePwbHistograms(e);
      
      //printf("PrintFeamChannels!\n");
      //PrintFeamChannels(e->hits);

      bool need_shift = false;
      need_shift |= (runinfo->fRunNo >= 1865 && runinfo->fRunNo <= 9999);

      if (need_shift) {
         CheckAndShiftFpn(e);
      }

      for (unsigned ii=0; ii<e->hits.size(); ii++) {
         FeamChannel* c = e->hits[ii];
         if (!c)
            continue;

         int imodule    = c->imodule;
         int pwb_column = c->pwb_column;
         int pwb_ring   = c->pwb_ring;

         int seqpwb = 0;

         if (pwb_column >= 0) {
            seqpwb = pwb_column*8 + pwb_ring;
         }

         PwbHistograms* hf = fHF[imodule];

         int nhitchan_feam = 0;

         bool fpn_is_ok = true;
         bool pad_is_ok = true;

         int isca = c->sca;

         assert(isca >= 0);
         assert(isca < 4);

         int seqpwbsca = seqpwb*4+isca;
         int ichan = c->sca_readout;

         assert(ichan > 0);
         assert(ichan < 80);

         //unsigned seqchan = ifeam*(aaa->nsca*aaa->nchan) + isca*aaa->nchan + ichan;
         //unsigned seqchan = 0;
         int seqsca = isca*80 + ichan;

         int nbins = c->adc_samples.size();

         int scachan = c->sca_chan;
         int col = c->pad_col; // TPC pad column
         int row = c->pad_row; // TPC pad row
         int seqpad = -1; // TPC sequential pad number col*4*72+row

         bool scachan_is_pad = PwbPadMap::chan_is_pad(scachan);
         bool scachan_is_fpn = PwbPadMap::chan_is_fpn(scachan);
         bool scachan_is_reset = PwbPadMap::chan_is_reset(scachan);

         if (scachan_is_pad) {
            assert(col>=0 && col<4);
            assert(row>=0 && row<MAX_FEAM_PAD_ROWS);
            seqpad = col*MAX_FEAM_PAD_ROWS + row;
         } else {
            row = scachan; // special channel
         }
         
         char xname[256];
         char xtitle[256];

         if (scachan_is_pad) {
            sprintf(xname, "pwb%02d_%03d_sca%d_chan%02d_scachan%02d_col%02d_row%02d", imodule, seqsca, isca, ichan, scachan, col, row);
            sprintf(xtitle, "pwb%02d, sca %d, readout chan %d, sca chan %d, col %d, row %d", imodule, isca, ichan, scachan, col, row);
         } else if (scachan_is_fpn) {
            sprintf(xname, "pwb%02d_%03d_sca%d_chan%02d_fpn%d", imodule, seqsca, isca, ichan, -scachan);
            sprintf(xtitle, "pwb%02d, sca %d, readout chan %d, fpn %d", imodule, isca, ichan, -scachan);
         } else if (scachan_is_reset) {
            sprintf(xname, "pwb%02d_%03d_sca%d_chan%02d_reset%d", imodule, seqsca, isca, ichan, -scachan-4);
            sprintf(xtitle, "pwb%02d, sca %d, readout chan %d, reset %d", imodule, isca, ichan, -scachan-4);
         } else {
            sprintf(xname, "pwb%02d_%03d_sca%d_chan%02d", imodule, seqsca, isca, ichan);
            sprintf(xtitle, "pwb%02d, sca %d, readout chan %d", imodule, isca, ichan);
         }
         
         // create per-channel data
         
         //if (seqchan >= fHC.size()) {
         //   for (unsigned i=fHC.size(); i<=seqchan; i++)
         //      fHC.push_back(NULL);
         //}
         
         //if (fHC[seqchan] == NULL) {
         //   fHC[seqchan] = new ChanHistograms(xname, xtitle, hdir_pads, nbins);
         //}
         
         // check for spikes
         
         bool spike = false;
         
         double spike_max = 0;
         int spike_num = 0;
         for (int i=1; i<nbins-1; i++) {
            double a0 = c->adc_samples[i-1];
            double a1 = c->adc_samples[i];
            double a2 = c->adc_samples[i+1];
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
            //h_spike_diff->Fill(da);
         }
         //h_spike_diff_max->Fill(spike_max);
         //h_spike_num->Fill(spike_num);
         
         if (spike_max > 500 && spike_num > 10) {
            spike = true;
         }

#if 0         
         if (spike) {
            spike = true;
            if (fHC[seqchan]->SaveBad(nbins, c->adc_samples)) {
               if (verbose)
                  printf("BBB pwb%02d, seqsca %d, spike %f %d\n", imodule, seqsca, spike_max, spike_num);
               
               for (int i=1; i<nbins-1; i++) {
                  double a0 = c->adc_samples[i-1];
                  double a1 = c->adc_samples[i];
                  double a2 = c->adc_samples[i+1];
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
#endif

         if (fFlags->fWfSuppress) { // compute data suppression
            hdir_wfsuppress->cd();

            if (fWfSuppressAdcAmp == NULL) {
               fWfSuppressAdcAmp = new TH1D("WfSuppress ADC amp", "WfSuppress ADC amp", 100, 0, 0xFFF + 2000);
            }

            if (fWfSuppressAdcAmpPos == NULL) {
               fWfSuppressAdcAmpPos = new TH1D("WfSuppress ADC amp pos", "WfSuppress ADC amp pos", 100, 0, 0xFFF + 200);
            }

            if (fWfSuppressAdcAmpNeg == NULL) {
               fWfSuppressAdcAmpNeg = new TH1D("WfSuppress ADC amp neg", "WfSuppress ADC amp neg", 100, -(0xFFF + 200), 0);
            }

            if (fWfSuppressAdcAmpCumulA == NULL) {
               fWfSuppressAdcAmpCumulA = new TH1D("WfSuppress cumul a", "WfSuppress cumul A", 100, 0, 0xFFF + 200);
            }

            if (fWfSuppressAdcAmpCumulB == NULL) {
               fWfSuppressAdcAmpCumulB = new TH1D("WfSuppress cumul b", "WfSuppress cumul B", 100, 0, 0xFFF + 200);
            }

            if (imodule >= (int) fWfSuppressAdcMin.size())
               fWfSuppressAdcMin.resize(imodule+1);

            if (fWfSuppressAdcMin[imodule] == NULL) {
               char name[256];
               sprintf(name,  "pwb%02d_adc_min", imodule);
               char title[256];
               sprintf(title, "pwb%02d adc min", imodule);
               fWfSuppressAdcMin[imodule] = new TH1D(name, title, 100, 0, 0xFFF);
            }
            
            if (imodule >= (int)fWfSuppress.size())
               fWfSuppress.resize(imodule+1);
            
            if (isca >= (int)fWfSuppress[imodule].size())
               fWfSuppress[imodule].resize(MAX_FEAM_SCA);
            
            if (ichan >= (int)fWfSuppress[imodule][isca].size())
               fWfSuppress[imodule][isca].resize(MAX_FEAM_READOUT);
            
            //printf("imodule %d, size %d\n", imodule, (int)fWfSuppress.size());
            //printf("isca %d, size %d\n", isca, (int)fWfSuppress[imodule].size());
            //printf("ichan %d, size %d\n", isca, (int)fWfSuppress[imodule][isca].size());
            
            WfSuppress *s = fWfSuppress[imodule][isca][ichan];
            if (!s) {
               s = new WfSuppress();
               fWfSuppress[imodule][isca][ichan] = s;
            }
            
            unsigned sfirst = 1;
            
            s->Init(c->adc_samples[sfirst], 10000);
            
            bool keep = false;
            for (unsigned i=sfirst+1; i<c->adc_samples.size(); i++) {
               bool k = s->Add(c->adc_samples[i]);
               uint16_t base = s->GetBase();
               int16_t amp = s->GetAmp();
               keep |= k;
               if (0) {
                  printf("pwb %02d, sca %d, chan %2d: bin %3d, adc %d, base %d, amp %4d, keep %d %d, state: ", imodule, isca, ichan, i, c->adc_samples[i], base, amp, k, keep);
                  s->Print();
                  printf("\n");
               }
            }
            
            //double xampmax = fabs(s->GetAmpMax());
            double xampmin = fabs(s->GetAmpMin());
            //double xamp = std::min(xampmax, xampmin);
            double xamp = xampmin;
            if (s->GetClipped())
               xamp = 0xFFF + 1;
            
            fWfSuppressAdcAmp->Fill(xamp);
            fWfSuppressAdcAmpPos->Fill(s->GetAmpMax());
            fWfSuppressAdcAmpNeg->Fill(s->GetAmpMin());

            fWfSuppressAdcAmpCumulA->Fill(xamp);
            for (int i=0; i<xamp; i++) {
               fWfSuppressAdcAmpCumulA->Fill(i);
            }

            fWfSuppressAdcAmpCumulB->Fill(xamp);
            for (int i=xamp+1; i<=0xFFF+10; i++) {
               fWfSuppressAdcAmpCumulB->Fill(i);
            }
            fWfSuppressAdcMin[imodule]->Fill(s->GetAdcMin());

            printf("pwb %02d, sca %d, chan %2d: wfsuppress: ", imodule, isca, ichan);
            s->Print();
            printf(", keep: %d, xamp %d\n", keep, (int)xamp);
         }

         //exit(1);

         // compute baseline
         
         double sum0 = 0;
         double sum1 = 0;
         double sum2 = 0;
         
         double bmin = c->adc_samples[ibaseline_start]; // baseline minimum
         double bmax = c->adc_samples[ibaseline_start]; // baseline maximum
         
         for (int i=ibaseline_start; i<ibaseline_end; i++) {
            double a = c->adc_samples[i];
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

#if 0
         if (bmean > 6000) {
            printf("bmean %f\n", bmean);
            printf("chan %3d: baseline %8.1f, rms %8.1f, min %8.1f, max %8.1f, sum0/1/2 %f/%f/%f\n", ichan, bmean, brms, bmin, bmax, sum0, sum1, sum2);
            abort();
         }
#endif
         
         // scan the whole waveform
         
         double wmin = c->adc_samples[0]; // waveform minimum
         double wmax = c->adc_samples[0]; // waveform maximum
         
         for (int i=0; i<nbins; i++) {
            double a = c->adc_samples[i];
            if (a < wmin)
               wmin = a;
            if (a > wmax)
               wmax = a;
         }
         
         // scan the drift time region of the waveform
         
         double dmin = c->adc_samples[idrift_start]; // waveform minimum
         double dmax = c->adc_samples[idrift_start]; // waveform maximum
         
         for (int i=idrift_start; i<idrift_end; i++) {
            double a = c->adc_samples[i];
            if (a < dmin)
               dmin = a;
            if (a > dmax)
               dmax = a;
         }
         
         // diagnostics
         
         if (scachan_is_fpn) {
            h_all_fpn_count->Fill(seqpwbsca, 1);
            h_all_fpn_mean_bis_prof->Fill(seqpwbsca, bmean);
            h_all_fpn_rms_bis_prof->Fill(seqpwbsca, brms);

            if (pwb_column >= 0 && pwb_column < (int)h_all_fpn_mean_per_col_prof.size()) {
               h_all_fpn_mean_per_col_prof[pwb_column]->Fill(pwb_ring*4+isca, bmean);
            }

            if (pwb_column >= 0 && pwb_column < (int)h_all_fpn_rms_per_col_prof.size()) {
               h_all_fpn_rms_per_col_prof[pwb_column]->Fill(pwb_ring*4+isca, brms);
            }

            //h_all_fpn_rms_bis->Fill(seqpwbsca, brms);

            if (brms > ADC_RMS_FPN_MIN && brms < ADC_RMS_FPN_MAX) {
            } else {
               printf("XXX bad fpn, pwb%02d, sca %d, readout %d, scachan %d, col %d, row %d, bmin %f, bmax %f, in hex 0x%04x, brms %f\n", imodule, isca, ichan, scachan, col, row, bmin, bmax, (uint16_t)bmin, brms);
               fpn_is_ok = false;
            }
         }

         if (scachan_is_pad) {
            if (brms > ADC_RMS_PAD_MIN) {
               //if (brms > 15.0) {
               //   printf("ZZZ bad pad, pwb%02d, sca %d, readout %d, scachan %d, col %d, row %d, bmin %f, bmax %f, in hex 0x%04x, brms %f\n", imodule, isca, ichan, scachan, col, row, bmin, bmax, (uint16_t)bmin, brms);
               //   if (imodule==10 && isca==3) {
               //      (*flags) |= TAFlag_DISPLAY;
               //   }
               //}
            } else {
               //printf("XXX bad pad, pwb%02d, sca %d, readout %d, scachan %d, col %d, row %d, bmin %f, bmax %f, in hex 0x%04x, brms %f\n", imodule, isca, ichan, scachan, col, row, bmin, bmax, (uint16_t)bmin, brms);
               pad_is_ok = false;
            }
         }
         
         // diagnostics
         
#if 0
         if (scachan_is_fpn) {
            printf("XXX fpn, pwb%02d, sca %d, readout %d, scachan %d, col %d, row %d, bmin %f, bmax %f, in hex 0x%04x, brms %f\n", imodule, isca, ichan, scachan, col, row, bmin, bmax, (uint16_t)bmin, brms);
         }
#endif

#if 1
         if (scachan_is_pad || scachan_is_fpn) {
            if (bmax-bmin == 0) {
               if (first_zero_range) {
                  first_zero_range = false;
                  printf("XXX zero baseline range, pwb%02d, sca %d, readout %d, scachan %d, col %d, row %d, bmin %f, bmax %f, in hex 0x%04x\n", imodule, isca, ichan, scachan, col, row, bmin, bmax, (uint16_t)bmin);
               }
            }
         }
#endif

         if (scachan_is_pad) {
            h_all_pad_baseline_count->Fill(seqpwbsca);
            h_all_pad_baseline_mean_prof->Fill(seqpwbsca, bmean);
            h_all_pad_baseline_rms_prof->Fill(seqpwbsca, brms);
         }
         
         // find pulses
         
         double wamp = bmean - wmin;
         
         if (wmin == ADC_MIN_ADC)
            wamp = ADC_OVERFLOW;
         else if (wmin < ADC_MIN_ADC+10)
            wamp = ADC_OVERFLOW;

         //if (wmin < ADC_MIN_ADC + 20) {
         //   printf("wmin %f\n", wmin);
         //   doPrint = true;
         //}
         
         if (0) {
            static double xwmin = 0;
            if (wmin < xwmin) {
               xwmin = wmin;
               printf("MMM --- mean %f, wmin %f, wamp %f\n", bmean, wmin, wamp);
            }
         }
         
         if (scachan_is_pad) {
            hf->h_amp->Fill(wamp);
         }
         
         //int wpos = find_pulse(c->adc_samples, nbins, bmean, -1.0, wamp/2.0);
         double wpos = find_pulse_time(c->adc_samples, nbins, bmean, -1.0, wamp/2.0);
         
         double wpos_offset_ns = 2350.0;
         
         double wpos_ns = wpos*16.0 - wpos_offset_ns + 1000.0;

         
         if (runinfo->fRunNo >= 2166) {
            wpos_ns += 200.0 + 150;
         } else if (runinfo->fRunNo >= 2028) {
            wpos_ns += 200.0;
         }
         
         //double damp = bmean - dmin;
         //int dpos = find_pulse(c->adc_samples, idrift_start, idrift_end, bmean, -1.0, damp/2.0);
         
         // decide if we have a hit
         
         bool hit_time = false;
         bool hit_amp = false;
         bool hit = false;
         
         if (scachan_is_pad && (wpos_ns > 800.0) && (wpos_ns < 5600.0)) {
            hit_time = true;
         }
         
         if (fPulser) {
            if ((wpos > ipulser_start) && (wpos < ipulser_end)) {
               hit_time = true;
            }
         }

         if (scachan_is_pad && (wamp > hit_amp_threshold)) {
            hit_amp = true;
         }
         
         hit = hit_time && hit_amp;
         
         if (hit_amp) {
            
            nhitchan++;
            nhitchan_feam++;
            
            hf->h_nhitchan_seqsca->Fill(seqsca);
            
            if (hit && col >= 0 && row >= 0) {
               assert(col >= 0 && col < MAX_FEAM_PAD_COL);
               assert(row >= 0 && row < MAX_FEAM_PAD_ROWS);
               
               AgPadHit h;
               h.imodule = imodule;
               h.seqsca = seqsca;
               if (c->pwb_column >= 0) {
                  h.tpc_col = c->pwb_column * MAX_FEAM_PAD_COL + col;
                  h.tpc_row = c->pwb_ring * MAX_FEAM_PAD_ROWS + row;
               } else {
                  h.tpc_col = -1;
                  h.tpc_row = -1;
               }
               h.time_ns = wpos_ns;
               h.amp  = wamp;
               hits->fPadHits.push_back(h);

               if (0) {
                  printf("hit: pwb%02d, c%dr%d, seqsca %3d, tpc col %2d, row %3d, time %4.0f, amp %4.0f\n", imodule, c->pwb_column, c->pwb_ring, seqsca, h.tpc_col, h.tpc_row, wpos_ns, wamp);
               }
            }
         }
         
         if (doPrint) {
            printf("chan %3d: baseline %8.1f, rms %8.1f, min %8.1f, max %8.1f, amp %8.1f, wpos %5.1f, hit %d\n", ichan, bmean, brms, wmin, wmax, wamp, wpos, hit);
            //exit(1);
         }
         
         // save first waveform
         
         //if (fHC[seqchan]->hwaveform_first->GetEntries() == 0) {
         //   if (doPrint)
         //      printf("saving first waveform %d\n", seqchan);
         //   for (int i=0; i<nbins; i++)
         //      fHC[seqchan]->hwaveform_first->SetBinContent(i+1, c->adc_samples[i]);
         //}
         
         // save biggest waveform
         
         //if (wamp > fHC[seqchan]->fMaxWamp) {
         //   fHC[seqchan]->fMaxWamp = wamp;
         //   if (doPrint)
         //      printf("saving biggest waveform %d\n", seqchan);
         //   for (int i=0; i<nbins; i++)
         //      fHC[seqchan]->hwaveform_max->SetBinContent(i+1, c->adc_samples[i]);
         //}
         
         // add to average waveform
         
         //for (int j=0; j< nbins; j++)
         //   fHC[seqchan]->hwaveform_avg->AddBinContent(j+1, c->adc_samples[j]);
         //fHC[seqchan]->nwf++;
         
         // save biggest drift region waveform

#if 0         
         if (dpos > idrift_cut){
            if(damp > fHC[seqchan]->fMaxWampDrift) {
               fHC[seqchan]->fMaxWampDrift = damp;
               if (doPrint)
                  printf("saving biggest drift waveform %d\n", seqchan);
               for (int i=0; i<nbins; i++)
                  fHC[seqchan]->hwaveform_max_drift->SetBinContent(i+1, c->adc_samples[i]);
            }
            
            // add to average waveform
            
            for (int j=0; j< nbins; j++)
               fHC[seqchan]->hwaveform_avg_drift->AddBinContent(j+1, c->adc_samples[j]);
            fHC[seqchan]->nwf_drift++;
            
         }
#endif

         if (scachan_is_pad || scachan_is_fpn) {
            hbmean_all->Fill(bmean);
            hbrms_all->Fill(brms);
            hbmean_pwb_prof->Fill(seqpwb, bmean);
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
            //h_adc_range_all->Fill(wmax-wmin);
            //h_adc_range_baseline->Fill(bmax-bmin);
            //h_adc_range_drift->Fill(dmax-dmin);
            
            //hf->hbmean_prof->Fill(seqsca, bmean);
            //hf->hbrms_prof->Fill(seqsca, brms);
            //hf->hbrange_prof->Fill(seqsca, bmax-bmin);

            if (brms < ADC_RMS_PAD_MAX) {
               h_all_pad_baseline_good_count->Fill(seqpwbsca);
               hf->hbmean_bis_prof->Fill(seqsca, bmean);
               hf->hbrms_bis_prof->Fill(seqsca, brms);
               hf->hbrange_bis_prof->Fill(seqsca, bmax-bmin);
            }
            
            if (scachan_is_pad) {
               hf->hbrms_pads->Fill(brms);
            }
            
            if (scachan_is_fpn) {
               hf->hbrms_fpn->Fill(brms);
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
         
         //if (dpos > idrift_cut) {
         //   hdrift_amp_all->Fill(damp);
         //   hdrift_amp_all_pedestal->Fill(damp);
         //   hdrift_amp_all_above_pedestal->Fill(damp);
         //   
         //   if (damp > hit_amp_threshold) {
         //      hdrift_led_all->Fill(dpos);
         //      hdrift_led2amp->Fill(dpos, damp);
         //      if (seqpad >= 0) {
         //         hf->hnhits_pad_drift->Fill(seqpad);
         //      }
         //   }
         //}
         
         if (hit) {
            //hnhits->Fill(seqchan);
            hled_hit->Fill(wpos);
            hamp_hit->Fill(wamp);
            
            if (fPulser) {
               h_pulser_led_hit->Fill(wpos);
               hf->h_pulser_hit_amp_seqpad_prof->Fill(-1, 0); // force plot to start from 0
               hf->h_pulser_hit_amp_seqpad_prof->Fill(seqpad, wamp);
               hf->h_pulser_hit_time_seqpad_prof->Fill(seqpad, wpos);
               hf->h_pulser_hit_time_seqsca_prof->Fill(seqsca, wpos);
               hf->h_pulser_hit_amp->Fill(wamp);
               hf->h_pulser_hit_time->Fill(wpos);
               hf->h_pulser_hit_time_zoom->Fill(wpos);
               //if (seqsca == 4)
               //   hf->h_pulser_hit_time_seqsca4_zoom->Fill(wpos);
            }
            
            hf->h_nhits_seqsca->Fill(seqsca);
            hf->h_hit_time_seqsca->Fill(seqsca, wpos);
            hf->h_hit_amp_seqsca->Fill(seqsca, wamp);
            hf->h_hit_amp_seqpad->Fill(seqpad, wamp);
            
            //if (wamp >= 10000 && wamp <= 40000) {
            //   hf->h_amp_seqsca_prof->Fill(seqsca, wamp);
            //   hf->h_amp_seqpad_prof->Fill(seqpad, wamp);
            //}
            
            //h_amp_hit_col->Fill((ifeam*4 + col)%(MAX_FEAM_PAD_COL*MAX_FEAM), wamp);
            
            if (seqpad >= 0) {
               hf->h_nhits_seqpad->Fill(seqpad);
               if (!spike) {
                  hf->hnhits_pad_nospike->Fill(seqpad);
               }
            }
         }
         
         if (spike) {
            hf->h_spike_seqsca->Fill(seqsca);
         }
   
         if (fpn_is_ok) {
            fCountGoodFpn ++;
            //printf("XXX good fpn count %d\n", fCountGoodFpn);
         } else {
            fCountBadFpn ++;
            //printf("XXX bad fpn count %d\n", fCountBadFpn);
         }
         
         if (pad_is_ok) {
            //fCountGoodPad ++;
            //printf("XXX good pad count %d\n", fCountGoodPad);
         } else {
            fCountBadPad ++;
            //printf("XXX bad pad count %d\n", fCountBadPad);
         }
         
         hf->h_spike_seqsca->Fill(1); // event counter marker
         hf->hnhitchan->Fill(nhitchan_feam);
      }

      hnhitchan->Fill(nhitchan);

      return flow;
   }

   void AnalyzeSpecialEvent(TARunInfo* runinfo, TMEvent* event)
   {
      if (fTrace)
         printf("AnalyzeSpecialEvent, run %d, event serno %d, id 0x%04x, data size %d\n", runinfo->fRunNo, event->serial_number, (int)event->event_id, event->data_size);
   }
};

class PwbModuleFactory: public TAFactory
{
public:
   PwbFlags fFlags;
   
public:
   void Usage()
   {
      printf("PwbModuleFactory flags:\n");
      printf("--plot1 <fPlotPad>\n");
      printf("--pwbwfsuppress -- enable waveform suppression code\n");
   }

   void Init(const std::vector<std::string> &args)
   {
      printf("PwbModuleFactory::Init!\n");

      for (unsigned i=0; i<args.size(); i++) {
         if (args[i] == "--plot1")
            fFlags.fPlotPad = atoi(args[i+1].c_str());
         if (args[i] == "--pwbwfsuppress")
            fFlags.fWfSuppress = true;
      }
   }

   void Finish()
   {
      printf("PwbModuleFactory::Finish!\n");
   }

   TARunObject* NewRunObject(TARunInfo* runinfo)
   {
      printf("PwbModuleFactory::NewRunObject, run %d, file %s\n", runinfo->fRunNo, runinfo->fFileName.c_str());
      return new PwbModule(runinfo, &fFlags);
   }
};

static TARegister tar(new PwbModuleFactory);

/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */
