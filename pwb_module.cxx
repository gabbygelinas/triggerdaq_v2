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

//#include "wfsuppress.h"
//#include "wfsuppress2.h"
#include "wfsuppress_pwb.h"

#include "ncfm.h"

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

   void CreateHistograms(TDirectory* hdir, int imodule, int icolumn, int iring, bool pulser, double pulser_start, double pulser_end)
   {
      char xname[128];
      char xtitle[128];

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

      sprintf(name,  "%s_nhitchan_map_seqsca", xname);
      sprintf(title, "%s hit channels vs seqsca; SCA*80 + readout index; number of hits", xtitle);
      h_nhitchan_seqsca = new TH1D(name, title, NUM_SEQSCA, 0.5, NUM_SEQSCA+0.5);

      sprintf(name,  "%s_hits_map_seqsca", xname);
      sprintf(title, "%s hits vs seqsca; SCA*80 + readout index; number of hits", xtitle);
      h_nhits_seqsca = new TH1D(name, title, NUM_SEQSCA, 0.5, NUM_SEQSCA+0.5);

      sprintf(name,  "%s_spike_seqsca", xname);
      sprintf(title, "%s spikes vs seqsca; SCA*80 + readout index; number of hits", xtitle);
      h_spike_seqsca = new TH1D(name, title, NUM_SEQSCA, 0.5, NUM_SEQSCA+0.5);

      sprintf(name,  "%s_hit_map_seqpad", xname);
      sprintf(title, "%s hits vs seqpad; col*4*18+row; number of hits", xtitle);
      h_nhits_seqpad  = new TH1D(name, title, MAX_FEAM_PAD_ROWS*MAX_FEAM_PAD_COL, -0.5, MAX_FEAM_PAD_ROWS*MAX_FEAM_PAD_COL-0.5);

      sprintf(name,  "%s_hit_map_seqpad_nospike", xname);
      sprintf(title, "%s hits with spikes removed vs seqpad; col*4*18+row; number of hits", xtitle);
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
         sprintf(title, "%s pulser hit amp; hit amp, adc counts", xtitle);
         h_pulser_hit_amp = new TH1D(name, title, ADC_BINS_PULSER, 0, ADC_RANGE_PULSER);

         sprintf(name,  "%s_pulser_hit_time", xname);
         sprintf(title, "%s pulser hit time; sca time bins", xtitle);
         h_pulser_hit_time = new TH1D(name, title, NUM_TIME_BINS, 0, NUM_TIME_BINS);

         sprintf(name,  "%s_pulser_hit_time_zoom", xname);
         sprintf(title, "%s pulser hit time zoom; hit amp, adc counts", xtitle);
         h_pulser_hit_time_zoom = new TH1D(name, title, 100, pulser_start, pulser_end);

         sprintf(name,  "%s_pulser_hit_amp_seqpad_prof", xname);
         sprintf(title, "%s pulser hit p.h. vs seqpad; seqpad = col*4*18+row; hit amp, adc counts", xtitle);
         h_pulser_hit_amp_seqpad_prof = new TProfile(name, title, MAX_FEAM_PAD_ROWS*MAX_FEAM_PAD_COL+1, -1.5, MAX_FEAM_PAD_ROWS*MAX_FEAM_PAD_COL-0.5);

         sprintf(name,  "%s_pulser_hit_time_seqpad_prof", xname);
         sprintf(title, "%s pulser hit time vs seqpad; seqpad = col*4*18+row; hit time, sca time bins", xtitle);
         h_pulser_hit_time_seqpad_prof = new TProfile(name, title, MAX_FEAM_PAD_ROWS*MAX_FEAM_PAD_COL, -0.5, MAX_FEAM_PAD_ROWS*MAX_FEAM_PAD_COL-0.5);
         h_pulser_hit_time_seqpad_prof->SetMinimum(pulser_start);
         h_pulser_hit_time_seqpad_prof->SetMaximum(pulser_end);

         sprintf(name,  "%s_pulser_hit_time_seqsca_prof", xname);
         sprintf(title, "%s pulser hit time vs seqsca; seqsca = SCA*80 + readout index; hit time, sca time bins", xtitle);
         h_pulser_hit_time_seqsca_prof = new TProfile(name, title, NUM_SEQSCA, -0.5, NUM_SEQSCA-0.5, pulser_start, pulser_end);
         h_pulser_hit_time_seqsca_prof->SetMinimum(pulser_start);
         h_pulser_hit_time_seqsca_prof->SetMaximum(pulser_end);
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

TH1D* WfToTH1D(const char* name, const char* title, const std::vector<int> &adc)
{
   size_t nbins = adc.size();
   TH1D* h = new TH1D(name, title, nbins, -0.5, nbins-0.5);
   
   for (size_t i=0; i<nbins; i++)
      h->SetBinContent(i+1, adc[i]);
   
   return h;
}

class PwbFlags
{
public:
   //int  fPlotPad = -1;
   //TCanvas* fPlotPadCanvas = NULL;
   bool fWfSuppress = false;
   int  fWfThreshold = 0;
   bool fWfSaveBad = false;

public:
   PwbFlags() // ctor
   {
   }

   ~PwbFlags() // dtor
   {
      //DELETE(fPlotPadCanvas);
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
   Ncfm* fCfm = NULL;
   NcfmParser* fCfmCuts = NULL;

   //std::vector<std::vector<std::vector<WfSuppress*>>> fWfSuppress;
   //std::vector<std::vector<std::vector<WfSuppress2*>>> fWfSuppress;
   std::vector<std::vector<std::vector<WfSuppressPwb*>>> fWfSuppress;
   std::vector<std::vector<std::vector<WfSuppressPwb*>>> fWfSuppressPwb;
   TH1D* fWfSuppressAdcAmp = NULL;
   TH1D* fWfSuppressAdcAmpPos = NULL;
   TH1D* fWfSuppressAdcAmpNeg = NULL;
   TH1D* fWfSuppressAdcAmpCumulKeepAll = NULL;
   TH1D* fWfSuppressAdcAmpCumulDropAll = NULL;
   TH2D* fWfSuppressAdcAmpCumulKeepMap = NULL;
   TH2D* fWfSuppressAdcAmpCumulDropMap = NULL;
   std::vector<TH1D*> fWfSuppressAdcAmpCumulKeep;
   std::vector<TH1D*> fWfSuppressAdcAmpCumulDrop;
   TH2D* fWfSuppressAdcMinMap = NULL;
   std::vector<TH1D*> fWfSuppressAdcMin;

   TH1D* h_all_fpn_count = NULL;

   //TProfile* h_all_fpn_rms_prof = NULL;
   TProfile* h_all_fpn_mean_bis_prof = NULL;
   TProfile* h_all_fpn_rms_bis_prof = NULL;
   TProfile* h_all_fpn_wrange_bis_prof = NULL;

#if 0
   std::vector<TProfile*> h_all_fpn_mean_per_col_prof;
   std::vector<TProfile*> h_all_fpn_rms_per_col_prof;
#endif

   std::vector<int> fCountWfSave;

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

   TH1D* h_fpn_wrange = NULL;

   TH1D* hpad_ph = NULL;
   TH1D* hpad_ph_zoom_pedestal = NULL;
   TH1D* hpad_ph_above_pedestal = NULL;
   TH2D* hpad_time_ph = NULL;
   TH1D* hpad_time_cut_ph = NULL;
   TH1D* hpad_time_cut_ph_ns = NULL;
   TH1D* hpad_ph_cut_time = NULL;

   //TH1D* hdrift_amp_all;
   //TH1D* hdrift_amp_all_pedestal;
   //TH1D* hdrift_amp_all_above_pedestal;
   //TH1D* hdrift_led_all;
   //TH2D* hdrift_led2amp;

   //TH1D* hnhits;
   TH1D* hhit_time_ns = NULL;
   TH1D* hhit_ph = NULL;
   //TH2D* h_amp_hit_col = NULL;

   bool  fPulser = true;
   TH1D* h_pulser_hit_amp  = NULL;
   TH1D* h_pulser_hit_time = NULL;
   TProfile* h_pulser_hit_amp_seqpwbsca_prof  = NULL;
   TProfile* h_pulser_hit_time_seqpwbsca_prof = NULL;

   TH1D* hnhitchan = NULL;

   //TH2D* hpadmap;

   TDirectory* hdir_summary = NULL;
   TDirectory* hdir_pulser  = NULL;
   TDirectory* hdir_wfsuppress = NULL;
   TDirectory* hdir_waveforms  = NULL;
   TDirectory* hdir_pwb  = NULL;
   TDirectory* hdir_pads = NULL;
   TDirectory* hdir_pwb_hit_map_pads = NULL;
   std::vector<PwbHistograms*> fHF;
   std::vector<TH1D*> fHPwbHitMapPads;

   int fCountTestScaEvents = 0;
   int fCountBadScaEvents = 0;
   int fCountBadSca = 0;
   int fCountGoodFpn = 0;
   int fCountBadFpn = 0;

   int fCountBadPad = 0;

   bool fEnableTestMode = false;
   int  fTestMode = 0; // see PWB manual

   std::vector<int> fTestPatternErrors;

   bool fTrace = false;

   int fPulserStart = 0;
   int fPulserEnd   = 0;

   PwbModule(TARunInfo* runinfo, PwbFlags* f)
      : TARunObject(runinfo)
   {
      if (fTrace)
         printf("PwbModule::ctor!\n");

      fModuleName = "pwb_module";

      fFlags = f;
      std::string agcfmdb_path="agcfmdb";
      if (getenv("AGRELEASE"))
      {
         agcfmdb_path=getenv("AGRELEASE");
         agcfmdb_path+="/agana/agcfmdb";
      }
      fCfm     = new Ncfm(agcfmdb_path.c_str());
      fCfmCuts = fCfm->ParseFile("pwb", "cuts", runinfo->fRunNo);

      fPulserStart = fCfmCuts->GetInt("sca_bin_pulser_start", 200);
      fPulserEnd   = fCfmCuts->GetInt("sca_bin_pulser_end",   240);

      runinfo->fRoot->fOutputFile->cd();
      hdir_pads = gDirectory->mkdir("pads");
      hdir_pads->cd(); // select correct ROOT directory

      hdir_summary = hdir_pads->mkdir("summary");
      hdir_pulser  = hdir_pads->mkdir("pulser");
      hdir_wfsuppress = hdir_pads->mkdir("wfsuppress");
      hdir_waveforms  = hdir_pads->mkdir("waveforms");
      hdir_pwb = hdir_pads->mkdir("pwb");
      hdir_pwb_hit_map_pads = hdir_pads->mkdir("pwb_hit_map_seqpad");

      if (fPulser) {
         hdir_pulser->cd();

         h_pulser_hit_amp = new TH1D("pulser_hit_amp", "pulser hit amp; hit amp , adc counts", ADC_BINS_PULSER, 0, ADC_RANGE_PULSER);
         h_pulser_hit_time = new TH1D("pulser_hit_time", "pulser hit time; sca time bins", 100, fPulserStart, fPulserEnd);

         h_pulser_hit_amp_seqpwbsca_prof  = new TProfile("pulser_hit_amp_seqpwbsca_prof", "pulser hit amp vs seqpwbsca; seqpwb*4 + isca; hit amp , adc counts", NUM_PWB*4, -0.5, NUM_PWB*4-0.5);
         h_pulser_hit_amp_seqpwbsca_prof->SetMinimum(0);
         h_pulser_hit_amp_seqpwbsca_prof->SetMaximum(ADC_RANGE_PULSER);

         h_pulser_hit_time_seqpwbsca_prof = new TProfile("pulser_hit_time_seqpwbsca_prof", "pulser hit time vs seqpwbsca; seqpwb*4 + isca; sca time bins", NUM_PWB*4, -0.5, NUM_PWB*4-0.5);
         h_pulser_hit_time_seqpwbsca_prof->SetMinimum(fPulserStart);
         h_pulser_hit_time_seqpwbsca_prof->SetMaximum(fPulserEnd);
      }

      fCountWfSave.resize(PWB_MODULE_LAST+1);
   }

   ~PwbModule()
   {
      if (fTrace)
         printf("PwbModule::dtor!\n");
      for (unsigned i=0; i<fHF.size(); i++) {
         DELETE(fHF[i]);
      }

      //for (unsigned i=0; i<h_all_fpn_mean_per_col.size(); i++) {
      //   DELETE(h_all_fpn_mean_per_col[i]);
      //}

      //for (unsigned i=0; i<h_all_fpn_rms_per_col.size(); i++) {
      //   DELETE(h_all_fpn_rms_per_col[i]);
      //}

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
         printf("BeginRun, run %d, file %s\n", runinfo->fRunNo, runinfo->fFileName.c_str());
      //time_t run_start_time = runinfo->fOdb->odbReadUint32("/Runinfo/Start time binary", 0, 0);
      //printf("ODB Run start time: %d: %s", (int)run_start_time, ctime(&run_start_time));


      runinfo->fRoot->fOutputFile->cd(); // select correct ROOT directory

      runinfo->fOdb->RB("Equipment/Ctrl/Settings/PWB/enable_test_mode", &fEnableTestMode);
      runinfo->fOdb->RI("Equipment/Ctrl/Settings/PWB/test_mode", &fTestMode);
      printf("test mode: enabled: %d, mode %d\n", fEnableTestMode, fTestMode);

      runinfo->fOdb->RB("Equipment/Ctrl/Settings/FwPulserEnable", &fPulser);
      printf("pulser mode: enabled: %d\n", fPulser);
   }

   void CreateHistograms(TARunInfo* /*runinfo*/)
   {
      if (hbmean_all) // already created
         return;

      hdir_summary->cd();

      //int max_fpn = NUM_PWB*MAX_FEAM_SCA*4;
      int max_pad = NUM_PWB*MAX_FEAM_SCA;

      h_all_fpn_count = new TH1D("all_fpn_count", "count of all fpn channels; sca+4*(ring+8*column); fpn count", max_pad, -0.5, max_pad-0.5);

      //h_all_fpn_rms_prof  = new TProfile("all_fpn_rms_prof", "rms of all fpn channels; fpn+4*(sca+4*(ring+8*column)); rms, adc counts", max_fpn, -0.5, max_fpn-0.5);

      h_all_fpn_mean_bis_prof  = new TProfile("all_fpn_mean_bis_prof", "mean of all fpn channels; sca+4*(ring+8*column)", max_pad, -0.5, max_pad-0.5);
      h_all_fpn_rms_bis_prof  = new TProfile("all_fpn_rms_bis_prof", "rms of all fpn channels; sca+4*(ring+8*column)", max_pad, -0.5, max_pad-0.5);

      h_all_fpn_wrange_bis_prof  = new TProfile("all_fpn_wrange_bis_prof", "waveform range of all fpn channels; sca+4*(ring+8*column)", max_pad, -0.5, max_pad-0.5);

#if 0
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
#endif

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

      h_fpn_wrange           = new TH1D("fpn_wrange",             "fpn waveform range (max-min); adc counts", 100, 0, ADC_RANGE);

      hpad_ph                = new TH1D("hpad_ph",                "pad waveform pulse height; adc counts", 100, 0, ADC_RANGE);
      hpad_ph_zoom_pedestal  = new TH1D("hpad_ph_zoom_pedestal",  "pad waveform pulse height, zoom on pedestal area; adc counts", 100, 0, ADC_RANGE_PED);
      hpad_ph_above_pedestal = new TH1D("hpad_ph_above_pedestal", "pad waveform pulse height, away from pedestal area; adc counts", 100, ADC_RANGE_PED, ADC_RANGE);

      hpad_time_ph        = new TH2D("hpad_time_ph", "pad hit p.h. vs time; sca time bins; adc counts", 100, 0, MAX_TIME_BINS, 100, 0, ADC_RANGE);

      hpad_time_cut_ph    = new TH1D("hpad_time_cut_ph",    "pad hit time with p.h. cut; sca time bins", 100, 0, MAX_TIME_BINS);
      hpad_time_cut_ph_ns = new TH1D("hpad_time_cut_ph_ns", "pad hit time with p.h. cut; time, ns", 100, 0, MAX_TIME_NS);

      hpad_ph_cut_time    = new TH1D("hpad_ph_cut_time",    "pad p.h. with time cut; adc counts", 100, 0, ADC_RANGE);

      //hdrift_amp_all = new TH1D("drift_amp", "drift region pulse height", 100, 0, ADC_RANGE);
      //hdrift_amp_all_pedestal = new TH1D("drift_amp_pedestal", "drift region pulse height, zoom on pedestal area", 100, 0, ADC_RANGE_PED);
      //hdrift_amp_all_above_pedestal = new TH1D("drift_amp_above_pedestal", "drift region pulse height, away from pedestal area", 100, ADC_RANGE_PED, ADC_RANGE);
      //hdrift_led_all = new TH1D("drift_led", "drift region pulse leading edge, sca time bins, above pedestal", 100, 0, NUM_TIME_BINS);
      //hdrift_led2amp = new TH2D("drift_led2amp", "drift region pulse amp vs time, sca time bins, above pedestal", 100, 0, NUM_TIME_BINS, 100, 0, ADC_RANGE);

      //hnhits = new TH1D("hnhits", "hits per channel", nchan, -0.5, nchan-0.5);
      hhit_time_ns = new TH1D("hhit_time_ns", "pad hit time; time, ns", 100, 0, MAX_TIME_NS);
      hhit_ph      = new TH1D("hhit_ph",   "pad hit pulse height; adc counts", 100, 0, ADC_RANGE);

      hnhitchan = new TH1D("hnhitchan", "number of hit channels per event", 100, 0, 1000);

      //hpadmap = new TH2D("hpadmap", "map from TPC pad number (col*4*18+row) to SCA readout channel (sca*80+chan)", 4*4*18, -0.5, 4*4*18-0.5, NUM_SEQSCA, 0.5, NUM_SEQSCA+0.5);

      if (1) {
         hdir_pwb_hit_map_pads->cd();
         for (int i=0; i<=PWB_MODULE_LAST; i++) {
            char xname[50];
            char xtitle[50];
            char name[100];
            char title[100];
            sprintf(xname, "pwb%02d", i);
            sprintf(xtitle, "pwb%02d", i);
            sprintf(name,  "%s_hit_map_seqpad", xname);
            sprintf(title, "%s hits vs seqpad; col*4*18+row", xtitle);
            TH1D* h = new TH1D(name, title, MAX_FEAM_PAD_ROWS*MAX_FEAM_PAD_COL, -0.5, MAX_FEAM_PAD_ROWS*MAX_FEAM_PAD_COL-0.5);
            h->SetMinimum(0);
            fHPwbHitMapPads.push_back(h);
         }
         hdir_summary->cd();
      }

      if (fFlags->fWfSuppress) {
         hdir_wfsuppress->cd();

         if (fWfSuppressAdcAmp == NULL) {
            int min = 0;
            int max = 4100;
            fWfSuppressAdcAmp = new TH1D("WfSuppress ADC amp", "WfSuppress ADC amp; adc counts", max-min, min, max);
         }
         
         if (fWfSuppressAdcAmpPos == NULL) {
            int min = 0;
            int max = 4100;
            fWfSuppressAdcAmpPos = new TH1D("WfSuppress ADC amp pos", "WfSuppress ADC amp pos; adc counts", max-min, min, max);
         }
         
         if (fWfSuppressAdcAmpNeg == NULL) {
            int min = -4100;
            int max = 0;
            fWfSuppressAdcAmpNeg = new TH1D("WfSuppress ADC amp neg", "WfSuppress ADC amp neg; adc counts", max-min, min, max);
         }
         
         if (fWfSuppressAdcAmpCumulKeepAll == NULL) {
            int min = -1;
            int max = 4200;
            fWfSuppressAdcAmpCumulKeepAll = new TH1D("WfSuppress cumul keep", "WfSuppress cumulative kept channels; ch_threshold, adc counts", max-min+1, min-0.5, max+0.5);
         }
         
         if (fWfSuppressAdcAmpCumulDropAll == NULL) {
            int min = -1;
            int max = 4200;
            fWfSuppressAdcAmpCumulDropAll = new TH1D("WfSuppress cumul drop", "WfSuppress cumulative dropped channels; ch_threshold, adc_counts", max-min+1, min-0.5, max+0.5);
         }

         if (fWfSuppressAdcAmpCumulKeepMap == NULL) {
            int min = -1;
            int max = 4200;
            fWfSuppressAdcAmpCumulKeepMap = new TH2D("WfSuppress cumul keep map", "WfSuppress cumulative kept channels; pwbNN; ch_threshold, adc counts", PWB_MODULE_LAST+1, 0-0.5, PWB_MODULE_LAST+1-0.5, max-min+1, min-0.5, max+0.5);
         }
         
         if (fWfSuppressAdcAmpCumulDropMap == NULL) {
            int min = -1;
            int max = 4200;
            fWfSuppressAdcAmpCumulDropMap = new TH2D("WfSuppress cumul drop map", "WfSuppress cumulative dropped channels; pwbNN; ch_threshold, adc_counts", PWB_MODULE_LAST+1, 0-0.5, PWB_MODULE_LAST+1-0.5, max-min+1, min-0.5, max+0.5);
         }

         for (int i=0; i<=PWB_MODULE_LAST; i++) {
            char name[100];
            char title[100];

            sprintf(name,  "pwb%02d_keep", i);
            sprintf(title, "pwb%02d cumulative kept channels; ch_threshold, adc counts", i);

            int min = -1;
            int max = 4200;
            TH1D* h = new TH1D(name, title, max-min+1, min-0.5, max+0.5);

            fWfSuppressAdcAmpCumulKeep.push_back(h);
         }
         
         for (int i=0; i<=PWB_MODULE_LAST; i++) {
            char name[100];
            char title[100];

            sprintf(name,  "pwb%02d_drop", i);
            sprintf(title, "pwb%02d cumulative dropped channels; ch_threshold, adc counts", i);

            int min = -1;
            int max = 4200;
            TH1D* h = new TH1D(name, title, max-min+1, min-0.5, max+0.5);
            fWfSuppressAdcAmpCumulDrop.push_back(h);
         }

         if (fWfSuppressAdcMinMap == NULL) {
            int min = -2050;
            int max = 2050;
            fWfSuppressAdcMinMap = new TH2D("WfSuppress_adc_min_map", "WfSuppress adc_min for each PWB; pwbNN; adc counts", PWB_MODULE_LAST+1, 0-0.5, PWB_MODULE_LAST+1-0.5, max-min, min, max);
         }

         for (int i=0; i<=PWB_MODULE_LAST; i++) {
            char name[100];
            char title[100];
            sprintf(name,  "pwb%02d_adc_min", i);
            sprintf(title, "pwb%02d adc_min for channel suppression; adc counts", i);
            int min = -2050;
            int max = 2050;
            TH1D* h = new TH1D(name, title, (max-min)/10.0, min, max);
            fWfSuppressAdcMin.push_back(h);
         }
         hdir_summary->cd();
      }
   }

   void EndRun(TARunInfo* runinfo)
   {
      if (fTrace)
         printf("PwbModule::EndRun, run %d\n", runinfo->fRunNo);
      //time_t run_stop_time = runinfo->fOdb->odbReadUint32("/Runinfo/Stop time binary", 0, 0);
      //printf("ODB Run stop time: %d: %s", (int)run_stop_time, ctime(&run_stop_time));

      printf("PwbModule::EndRun: test for bad SCA: total events %d, bad events %d, bad sca %d, bad fpn %d, good fpn %d, bad pad %d\n", fCountTestScaEvents, fCountBadScaEvents, fCountBadSca, fCountBadFpn, fCountGoodFpn, fCountBadPad);

      if (fTestPatternErrors.size() > 0) {
         for (unsigned imodule = 0; imodule < fTestPatternErrors.size(); imodule++) {
            int count = fTestPatternErrors[imodule];
            if (count) {
               printf("PwbModule::EndRun: pwb%02d: %d test pattern mismatches\n", imodule, count);
            }
         }
      }
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

         int imodule    = c->imodule;
         int isca       = c->sca;
         int iri        = c->sca_readout;

         assert(isca>=0 && isca<4);
         assert(iri>=0 && iri<80);

         while (imodule >= (int)pwb_ptr.size()) {
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
            fHF[imodule]->CreateHistograms(hdir_pwb, imodule, e->hits[i]->pwb_column, e->hits[i]->pwb_ring, fPulser, fPulserStart, fPulserEnd);
         }
      }
   }

   TAFlowEvent* AnalyzeFlowEvent(TARunInfo* runinfo, TAFlags* flags, TAFlowEvent* flow)
   {
      //bool verbose = false;
      
      //printf("Analyze, run %d, event serno %d, id 0x%04x, data size %d\n", runinfo->fRunNo, event->serial_number, (int)event->event_id, event->data_size);

      AgEventFlow *ef = flow->Find<AgEventFlow>();

      if (!ef || !ef->fEvent)
      {
         *flags|=TAFlag_SKIP_PROFILE;
         return flow;
      }

      FeamEvent* e = ef->fEvent->feam;

      if (!e)
      {
         *flags|=TAFlag_SKIP_PROFILE;
         return flow;
      }


      if (0) {
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

      int ibaseline_start = fCfmCuts->GetInt("sca_bin_baseline_start", 10);
      int ibaseline_end   = fCfmCuts->GetInt("sca_bin_baseline_end",  100);

      // unused
      //int iwire_start   = fCfmCuts->GetInt("sca_bin_pc_start",  130);
      double iwire_middle = fCfmCuts->GetDouble("sca_bin_pc_middle", 145);
      int iwire_end     = fCfmCuts->GetInt("sca_bin_pc_end",    160);

      int idrift_start  = fCfmCuts->GetInt("sca_bin_drift_start", iwire_end);
      // unused
      //int idrift_cut    = fCfmCuts->GetInt("sca_bin_drift_cut",   iwire_end);
      int idrift_end    = fCfmCuts->GetInt("sca_bin_drift_end", 410);

      double wpos_min_ns = fCfmCuts->GetDouble("pad_hit_time_min_ns", 800.0);
      double wpos_max_ns = fCfmCuts->GetDouble("pad_hit_time_max_ns", 5600.0);

      double hit_amp_threshold = fCfmCuts->GetDouble("pad_hit_thr", 100);

      int nhitchan = 0;

      CreatePwbHistograms(e);
      
      //printf("PrintFeamChannels!\n");
      //PrintFeamChannels(e->hits);

      bool need_shift = false;
      need_shift |= (runinfo->fRunNo >= 1865 && runinfo->fRunNo <= 9999);

      if (need_shift) {
         CheckAndShiftFpn(e);
      }

      std::vector<bool> test_pattern_ok;

      WfSuppressPwb ch_supp;
      int ch_threshold = 0;

      if (fEnableTestMode) {
         runinfo->fOdb->RI("Equipment/Ctrl/Settings/PWB/ch_threshold", &ch_threshold);
         printf("PWB ch_threshold %d\n", ch_threshold);
         ch_supp.fThreshold = ch_threshold;
      }

      for (unsigned ii=0; ii<e->hits.size(); ii++) {
         FeamChannel* c = e->hits[ii];
         if (!c)
            continue;

         bool bad_wf = false;

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

         if (fEnableTestMode) {
            bool ok = true;
            //printf("imodule %d, sca %d, ri %d\n", c->imodule, c->sca, c->sca_readout);
            //if (c->imodule != 78)
            //   continue;
            //if (c->sca_readout != 2)
            //   continue;

            int ri = c->sca_readout;

            switch (fTestMode) {
               // see https://daqstore.triumf.ca/AgWiki/index.php/PWB#ESPER_Variables
            default: break;
            case 0: { // fixed pattern 0xa5a
               const int pattern = 0xa5a;
               for (unsigned i=0; i<c->adc_samples.size(); i++) {
                  int a = c->adc_samples[i];
                  if ((a&0xFFF) != pattern) {
                     ok = false;
                     printf("BBB0 imodule %02d, sca %d, ri %2d, bin %d: sample %d 0x%03x expected 0x%03x\n", c->imodule, c->sca, ri, i, a, a&0xFFF, pattern&0xFFF);
                  }
               }
               break;
            }
            case 1: { // time bin counter
               for (unsigned i=0; i<c->adc_samples.size(); i++) {
                  int a = c->adc_samples[i];
                  int exp = i;
                  if ((a&0xFFF) != (exp&0xFFF)) {
                     ok = false;
                     printf("BBB1 imodule %02d, sca %d, ri %2d, bin %d: sample %d 0x%03x expected 0x%03x\n", c->imodule, c->sca, ri, i, a, a&0xFFF, exp&0xFFF);
                  }
               }
               break;
            }
            case 2: { // time bin counter with channel number
               for (unsigned i=0; i<c->adc_samples.size(); i++) {
                  int a = c->adc_samples[i];
                  int xa = a & 0x1FF; // 9 bits of bin number
                  int xri = (a>>9)&0x7; // 3 bits of readout index
                  int exp_xri = (ri-1)&0x7;
                  
                  //printf("bin %d: sample %d 0x%04x, ri %d, xri %d\n", i, a, a, exp_xri, xri);
                  
                  int exp_xa = i;
                  //if (i==0) {
                  //   exp_xa = 509;
                  //   exp_xri = (exp_xri-1)&0x7;
                  //}
                  if ((c->sca_readout != 1) && (xa != exp_xa || xri != exp_xri)) {
                     ok = false;
                     printf("BBB2 imodule %02d, sca %d, ri %2d, bin %d: sample %d 0x%04x, xa %d expected %d, xri %d expected %d\n", c->imodule, c->sca, ri, i, a, a, xa, exp_xa, xri, exp_xri);
                  }
               }
               break;
            }
            case 3: { // sequential adc sample counter
               for (unsigned i=0; i<c->adc_samples.size(); i++) {
                  int a = c->adc_samples[i];
                  int xa = a & 0x1FF; // 9 bits of bin number
                  
                  int exp_xa = ((ri-1) + (i-0)*79) & 0x1FF;
                  
                  //printf("bin %d: sample %d 0x%04x, ri %d, xa %d expected %d\n", i, a, a, ri, xa, exp_xa);
                  
                  if (ri != 1) {
                     if (xa != exp_xa) {
                        ok = false;
                        printf("BBB3 imodule %02d, sca %d, ri %2d, bin %d: sample %d 0x%04x, xa %d expected %d\n", c->imodule, c->sca, c->sca_readout, i, a, a, xa, exp_xa);
                     }
                  }
               }
               break;
            }
            case 4: { // channel suppression test {ch_crossed_out,trig_pos,trig_neg,adc[8:0]}
               break;
            }
            case 5: { // channel suppression test {trig,adc[10:0]}
               //if (c->imodule != 78)
               //   break;
               //if (c->sca_readout != 5)
               //   break;
               ch_supp.Reset();
               bool trigprev = false;
               bool crossed = false;
               bool tcrossed = false;
               std::string buf;
               for (unsigned i=0; i<c->adc_samples.size(); i++) {
                  int a = c->adc_samples[i];
                  int t = (a & 0x800) != 0;
                  int xa = (a&0x7FF);
                  bool trig = ch_supp.Add(xa);
                  crossed |= trigprev;
                  tcrossed |= t;

#if 0
                  char xbuf[1024];
                  sprintf(xbuf, "BBB5 imodule %02d, sca %d, ri %2d, bin %d: sample %d 0x%03x 0x%03x, trig %d, t %d/%d c %d/%d, supp: ", c->imodule, c->sca, c->sca_readout, i, a, a&0xFFF, xa, trig, t, trigprev, tcrossed, crossed);
                  buf += xbuf;
                  buf += ch_supp.PrintToString();
                  buf += "\n";
#endif

                  if ((t != trigprev) || (tcrossed != crossed)) {
                     ok = false;
                  }
                  trigprev = trig;
               }
               //if (!ok)
               //   exit(1);
               if (!ok) {
                  printf("%s", buf.c_str());
               }
               break;
            }
            case 6: { // channel suppression test {ch_crossed_min,adc[10:0]}
               //if (c->imodule != 78)
               //   break;
               //if (c->sca_readout != 5)
               //   break;
               std::string buf;
               bool trigprev = false;
               bool trigsum = false;
               for (unsigned i=0; i<c->adc_samples.size(); i++) {
                  int a = c->adc_samples[i];
                  bool t = (a & 0x800) != 0;
                  int xa = (a&0x7FF);
                  bool trig = (xa <= ch_threshold);
                  if (i<=32)
                     trig = false;

                  trigsum |= trigprev;

                  char xbuf[1024];
                  sprintf(xbuf, "BBB6 imodule %02d, sca %d, ri %2d, bin %d: sample %d 0x%03x 0x%03x, t %d, cmp %d/%d", c->imodule, c->sca, c->sca_readout, i, a, a&0xFFF, xa, trig, trigsum, t);
                  buf += xbuf;
                  buf += "\n";

                  if (trigsum != t) {
                     ok = false;
                  }

                  trigprev = trig;
               }
               if (!ok) {
                  printf("%s", buf.c_str());
               }
               break;
            }
            } // switch

            if (!ok) {
               printf("BBBB imodule %02d, sca %d, ri %2d, failed test pattern test\n", c->imodule, c->sca, c->sca_readout);
            }
            
            while ((int)test_pattern_ok.size() <= c->imodule) {
               test_pattern_ok.push_back(true);
            }

            assert(c->imodule < (int)test_pattern_ok.size());

            test_pattern_ok[c->imodule] = test_pattern_ok[c->imodule] && ok;

            if (!ok) {
               while ((int)fTestPatternErrors.size() <= c->imodule) {
                  fTestPatternErrors.push_back(0);
               }
               
               assert(c->imodule < (int)fTestPatternErrors.size());
               
               fTestPatternErrors[c->imodule] = fTestPatternErrors[c->imodule] + 1;
            }

            continue; // ADC test patterns do not go through the rest of normal analysis
         }

         if (fFlags->fWfSuppress) { // compute data suppression

            if (imodule >= (int)fWfSuppress.size())
               fWfSuppress.resize(imodule+1);
            
            if (isca >= (int)fWfSuppress[imodule].size())
               fWfSuppress[imodule].resize(MAX_FEAM_SCA);
            
            if (ichan >= (int)fWfSuppress[imodule][isca].size())
               fWfSuppress[imodule][isca].resize(MAX_FEAM_READOUT);
            
            //printf("imodule %d, size %d\n", imodule, (int)fWfSuppress.size());
            //printf("isca %d, size %d\n", isca, (int)fWfSuppress[imodule].size());
            //printf("ichan %d, size %d\n", isca, (int)fWfSuppress[imodule][isca].size());

#if 0            
            WfSuppress *s = fWfSuppress[imodule][isca][ichan];
            if (!s) {
               s = new WfSuppress();
               fWfSuppress[imodule][isca][ichan] = s;
            }
#endif

#if 0            
            WfSuppress2 *s = fWfSuppress[imodule][isca][ichan];
            if (!s) {
               s = new WfSuppress2();
               fWfSuppress[imodule][isca][ichan] = s;
            }
#endif

            WfSuppressPwb *s = fWfSuppress[imodule][isca][ichan];
            if (!s) {
               s = new WfSuppressPwb();
               fWfSuppress[imodule][isca][ichan] = s;
            }
            
            //s->Init(c->adc_samples[sfirst], fFlags->fWfThreshold);
            s->Reset();
            s->fThreshold = fFlags->fWfThreshold;

            //std::vector<int> samples_base;
            std::vector<int> samples_amp;
            
            bool keep = false;
            int ampmin = 0;
            int ampmax = 0;
            int adcmin = c->adc_samples[0];
            for (unsigned i=0; i<c->adc_samples.size(); i++) {
               if (c->adc_samples[i] < adcmin)
                  adcmin = c->adc_samples[i];
               bool k = s->Add(c->adc_samples[i]);
               //uint16_t base = s->GetBase();
               //int16_t amp = s->GetAmp();
               //samples_base.push_back(base);
               //samples_amp.push_back(amp);
               int base = s->fBaseline;
               int amp = s->fAdcValue;
               if (amp > ampmax)
                  ampmax = amp;
               if (amp < ampmin)
                  ampmin = amp;
               samples_amp.push_back(amp);
               keep |= k;
               if (0) {
                  printf("pwb %02d, sca %d, chan %2d: bin %3d, adc %d, base %d, amp %4d, keep %d %d, state: %s\n", imodule, isca, ichan, i, c->adc_samples[i], base, amp, k, keep, s->PrintToString().c_str());
               }
            }
            
            //double xampmax = fabs(s->GetAmpMax());
            double xampmin = fabs(ampmin);
            //double xamp = std::min(xampmax, xampmin);
            double xamp = xampmin;
            //if (s->GetClipped())
            //   xamp = 0xFFF + 1;
            
            fWfSuppressAdcAmp->Fill(xamp);
            fWfSuppressAdcAmpPos->Fill(ampmax);
            fWfSuppressAdcAmpNeg->Fill(ampmin);

            //fWfSuppressAdcAmpCumulKeep->Fill(xamp);
            for (int i=0; i<xamp; i++) {
               fWfSuppressAdcAmpCumulKeepAll->Fill(i);
               fWfSuppressAdcAmpCumulKeep[imodule]->Fill(i);
               fWfSuppressAdcAmpCumulKeepMap->Fill(imodule, i);
            }

            //fWfSuppressAdcAmpCumulDrop->Fill(xamp);
            for (int i=xamp; i<4200; i++) {
               fWfSuppressAdcAmpCumulDropAll->Fill(i);
               fWfSuppressAdcAmpCumulDrop[imodule]->Fill(i);
               fWfSuppressAdcAmpCumulDropMap->Fill(imodule, i);
            }
            fWfSuppressAdcMinMap->Fill(imodule, adcmin);
            fWfSuppressAdcMin[imodule]->Fill(adcmin);

            int range = ampmax - ampmin;

            //printf("pwb %02d, sca %d, chan %2d: wfsuppress: %s, keep: %d, xamp %d\n", imodule, isca, ichan, s->PrintToString().c_str(), keep, (int)xamp);

            if (!keep) {
               printf("TTTT: ");
               printf("pwb %02d, sca %d, chan %2d: wfsuppress: %s, , keep: %d, xamp %d, mismatch!\n", imodule, isca, ichan, s->PrintToString().c_str(), keep, (int)xamp);

               static int count = 0;

               if (count < 100) {
                  char name[256];
                  char title[256];
                  
                  hdir_wfsuppress->cd();
                  
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

#if 0
            if (pwb_column >= 0 && pwb_column < (int)h_all_fpn_mean_per_col_prof.size()) {
               h_all_fpn_mean_per_col_prof[pwb_column]->Fill(pwb_ring*4+isca, bmean);
            }

            if (pwb_column >= 0 && pwb_column < (int)h_all_fpn_rms_per_col_prof.size()) {
               h_all_fpn_rms_per_col_prof[pwb_column]->Fill(pwb_ring*4+isca, brms);
            }
#endif

            //h_all_fpn_rms_bis->Fill(seqpwbsca, brms);

            double wrange = wmax - wmin;

            h_fpn_wrange->Fill(wrange);
            h_all_fpn_wrange_bis_prof->Fill(seqpwbsca, wrange);

            if (wrange > 10) {
               printf("XXX bad fpn, pwb%02d, sca %d, readout %d, scachan %d, col %d, row %d, wmin %f, wmax %f, wrange %f\n", imodule, isca, ichan, scachan, col, row, wmin, wmax, wrange);
               fpn_is_ok = false;
               bad_wf = true;
            } else if (brms > ADC_RMS_FPN_MIN && brms < ADC_RMS_FPN_MAX) {
            } else {
               printf("XXX bad fpn, pwb%02d, sca %d, readout %d, scachan %d, col %d, row %d, bmin %f, bmax %f, in hex 0x%04x, brms %f\n", imodule, isca, ichan, scachan, col, row, bmin, bmax, (uint16_t)bmin, brms);
               fpn_is_ok = false;
               bad_wf = true;
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
         
#if 0
         if (fFlags->fWfSuppress) {
            // this is special code, only use on pulser data
            // where no real hits are expected. K.O.
            if (fabs(wmax - wmin) > 1000) {
               printf("XXX bad suppress, pwb%02d, sca %d, readout %d, scachan %d, col %d, row %d, wmin %7.1f, wmax %7.1f, diff %7.1f\n", imodule, isca, ichan, scachan, col, row, wmin, wmax, wmax - wmin);
               bad_wf = true;
            } else {
               bad_wf = false;
            }
         }
#endif

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
                  bad_wf = true;
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

         double cfd_thr_pct = fCfmCuts->GetDouble("pad_cfd_thr_pct", 0.5);
         double cfd_thr = cfd_thr_pct*wamp;
         
         //int wpos = find_pulse(c->adc_samples, nbins, bmean, -1.0, cfd_thr);
         double wpos = find_pulse_time(c->adc_samples, nbins, bmean, -1.0, cfd_thr);

         double time_bin = 1000.0/62.5; // 62.5 MHz SCA sampling (write) clock
         
         double wpos_ns = (wpos - iwire_middle)*time_bin + 1000.0;
         
         //printf("ZZZ wpos %f, middle %f, time bin %f, wpos_ns %f %f\n", wpos, iwire_middle, time_bin, wpos_ns, (wpos - iwire_middle)*time_bin + 1000.0);
         
         //double damp = bmean - dmin;
         //int dpos = find_pulse(c->adc_samples, idrift_start, idrift_end, bmean, -1.0, damp/2.0);
         
         // decide if we have a hit
         
         bool hit_time = false;
         bool hit_amp = false;
         bool hit = false;
         
         if (scachan_is_pad && (wpos_ns > wpos_min_ns) && (wpos_ns < wpos_max_ns)) {
            hit_time = true;
         }
         
         if (fPulser) {
            if ((wpos > fPulserStart) && (wpos < fPulserEnd)) {
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

         // save bad waveform

         if (fFlags->fWfSaveBad && bad_wf) {
            if (fCountWfSave[imodule] < 10) {
               fCountWfSave[imodule]++;

               char name[256];
               char title[256];
               
               hdir_waveforms->cd();
               
               if (scachan_is_reset) {
                  sprintf(name, "pwb_%02d_sca_%d_chan_%02d_bad_reset_event_%d", imodule, isca, ichan, ef->fEvent->counter);
                  sprintf(title, "bad reset waveform pwb %02d, sca %d, chan %2d, event %d", imodule, isca, ichan, ef->fEvent->counter);
               } else if (scachan_is_fpn) {
                  sprintf(name, "pwb_%02d_sca_%d_chan_%02d_bad_fpn_event_%d", imodule, isca, ichan, ef->fEvent->counter);
                  sprintf(title, "bad FPN waveform pwb %02d, sca %d, chan %2d, event %d", imodule, isca, ichan, ef->fEvent->counter);
               } else {
                  sprintf(name, "pwb_%02d_sca_%d_chan_%02d_bad_pad_event_%d", imodule, isca, ichan, ef->fEvent->counter);
                  sprintf(title, "bad pad waveform pwb %02d, sca %d, chan %2d, event %d", imodule, isca, ichan, ef->fEvent->counter);
               }

               printf("XXX saving waveform %s %s\n", name, title);

               WfToTH1D(name, title, c->adc_samples);
            } else {
               //printf("XXX pwb%02d waveform not saved, limit %d\n", imodule, fCountWfSave[imodule]);
            }
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
            hpad_ph->Fill(wamp);
            hpad_ph_zoom_pedestal->Fill(wamp);
            hpad_ph_above_pedestal->Fill(wamp);
            hpad_time_ph->Fill(wpos, wamp);
            if (hit_amp) {
               hpad_time_cut_ph->Fill(wpos);
               hpad_time_cut_ph_ns->Fill(wpos_ns);
            }

            if (hit_time) {
               hpad_ph_cut_time->Fill(wamp);
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
            hhit_time_ns->Fill(wpos_ns);
            hhit_ph->Fill(wamp);
            
            if (fPulser) {
               h_pulser_hit_amp->Fill(wamp);
               h_pulser_hit_time->Fill(wpos);
               h_pulser_hit_amp_seqpwbsca_prof->Fill(seqpwbsca, wamp);
               h_pulser_hit_time_seqpwbsca_prof->Fill(seqpwbsca, wpos);
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
               fHPwbHitMapPads[c->imodule]->Fill(seqpad);
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

      if (test_pattern_ok.size() > 0) {
         for (unsigned imodule = 0; imodule < test_pattern_ok.size(); imodule++) {
            bool ok = test_pattern_ok[imodule];
            if (!ok) {
               printf("BBBBB pwb%02d: test pattern mismatch error\n", imodule);
            }
         }
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
      //printf("--plot1 <fPlotPad>\n");
      printf("--pwb-wf-suppress -- enable waveform suppression code\n");
      printf("--pwb-wf-threshold -- set the waveform suppression threshold\n");
      printf("--pwb-wf-save-bad  -- write bad waveforms to the root output file\n");
   }

   void Init(const std::vector<std::string> &args)
   {
      printf("PwbModuleFactory::Init!\n");

      for (unsigned i=0; i<args.size(); i++) {
         //if (args[i] == "--plot1")
         //   fFlags.fPlotPad = atoi(args[i+1].c_str());
         if (args[i] == "--pwb-wf-suppress")
            fFlags.fWfSuppress = true;
         if (args[i] == "--pwb-wf-threshold")
            fFlags.fWfThreshold = atoi(args[i+1].c_str());
         if (args[i] == "--pwb-wf-save-bad")
            fFlags.fWfSaveBad = true;
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
