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
#define ADC_RANGE 65000
#define ADC_RANGE_RMS 500

#define NUM_SEQSCA (3*80+79)

class FeamModule: public TAModuleInterface
{
public:
   void Init(const std::vector<std::string> &args);
   void Finish();
   TARunInterface* NewRun(TARunInfo* runinfo);

   bool fDoPads;
   int  fPlotPad;
   TCanvas* fPlotPadCanvas;
};

class FeamHistograms
{
public:
   TProfile* hbmean_prof;
   TProfile* hbrms_prof;
   TH1D* hnhits;
   TH1D* hnhits_pad;
   TH2D* htime;
   TH2D* hamp;

public:
   FeamHistograms()
   {
      hbmean_prof = NULL;
      hbrms_prof = NULL;
      hnhits = NULL;
      hnhits_pad = NULL;
      htime = NULL;
      hamp = NULL;
   };

   void CreateHistograms(int position, int nbins)
   {
      char name[256];
      char title[256];

      sprintf(name,  "pos%02d_baseline_mean_prof", position);
      sprintf(title, "feam pos %2d baseline mean vs (SCA*80 + readout index)", position);
      hbmean_prof = new TProfile(name, title, NUM_SEQSCA, 0.5, NUM_SEQSCA+0.5);

      sprintf(name,  "pos%02d_baseline_rms_prof", position);
      sprintf(title, "feam pos %2d baseline rms vs (SCA*80 +  readout index)", position);
      hbrms_prof  = new TProfile(name, title,  NUM_SEQSCA, 0.5, NUM_SEQSCA+0.5);

      sprintf(name,  "pos%02d_hit_map", position);
      sprintf(title, "feam pos %2d hits vs (SCA*80 + readout index)", position);
      hnhits      = new TH1D(name, title, NUM_SEQSCA, 0.5, NUM_SEQSCA+0.5);

      sprintf(name,  "pos%02d_hit_map_pads", position);
      sprintf(title, "feam pos %2d hits vs TPC seq.pad (col*4*18+row)", position);
      hnhits_pad  = new TH1D(name, title, MAX_FEAM_PAD_ROWS*MAX_FEAM_PAD_COL, -0.5, MAX_FEAM_PAD_ROWS*MAX_FEAM_PAD_COL-0.5);

      sprintf(name,  "pos%02d_hit_time", position);
      sprintf(title, "feam pos %2d hit time vs (SCA*80 + readout index)", position);
      htime       = new TH2D(name, title, NUM_SEQSCA, 0.5, NUM_SEQSCA+0.5, 50, 0, 500);

      sprintf(name,  "pos%02d_hit_amp", position);
      sprintf(title, "feam pos %2d hit p.h. vs (SCA*80 + readout index)", position);
      hamp        = new TH2D(name, title, NUM_SEQSCA, 0.5, NUM_SEQSCA+0.5, 50, 0, ADC_RANGE);
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

class FeamRun: public TARunInterface
{
private:
   const padMap padMapper;
public:
   FeamModule* fModule;

   FILE *fin;
   TCanvas* fC;

   TH1D** hbmean;
   TH1D** hbrms;
   std::vector<TH1D*> hwaveform_first;
   std::vector<TH1D*> hwaveform_max;
   double *fMaxWamp;

   TH1D** hamp;
   TH1D** hled;


   TH1D* hbmean_all;
   TH1D* hbrms_all;

   TH1D* hamp_all;
   TH1D* hled_all;

   TH1D* hled_all_cut;
   TH1D* hamp_all_cut;

   TProfile* hbmean_prof;
   TProfile* hbrms_prof;

   TH2D* h2led2amp;

   TH1D* hnhits;
   TH1D* hled_hit;
   TH1D* hamp_hit;
   TH1D* hamp_hit_pedestal;

   TDirectory* hdir_waveform_first;


   FeamHistograms fHF[MAX_FEAM];

   FeamRun(TARunInfo* runinfo, FeamModule* m)
      : TARunInterface(runinfo)
   {
      printf("FeamRun::ctor!\n");
      fModule = m;

      hbmean_all = NULL;

      fC = NULL;
      fin = NULL;

      if (!m->fDoPads)
         return;

      fC = new TCanvas();

      fin = NULL;

      if (0) {
         //fin = fopen("/pool8tb/agdaq/pads/yair1485457343.txt", "r");
         //fin = fopen("/pool8tb/agdaq/pads/yair1485479547.txt", "r");
         //fin = fopen("/pool8tb/agdaq/pads/yair1485479547.txt", "r"); // 400 events
         //fin = fopen("/home/agdaq/online/src/yair1485563694.txt", "r");
         //fin = fopen("/home/agdaq/online/src/yair1485564028.txt", "r");
         //fin = fopen("/home/agdaq/online/src/yair1485564199.txt", "r");
         //fin = fopen("/pool8tb/agdaq/pads/tpc01.1485570050.txt", "r");
         //fin = fopen("/pool8tb/agdaq/pads/tpc01.1485570419.txt", "r");
         //fin = fopen("/pool8tb/agdaq/pads/yair1485564199.txt", "r"); // 84 events
         //fin = fopen("/pool8tb/agdaq/pads/tpc01.1485571153.txt", "r"); // ??? events
         //fin = fopen("/pool8tb/agdaq/pads/tpc02.1485571169.txt", "r"); // ??? events
         //fin = fopen("/pool8tb/agdaq/pads/tpc01.1485577866.txt", "r"); // 2700 events
         //fin = fopen("/pool8tb/agdaq/pads/tpc02.1485577870.txt", "r"); // 2700 events
         //fin = fopen("/home/agdaq/online/src/tpc01.1485825981.txt", "r");
         //fin = fopen("/home/agdaq/online/src/tpc01.1485987366.txt", "r"); // delay 450, pulse at bin 60
         fin = fopen("/home/agdaq/online/src/tpc01.1485989371.txt", "r"); // delay 350, pulse at bin 160
         assert(fin);
      }
   }

   ~FeamRun()
   {
      printf("FeamRun::dtor!\n");
      DELETE(fC);
   }

   void BeginRun(TARunInfo* runinfo)
   {
      printf("BeginRun, run %d, file %s\n", runinfo->fRunNo, runinfo->fFileName.c_str());
      time_t run_start_time = runinfo->fOdb->odbReadUint32("/Runinfo/Start time binary", 0, 0);
      printf("ODB Run start time: %d: %s", (int)run_start_time, ctime(&run_start_time));
      runinfo->fRoot->fOutputFile->cd(); // select correct ROOT directory

      if (!fModule->fDoPads)
         return;
   }

   void CreateHistograms(TARunInfo* runinfo, int nfeam, int nchan_feam, int nchan, int nbins)
   {
      if (hbmean_all) // already created
         return;

      runinfo->fRoot->fOutputFile->cd();
      TDirectory* pads = gDirectory->mkdir("pads");
      pads->cd(); // select correct ROOT directory

      pads->mkdir("summary")->cd();

      hbmean_all = new TH1D("hbmean", "baseline mean", 100, ADC_MIN, ADC_MAX);
      hbrms_all  = new TH1D("hbrms",  "baseline rms",  100, 0, ADC_RANGE_RMS);

      hbmean_prof = new TProfile("hbmean_prof", "baseline mean vs channel", nchan, -0.5, nchan-0.5);
      hbrms_prof  = new TProfile("hbrms_prof",  "baseline rms vs channel",  nchan, -0.5, nchan-0.5);

      hamp_all   = new TH1D("hamp",   "pulse height", 100, 0, ADC_RANGE);
      hled_all   = new TH1D("hled",   "pulse leading edge, adc time bins", 100, 0, 900);

      h2led2amp  = new TH2D("h2led2amp", "pulse amp vs time, adc time bins", 100, 0, 900, 100, 0, ADC_RANGE);

      hled_all_cut = new TH1D("hled_cut",   "pulse leading edge, adc time bins, with cuts", 100, 0, 900);
      hamp_all_cut = new TH1D("hamp_cut",   "pulse height, with cuts", 100, 0, ADC_RANGE);

      hnhits = new TH1D("hnhits", "hits per channel", nchan, -0.5, nchan-0.5);
      hled_hit = new TH1D("hled_hit", "hit time, adc time bins", 100, 0, 900);
      hamp_hit = new TH1D("hamp_hit", "hit pulse height", 100, 0, ADC_RANGE);
      hamp_hit_pedestal = new TH1D("hamp_hit_pedestal", "hit pulse height, zoom on pedestal", 100, 0, 300);

      for (int i=0; i<nfeam; i++) {
         fHF[i].CreateHistograms(i, nbins);
      }

      // FIXME: who deletes this?
      hbmean = new TH1D*[nchan];
      hbrms  = new TH1D*[nchan];
      hamp = new TH1D*[nchan];
      hled = new TH1D*[nchan];

      fMaxWamp = new double[nchan];
      for (int i=0; i<nchan; i++) {
         fMaxWamp[i] = 0;
      }

      pads->mkdir("baseline_mean")->cd();

      for (int i=0; i<nchan; i++) {
         char name[256];
         char title[256];
         sprintf(name, "hbmean%04d", i);
         sprintf(title, "chan %04d baseline mean", i);
         hbmean[i] = new TH1D(name, title, 100, ADC_MIN, ADC_MAX);
      }

      pads->mkdir("baseline_rms")->cd();

      for (int i=0; i<nchan; i++) {
         char name[256];
         char title[256];
         sprintf(name, "hbrms%04d", i);
         sprintf(title, "chan %04d baseline rms", i);
         hbrms[i] = new TH1D(name, title, 100, 0, ADC_RANGE_RMS);
      }

      hdir_waveform_first = pads->mkdir("chan_waveform_first");

      pads->mkdir("chan_waveform_max")->cd();

      for (int i=0; i<nchan; i++) {
         char name[256];
         char title[256];
         sprintf(name, "hwaveform%04d_max", i);
         sprintf(title, "chan %04d biggest waveform", i);
         hwaveform_max.push_back(new TH1D(name, title, nbins, -0.5, nbins-0.5));
      }

      pads->mkdir("chan_amp")->cd();

      for (int i=0; i<nchan; i++) {
         char name[256];
         char title[256];
         sprintf(name, "hamp%04d", i);
         sprintf(title, "chan %04d pulse height", i);
         hamp[i] = new TH1D(name, title, 100, 0, ADC_RANGE);
      }

      pads->mkdir("chan_led")->cd();

      for (int i=0; i<nchan; i++) {
         char name[256];
         char title[256];
         sprintf(name, "hled%04d", i);
         sprintf(title, "chan %04d pulse leading edge, adc bins", i);
         hled[i] = new TH1D(name, title, 100, 0, 900);
      }
   }

   void EndRun(TARunInfo* runinfo)
   {
      printf("EndRun, run %d\n", runinfo->fRunNo);
      time_t run_stop_time = runinfo->fOdb->odbReadUint32("/Runinfo/Stop time binary", 0, 0);
      printf("ODB Run stop time: %d: %s", (int)run_stop_time, ctime(&run_stop_time));
   }

   void PauseRun(TARunInfo* runinfo)
   {
      printf("PauseRun, run %d\n", runinfo->fRunNo);
   }

   void ResumeRun(TARunInfo* runinfo)
   {
      printf("ResumeRun, run %d\n", runinfo->fRunNo);
   }

   TAFlowEvent* Analyze(TARunInfo* runinfo, TMEvent* event, TAFlags* flags, TAFlowEvent* flow)
   {
      //printf("Analyze, run %d, event serno %d, id 0x%04x, data size %d\n", runinfo->fRunNo, event->serial_number, (int)event->event_id, event->data_size);

      if (!fModule->fDoPads)
         return flow;

      AgEventFlow *ef = flow->Find<AgEventFlow>();

      if (!ef || !ef->fEvent)
         return flow;

      FeamEvent* e = ef->fEvent->feam;

      if (!e) {
         return flow;
      }

      int force_plot = false;

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
            printf("ZZZ Processing FEAM event: ");
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
               printf("position %2d, cnt %4d, ts_trig: 0x%08x %14.3f usec, ts_incr %14.3f usec\n", m->fPosition, m->cnt, m->ts_trig, m->fTime*1e6, m->fTimeIncr*1e6);
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

      // loop over all waveforms
      
      //int iplot = 0;
      double zmax = 0;

      int ibaseline_start = 10;
      int ibaseline_end = 100;

      for (unsigned ifeam=0; ifeam<e->adcs.size(); ifeam++) {
         FeamAdcData* aaa = e->adcs[ifeam];
         if (!aaa)
            continue;
         for (int isca=0; isca<aaa->nsca; isca++) {
            for (int ichan=1; ichan<=aaa->nchan; ichan++) {
               const int* aptr = &aaa->adc[isca][ichan][0];

               unsigned seqchan = ifeam*(aaa->nsca*aaa->nchan) + isca*aaa->nchan + ichan;
               int seqsca = isca*80 + ichan;

               // consult the pad map

               int scachan = padMapper.channel[ichan];
               int col = -1;
               int row = -1;

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
               }
                
               if (scachan > 0) {
                  col = padMapper.padcol[isca][scachan];
                  row = padMapper.padrow[isca][scachan];
                  //printf("isca %d, ichan %d, scachan %d, col %d, row %d\n", isca, ichan, scachan, col, row);
                  assert(col>=0 && col<4);
                  assert(row>=0 && row<4*72);
               } else {
                  row = scachan; // special channel
               }
               
               char xname[256];
               char xtitle[256];
               sprintf(xname, "pos%02d_sca%d_chan%02d_scachan%02d_col%02d_row%02d", ifeam, isca, ichan, scachan, col, row);
               sprintf(xtitle, "FEAM pos %d, sca %d, readout chan %d, sca chan %d, col %d, row %d", ifeam, isca, ichan, scachan, col, row);

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

               // find pulses

               double wamp = bmean - wmin;

               int xpos = find_pulse(aptr, nbins, bmean, -1.0, wamp/2.0);

               // decide if we have a hit

               bool hit = false;
               
               if ((xpos > 150) && (xpos < 450) && (wamp > 600)) {
                  hit = true;
               }

               if (hit) {
                  AgPadHit h;
                  h.chan = ichan;
                  h.time = xpos;
                  h.amp  = wamp;
                  hits->fPadHits.push_back(h);
               }

               if (doPrint) {
                  printf("chan %3d: baseline %8.1f, rms %8.1f, min %8.1f, max %8.1f, amp %8.1f, xpos %3d, hit %d\n", ichan, bmean, brms, wmin, wmax, wamp, xpos, hit);
               }
               
               if (1 || (xpos > 0 && xpos < 4000 && wamp > 1000)) {
                  if (wamp > zmax) {
                     if (doPrint)
                        printf("plot this one.\n");
                     //iplot = ichan;
                     zmax = wamp;
                  }
               }

               // create plots and histograms

               //hbmean[seqchan]->Fill(bmean);
               //hbrms[seqchan]->Fill(brms);

               if (seqchan >= hwaveform_first.size()) {
                  for (unsigned i=hwaveform_first.size(); i<=seqchan; i++)
                     hwaveform_first.push_back(NULL);
               }

               if (hwaveform_first[seqchan] == NULL) {
                  char name[256];
                  char title[256];
                  sprintf(name, "hwaveform%04d_first_%s", seqchan, xname);
                  sprintf(title, "%s first waveform", xtitle);
                  hdir_waveform_first->cd();
                  hwaveform_first[seqchan] = new TH1D(name, title, nbins, -0.5, nbins-0.5);
                  //printf("seqchan %d, size %d, ptr %p, name %s, title %s\n", seqchan, hwaveform_first.size(), hwaveform_first[seqchan], name, title);
               }

               if (hwaveform_first[seqchan]->GetEntries() == 0) {
                  if (doPrint)
                     printf("saving first waveform %d\n", seqchan);
                  for (int i=0; i<nbins; i++)
                     hwaveform_first[seqchan]->SetBinContent(i+1, aptr[i]);
               }

#if 0
               if (wamp > fMaxWamp[seqchan]) {
                  fMaxWamp[seqchan] = wamp;
                  if (doPrint)
                     printf("saving biggest waveform %d\n", seqchan);
                  for (int i=0; i<nbins; i++)
                     hwaveform_max[seqchan]->SetBinContent(i+1, aptr[i]);
               }
#endif

               //hamp[seqchan]->Fill(wamp);
               //hled[seqchan]->Fill(xpos);
               
               hbmean_all->Fill(bmean);
               hbrms_all->Fill(brms);
               hamp_all->Fill(wamp);
               hled_all->Fill(xpos);
               
               //hbmean_prof->Fill(seqchan, bmean);
               //hbrms_prof->Fill(seqchan, brms);

               fHF[ifeam].hbmean_prof->Fill(seqsca, bmean);
               fHF[ifeam].hbrms_prof->Fill(seqsca, brms);
               
               h2led2amp->Fill(xpos, wamp);
               
               if (wamp > 1000) {
                  hled_all_cut->Fill(xpos);
               }
               
               if (xpos > 100 && xpos < 500) {
                  hamp_all_cut->Fill(wamp);
               }

               if (hit) {
                  hnhits->Fill(seqchan);
                  hled_hit->Fill(xpos);
                  hamp_hit->Fill(wamp);
                  hamp_hit_pedestal->Fill(wamp);
                  
                  fHF[ifeam].hnhits->Fill(seqsca);
                  fHF[ifeam].htime->Fill(seqsca, xpos);
                  fHF[ifeam].hamp->Fill(seqsca, wamp);

                  if (scachan > 0) {
                     fHF[ifeam].hnhits_pad->Fill(col*MAX_FEAM_PAD_ROWS + row);
                  }
               }
            }
         }
      }

      bool do_plot = (runinfo->fRoot->fgApp != NULL);

      if (do_plot) {
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

      if (fModule->fPlotPad >= 0) {
         if (!fModule->fPlotPadCanvas)
            fModule->fPlotPadCanvas = new TCanvas("FEAM PAD", "FEAM PAD", 900, 650);

         TCanvas*c = fModule->fPlotPadCanvas;

         c->cd();

#if 0
         int nbins = ww[fModule->fPlotPad]->nsamples;
         TH1D* h = new TH1D("h", "h", nbins, 0, nbins);
         for (int ibin=0; ibin<nbins; ibin++) {
            h->SetBinContent(ibin+1, ww[fModule->fPlotPad]->samples[ibin]);
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

      *flags |= TAFlag_DISPLAY;

      return flow;
   }

   void AnalyzeSpecialEvent(TARunInfo* runinfo, TMEvent* event)
   {
      printf("AnalyzeSpecialEvent, run %d, event serno %d, id 0x%04x, data size %d\n", runinfo->fRunNo, event->serial_number, (int)event->event_id, event->data_size);
   }
};

void FeamModule::Init(const std::vector<std::string> &args)
{
   printf("FeamModule::Init!\n");

   fDoPads = true;
   fPlotPad = -1;
   fPlotPadCanvas = NULL;

   for (unsigned i=0; i<args.size(); i++) {
      if (args[i] == "--nopads")
         fDoPads = false;
      if (args[i] == "--plot1")
         fPlotPad = atoi(args[i+1].c_str());
   }
}

void FeamModule::Finish()
{
   printf("FeamModule::Finish!\n");

   DELETE(fPlotPadCanvas);
}

TARunInterface* FeamModule::NewRun(TARunInfo* runinfo)
{
   printf("FeamModule::NewRun, run %d, file %s\n", runinfo->fRunNo, runinfo->fFileName.c_str());
   return new FeamRun(runinfo, this);
}

static TARegisterModule tarm(new FeamModule);

/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */
