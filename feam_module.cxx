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

#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"

#include "Waveform.h"
#include "FeamEVB.h"
#include "Unpack.h"
#include "AgFlow.h"

#define DELETE(x) if (x) { delete (x); (x) = NULL; }

#define MEMZERO(p) memset((p), 0, sizeof(p))

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
   TH2D* hpadhits;

public:
   FeamHistograms()
   {
      hbmean_prof = NULL;
      hbrms_prof = NULL;
      hnhits = NULL;
      hnhits_pad = NULL;
      htime = NULL;
      hamp = NULL;
      hpadhits = NULL;
   };

   void CreateHistograms(int position, int nchan, int nbins)
   {
      char name[256];
      char title[256];

      sprintf(name,  "pos%02d_bmean_prof", position);
      sprintf(title, "feam pos %2d baseline mean vs chan", position);
      hbmean_prof = new TProfile(name, title, nchan, -0.5, nchan-0.5);

      sprintf(name,  "pos%02d_brms_prof", position);
      sprintf(title, "feam pos %2d baseline rms vs chan", position);
      hbrms_prof  = new TProfile(name, title,  nchan, -0.5, nchan-0.5);

      sprintf(name,  "pos%02d_hit_map", position);
      sprintf(title, "feam pos %2d hits vs SCA readout", position);
      hnhits      = new TH1D(name, title, nchan, -0.5, nchan-0.5);

      sprintf(name,  "pos%02d_hit_map_pads", position);
      sprintf(title, "feam pos %2d hits vs TPC pad", position);
      hnhits_pad  = new TH1D(name, title, MAX_FEAM_PAD_ROWS*MAX_FEAM_PAD_COL, -0.5, MAX_FEAM_PAD_ROWS*MAX_FEAM_PAD_COL-0.5);

      sprintf(name,  "pos%02d_hit_time", position);
      sprintf(title, "feam pos %2d hit time vs chan", position);
      htime       = new TH2D(name, title, nchan, -0.5, nchan-0.5, 50, 0, 500);

      sprintf(name,  "pos%02d_hit_amp", position);
      sprintf(title, "feam pos %2d hit p.h. vs chan", position);
      hamp        = new TH2D(name, title, nchan, -0.5, nchan-0.5, 50, 0, 17000);

      sprintf(name,  "pos%02d_hit_pad", position);
      sprintf(title, "feam pos %2d hit column vs pad", position);
      hpadhits    = new TH2D(name, title, MAX_FEAM_PAD_ROWS, -0.5, MAX_FEAM_PAD_ROWS-0.5, MAX_FEAM_PAD_COL, 0, MAX_FEAM_PAD_COL);
   }
};

class FeamRun: public TARunInterface
{
public:
   FeamModule* fModule;

   FILE *fin;
   TCanvas* fC;

   TH1D** hbmean;
   TH1D** hbrms;
   TH1D** hwaveform;
   TH1D** hwaveform1;
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

      hbmean_all = new TH1D("hbmean", "baseline mean", 100, 0, 17000);
      hbrms_all  = new TH1D("hbrms",  "baseline rms",  100, 0, 200);

      hbmean_prof = new TProfile("hbmean_prof", "baseline mean vs channel", nchan, -0.5, nchan-0.5);
      hbrms_prof  = new TProfile("hbrms_prof",  "baseline rms vs channel",  nchan, -0.5, nchan-0.5);

      hamp_all   = new TH1D("hamp",   "pulse height", 100, 0, 17000);
      hled_all   = new TH1D("hled",   "pulse leading edge, adc time bins", 100, 0, 900);

      h2led2amp  = new TH2D("h2led2amp", "pulse amp vs time, adc time bins", 100, 0, 900, 100, 0, 17000);

      hled_all_cut = new TH1D("hled_cut",   "pulse leading edge, adc time bins, with cuts", 100, 0, 900);
      hamp_all_cut = new TH1D("hamp_cut",   "pulse height, with cuts", 100, 0, 17000);

      hnhits = new TH1D("hnhits", "hits per channel", nchan, -0.5, nchan-0.5);
      hled_hit = new TH1D("hled_hit", "hit time, adc time bins", 100, 0, 900);
      hamp_hit = new TH1D("hamp_hit", "hit pulse height", 100, 0, 17000);
      hamp_hit_pedestal = new TH1D("hamp_hit_pedestal", "hit pulse height, zoom on pedestal", 100, 0, 300);

      for (int i=0; i<nfeam; i++) {
         fHF[i].CreateHistograms(i, nchan_feam, nbins);
      }

      // FIXME: who deletes this?
      hbmean = new TH1D*[nchan];
      hbrms  = new TH1D*[nchan];
      hwaveform = new TH1D*[nchan];
      hwaveform1 = new TH1D*[nchan];
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
         hbmean[i] = new TH1D(name, title, 100, 0, 17000);
      }

      pads->mkdir("baseline_rms")->cd();

      for (int i=0; i<nchan; i++) {
         char name[256];
         char title[256];
         sprintf(name, "hbrms%04d", i);
         sprintf(title, "chan %04d baseline rms", i);
         hbrms[i] = new TH1D(name, title, 100, 0, 200);
      }

      pads->mkdir("chan_waveform")->cd();

      for (int i=0; i<nchan; i++) {
         char name[256];
         char title[256];
         sprintf(name, "hwaveform%04d", i);
         sprintf(title, "chan %04d waveform", i);
         hwaveform[i] = new TH1D(name, title, nbins, -0.5, nbins-0.5);
      }

      pads->mkdir("chan_waveform_max")->cd();

      for (int i=0; i<nchan; i++) {
         char name[256];
         char title[256];
         sprintf(name, "hwaveform%04d_max", i);
         sprintf(title, "chan %04d waveform max p.h.", i);
         hwaveform1[i] = new TH1D(name, title, nbins, -0.5, nbins-0.5);
      }

      pads->mkdir("chan_amp")->cd();

      for (int i=0; i<nchan; i++) {
         char name[256];
         char title[256];
         sprintf(name, "hamp%04d", i);
         sprintf(title, "chan %04d pulse height", i);
         hamp[i] = new TH1D(name, title, 100, 0, 17000);
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

      int nchan = nfeam * nchan_feam;

      if (nbins == 0 || nchan == 0)
         return flow;

      Waveform** ww = new Waveform*[nchan];

      for (int i=0; i<nchan; i++)
         ww[i] = NULL;

      for (unsigned ifeam=0; ifeam<e->adcs.size(); ifeam++) {
         FeamAdcData* aaa = e->adcs[ifeam];
         if (!aaa)
            continue;
         for (int isca=0; isca<aaa->nsca; isca++) {
            for (int ichan=0; ichan<aaa->nchan; ichan++) {
               int xchan = ifeam*(aaa->nsca*aaa->nchan) + isca*aaa->nchan + ichan;
               ww[xchan] = new Waveform(aaa->nbins);
               for (int ibin=0; ibin<aaa->nbins; ibin++) {
                  ww[xchan]->samples[ibin] = (aaa->adc[isca][ichan][ibin])/4;
               }
            }
         }
      }

      // create pad hits flow event

      AgPadHitsFlow* hits = new AgPadHitsFlow(flow);
      flow = hits;

      // create histograms

      CreateHistograms(runinfo, nfeam, nchan_feam, nchan, nbins);

      int iplot = 0;
      double zmax = 0;

      for (int ichan=0; ichan<nchan; ichan++) {
         if (!ww[ichan])
            continue;

         double r;
         double b = baseline(ww[ichan], 10, 100, NULL, &r);

         if (b==0 && r==0)
            continue;

         double wmin = min(ww[ichan]);
         double wmax = max(ww[ichan]);
         double wamp = b - wmin;

         int xpos = led(ww[ichan], b, -1, wamp/2.0);

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
            printf("chan %3d: baseline %8.1f, rms %8.1f, min %8.1f, max %8.1f, amp %8.1f, xpos %3d, hit %d\n", ichan, b, r, wmin, wmax, wamp, xpos, hit);
         }

         if (1 || (xpos > 0 && xpos < 4000 && wamp > 1000)) {
            if (wamp > zmax) {
               if (doPrint)
                  printf("plot this one.\n");
               iplot = ichan;
               zmax = wamp;
            }
         }

         hbmean[ichan]->Fill(b);
         hbrms[ichan]->Fill(r);

         if (hwaveform[ichan]->GetEntries() == 0) {
            if (doPrint)
               printf("saving first waveform %d\n", ichan);
            for (int i=0; i<ww[ichan]->nsamples; i++)
               hwaveform[ichan]->SetBinContent(i+1, ww[ichan]->samples[i]);
         }

         if (wamp > fMaxWamp[ichan]) {
            fMaxWamp[ichan] = wamp;
            if (doPrint)
               printf("saving biggest waveform %d\n", ichan);
            for (int i=0; i<ww[ichan]->nsamples; i++)
               hwaveform1[ichan]->SetBinContent(i+1, ww[ichan]->samples[i]);
         }

         hamp[ichan]->Fill(wamp);
         hled[ichan]->Fill(xpos);

         hbmean_all->Fill(b);
         hbrms_all->Fill(r);
         hamp_all->Fill(wamp);
         hled_all->Fill(xpos);

         hbmean_prof->Fill(ichan, b);
         hbrms_prof->Fill(ichan, r);

         int ifeam = ichan/nchan_feam;
         int ichan_feam = ichan%nchan_feam;
         int nc_after = e->adcs[ifeam]->nchan;

         fHF[ifeam].hbmean_prof->Fill(ichan_feam, b);
         fHF[ifeam].hbrms_prof->Fill(ichan_feam, r);

         h2led2amp->Fill(xpos, wamp);

         if (wamp > 1000) {
            hled_all_cut->Fill(xpos);
         }

         if (xpos > 100 && xpos < 500) {
            hamp_all_cut->Fill(wamp);
         }

         if (hit) {
            hnhits->Fill(ichan);
            hled_hit->Fill(xpos);
            hamp_hit->Fill(wamp);
            hamp_hit_pedestal->Fill(wamp);

            fHF[ifeam].hnhits->Fill(ichan_feam);
            fHF[ifeam].htime->Fill(ichan_feam, xpos);
            fHF[ifeam].hamp->Fill(ichan_feam, wamp);
            auto pad = getPad(ichan_feam/nc_after, ichan_feam%nc_after);
            int col = pad.first;
            int row = pad.second;
            fHF[ifeam].hpadhits->Fill(row, col);
            fHF[ifeam].hnhits_pad->Fill(col*MAX_FEAM_PAD_ROWS + row);
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

         fC->Modified();
         fC->Draw();
         fC->Update();
      }

      if (fModule->fPlotPad >= 0) {
         if (!fModule->fPlotPadCanvas)
            fModule->fPlotPadCanvas = new TCanvas("FEAM PAD", "FEAM PAD", 900, 650);

         TCanvas*c = fModule->fPlotPadCanvas;

         c->cd();

         int nbins = ww[fModule->fPlotPad]->nsamples;
         TH1D* h = new TH1D("h", "h", nbins, 0, nbins);
         for (int ibin=0; ibin<nbins; ibin++) {
            h->SetBinContent(ibin+1, ww[fModule->fPlotPad]->samples[ibin]);
         }

         h->SetMinimum(-9000);
         h->SetMaximum(+9000);
         h->Draw();

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

      for (int i=0; i<nchan; i++) {
         if (ww[i]) {
            delete ww[i];
            ww[i] = NULL;
         }
      }

      delete ww;

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
