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

//#include "Unpack.h"

#include "Waveform.h"

#define DELETE(x) if (x) { delete (x); (x) = NULL; }

class FeamModule: public TAModuleInterface
{
public:
   void Init(const std::vector<std::string> &args);
   void Finish();
   TARunInterface* NewRun(TARunInfo* runinfo);

   bool fDoPads;
};

#define MAX_CHAN 80

struct FeamRun: public TARunInterface
{
   FeamModule* fModule;

   FILE *fin;
   TCanvas* fC;

   TH1D* hbmean[MAX_CHAN];
   TH1D* hbrms[MAX_CHAN];
   TH1D* hamp[MAX_CHAN];
   TH1D* hled[MAX_CHAN];
   
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
   
   FeamRun(TARunInfo* runinfo, FeamModule* m)
      : TARunInterface(runinfo)
   {
      printf("FeamRun::ctor!\n");
      fModule = m;

      fC = NULL;
      fin = NULL;

      if (!m->fDoPads)
         return;
      
      fC = new TCanvas();
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
      fin = fopen("/pool8tb/agdaq/pads/tpc02.1485571169.txt", "r"); // ??? events
      assert(fin);
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
      
      hbmean_all = new TH1D("hbmean", "baseline mean", 100, 0, 17000);
      hbrms_all  = new TH1D("hbrms",  "baseline rms",  100, 0, 200);

      hbmean_prof = new TProfile("hbmean_prof", "baseline mean vs channel", MAX_CHAN, -0.5, MAX_CHAN-0.5);
      hbrms_prof  = new TProfile("hbrms_prof",  "baseline rms vs channel",  MAX_CHAN, -0.5, MAX_CHAN-0.5);

      hamp_all   = new TH1D("hamp",   "pulse height", 100, 0, 17000);
      hled_all   = new TH1D("hled",   "pulse leading edge, adc time bins", 100, 0, 900);

      h2led2amp  = new TH2D("h2led2amp", "pulse amp vs time, adc time bins", 100, 0, 900, 100, 0, 17000);

      hled_all_cut = new TH1D("hled_cut",   "pulse leading edge, adc time bins, with cuts", 100, 0, 900);
      hamp_all_cut = new TH1D("hamp_cut",   "pulse height, with cuts", 100, 0, 17000);

      hnhits = new TH1D("hnhits", "hits per channel", MAX_CHAN, -0.5, MAX_CHAN-0.5);
      hled_hit = new TH1D("hled_hit", "hit time, adc time bins", 100, 0, 900);
      hamp_hit = new TH1D("hamp_hit", "hit pulse height", 100, 0, 17000);
      hamp_hit_pedestal = new TH1D("hamp_hit_pedestal", "hit pulse height, zoom on pedestal", 100, 0, 300);
   
      for (int i=0; i<MAX_CHAN; i++) {
         char name[256];
         char title[256];
         sprintf(name, "hbmean%02d", i);
         sprintf(title, "chan %02d baseline mean", i);
         hbmean[i] = new TH1D(name, title, 100, 0, 17000);
      }

      for (int i=0; i<MAX_CHAN; i++) {
         char name[256];
         char title[256];
         sprintf(name, "hbrms%02d", i);
         sprintf(title, "chan %02d baseline rms", i);
         hbrms[i] = new TH1D(name, title, 100, 0, 200);
      }

      for (int i=0; i<MAX_CHAN; i++) {
         char name[256];
         char title[256];
         sprintf(name, "hamp%02d", i);
         sprintf(title, "chan %02d pulse height", i);
         hamp[i] = new TH1D(name, title, 100, 0, 17000);
      }

      for (int i=0; i<MAX_CHAN; i++) {
         char name[256];
         char title[256];
         sprintf(name, "hled%02d", i);
         sprintf(title, "chan %02d pulse leading edge, adc bins", i);
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
      printf("Analyze, run %d, event serno %d, id 0x%04x, data size %d\n", runinfo->fRunNo, event->serial_number, (int)event->event_id, event->data_size);

      if (!fModule->fDoPads)
         return flow;
      
      if (event->event_id != 1)
         return flow;

      int force_plot = false;

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

      const int xbins = 829;
      const int xchan = 79;

      int adc[80][5120];

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

      // got all the data here

      Waveform** ww;

      ww = new Waveform*[xchan];

      for (int ichan=0; ichan<xchan; ichan++) {
         ww[ichan] = new Waveform(xbins);
         for (int ibin=0; ibin<xbins; ibin++)
            ww[ichan]->samples[ibin] = adc[ichan][ibin]/4;
      }

      int iplot = 0;
      double zmax = 0;

      for (int ichan=0; ichan<xchan; ichan++) {
         double r;
         double b = baseline(ww[ichan], 10, 60, NULL, &r);
         double wmin = min(ww[ichan]);
         double wmax = max(ww[ichan]);
         double wamp = b - wmin;

         int xpos = led(ww[ichan], b, -1, wamp/2.0);

         bool hit = false;

         if ((xpos > 0) && (xpos < 500) && (wamp > 200)) {
            hit = true;
         }

         printf("chan %3d: baseline %8.1f, rms %8.1f, min %8.1f, max %8.1f, amp %8.1f, xpos %3d, hit %d\n", ichan, b, r, wmin, wmax, wamp, xpos, hit);

         if (1 || (xpos > 0 && xpos < 4000 && wamp > 1000)) {
            if (wamp > zmax) {
               printf("plot this one.\n");
               iplot = ichan;
               zmax = wamp;
            }
         }

         hbmean[ichan]->Fill(b);
         hbrms[ichan]->Fill(r);
         hamp[ichan]->Fill(wamp);
         hled[ichan]->Fill(xpos);

         hbmean_all->Fill(b);
         hbrms_all->Fill(r);
         hamp_all->Fill(wamp);
         hled_all->Fill(xpos);

         hbmean_prof->Fill(ichan, b);
         hbrms_prof->Fill(ichan, r);

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
         }
      }

      bool do_plot = (runinfo->fRoot->fgApp != NULL);

      if (do_plot) {
         // plot waveforms

         fC->Clear();
         fC->Divide(2,3);
         
         if (1) {
            fC->cd(1);
            TH1D* hh = new TH1D("hh", "hh", xbins, 0, xbins);
            for (int ibin=0; ibin<xbins; ibin++) {
               hh->SetBinContent(ibin+1, adc[0][ibin]);
            }
            hh->Draw();
         }
         
         if (1) {
            fC->cd(2);
            TH1D* hhh = new TH1D("hhh", "hhh", xbins, 0, xbins);
            for (int ibin=0; ibin<xbins; ibin++) {
               hhh->SetBinContent(ibin+1, adc[0][ibin]);
            }
            hhh->SetMinimum(0);
            hhh->SetMaximum(66000);
            hhh->Draw();
         }
         
         if (1) {
            fC->cd(3);
            TH1D* hhh = new TH1D("hhhh", "hhhh", xbins, 0, xbins);
            for (int ibin=0; ibin<xbins; ibin++) {
               hhh->SetBinContent(ibin+1, ww[39]->samples[ibin]);
            }
            hhh->SetMinimum(0);
            hhh->SetMaximum(17000);
            hhh->Draw();
         }
         
         if (1) {
            fC->cd(4);
            TH1D* hhh = new TH1D("hhhhh", "hhhhh", xbins, 0, xbins);
            for (int ibin=0; ibin<xbins; ibin++) {
               hhh->SetBinContent(ibin+1, ww[iplot]->samples[ibin]);
            }
            hhh->SetMinimum(0);
            hhh->SetMaximum(17000);
            hhh->Draw();
         }
         
         if (1) {
            fC->cd(5);
            TH1D* h33 = new TH1D("h33", "h33", xbins, 0, xbins);
            for (int ibin=0; ibin<xbins; ibin++) {
               h33->SetBinContent(ibin+1, ww[0]->samples[ibin]);
            }
            h33->SetMinimum(0);
            h33->SetMaximum(17000);
            h33->Draw();
         }
         
         if (1) {
            fC->cd(6);
            TH1D* h34 = new TH1D("h34", "h34", xbins, 0, xbins);
            for (int ibin=0; ibin<xbins; ibin++) {
               h34->SetBinContent(ibin+1, ww[1]->samples[ibin]);
            }
            h34->SetMinimum(0);
            h34->SetMaximum(17000);
            h34->Draw();
         }
         
         fC->Modified();
         fC->Draw();
         fC->Update();
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

      return flow;
   }

   void AnalyzeSpecialEvent(TARunInfo* runinfo, TMEvent* event)
   {
      printf("AnalyzeSpecialEvent, run %d, event serno %d, id 0x%04x, data size %d\n", runinfo->fRunNo, event->serial_number, (int)event->event_id, event->data_size);
   }
};

void FeamModule::Init(const std::vector<std::string> &args)
{
   printf("Init!\n");

   fDoPads = true;

   for (unsigned i=0; i<args.size(); i++) {
      if (args[i] == "--nopads")
         fDoPads = false;
   }
}
   
void FeamModule::Finish()
{
   printf("Finish!\n");
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
