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

//#include "Unpack.h"

#include "Waveform.h"

#define DELETE(x) if (x) { delete (x); (x) = NULL; }

#define MEMZERO(p) memset((p), 0, sizeof(p))

static uint8_t getUint8(const void* ptr, int offset)
{
   return *(uint8_t*)(((char*)ptr)+offset);
}

static uint16_t getUint16be(const void* ptr, int offset)
{
   uint8_t *ptr8 = (uint8_t*)(((char*)ptr)+offset);
   return ((ptr8[0]<<8) | ptr8[1]);
}

static uint32_t getUint32be(const void* ptr, int offset)
{
   uint8_t *ptr8 = (uint8_t*)(((char*)ptr)+offset);
   return (ptr8[0]<<24) | (ptr8[1]<<16) | (ptr8[2]<<8) | ptr8[3];
}

static uint16_t getUint16le(const void* ptr, int offset)
{
   uint8_t *ptr8 = (uint8_t*)(((char*)ptr)+offset);
   return ((ptr8[1]<<8) | ptr8[0]);
}

static uint32_t getUint32le(const void* ptr, int offset)
{
   uint8_t *ptr8 = (uint8_t*)(((char*)ptr)+offset);
   return (ptr8[3]<<24) | (ptr8[2]<<16) | (ptr8[1]<<8) | ptr8[0];
}

class FeamPacket
{
public:
   uint32_t cnt;
   uint16_t n;
   uint16_t x511;
   uint16_t buf_len;
   uint32_t ts_start;
   uint32_t ts_trig;
   int off;
   bool error;

public:
   FeamPacket(); // ctor
   void Unpack(const char* data, int size);
   void Print() const;
};

FeamPacket::FeamPacket()
{
   error = true;
}

void FeamPacket::Unpack(const char* data, int size)
{
   error = true;
   
   off = 0;
   cnt = getUint32le(data, off); off += 4;
   n   = getUint16le(data, off); off += 2;
   x511 = getUint16le(data, off); off += 2;
   buf_len = getUint16le(data, off); off += 2;
   if (n == 0) {
      ts_start = getUint32le(data, off); off += 8;
      ts_trig  = getUint32le(data, off); off += 8;
   } else {
      ts_start = 0;
      ts_trig  = 0;
   }

   error = false;
}

void FeamPacket::Print() const
{
   printf("decoded %2d bytes, ", off);
   printf("cnt %6d, n %3d, x511 %3d, buf_len %4d, ts_start 0x%08x, ts_trig 0x%08x, ",
          cnt,
          n,
          x511,
          buf_len,
          ts_start,
          ts_trig);
   printf("error %d", error);
}

class FeamAdcData
{
public:
   uint32_t cnt;
   uint32_t ts_start;
   uint32_t ts_trig;

   bool error;

   uint32_t next_n;

   int fSize;
   char* fPtr;

public:
   FeamAdcData(const FeamPacket* p)
   {
      assert(p->n == 0);

      cnt = p->cnt;
      ts_start = p->ts_start;
      ts_trig = p->ts_trig;

      fSize = 0;
      fPtr = NULL;

      next_n = 0;

      error = false;
   }

   ~FeamAdcData() // dtor
   {
      if (fPtr)
         free(fPtr);
      fPtr = NULL;
      fSize = 0;
   }

   void AddData(const FeamPacket*p, const char* ptr, int size);
   void Print() const;
};

void FeamAdcData::Print() const
{
   printf("cnt %6d, ts_start 0x%08x, ts_trig 0x%08x, ",
          cnt,
          ts_start,
          ts_trig);
   printf("next_n %d, ", next_n);
   printf("size %d, ", fSize);
   printf("error %d", error);
}

void FeamAdcData::AddData(const FeamPacket*p, const char* ptr, int size)
{
   //printf("add %d size %d\n", p->n, size);
   if (p->n != next_n) {
      printf("wrong packet sequence: expected %d, got %d!\n", next_n, p->n);
      error = true;
      return;
   }

   assert(size >= 0);
   assert(size < 12000); // UDP packet size is 1500 bytes, jumbo frame up to 9000 bytes

   int new_size = fSize + size;
   char* new_ptr = (char*)realloc(fPtr, new_size);

   if (!new_ptr) {
      printf("cannot reallocate ADC buffer from %d to %d bytes!\n", fSize, new_size);
      error = true;
      return;
   }

   memcpy(new_ptr + fSize, ptr, size);

   fPtr = new_ptr;
   fSize = new_size;

   next_n = p->n + 1;
}

class FeamEVB
{
public:
   std::vector<FeamAdcData*> data;

   std::deque<FeamAdcData*> buf;

   FeamEVB()
   {
      data.push_back(NULL);
      data.push_back(NULL);
   }

   void AddPacket(int ifeam, const FeamPacket* p, const char* ptr, int size)
   {
      if (p->n == 0) {
         // 1st packet

         if (data[ifeam]) {
            printf("Complete event: FEAM %d: ", ifeam);
            data[ifeam]->Print();
            printf("\n");

            FeamAdcData* a = data[ifeam];
            data[ifeam] = NULL;

            // xxx

            buf.push_back(a);
         }

         printf("Start ew event: FEAM %d: ", ifeam);
         p->Print();
         printf("\n");

         data[ifeam] = new FeamAdcData(p);
      }

      FeamAdcData* a = data[ifeam];

      if (a == NULL) {
         // did not see the first event yet, cannot unpack
         delete p;
         return;
      }

      a->AddData(p, ptr, size);

      //a->Print();
      //printf("\n");

      delete p;
   }

   void Print() const
   {
      printf("FEAM evb status: %p %p, buffered %d\n", data[0], data[1], (int)buf.size());
   }

   FeamAdcData* Get()
   {
      if (buf.size() < 1)
         return NULL;

      FeamAdcData* a = buf.front();
      buf.pop_front();
      return a;
   }
};

static FeamEVB* xevb = NULL;

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
      //printf("Analyze, run %d, event serno %d, id 0x%04x, data size %d\n", runinfo->fRunNo, event->serial_number, (int)event->event_id, event->data_size);

      if (!fModule->fDoPads)
         return flow;
      
      if (event->event_id != 1 && event->event_id != 2)
         return flow;

      int force_plot = false;

      int adc[80][5120];

      const int xbins = 829;
      const int xchan = 79;

      //printf("event id %d\n", event->event_id);

      if (event->event_id == 1) {
         //*flags |= TAFlag_SKIP; // enable this to skip GRIF16 events

         const char* banks[] = { "BB01", "BB02", NULL };
         char *data = NULL;
         
         for (int i=0; banks[i]; i++) {
            TMBank* b = event->FindBank(banks[i]);
            if (b) {
               data = event->GetBankData(b);
               if (data) {
                  //printf("Have bank %s\n", banks[i]);
                  //HandleFeam(i, data, b->data_size);

                  if (b->data_size < 26) {
                     printf("bad FEAM %d packet length %d\n", i, b->data_size);
                     continue;
                  }
                  
                  FeamPacket* p = new FeamPacket();
                  
                  p->Unpack(data, b->data_size);
                  
                  assert(!p->error);
                  
                  if (!xevb)
                     xevb = new FeamEVB();
                  
                  xevb->AddPacket(i, p, data + p->off, p->buf_len);
               }
            }
         }

         //xevb->Print();

         FeamAdcData *a = xevb->Get();

         if (!a) {
            return flow;
         }

         printf("ZZZZZZZZZZZZZZZZZZ Processing FEAM event: ");
         a->Print();
         printf("\n");

         if (a->error) {
            delete a;
            return flow;
         }

         assert(a->next_n == 256);
         assert(a->fSize == 310688);

         const unsigned char* ptr = (const unsigned char*)a->fPtr;

         MEMZERO(adc);

         int count = 0;
         for (int ibin = 0; ibin < 511; ibin++) {
            for (int ichan = 0; ichan < 76; ichan++) {
               for (int iasic = 0; iasic < 4; iasic++) {
                  int v = ptr[0] | ((ptr[1])<<8);
                  if (iasic == 0) {
                     adc[ichan][ibin] = v;
                  }
                  ptr += 2;
                  count += 2;
               }
            }
         }

         printf("count %d\n", count);

         for (int ibin = 511; ibin < xbins; ibin++) {
            for (int ichan = 0; ichan < 76; ichan++) {
               adc[ichan][ibin] = adc[ichan][ibin-511];
            }
         }

         delete a;

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
         
         int count = 0;
         for (int ibin=0; ibin<xbins; ibin++) {
            for (int ichan=0; ichan<xchan; ichan++) {
               adc[ichan][ibin] = samples[count];
               count++;
            }
         }
         printf("got %d samples\n", count);
         
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
      }
         
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
               h33->SetBinContent(ibin+1, ww[33]->samples[ibin]);
            }
            h33->SetMinimum(0);
            h33->SetMaximum(17000);
            h33->Draw();
         }
         
         if (1) {
            fC->cd(6);
            TH1D* h34 = new TH1D("h34", "h34", xbins, 0, xbins);
            for (int ibin=0; ibin<xbins; ibin++) {
               h34->SetBinContent(ibin+1, ww[34]->samples[ibin]);
            }
            h34->SetMinimum(0);
            h34->SetMaximum(17000);
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

         TH1D* h = new TH1D("h", "h", xbins, 0, xbins);
         for (int ibin=0; ibin<xbins; ibin++) {
            h->SetBinContent(ibin+1, ww[fModule->fPlotPad]->samples[ibin]);
         }

         h->SetMinimum(0);
         h->SetMaximum(17000);
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
   printf("Finish!\n");

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
