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
   ~FeamPacket(); // dtor
   void Unpack(const char* data, int size);
   void Print() const;
};

//static int x1count = 0;

FeamPacket::FeamPacket()
{
   //printf("FeamPacket: ctor!\n");
   //printf("FeamPacket: count %d!\n", x1count++);
   error = true;
}

FeamPacket::~FeamPacket()
{
   //printf("FeamPacket: dtor!\n");
   //x1count--;
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

//static int x2count = 0;

class FeamAdcData
{
public:
   int module;

   uint32_t cnt;
   uint32_t ts_start;
   uint32_t ts_trig;

   bool error;

   uint32_t next_n;

   int fSize;
   char* fPtr;

   double   fTime;
   double   fTsAbsNs;

   uint32_t fTsIncr;
   double   fTsIncrNs;

public:
   FeamAdcData(const FeamPacket* p, int xmodule)
   {
      //printf("FeamAdcData: ctor! %d\n", x2count++);
      assert(p->n == 0);

      module = xmodule;

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
      //printf("FeamAdcData: dtor!\n"); x2count--;
      if (fPtr)
         free(fPtr);
      fPtr = NULL;
      fSize = 0;
   }

   void AddData(const FeamPacket*p, int module, const char* ptr, int size);
   void Finalize();
   void Print() const;
};

void FeamAdcData::Finalize()
{
   if (error) {
      return;
   }

   if (next_n != 256) {
      error = true;
      return;
   }
   
   if (fSize != 310688) {
      error = true;
      return;
   }
}

/*
   ZZZ Run 421
   ZZZ Processing FEAM event: module  4, cnt     84, ts_start 0xc8f30f8a, ts_trig 0xc8f311fa, next_n 256, size 310688, error 0
   ZZZ Processing FEAM event: module  5, cnt     97, ts_start 0xc8f6d5d4, ts_trig 0xc8f6d844, next_n 256, size 310688, error 0
   ZZZ Processing FEAM event: module  6, cnt     93, ts_start 0xc8f91bb4, ts_trig 0xc8f91e24, next_n 256, size 310688, error 0
   ZZZ Processing FEAM event: module  7, cnt     83, ts_start 0xc8f5c6dc, ts_trig 0xc8f5c94c, next_n 256, size 310688, error 0
   ZZZ Processing FEAM event: module  0, cnt     80, ts_start 0xc8f46ecc, ts_trig 0xc8f4713c, next_n 256, size 310688, error 0
   ZZZ Processing FEAM event: module  1, cnt   1465, ts_start 0xc8f2ce84, ts_trig 0xc8f2d0f4, next_n 256, size 310688, error 0
   ZZZ Processing FEAM event: module  2, cnt     63, ts_start 0xc8f2892c, ts_trig 0xc8f28b9c, next_n 256, size 310688, error 0
   ZZZ Processing FEAM event: module  3, cnt     62, ts_start 0xc8ef4a56, ts_trig 0xc8ef4cc6, next_n 256, size 310688, error 0
*/

const unsigned xts[8] = { 0xc8f4713c, 0xc8f2d0f4, 0xc8f28b9c, 0xc8ef4cc6, 0xc8f311fa, 0xc8f6d844, 0xc8f91e24, 0xc8f5c94c };

void FeamAdcData::Print() const
{
   printf("module %2d, ", module);
   printf("cnt %6d, ts_start 0x%08x, ts_trig 0x%08x, ",
          cnt,
          ts_start,
          ts_trig);
   printf("next_n %d, ", next_n);
   printf("size %d, ", fSize);
   printf("relts 0x%08x, ", ts_trig-xts[module]);
   printf("ts_incr 0x%08x, ", fTsIncr);
   printf("error %d", error);
}

void FeamAdcData::AddData(const FeamPacket*p, int xmodule, const char* ptr, int size)
{
   assert(xmodule == module);

   //printf("add %d size %d\n", p->n, size);
   if (p->n != next_n) {
      printf("module %2d, cnt %6d, wrong packet sequence: expected %d, got %d!\n", module, cnt, next_n, p->n);
      next_n = p->n + 1;
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

struct TsSyncEntry
{
   uint32_t ts;
   int epoch;
   double time;

   TsSyncEntry(uint32_t xts, int xepoch, double xtime) // ctor
      : ts(xts), epoch(xepoch), time(xtime)
   {
   }
};

class TsSyncModule
{
public:
   double   fFreqHz;
   int      fEpoch;
   uint32_t fFirstTs;
   uint32_t fPrevTs;
   uint32_t fLastTs;
   double   fOffsetSec;
   double   fPrevTimeSec;
   double   fLastTimeSec;
   double   fEps;
   int      fSyncedWith;

   std::vector<TsSyncEntry> fBuf;

public:
   TsSyncModule() // ctor
   {
      fFreqHz   = 0;
      fEpoch    = 0;
      fFirstTs  = 0;
      fPrevTs   = 0;
      fLastTs   = 0;
      fOffsetSec   = 0;
      fPrevTimeSec = 0;
      fLastTimeSec = 0;
      fEps = 200*1e-9; // 200ns
      fSyncedWith  = -1;
   }

   void Print() const
   {
      printf("ts 0x%08x, prev 0x%08x, first 0x%08x, offset %f, time %f %f, diff %f, buf %d", fLastTs, fPrevTs, fFirstTs, fOffsetSec, fLastTimeSec, fPrevTimeSec, fLastTimeSec-fPrevTimeSec, (int)fBuf.size());
   }

   double GetTime(uint32_t ts, int epoch) const
   {
      // must do all computations in double precision to avoid trouble with 32-bit timestamp wraparound.
      return ts/fFreqHz - fFirstTs/fFreqHz + fOffsetSec + epoch*2.0*0x80000000/fFreqHz;
   }

   void Add(uint32_t ts)
   {
      if (fFirstTs == 0) {
         fFirstTs = ts;
      }
      fPrevTs = fLastTs;
      fPrevTimeSec = fLastTimeSec;
      fLastTs = ts;

      if (ts < fPrevTs) {
         fEpoch += 1;
      }

      // must do all computations in double precision to avoid trouble with 32-bit timestamp wraparound.
      fLastTimeSec = GetTime(fLastTs, fEpoch); // fLastTs/fFreqHz + fOffsetSec + fEpoch*2.0*0x80000000/fFreqHz;

      if (fBuf.size() < 100) {
         fBuf.push_back(TsSyncEntry(fLastTs, fEpoch, fLastTimeSec));
      }
   }

   void Retime()
   {
      for (unsigned i=0; i<fBuf.size(); i++) {
         fBuf[i].time = GetTime(fBuf[i].ts, fBuf[i].epoch);
      }
   }

   double GetDt(unsigned j)
   {
      //assert(j>0);
      //assert(j<fBuf.size());
      return fBuf[j].time - fBuf[j-1].time;
   }

   unsigned FindDt(double dt)
   {
      for (unsigned j=fBuf.size()-1; j>1; j--) {
         double jdt = GetDt(j);
         //printf("size %d, buf %d, jdt %f, dt %f\n", (int)fBuf.size(), j, jdt, dt);
         if (fabs(dt - jdt) < fEps) {
            //printf("found %d %f\n", j, jdt);
            return j;
         }
      }
      //printf("not found!\n");
      return 0;
   }
};

class TsSync
{
public:
   std::vector<TsSyncModule> fModules;
   bool fSyncOk;

public:
   TsSync() // ctor
   {
      fSyncOk = false;
   }

   ~TsSync() // dtor
   {
   }

   void Configure(unsigned i, double freq_hz)
   {
      // grow the array if needed
      TsSyncModule m;
      for (unsigned j=fModules.size(); j<=i; j++)
         fModules.push_back(m);

      fModules[i].fFreqHz = freq_hz;
   }

   void CheckSync(unsigned ii, unsigned i)
   {
      unsigned ntry = 3;
      unsigned jj = fModules[ii].fBuf.size()-1;

      double tt = fModules[ii].GetDt(jj);

      unsigned j = fModules[i].FindDt(tt);
      if (j == 0)
         return;

      // demand a few more good matches
      for (unsigned itry=1; itry<=ntry; itry++) {
         if (jj-itry == 0)
            return;
         if (j-itry == 0)
            return;
         double xtt = fModules[ii].GetDt(jj-itry);
         double xt  = fModules[i].GetDt(j-itry);

         if (fabs(xt - xtt) > fModules[ii].fEps) {
            return;
         }
      }

      fModules[ii].fSyncedWith = i;

      // check for sync loop
      if (fModules[i].fSyncedWith >= 0)
         fModules[ii].fSyncedWith = fModules[i].fSyncedWith;

      double off = fModules[i].fBuf[j].time - fModules[ii].fBuf[jj].time;
      fModules[ii].fOffsetSec = off;
      fModules[ii].Retime();
   }

   void Check(unsigned inew)
   {
      unsigned min = 0;
      unsigned max = 0;

      for (unsigned i=0; i<fModules.size(); i++) {
         unsigned s = fModules[i].fBuf.size();
         if (s > 0) {
            if (min == 0)
               min = s;
            if (s < min)
               min = s;
            if (s > max)
               max = s;
         }
      }

      //printf("min %d, max %d\n", min, max);

      if (min < 3)
         return;

      for (unsigned i=0; i<fModules.size(); i++) {
         if (inew != i && fModules[inew].fSyncedWith < 0) {
            CheckSync(inew, i);
         }
      }

      for (unsigned j=1; j<max; j++) {
         printf("buf %2d: ", j);
         for (unsigned i=0; i<fModules.size(); i++) {
            if (j<fModules[i].fBuf.size()) {
               double dt = fModules[i].fBuf[j].time - fModules[i].fBuf[j-1].time;
               printf(" %f", dt);
            } else {
               printf(" -");
            }
         }
         printf("\n");
      }

      printf("buf %2d: ", 99);
      for (unsigned i=0; i<fModules.size(); i++) {
         printf(" %d", fModules[i].fSyncedWith);
      }
      printf("\n");

      bool ok = true;
      for (unsigned i=0; i<fModules.size(); i++) {
         if (fModules[i].fSyncedWith < 0) {
            // FIXME: if there is no data from this module, mark sync ok.
            ok = false;
         }
      }
      fSyncOk = ok;
   }

   void Add(unsigned i, uint32_t ts)
   {
      fModules[i].Add(ts);

      if (0) {
         printf("module %2d: ", i);
         fModules[i].Print();
         printf("\n");
      }

      if (!fSyncOk) {
         Check(i);
      }
   }
};

//static const double TSNS = 16.0;
static const double TSNS = 8.0/.99999821102751183809;

struct FeamData
{
   FeamAdcData *adc;

   ~FeamData() // dtor
   {
      assert(adc);
      delete adc;
   }
};

struct FeamEvent
{
   bool   complete;
   bool   error;
   int    counter;
   double time;
   std::vector<FeamData*> modules;

   FeamEvent(); // ctor
   ~FeamEvent(); // dtor
   void Print() const;
};

FeamEvent::FeamEvent() // ctor
{
   complete = false;
   error = false;
   counter = 0;
   time = 0;
}

FeamEvent::~FeamEvent() // dtor
{
   for (unsigned i=0; i<modules.size(); i++)
      if (modules[i]) {
         delete modules[i];
         modules[i] = NULL;
      }
}

void FeamEvent::Print() const
{
   printf("FeamEvent %d, time %f, complete %d, error %d, modules: ", counter, time, complete, error);
   for (unsigned i=0; i<modules.size(); i++) {
      if (modules[i] == NULL) {
            printf(" null");
      } else {
         printf(" %d", modules[i]->adc->cnt);
      }
   }
}

class FeamEVB
{
public:
   unsigned fNumModules;
   std::vector<FeamAdcData*> fData;

   std::deque<FeamAdcData*> fBuf;

   TsSync fSync;

   std::deque<FeamEvent*> fEvents;

   int fCounter;

   FeamEVB(int num_modules)
   {
      fNumModules = num_modules;
      fCounter = 0;
      for (unsigned i=0; i<fNumModules; i++) {
         fData.push_back(NULL);
         fSync.Configure(i, 1.0/TSNS*1e9);
      }
   }

   FeamEvent* FindEvent(double t)
   {
      for (unsigned i=0; i<fEvents.size(); i++) {
         if (fabs(fEvents[i]->time - t) < 200.0/1e9) {
            //printf("Found event for time %f\n", t);
            return fEvents[i];
         }
      }

      FeamEvent* e = new FeamEvent();
      e->complete = false;
      e->error = false;
      e->counter = fCounter++;
      e->time = t;

      for (unsigned i=0; i<fNumModules; i++)
         e->modules.push_back(NULL);

      fEvents.push_back(e);

      //printf("New event for time %f\n", t);

      return e;
   }

   void CheckFeam(FeamEvent *e)
   {
      bool c = true;
      for (unsigned i=0; i<e->modules.size(); i++) {
         if (e->modules[i] == NULL) {
            c = false;
         } else {
            if (e->modules[i]->adc->error)
               e->error = true;
         }
      }
      e->complete = c;

      //PrintFeam(e);
   }
   
   void AddFeam(int ifeam, FeamAdcData *a)
   {
      FeamEvent* e = FindEvent(a->fTime);

      if (e->modules[ifeam]) {
         // FIXME: duplicate data
         printf("duplicate data!\n");
         delete a;
         return;
      }

      e->modules[ifeam] = new FeamData;
      e->modules[ifeam]->adc = a;

      CheckFeam(e);
   }

   void Build()
   {
      while (fBuf.size() > 0) {
         FeamAdcData* a = fBuf.front();
         fBuf.pop_front();
         AddFeam(a->module, a);
      }
   }

   void AddPacket(int ifeam, const FeamPacket* p, const char* ptr, int size)
   {
      if (p->n == 0) {
         // 1st packet

         if (fData[ifeam]) {
            //printf("Complete event: FEAM %d: ", ifeam);
            //data[ifeam]->Print();
            //printf("\n");

            FeamAdcData* a = fData[ifeam];
            fData[ifeam] = NULL;

            fSync.Add(ifeam, a->ts_trig);

            a->fTime = fSync.fModules[ifeam].fLastTimeSec;

            a->fTsAbsNs = fSync.fModules[ifeam].fLastTimeSec*1e9; // a->ts_trig*TSNS - xts[ifeam]*TSNS + fTsEpoch[ifeam]*(2.0*TSNS*0x80000000);
            a->fTsIncrNs = (fSync.fModules[ifeam].fLastTimeSec - fSync.fModules[ifeam].fPrevTimeSec)*1e9; // a->ts_trig*TSNS - xts[ifeam]*TSNS + fTsEpoch[ifeam]*(2.0*TSNS*0x80000000);

            a->Finalize();

            // xxx

            fBuf.push_back(a);
         }

         //printf("Start ew event: FEAM %d: ", ifeam);
         //p->Print();
         //printf("\n");

         fData[ifeam] = new FeamAdcData(p, ifeam);
      }

      FeamAdcData* a = fData[ifeam];

      if (a == NULL) {
         // did not see the first event yet, cannot unpack
         delete p;
         return;
      }

      a->AddData(p, ifeam, ptr, size);

      //a->Print();
      //printf("\n");

      delete p;
   }

   void Print() const
   {
      //printf("FEAM evb status: %p %p, buffered %d\n", data[0], data[1], (int)buf.size());
   }

   FeamEvent* Get()
   {
      if (fBuf.size() < 1)
         return NULL;

      if (fSync.fSyncOk)
         Build();

      if (fEvents.size() < 1)
         return NULL;

      // check if the oldest event is complete
      if (!fEvents.front()->complete) {
         // oldest event is incomplete,
         // check if any newer events are completed,
         // if they are, pop this incomplete event
         bool c = false;
         for (unsigned i=0; i<fEvents.size(); i++) {
            if (fEvents[i]->complete) {
               c = true;
               break;
            }
         }
         // if there are too many buffered events, all incomplete,
         // something is wrong, push them out anyway
         if (!c && fEvents.size() < 10)
            return NULL;

         printf("popping in incomplete event! have %d buffered events\n", (int)fEvents.size());
      }

      FeamEvent* e = fEvents.front();
      fEvents.pop_front();
      return e;
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

#define MAX_FEAM   8
#define MAX_SCA    4
#define MAX_CHAN  80
#define MAX_BINS 900

static int x3count = 0;

struct AdcData
{
   int nfeam;
   int nsca;
   int nchan;
   int nbins;

   int adc[MAX_FEAM][MAX_SCA][MAX_CHAN][MAX_BINS];

   AdcData() // ctor
   {
      printf("AdcData::ctor, count %d\n", x3count++);

      nfeam = 0;
      nsca  = 0;
      nchan = 0;
      nbins = 0;
   }

   ~AdcData() // dtor
   {
      printf("AdcData::dtor!\n"); x3count--;
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

   void CreateHistograms(TARunInfo* runinfo, int nchan, int nbins)
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

      // FIXME: who deletes this?
      hbmean = new TH1D*[nchan];
      hbrms  = new TH1D*[nchan];
      hwaveform = new TH1D*[nchan];
      hamp = new TH1D*[nchan];
      hled = new TH1D*[nchan];
   
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
      
      if (event->event_id != 1 && event->event_id != 2)
         return flow;

      int force_plot = false;

      AdcData *aaa = NULL;

#if 0
      int adc[80][5120];

      const int xbins = 829;
      const int xchan = 79;
#endif

      //printf("event id %d\n", event->event_id);

      if (event->event_id == 1) {
         //*flags |= TAFlag_SKIP; // enable this to skip GRIF16 events

         const char* banks[] = { "BB01", "BB02", "BB03", "BB04", "BB05", "BB06", "BB07", "BB08", NULL };
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
                     xevb = new FeamEVB(MAX_FEAM);
                  
                  xevb->AddPacket(i, p, data + p->off, p->buf_len);
               }
            }
         }

         if (!xevb) {
            return flow;
         }

         //xevb->Print();

         FeamEvent *e = xevb->Get();

         if (!e) {
            return flow;
         }

         if (1) {
            printf("ZZZ Processing FEAM event: ");
            e->Print();
            printf("\n");
         }

         if (1) {
            for (unsigned i=0; i<e->modules.size(); i++) {
               if (!e->modules[i])
                  break;
               FeamAdcData* a = e->modules[i]->adc;
               printf("module %2d, cnt %4d, ts_trig: 0x%08x %14.3f usec, ts_incr %14.3f usec\n", a->module, a->cnt, a->ts_trig, a->fTsAbsNs/1e3, a->fTsIncrNs/1e3);
               //a->Print();
               //printf("\n");
            }
         }

         if (e->error) {
            delete e;
            return flow;
         }

         if (1) {
            delete e;
            return flow;
         }

         //assert(a->next_n == 256);
         //assert(a->fSize == 310688);

         if (!aaa)
            aaa = new AdcData;

         //MEMZERO(adc);
         MEMZERO(aaa->adc);

         aaa->nfeam = e->modules.size();
         aaa->nsca  = 4;
         aaa->nchan = 76;
         aaa->nbins = 511;

         for (unsigned ifeam=0; ifeam<e->modules.size(); ifeam++) {
            FeamAdcData* a = e->modules[ifeam]->adc;
            const unsigned char* ptr = (const unsigned char*)a->fPtr;
            int count = 0;

            for (int ibin = 0; ibin < 511; ibin++) {
               for (int ichan = 0; ichan < 76; ichan++) {
                  for (int isca = 0; isca < 4; isca++) {
                     unsigned v = ptr[0] | ((ptr[1])<<8);
                     // manual sign extension
                     if (v & 0x8000)
                        v |= 0xffff0000;
                     //if (isca == 0) {
                     //   adc[ichan][ibin] = v;
                     //}
                     aaa->adc[ifeam][isca][ichan][ibin] = v;
                     ptr += 2;
                     count += 2;
                  }
               }
            }
            
            printf("count %d\n", count);
         }

         //for (int ibin = 511; ibin < xbins; ibin++) {
         //   for (int ichan = 0; ichan < 76; ichan++) {
         //      adc[ichan][ibin] = adc[ichan][ibin-511];
         //   }
         //}

         delete e;

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
         
      // got all the data here

      int nchan = aaa->nfeam * aaa->nsca * aaa->nchan;

      printf("nchan %d\n", nchan);

      Waveform** ww = new Waveform*[nchan];

      for (int ifeam=0; ifeam<aaa->nfeam; ifeam++) {
         for (int isca=0; isca<aaa->nsca; isca++) {
            for (int ichan=0; ichan<aaa->nchan; ichan++) {
               int xchan = ifeam*(aaa->nsca*aaa->nchan) + isca*aaa->nchan + ichan;
               ww[xchan] = new Waveform(aaa->nbins);
               for (int ibin=0; ibin<aaa->nbins; ibin++) {
                  ww[xchan]->samples[ibin] = (aaa->adc[ifeam][isca][ichan][ibin])/4;
               }
            }
         }
      }

      // create histograms

      CreateHistograms(runinfo, nchan, aaa->nbins);

      int iplot = 0;
      double zmax = 0;

      for (int ichan=0; ichan<nchan; ichan++) {
         double r;
         double b = baseline(ww[ichan], 10, 60, NULL, &r);
         
         if (b==0 && r==0)
            continue;

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

         if (hwaveform[ichan]->GetEntries() == 0) {
            printf("saving waveform %d\n", ichan);
            for (int i=0; i<aaa->nbins; i++)
               hwaveform[ichan]->SetBinContent(i+1, ww[ichan]->samples[i]);
         }

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
            TH1D* hh = new TH1D("hh", "hh", aaa->nbins, 0, aaa->nbins);
            for (int ibin=0; ibin<aaa->nbins; ibin++) {
               hh->SetBinContent(ibin+1, aaa->adc[0][0][0][ibin]);
            }
            hh->Draw();
         }
         
         if (1) {
            fC->cd(2);
            TH1D* hhh = new TH1D("hhh", "hhh", aaa->nbins, 0, aaa->nbins);
            for (int ibin=0; ibin<aaa->nbins; ibin++) {
               hhh->SetBinContent(ibin+1, aaa->adc[0][0][0][ibin]);
            }
            hhh->SetMinimum(-33000);
            hhh->SetMaximum(+33000);
            hhh->Draw();
         }
         
         if (1) {
            fC->cd(3);
            TH1D* hhh = new TH1D("hhhh", "hhhh", aaa->nbins, 0, aaa->nbins);
            for (int ibin=0; ibin<aaa->nbins; ibin++) {
               hhh->SetBinContent(ibin+1, ww[39]->samples[ibin]);
            }
            hhh->SetMinimum(-9000);
            hhh->SetMaximum(+9000);
            hhh->Draw();
         }
         
         if (1) {
            fC->cd(4);
            TH1D* hhh = new TH1D("hhhhh", "hhhhh", aaa->nbins, 0, aaa->nbins);
            for (int ibin=0; ibin<aaa->nbins; ibin++) {
               hhh->SetBinContent(ibin+1, ww[iplot]->samples[ibin]);
            }
            hhh->SetMinimum(-9000);
            hhh->SetMaximum(+9000);
            hhh->Draw();
         }
         
         if (1) {
            fC->cd(5);
            TH1D* h33 = new TH1D("h33", "h33", aaa->nbins, 0, aaa->nbins);
            for (int ibin=0; ibin<aaa->nbins; ibin++) {
               h33->SetBinContent(ibin+1, ww[33]->samples[ibin]);
            }
            h33->SetMinimum(-9000);
            h33->SetMaximum(+9000);
            h33->Draw();
         }
         
         if (1) {
            fC->cd(6);
            TH1D* h34 = new TH1D("h34", "h34", aaa->nbins, 0, aaa->nbins);
            for (int ibin=0; ibin<aaa->nbins; ibin++) {
               h34->SetBinContent(ibin+1, ww[531]->samples[ibin]);
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

         TH1D* h = new TH1D("h", "h", aaa->nbins, 0, aaa->nbins);
         for (int ibin=0; ibin<aaa->nbins; ibin++) {
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
