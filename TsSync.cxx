//
// TsSync.cxx
// Timestamp synchronization
// K.Olchanski
//

#include <stdio.h>
#include <math.h>
#include <assert.h> // assert()

#include "TsSync.h"

TsSyncEntry::TsSyncEntry(uint32_t xts, int xepoch, double xtime) // ctor
      : ts(xts), epoch(xepoch), time(xtime)
{
}

TsSyncModule::TsSyncModule() // ctor
{
   fFreqHz   = 0;
   fEpoch    = 0;
   fFirstTs  = 0;
   fPrevTs   = 0;
   fLastTs   = 0;
   fOffsetSec   = 0;
   fPrevTimeSec = 0;
   fLastTimeSec = 0;
   fEps = 2000*1e-9; // in ns
   fSyncedWith  = -1;
   fOverflow = false;
   fDead = false;
   fBufMax = 100;
}

void TsSyncModule::Print() const
{
   printf("ts 0x%08x, prev 0x%08x, first 0x%08x, offset %f, time %f %f, diff %f, buf %d", fLastTs, fPrevTs, fFirstTs, fOffsetSec, fLastTimeSec, fPrevTimeSec, fLastTimeSec-fPrevTimeSec, (int)fBuf.size());
}

double TsSyncModule::GetTime(uint32_t ts, int epoch) const
{
   // must do all computations in double precision to avoid trouble with 32-bit timestamp wraparound.
   return ts/fFreqHz - fFirstTs/fFreqHz + fOffsetSec + epoch*2.0*0x80000000/fFreqHz;
}

void TsSyncModule::Add(uint32_t ts)
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
   
   if (fBuf.size() > fBufMax) {
      fOverflow = true;
      return;
   }

   fBuf.push_back(TsSyncEntry(fLastTs, fEpoch, fLastTimeSec));
}

void TsSyncModule::Retime()
{
   for (unsigned i=0; i<fBuf.size(); i++) {
      fBuf[i].time = GetTime(fBuf[i].ts, fBuf[i].epoch);
   }
}

double TsSyncModule::GetDt(unsigned j)
{
   assert(j>0);
   assert(j<fBuf.size());
   return fBuf[j].time - fBuf[j-1].time;
}

unsigned TsSyncModule::FindDt(double dt)
{
   //printf("FindDt: fBuf.size %d\n", fBuf.size());
   assert(fBuf.size() > 0);
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

TsSync::TsSync() // ctor
{
   fSyncOk = false;
   fTrace = false;
   fOverflow = false;
   fDeadMin = 0;
}

TsSync::~TsSync() // dtor
{
}

void TsSync::SetDeadMin(int dead_min)
{
   fDeadMin = dead_min;
}

void TsSync::Configure(unsigned i, double freq_hz, int buf_max)
{
   // grow the array if needed
   TsSyncModule m;
   for (unsigned j=fModules.size(); j<=i; j++)
      fModules.push_back(m);
   
   fModules[i].fFreqHz = freq_hz;
   fModules[i].fBufMax = buf_max;
}

void TsSync::CheckSync(unsigned ii, unsigned i)
{
   unsigned ntry = 3;
   unsigned jj = fModules[ii].fBuf.size()-1;
   
   double tt = fModules[ii].GetDt(jj);
   
   unsigned j = fModules[i].FindDt(tt);

   //printf("TsSync::CheckSync: module %d buf %d, dt %f with module %d buf %d\n", ii, jj, tt, i, j);
   
   if (j == 0) {
      if (jj < 2) {
         return;
      }

      jj -= 1;
      tt = fModules[ii].GetDt(jj);
      j = fModules[i].FindDt(tt);

      //printf("TsSync::CheckSync: module %d buf %d, dt %f with module %d buf %d (2nd try)\n", ii, jj, tt, i, j);
      
      if (j == 0) {
         return;
      }
   }
   
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

   printf("TsSync: module %d buf %d synced with module %d buf %d, offset %f\n", ii, jj, i, j, off);

   Dump();
}

void TsSync::Check(unsigned inew)
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

   if (fTrace)
      printf("TsSync::Check: min %d, max %d\n", min, max);
   
   if (min < 3)
      return;
   
   for (unsigned i=0; i<fModules.size(); i++) {
      if (fModules[i].fBuf.size() < 1)
         continue;
      if (inew != i && fModules[inew].fSyncedWith < 0) {
         if (i==1) { // kludge: only sync with module 1
            CheckSync(inew, i);
         }
      }
   }

   if (fTrace) {
      Dump();
   }

   int modules_with_data = 0;
   int no_sync = 0;
   for (unsigned i=0; i<fModules.size(); i++) {
      if (fModules[i].fBuf.size() > 0) {
         modules_with_data++;
         if (fModules[i].fSyncedWith < 0) {
            no_sync += 1;
         }
      }
   }

   if (fTrace)
      printf("modules: %d, with data: %d, no_sync: %d\n", (int)fModules.size(), modules_with_data, no_sync);

   if (modules_with_data > 1 && no_sync <= 1) {
      // at least one module has data and
      // only one unsynced module (all other modules synced to it)
      fSyncOk = true;
   } else if (min > fDeadMin && fModules.size() == 2 && modules_with_data == 1 && no_sync == 1) {
      // total 2 modules, one of them has data, the other one is dead
      fSyncOk = true;
   }

   if (fSyncOk) {
      int count_dead = 0;
      for (unsigned i=0; i<fModules.size(); i++) {
         if ((fModules[i].fSyncedWith < 0) && (fModules[i].fBuf.size() == 0)) {
            fModules[i].fDead = true;
            count_dead += 1;
         }
      }

      printf("TsSync: synchronization completed, %d dead modules.\n", count_dead);
   }
}

void TsSync::Add(unsigned i, uint32_t ts)
{
   if (fOverflow)
      return;

   if (0 && fTrace) {
      printf("Add %d, ts 0x%08x\n", i, ts);
   }

   fModules[i].Add(ts);

   if (!fSyncOk && fModules[i].fOverflow) {
      fOverflow = true;
      printf("TsSync: module %d buffer overflow, synchronization failed.\n", i);
      Dump();
      return;
   }
   
   if (0 && fTrace) {
      printf("module %2d: ", i);
      fModules[i].Print();
      printf("\n");
   }
   
   if (!fSyncOk) {
      Check(i);
   }
}

void TsSync::Dump() const
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
}

/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */
