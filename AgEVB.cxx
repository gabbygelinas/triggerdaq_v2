//
// AgEVB.cxx
//
// ALPHA-g event builder
// K.Olchanski
//

#include "AgEVB.h"

#include <stdio.h> // NULL, printf()
#include <math.h> // fabs()
#include <assert.h> // assert()

AgEVB::AgEVB(double a16_ts_freq, double feam_ts_freq, double eps_sec, int max_skew, int max_dead, bool clock_drift)
{
   fMaxSkew = max_skew;
   fMaxDead = max_dead;
   fEpsSec = eps_sec;
   fClockDrift = clock_drift;

   fCounter = 0;
   fSync.SetDeadMin(fMaxDead);
   fSync.Configure(0, a16_ts_freq,  0, 10.0*1e-6, fMaxSkew);
   fSync.Configure(1, feam_ts_freq, 0, 10.0*1e-6, fMaxSkew);
   fLastA16Time = 0;
   fLastFeamTime = 0;
   fMaxDt = 0;
   fMinDt = 0;
}

AgEVB::~AgEVB()
{
   printf("AgEVB: max dt: %.0f ns, min dt: %.0f ns\n", fMaxDt*1e9, fMinDt*1e9);
}

AgEvent* AgEVB::FindEvent(double t)
{
   double amin = 0;
   for (unsigned i=0; i<fEvents.size(); i++) {
      //printf("find event for time %f: event %d, %f, diff %f\n", t, i, fEvents[i]->time, fEvents[i]->time - t);

      double dt = fEvents[i]->time - t;
      double adt = fabs(dt);

      if (adt < fEpsSec) {
         if (adt > fMaxDt) {
            //printf("AgEVB: for time %f found event at time %f, new max dt %.0f ns, old max dt %.0f ns\n", t, fEvents[i]->time, adt*1e9, fMaxDt*1e9);
            fMaxDt = adt;
         }
         //printf("Found event for time %f\n", t);
         //printf("AgEVB: Found event for time %f: event %d, %f, diff %f %.0f ns\n", t, i, fEvents[i]->time, dt, dt*1e9);
         return fEvents[i];
      }

      if (amin == 0)
         amin = adt;
      if (adt < amin)
         amin = adt;
   }
   
   if (fMinDt == 0)
      fMinDt = amin;

   if (amin < fMinDt)
      fMinDt = amin;
   
   AgEvent* e = new AgEvent();
   e->complete = false;
   e->error = false;
   e->counter = (++fCounter);
   e->time = t;
   e->timeIncr = t - fLastEventTime;

   fLastEventTime = e->time;

   fEvents.push_back(e);
   
   //printf("New event for time %f\n", t);
   
   return e;
}

void AgEVB::CheckEvent(AgEvent *e)
{
   e->complete = true;

   if (!e->a16 && !fSync.fModules[0].fDead)
      e->complete = false;

   if (!e->feam && !fSync.fModules[1].fDead)
      e->complete = false;

   e->error = false;

   if (e->a16 && e->a16->error)
      e->error = true;

   if (e->feam && e->feam->error)
      e->error = true;

   //e->Print();
}
   
void AgEVB::Build(int index, AgEvbBuf *m)
{
   m->time = fSync.fModules[index].GetTime(m->ts, m->epoch);

   AgEvent* e = FindEvent(m->time);

   if (0 && index == 1) {
      printf("offset: %f %f, index %d, ts 0x%08x, epoch %d, feam time %f\n", e->time, m->time, index, m->ts, m->epoch, m->feam->time);
   }

   if (fClockDrift) { // adjust offset for clock drift
      double off = e->time - m->time;
      //printf("offset: %f %f, diff %f, index %d\n", e->time, m->time, off, index);
      fSync.fModules[index].fOffsetSec += off/2.0;
   }

   if (m->a16) {
      if (e->a16) {
         // FIXME: duplicate data
         printf("AgEVB: A16 duplicate data!\n");
         // FIXME memory leak
         delete m;
         return;
      }

      e->a16 = m->a16;
      m->a16 = NULL;
   }

   if (m->feam) {
      if (e->feam) {
         // FIXME: duplicate data
         printf("AgEVB: FEAM duplicate data!\n");
         // FIXME memory leak
         delete m;
         return;
      }

      e->feam = m->feam;
      m->feam = NULL;
   }

   delete m;

   CheckEvent(e);
}

void AgEVB::Build()
{
   while (fBuf[0].size() > 0) {
      AgEvbBuf* m = fBuf[0].front();
      fBuf[0].pop_front();
      Build(0, m);
   }

   while (fBuf[1].size() > 0) {
      AgEvbBuf* m = fBuf[1].front();
      fBuf[1].pop_front();
      Build(1, m);
   }
}

void AgEVB::AddAlpha16Event(Alpha16Event* e)
{
   fCountA16++;

   if (fCountA16 > 1 && e->time <= fLastA16Time) {
      fCountRejectedA16++;
      printf("AgEVB::AddA16Event: Alpha16 event time did not increase: new time %f, last seen time %f\n", e->time, fLastA16Time);
      delete e;
      return;
   }

   if (e->error)
      fCountErrorA16++;

   if (e->complete)
      fCountCompleteA16++;

   fLastA16Time = e->time;
   uint32_t ts = e->time*fSync.fModules[0].fFreqHz;
   //printf("Alpha16Event: t %f, ts 0x%08x", e->eventTime, ts);
   //e->Print();
   //printf("\n");
   fSync.Add(0, ts);
   AgEvbBuf* m = new AgEvbBuf;
   m->a16 = e;
   m->feam = NULL;
   m->ts = fSync.fModules[0].fLastTs;
   m->epoch = fSync.fModules[0].fEpoch;
   m->time = 0;
   m->timeIncr = fSync.fModules[0].fLastTimeSec - fSync.fModules[0].fPrevTimeSec;
   fBuf[0].push_back(m);
}

void AgEVB::AddFeamEvent(FeamEvent* e)
{
   fCountFeam++;

   if (fCountFeam > 1 && e->time <= fLastFeamTime) {
      fCountRejectedFeam++;
      printf("AgEVB::AddFeamEvent: FEAM event time did not increase: new time %f, last seen time %f\n", e->time, fLastFeamTime);
      delete e;
      return;
   }

   if (e->error)
      fCountErrorFeam++;

   if (e->complete)
      fCountCompleteFeam++;

   fLastFeamTime = e->time;
   uint32_t ts = e->time*fSync.fModules[1].fFreqHz;
   //printf("FeamEvent: t %f, ts 0x%08x", e->time, ts);
   //printf("\n");
   fSync.Add(1, ts);
   AgEvbBuf* m = new AgEvbBuf;
   m->a16 = NULL;
   m->feam = e;
   m->ts = fSync.fModules[1].fLastTs;
   m->epoch = fSync.fModules[1].fEpoch;
   m->time = 0;
   m->timeIncr = fSync.fModules[1].fLastTimeSec - fSync.fModules[1].fPrevTimeSec;
   fBuf[1].push_back(m);
}

void AgEVB::Print() const
{
   printf("AgEVB status:\n");
   printf("  Sync: "); fSync.Print(); printf("\n");
   printf("  A16 events: in %d, rejected %d, complete %d, error %d\n", fCountA16, fCountRejectedA16, fCountCompleteA16, fCountErrorA16);
   printf("  Feam events: in %d, rejected %d, complete %d, error %d\n", fCountFeam, fCountRejectedFeam, fCountCompleteFeam, fCountErrorFeam);
   printf("  Buffered A16:  %d\n", (int)fBuf[0].size());
   printf("  Buffered FEAM: %d\n", (int)fBuf[1].size());
   printf("  Buffered output: %d\n", (int)fEvents.size());
   printf("  Output %d events: %d complete, %d with errors, %d incomplete (%d no A16, %d no FEAM)\n", fCount, fCountComplete, fCountError, fCountIncomplete, fCountIncompleteA16, fCountIncompleteFeam);
}

void AgEVB::UpdateCounters(const AgEvent* e)
{
   fCount++;
   if (e->error) {
      fCountError++;
   }

   if (e->complete) {
      fCountComplete++;
   } else {
      fCountIncomplete++;

      if (!e->a16 && e->feam) {
         fCountIncompleteA16++;
      } else if (e->a16 && !e->feam) {
         fCountIncompleteFeam++;
      } else if (!e->a16 && !e->feam) {
         assert(!"This cannot happen, both A16 and FEAM missing!");
      } else {
         assert(!"This cannot happen!");
      }
   }
}

AgEvent* AgEVB::Get()
{
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
      if (!c && fEvents.size() < fMaxSkew)
         return NULL;
      
      printf("AgEVB:: popping an incomplete event! have %d buffered events, have complete %d\n", (int)fEvents.size(), c);
   }
   
   AgEvent* e = fEvents.front();
   fEvents.pop_front();
   UpdateCounters(e);
   return e;
}

AgEvent* AgEVB::GetLastEvent()
{
   Build();
   
   if (fEvents.size() < 1)
      return NULL;
   
   AgEvent* e = fEvents.front();
   fEvents.pop_front();
   UpdateCounters(e);
   return e;
}

/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */


