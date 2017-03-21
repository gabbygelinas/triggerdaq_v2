//
// AgEVB.cxx
//
// ALPHA-g event builder
// K.Olchanski
//

#include "AgEVB.h"

#include <stdio.h> // NULL, printf()
#include <math.h> // fabs()

AgEVB::AgEVB(double a16_ts_freq, double feam_ts_freq)
{
   fCounter = 0;
   fSync.Configure(0, a16_ts_freq, 100);
   fSync.Configure(1, feam_ts_freq, 100);
}

AgEvent* AgEVB::FindEvent(double t)
{
   for (unsigned i=0; i<fEvents.size(); i++) {
      //printf("find event for time %f: event %d, %f, diff %f\n", t, i, fEvents[i]->time, fEvents[i]->time - t);

      if (fabs(fEvents[i]->time - t) < 5000.0/1e9) {
         //printf("Found event for time %f\n", t);
         printf("find event for time %f: event %d, %f, diff %f\n", t, i, fEvents[i]->time, fEvents[i]->time - t);
         return fEvents[i];
      }
   }
   
   AgEvent* e = new AgEvent();
   e->complete = false;
   e->error = false;
   e->counter = fCounter++;
   e->time = t;
   
   fEvents.push_back(e);
   
   //printf("New event for time %f\n", t);
   
   return e;
}

void AgEVB::CheckEvent(AgEvent *e)
{
   e->complete = true;

   if (!e->a16)
      e->complete = false;

   if (!e->feam)
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

   if (1) { // adjust offset for clock drift
      double off = e->time - m->time;
      printf("offset: %f %f, diff %f\n", e->time, m->time, off);
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
   uint32_t ts = e->eventTime*fSync.fModules[0].fFreqHz/1e9;
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
   //printf("FEAM evb status: %p %p, buffered %d\n", data[0], data[1], (int)buf.size());
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
      if (!c && fEvents.size() < 200)
         return NULL;
      
      printf("AgEVB:: popping in incomplete event! have %d buffered events, have complete %d\n", (int)fEvents.size(), c);
   }
   
   AgEvent* e = fEvents.front();
   fEvents.pop_front();
   return e;
}

/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */


