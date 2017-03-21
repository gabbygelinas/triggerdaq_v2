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
   fSync.Configure(0, a16_ts_freq);
   fSync.Configure(1, feam_ts_freq);
}

AgEvent* AgEVB::FindEvent(double t)
{
   for (unsigned i=0; i<fEvents.size(); i++) {
      if (fabs(fEvents[i]->time - t) < 200.0/1e9) {
         //printf("Found event for time %f\n", t);
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
   
void AgEVB::BuildAlpha16(Alpha16Event *m)
{
   AgEvent* e = FindEvent(m->eventTime);
   
   if (e->a16) {
      // FIXME: duplicate data
      printf("duplicate data!\n");
      delete m;
      return;
   }
   
   e->a16 = m;

   CheckEvent(e);
}

void AgEVB::BuildFeam(FeamEvent *m)
{
   AgEvent* e = FindEvent(m->time);
   
   if (e->feam) {
      // FIXME: duplicate data
      printf("duplicate data!\n");
      delete m;
      return;
   }
   
   e->feam = m;

   CheckEvent(e);
}

void AgEVB::Build()
{
   while (fBuf0.size() > 0) {
      Alpha16Event* m = fBuf0.front();
      fBuf0.pop_front();
      BuildAlpha16(m);
   }

   while (fBuf1.size() > 0) {
      FeamEvent* m = fBuf1.front();
      fBuf1.pop_front();
      BuildFeam(m);
   }
}

void AgEVB::AddAlpha16Event(Alpha16Event* e)
{
   uint32_t ts = e->eventTime*fSync.fModules[0].fFreqHz/1e9;
   //printf("Alpha16Event: t %f, ts 0x%08x", e->eventTime, ts);
   //e->Print();
   //printf("\n");
   fSync.Add(0, ts);
   fBuf0.push_back(e);
}

void AgEVB::AddFeamEvent(FeamEvent* e)
{
   uint32_t ts = e->time*fSync.fModules[1].fFreqHz;
   //printf("FeamEvent: t %f, ts 0x%08x", e->time, ts);
   //printf("\n");
   fSync.Add(1, ts);
   fBuf1.push_back(e);
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
      if (!c && fEvents.size() < 10)
         return NULL;
      
      printf("AgEVB:: popping in incomplete event! have %d buffered events\n", (int)fEvents.size());
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


