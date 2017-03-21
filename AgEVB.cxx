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
   
void AgEVB::AddAlpha16(Alpha16Event *m)
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

void AgEVB::AddFeam(FeamEvent *m)
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
   while (fBuf.size() > 0) {
      FeamModuleData* m = fBuf.front();
      fBuf.pop_front();
      AddFeam(m->module, m);
   }
}

void FeamEVB::AddPacket(int ifeam, const FeamPacket* p, const char* ptr, int size)
{
   if (p->n == 0) {
      // 1st packet
      
      if (fData[ifeam]) {
         //printf("Complete event: FEAM %d: ", ifeam);
         //data[ifeam]->Print();
         //printf("\n");
         
         FeamModuleData* m = fData[ifeam];
         fData[ifeam] = NULL;
         
         fSync.Add(ifeam, m->ts_trig);
         
         m->fTime = fSync.fModules[ifeam].fLastTimeSec;
         m->fTimeIncr = fSync.fModules[ifeam].fLastTimeSec - fSync.fModules[ifeam].fPrevTimeSec;
         
         m->Finalize();
         
         // xxx
         
         fBuf.push_back(m);
      }
      
      //printf("Start ew event: FEAM %d: ", ifeam);
      //p->Print();
      //printf("\n");
      
      fData[ifeam] = new FeamModuleData(p, ifeam);
   }
   
   FeamModuleData* m = fData[ifeam];
   
   if (m == NULL) {
      // did not see the first event yet, cannot unpack
      printf("dropped packet!\n");
      delete p;
      return;
   }
   
   m->AddData(p, ifeam, ptr, size);
   
   //a->Print();
   //printf("\n");
   
   delete p;
}

void AgEVB::Print() const
{
   //printf("FEAM evb status: %p %p, buffered %d\n", data[0], data[1], (int)buf.size());
}

AgEvent* AgEVB::Get()
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


