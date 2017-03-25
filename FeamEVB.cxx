//
// Unpacking FEAM data
// K.Olchanski
//

#include "FeamEVB.h"

#include <stdio.h> // NULL, printf()
#include <math.h> // fabs()
#include <assert.h> // assert()

FeamEVB::FeamEVB(int num_modules, double ts_freq)
{
   fNumModules = num_modules;
   fCounter = 0;
   fSync.SetDeadMin(10);
   for (unsigned i=0; i<fNumModules; i++) {
      fData.push_back(NULL);
      fSync.Configure(i, ts_freq, 1000.0*1e-9, 0, 50);
   }
}

FeamEvent* FeamEVB::FindEvent(double t)
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
   
   for (unsigned i=0; i<fNumModules; i++) {
      e->modules.push_back(NULL);
      e->adcs.push_back(NULL);
   }
   
   fEvents.push_back(e);
   
   //printf("New event for time %f\n", t);
   
   return e;
}

void FeamEVB::CheckFeam(FeamEvent *e)
{
   bool c = true;
   
   for (unsigned i=0; i<e->modules.size(); i++) {
      if (!fSync.fModules[i].fDead) {
         if (e->modules[i] == NULL) {
            c = false;
         } else {
            if (e->modules[i]->error)
               e->error = true;
         }
      }
   }
   
   if (e->adcs.size() != e->modules.size()) {
      e->error = true;
   }
   
   for (unsigned i=0; i<e->adcs.size(); i++) {
      if (!fSync.fModules[i].fDead) {
         if (e->adcs[i] == NULL) {
            c = false;
         }
      }
   }
   
   e->complete = c;
   
   //PrintFeam(e);
}
   
void FeamEVB::AddFeam(int ifeam, FeamModuleData *m)
{
   m->fTime = fSync.fModules[ifeam].GetTime(m->fTs, m->fTsEpoch);

   //printf("FeamEVB::AddFeam: module %d, ts 0x%08x, epoch %d, time %f\n", ifeam, m->fTs, m->fTsEpoch, m->fTime);

   FeamEvent* e = FindEvent(m->fTime);
   
   if (e->modules[ifeam]) {
      // FIXME: duplicate data
      printf("duplicate data!\n");
      delete m;
      return;
   }

   if (e->timeIncr == 0)
      e->timeIncr = m->fTimeIncr;
   // use smaller time, in case of missing packets or events
   if (m->fTimeIncr < e->timeIncr)
      e->timeIncr = m->fTimeIncr;
   
   e->modules[ifeam] = m;
   e->adcs[ifeam] = new FeamAdcData;
   
   Unpack(e->adcs[ifeam], m);
   
   CheckFeam(e);
}

void FeamEVB::Build()
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
         
         m->fTs = fSync.fModules[ifeam].fLastTs;
         m->fTsEpoch = fSync.fModules[ifeam].fEpoch;
         m->fTime = 0; // fSync.fModules[ifeam].fLastTimeSec;
         m->fTimeIncr = fSync.fModules[ifeam].fLastTimeSec - fSync.fModules[ifeam].fPrevTimeSec;
         
         m->Finalize();
         
         // xxx
         
         fBuf.push_back(m);
      } else {
         printf("FeamEVB: Received first data from FEAM %d\n", ifeam);
      }
      
      //printf("Start ew event: FEAM %d: ", ifeam);
      //p->Print();
      //printf("\n");

      assert(fData[ifeam] == NULL);
      
      fData[ifeam] = new FeamModuleData(p, ifeam);
   }
   
   FeamModuleData* m = fData[ifeam];
   
   if (m == NULL) {
      // did not see the first event yet, cannot unpack
      printf("FeamEVB: dropped packet, feam %d: ", ifeam);
      p->Print();
      printf("\n");
      delete p;
      return;
   }
   
   m->AddData(p, ifeam, ptr, size);
   
   //a->Print();
   //printf("\n");
   
   delete p;
}

void FeamEVB::Print() const
{
   //printf("FEAM evb status: %p %p, buffered %d\n", data[0], data[1], (int)buf.size());
}

FeamEvent* FeamEVB::Get()
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
      
      printf("FeamEVB: popping in incomplete event! have %d buffered events, have complete %d\n", (int)fEvents.size(), c);

      if (c == 0) {
         printf("First incomplete event: ");
         fEvents.front()->Print();
         printf("\n");
      }
   }
   
   FeamEvent* e = fEvents.front();
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


