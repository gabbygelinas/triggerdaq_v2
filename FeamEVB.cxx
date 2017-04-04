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
   fMaxDt = 0;
   fMinDt = 0;
   fCountComplete = 0;
   fCountIncomplete = 0;
   fCountDuplicate = 0;
   fCountError = 0;
   fCountDropped = 0;
}

FeamEVB::~FeamEVB()
{
}

FeamEvent* FeamEVB::FindEvent(double t)
{
   double amin = 0;
   for (unsigned i=0; i<fEvents.size(); i++) {
      double dt = fEvents[i]->time - t;
      double adt = fabs(dt);
      if (adt < 200.0/1e9) {
         if (adt > fMaxDt) {
            printf("FeamEVB: for time %f found event at time %f, new max dt %.0f ns, old max dt %.0f ns\n", t, fEvents[i]->time, adt*1e9, fMaxDt*1e9);
            fMaxDt = adt;
         }
         //printf("Found event for time %f\n", t);
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
   
void FeamEVB::AddFeam(int position, FeamModuleData *m)
{
   m->fTime = fSync.fModules[position].GetTime(m->fTs, m->fTsEpoch);

   //printf("FeamEVB::AddFeam: module %d, ts 0x%08x, epoch %d, time %f\n", position, m->fTs, m->fTsEpoch, m->fTime);

   FeamEvent* e = FindEvent(m->fTime);
   
   if (e->modules[position]) {
      fCountDuplicate++;
      printf("FeamEVB::AddFeam: Error: duplicate data for time %f position %d: ", m->fTime, position);
      m->Print();
      delete m;
      return;
   }

   if (e->timeIncr == 0)
      e->timeIncr = m->fTimeIncr;
   // use smaller time, in case of missing packets or events
   if (m->fTimeIncr < e->timeIncr)
      e->timeIncr = m->fTimeIncr;
   
   e->modules[position] = m;
   e->adcs[position] = new FeamAdcData;
   
   Unpack(e->adcs[position], m);
   
   CheckFeam(e);
}

void FeamEVB::BuildLastEvent()
{
   for (unsigned i=0; i<fData.size(); i++) {
      if (fData[i]) {
         // complete the buffered
         Finalize(i);
      }
   }

   Build();
}

void FeamEVB::Build()
{
   while (fBuf.size() > 0) {
      FeamModuleData* m = fBuf.front();
      fBuf.pop_front();
      AddFeam(m->fPosition, m);
   }
}

void FeamEVB::Finalize(int position)
{
   FeamModuleData* m = fData[position];
   fData[position] = NULL;

   assert(m != NULL); // ensured by caller
   
   fSync.Add(position, m->ts_trig);
   
   m->fTs = fSync.fModules[position].fLastTs;
   m->fTsEpoch = fSync.fModules[position].fEpoch;
   m->fTime = 0; // fSync.fModules[position].fLastTimeSec;
   m->fTimeIncr = fSync.fModules[position].fLastTimeSec - fSync.fModules[position].fPrevTimeSec;
   
   m->Finalize();
   
   fBuf.push_back(m);
}

void FeamEVB::AddPacket(const char* bank, int position, const FeamPacket* p, const char* ptr, int size)
{
   if (p->n == 0) {
      // 1st packet
      
      if (fData[position]) {
         // complete the previous event
         Finalize(position);
      } else {
         printf("FeamEVB: Received first data from FEAM %d\n", position);
      }
      
      //printf("Start ew event: FEAM %d: ", position);
      //p->Print();
      //printf("\n");

      assert(fData[position] == NULL);
      
      fData[position] = new FeamModuleData(p, bank, position);
   }
   
   FeamModuleData* m = fData[position];
   
   if (m == NULL) {
      // did not see the first event yet, cannot unpack
      fCountDropped++;
      printf("FeamEVB::AddPacket: Error: no 1st packet for position %d, dropping packet: ", position);
      p->Print();
      printf("\n");
      delete p;
      return;
   }
   
   m->AddData(p, position, ptr, size);
   
   //a->Print();
   //printf("\n");
   
   delete p;
}

void FeamEVB::Print() const
{
   printf("FeamEVB::Print: FEAM event builder status:\n");

   printf("Sync status: ");
   fSync.Print();
   printf("\n");

   printf("event counter: %d\n", fCounter);
   printf("max dt: %.0f ns (time between modules in one event)\n", fMaxDt*1e9);
   printf("min dt: %.0f ns (time between events)\n", fMinDt*1e9);

   printf("complete events: %d\n", fCountComplete);
   printf("incomplete events: %d\n", fCountIncomplete);
   printf("events with error: %d\n", fCountError);
   printf("duplicate data error: %d\n", fCountDuplicate);
   printf("dropped packet error: %d\n", fCountDropped);
   
   printf("Data buffer: %d entries\n", (int)fData.size());
   for (unsigned i=0; i<fData.size(); i++) {
      printf("position %d: ", i);
      if (fData[i]) {
         fData[i]->Print();
      } else {
         printf("null");
      }
      printf("\n");
   }

   printf("Incomplete events buffer: %d entries\n", (int)fBuf.size());
   for (unsigned i=0; i<fBuf.size(); i++) {
      printf("entry %d: ", i);
      fBuf[i]->Print();
      printf("\n");
   }

   printf("Completed events buffer: %d entries\n", (int)fEvents.size());
   for (unsigned i=0; i<fEvents.size(); i++) {
      printf("entry %d: ", i);
      fEvents[i]->Print();
      printf("\n");
   }
   printf("FeamEVB::Print: done.\n");
}

FeamEvent* FeamEVB::Get()
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
      
      printf("FeamEVB: popping in incomplete event! have %d buffered events, have complete %d\n", (int)fEvents.size(), c);

      if (c == 0) {
         printf("FeamEVB: First incomplete event: ");
         fEvents.front()->Print();
         printf("\n");
      }
   }
   
   FeamEvent* e = fEvents.front();
   fEvents.pop_front();

   if (e->error)
      fCountError++;

   if (!e->complete)
      fCountIncomplete++;
   else
      fCountComplete++;
   
   return e;
}

FeamEvent* FeamEVB::GetLastEvent()
{
   BuildLastEvent();
   
   if (fEvents.size() < 1)
      return NULL;
   
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


