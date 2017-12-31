//
// Unpacking FEAM data
// K.Olchanski
//

#include "FeamEVB.h"

#include <stdio.h> // NULL, printf()
#include <math.h> // fabs()
#include <assert.h> // assert()

FeamEVB::FeamEVB(int num_modules, double ts_freq, double eps_sec)
{
   fNumModules = num_modules;
   fEpsSec  = eps_sec;
   fSync.SetDeadMin(10);
   for (unsigned i=0; i<fNumModules; i++) {
      fAsm.push_back(new FeamAsm);
      fSync.Configure(i, ts_freq, eps_sec, 0, 50);
   }
}

FeamEVB::~FeamEVB()
{
   for (unsigned i=0; i<fAsm.size(); i++) {
      if (fAsm[i]) {
         delete fAsm[i];
         fAsm[i] = NULL;
      }
   }

   for (unsigned i=0; i<fBuf.size(); i++) {
      if (fBuf[i]) {
         delete fBuf[i];
         fBuf[i] = NULL;
      }
   }

   for (unsigned i=0; i<fEvents.size(); i++) {
      if (fEvents[i]) {
         delete fEvents[i];
         fEvents[i] = NULL;
      }
   }
}

FeamEvent* FeamEVB::FindEvent(double t)
{
   double amin = 0;
   for (unsigned i=0; i<fEvents.size(); i++) {
      double dt = fEvents[i]->time - t;
      double adt = fabs(dt);
      if (adt < fEpsSec) {
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
   e->counter = (++fCounter);
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
      printf("\n");
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
   for (unsigned i=0; i<fAsm.size(); i++) {
      fAsm[i]->Finalize();
      Flush(i);
   }

   Build(true);
}

void FeamEVB::Build(bool force_build)
{
   bool ok_to_build = fSync.fSyncOk || fSync.fSyncFailed || force_build;

   if (!ok_to_build)
      return;

   while (fBuf.size() > 0) {
      FeamModuleData* m = fBuf.front();
      fBuf.pop_front();
      AddFeam(m->fPosition, m);
   }
}

void FeamEVB::Finalize(int position, FeamModuleData* m)
{
   fSync.Add(position, m->ts_trig);
   
   m->fTs = fSync.fModules[position].fLastTs;
   m->fTsEpoch = fSync.fModules[position].fEpoch;
   m->fTime = 0; // fSync.fModules[position].fLastTimeSec;
   m->fTimeIncr = fSync.fModules[position].fLastTimeSec - fSync.fModules[position].fPrevTimeSec;
   
   m->Finalize();
   
   fBuf.push_back(m);
}

void FeamEVB::Flush(int position)
{
   bool flag = false;
   while (fAsm[position]->fBuffer.size() > 0) {
      FeamModuleData* m = fAsm[position]->fBuffer.front();
      fAsm[position]->fBuffer.pop_front();
      Finalize(position, m);
      flag = true;
   }

   if (flag) {
      Build();
   }
}

void FeamEVB::AddPacket(const char* bank, int position, int format, const FeamPacket* p, const char* ptr, int size)
{
   fAsm[position]->AddPacket(p, bank, position, format, ptr, size);

   Flush(position);

   delete p;
}

void FeamEVB::Print() const
{
   printf("FeamEVB::Print: FEAM event builder status:\n");

   printf("  Sync status: ");
   fSync.Print();
   printf("\n");

   printf("  event assembler time threshold: %.0f ns\n", fEpsSec*1e9);
   printf("  event counter: %d\n", fCounter);
   printf("  max dt: %.0f ns (time between modules in one event)\n", fMaxDt*1e9);
   printf("  min dt: %.0f ns (time between events)\n", fMinDt*1e9);

   printf("  complete events: %d\n", fCountComplete);
   printf("  incomplete events: %d\n", fCountIncomplete);
   printf("  events with error: %d\n", fCountError);
   printf("  duplicate data error: %d\n", fCountDuplicate);
   
   printf("  Assembler: %d entries\n", (int)fAsm.size());
   for (unsigned i=0; i<fAsm.size(); i++) {
      printf("    position %d: ", i);
      fAsm[i]->Print();
      printf("\n");
   }

   printf("  Incomplete events buffer: %d entries\n", (int)fBuf.size());
   for (unsigned i=0; i<fBuf.size(); i++) {
      printf("    entry %d: ", i);
      fBuf[i]->Print();
      printf("\n");
   }

   printf("  Completed events buffer: %d entries\n", (int)fEvents.size());
   for (unsigned i=0; i<fEvents.size(); i++) {
      printf("    entry %d: ", i);
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
      int c = -1;
      for (unsigned i=0; i<fEvents.size(); i++) {
         if (fEvents[i]->complete) {
            c = i;
            break;
         }
      }
      // if there are too many buffered events, all incomplete,
      // something is wrong, push them out anyway
      if ((c < 0) && (fEvents.size() < 20))
         return NULL;
      
      printf("FeamEVB: popping an incomplete event! have %d buffered events, have complete at %d\n", (int)fEvents.size(), c);
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

   if (e->error)
      fCountError++;

   if (!e->complete)
      fCountIncomplete++;
   else
      fCountComplete++;

   return e;
}

/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */


