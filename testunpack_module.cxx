//
// testunpack_module.cxx
//
// K.Olchanski
//

#include <stdio.h>
#include <assert.h>
#include <unistd.h>

#include "manalyzer.h"
#include "midasio.h"

#include "AgAsm.h"
#include "AgFlow.h"

#include "ncfm.h"

#ifdef _TIME_ANALYSIS_
#include "AnalysisTimer.h"
#endif

class TestUnpackModule: public TARunObject
{
public:
   AgAsm*       fAgAsm = NULL;

   TestUnpackModule(TARunInfo* runinfo)
      : TARunObject(runinfo)
   {
      printf("TestUnpackModule::ctor!\n");
      fModuleName = "testunpack";
   }

   ~TestUnpackModule()
   {
      printf("TestUnpackModule::dtor!\n");
      if (fAgAsm) {
         delete fAgAsm;
         fAgAsm = NULL;
      }
   }

   void BeginRun(TARunInfo* runinfo)
   {
      printf("TestUnpackModule::BeginRun, run %d, file %s\n", runinfo->fRunNo, runinfo->fFileName.c_str());
      fAgAsm = new AgAsm();
      fAgAsm->BeginRun(runinfo->fRunNo);
   }

   void PreEndRun(TARunInfo* runinfo)
   {
      printf("TestUnpackModule::PreEndRun, run %d\n", runinfo->fRunNo);
   }
   
   void EndRun(TARunInfo* runinfo)
   {
      printf("TestUnpackModule::EndRun, run %d\n", runinfo->fRunNo);
   }
   
   void PauseRun(TARunInfo* runinfo)
   {
      printf("TestUnpackModule::PauseRun, run %d\n", runinfo->fRunNo);
   }

   void ResumeRun(TARunInfo* runinfo)
   {
      printf("TestUnpackModule::ResumeRun, run %d\n", runinfo->fRunNo);
   }

   TAFlowEvent* Analyze(TARunInfo* runinfo, TMEvent* event, TAFlags* flags, TAFlowEvent* flow)
   {
      //printf("Analyze, run %d, event serno %d, id 0x%04x, data size %d\n", runinfo->fRunNo, event->serial_number, (int)event->event_id, event->data_size);

      if (event->event_id != 1)
         return flow;

      assert(fAgAsm != NULL);

      AgEvent* e = fAgAsm->UnpackEvent(event);

      printf("Unpacked AgEvent:   ");
      e->Print();
      printf("\n");

      if( (e->counter%1000)==0 )
         printf("Unpack Event %d\n",e->counter);

      delete e;
      e = NULL;

      return flow;
   }

   void AnalyzeSpecialEvent(TARunInfo* runinfo, TMEvent* event)
   {
      printf("TestUnpackModule::AnalyzeSpecialEvent, run %d, event serno %d, id 0x%04x, data size %d\n", runinfo->fRunNo, event->serial_number, (int)event->event_id, event->data_size);
   }
};

class TestUnpackModuleFactory: public TAFactory
{
public:
   void Usage()
   {
   }

   void Init(const std::vector<std::string> &args)
   {
      printf("TestUnpackModuleFactory::Init!\n");
   }

   void Finish()
   {
      printf("TestUnpackModuleFactory::Finish!\n");
   }
   
   TARunObject* NewRunObject(TARunInfo* runinfo)
   {
      printf("TestUnpackModuleFactory::NewRunObject, run %d, file %s\n", runinfo->fRunNo, runinfo->fFileName.c_str());
      return new TestUnpackModule(runinfo);
   }
};

static TARegister tar(new TestUnpackModuleFactory);

/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */
