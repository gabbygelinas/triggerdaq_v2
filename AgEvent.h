//
// AgEvent.h
//
// Unpacking ALPHA-g event
// K.Olchanski
//

#ifndef AgEvent_H
#define AgEvent_H

#include "Trig.h"
#include "Alpha16.h"
#include "Feam.h"

struct AgEvent
{
   bool   complete = false;
   bool   error    = false;
   int    counter  = 0;
   double time     = 0;
   double timeIncr = 0;

   TrigEvent*    trig = NULL;
   Alpha16Event* a16  = NULL;
   FeamEvent*    feam = NULL;

   AgEvent(); // ctor
   ~AgEvent(); // dtor
   void Print() const;
};

#endif

/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */


