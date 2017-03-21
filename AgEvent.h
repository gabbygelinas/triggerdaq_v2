//
// AgEvent.h
//
// Unpacking ALPHA-g event
// K.Olchanski
//

#ifndef AgEvent_H
#define AgEvent_H

#include "Alpha16.h"
#include "Feam.h"

struct AgEvent
{
   bool   complete;
   bool   error;
   int    counter;
   double time;

   Alpha16Event* a16;
   FeamEvent*    feam;

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


