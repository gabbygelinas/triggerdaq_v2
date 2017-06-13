//
// AgEvent.cxx
//
// Unpacking ALPHA-g event
// K.Olchanski
//

#include "AgEvent.h"

#include <stdio.h> // NULL

AgEvent::AgEvent() // ctor
{
   complete = false;
   error = false;
   counter = 0;
   time = 0;
   timeIncr = 0;
   a16 = NULL;
   feam = NULL;
};

AgEvent::~AgEvent() // dtor
{
   if (a16) {
      delete a16;
      a16 = NULL;
   }

   if (feam) {
      delete feam;
      feam = NULL;
   }
}

void AgEvent::Print() const
{
   printf("AgEvent: complete %d, error %d, counter %d, time %f, incr %f", complete, error, counter, time, timeIncr);
   printf("\n");
   printf("  ");
   if (a16) {
      a16->Print();
   } else {
      printf("A16Event:  NULL");
   }
   printf("\n");
   printf("  ");
   if (feam) {
      feam->Print();
   } else {
      printf("FeamEvent: NULL");
   }
}

/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */


