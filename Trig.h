//
// Trigger data
// K.Olchanski
//

#ifndef Trig_H
#define Trig_H

#include <stdint.h>
#include <vector>

struct TrigEvent
{
   bool   complete; // event is complete
   bool   error;    // event has an error
   int    counter;  // event sequential counter
   double time;     // event time, sec
   double timeIncr; // time from previous event, sec

   std::vector<uint32_t> udpData;

   TrigEvent(); // ctor
   ~TrigEvent(); // dtor
   void Print(int level=0) const;
};

#endif

/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */
