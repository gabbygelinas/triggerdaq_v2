//
// Trg.cxx
// K.Olchanski
//

#include <stdio.h>
#include <assert.h> // assert()

#include "Trig.h"

TrigEvent::TrigEvent() // ctor
{
   // empty
}

TrigEvent::~TrigEvent() // dtor
{
   // empty
}

void TrigEvent::Print(int level) const
{
   printf("TrgEvent %d, time %f, incr %f, complete %d, error %d ", counter, time, timeIncr, complete, error);
   if (level > 0) {
      printf("\n");
      for (unsigned i=0; i<udpData.size(); i++) {
         printf("udpData[%d]: 0x%08x (%d)\n", i, udpData[i], udpData[i]);
      }
   }
}

/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */
