//
// Unpacking trigger data
// K.Olchanski
//

#ifndef TrgAsm_H
#define TrgAsm_H

#include "Trig.h"


class TrgAsm
{
public: // state
   int      fCounter        = 0; // our packet counter
   uint32_t fFirstUdpPacket = 0; // seqno of first udp packet
   double   fFirstTime = 0;
   double   fPrevTime  = 0;
   uint32_t fLastTs    = 0;
   int      fEpoch     = 0;

public:
   TrgAsm();
   ~TrgAsm();

public:
   void Reset();
   TrigEvent* UnpackBank(const char* ptr, int size);
};

#endif

/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */


