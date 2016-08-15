
#include <stdio.h>

#include "Unpack.h"

Alpha16Event* UnpackAlpha16Event(Alpha16EVB* evb, const TMidasEvent* me)
{
   Alpha16Event* e = evb->NewEvent();

   for (int imodule = 0; imodule < 9; imodule++) {
      for (int i=0; i<NUM_CHAN_ALPHA16; i++) {
         void *ptr;
         int bklen, bktype;
         int status;
         
         char bname[5];
         sprintf(bname, "A%01d%02d", 1+imodule, i);
         
         status = me->FindBank(bname, &bklen, &bktype, &ptr);
         if (status == 1) {
            //printf("ALPHA16 bname %s, pointer: %p, status %d, len %d, type %d\n", bname, ptr, status, bklen, bktype);
            
            // print header
            
            int packetType = Alpha16Packet::PacketType(ptr, bklen);
            int packetVersion = Alpha16Packet::PacketVersion(ptr, bklen);
            
            if (0) {
               printf("Header:\n");
               printf("  packet type:    0x%02x (%d)\n", packetType, packetType);
               printf("  packet version: 0x%02x (%d)\n", packetVersion, packetVersion);
            }
            
            if (packetType == 1 && packetVersion == 1) {
               evb->AddBank(e, imodule, ptr, bklen);
            } else {
               printf("unknown packet type %d, version %d\n", packetType, packetVersion);
            }
         }
      }
   }

   evb->CheckEvent(e);

   return e;
};

/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */
