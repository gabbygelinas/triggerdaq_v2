//
// Unpack.cxx
// K.Olchanski
//

#include <stdio.h>
#include <assert.h> // assert()

#include "Unpack.h"

static int bank_name_to_module(const char* s)
{
   if (s[0] == 'B' || s[0] == 'C') {
      int c1 = s[1]-'0';
      int c2 = s[2]-'0';
      return c1*10 + c2;
   } else {
      fprintf(stderr, "Cannot extract module name from bank name [%s]\n", s);
      abort();
   }
}

Alpha16Event* UnpackAlpha16Event(Alpha16EVB* evb, TMEvent* me)
{
   me->FindAllBanks();
   
   Alpha16Event* e = NULL;

   for (unsigned i=0; i<me->banks.size(); i++) {
      const TMBank* b = &me->banks[i];
      if (b->name[0] == 'A') {
         // UDP packet bank from feudp
      } else if (b->name[0] == 'B' && b->name[1] == 'B') {
         // obsolete bank "BBnn" from FEAM rev0 board
      } else if (b->name[0] == 'B' || b->name[0] == 'C') {
         // ALPHA16 bank from feevb
         int module = bank_name_to_module(b->name.c_str());

         if (module < 1 || module > 20) {
            fprintf(stderr, "UnpackAlpha16Event: bank name [%s] has invalid module number %d\n", b->name.c_str(), module);
            abort();
         }

         const void* bkptr = me->GetBankData(b);
         int bklen = b->data_size;

         int packetType = Alpha16Packet::PacketType(bkptr, bklen);
         int packetVersion = Alpha16Packet::PacketVersion(bkptr, bklen);

         if (0) {
            printf("ALPHA16 bname %s, bank %p, len %d\n", b->name.c_str(), b, bklen);
            printf("Header:\n");
            printf("  packet type:    0x%02x (%d)\n", packetType, packetType);
            printf("  packet version: 0x%02x (%d)\n", packetVersion, packetVersion);
            
            Alpha16Packet *p = Alpha16Packet::Unpack(bkptr, bklen);
            p->Print();
            delete p;
         }
            
         if (packetType == 1 && packetVersion == 1) {
            Alpha16Packet* p = Alpha16Packet::Unpack(bkptr, bklen);
            Alpha16Channel* c = Unpack(b->name.c_str(), module, p, bkptr, bklen);

            //p->Print();
            //printf("\n");
            
            //c->Print();
            //printf("\n");
            
            if (!e) {
               e = evb->NewEvent();
            }
            evb->AddBank(e, p, c);
         } else {
            printf("unknown packet type %d, version %d\n", packetType, packetVersion);
         }
      }
   }

   if (e) {
      evb->CheckEvent(e);
      e->Print();
      printf("\n");
   }

   return e;
};

FeamEvent* UnpackFeamEventNoEvb(FeamEVB* evb, TMEvent* event, const std::vector<std::string> &banks)
{
   //printf("event id %d\n", event->event_id);
   
   if (event->event_id != 1)
      return NULL;
      
   char *data = NULL;
   
   for (unsigned i=0; i<banks.size(); i++) {
      TMBank* b = event->FindBank(banks[i].c_str());
      if (b) {
         data = event->GetBankData(b);
         if (data) {
            //printf("Have bank %s\n", banks[i]);
            //HandleFeam(i, data, b->data_size);
            
            if (b->data_size < 26) {
               printf("bad FEAM %d packet length %d\n", i, b->data_size);
               continue;
            }
            
            FeamPacket* p = new FeamPacket();
            
            p->Unpack(data, b->data_size);
               
            assert(!p->error);
            
            int f = 0;
            if (b->name[0] == 'P' && b->name[1] == 'A')
               f = 1;
            else if (b->name[0] == 'P' && b->name[1] == 'B')
               f = 2;
            else if (b->name[0] == 'B' && b->name[1] == 'B')
               f = 1;
            else
               assert(!"invalid PWB bank name");
            
            evb->AddPacket(b->name.c_str(), i, f, p, data + p->off, p->buf_len);
         }
      }
   }
   
   return evb->Get();
}

FeamEvent* UnpackFeamEvent(FeamEVB* evb, TMEvent* event, const std::vector<std::string> &banks)
{
   //printf("event id %d\n", event->event_id);
   
   if (event->event_id != 1)
      return NULL;
      
   for (unsigned k=0; k<event->banks.size(); k++) {
      for (unsigned i=0; i<banks.size(); i++) {
         const TMBank* b = &event->banks[k];
         if (banks[i] == b->name) {
            char *data = event->GetBankData(b);
            if (data) {
               //printf("Have bank %s\n", banks[i]);
               //HandleFeam(i, data, b->data_size);
               
               if (b->data_size < 26) {
                  printf("bad FEAM %d packet length %d\n", i, b->data_size);
                  continue;
               }
               
               FeamPacket* p = new FeamPacket();
               
               p->Unpack(data, b->data_size);
               
               assert(!p->error);

               int f = 0;
               if (b->name[0] == 'P' && b->name[1] == 'A')
                  f = 1;
               else if (b->name[0] == 'P' && b->name[1] == 'B')
                  f = 2;
               else if (b->name[0] == 'B' && b->name[1] == 'B')
                  f = 1;
               else
                  assert(!"invalid PWB bank name");
            
               evb->AddPacket(b->name.c_str(), i, f, p, data + p->off, p->buf_len);
            }
            break;
         }
      }
   }
   
   return evb->Get();
}

/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */
