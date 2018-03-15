//
// Unpack.cxx
// K.Olchanski
//

#include <stdio.h>
#include <assert.h> // assert()

#include "Unpack.h"

static int adc_bank_name_to_module(const char* s)
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

static int pwb_bank_name_to_module(const char* s)
{
   if (s[0] == 'P' && (s[1] == 'A' || s[1] == 'B')) {
      int c1 = s[2]-'0';
      int c2 = s[3]-'0';
      return c1*10 + c2;
   } else {
      fprintf(stderr, "Cannot extract module name from bank name [%s]\n", s);
      abort();
   }
}

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

TrigEvent* UnpackTrigEvent(TMEvent* event, TMBank* atat_bank)
{
   TrigEvent* e = new TrigEvent;
   
   const char* p8 = event->GetBankData(atat_bank);
   const uint32_t *p32 = (const uint32_t*)p8;
   for (unsigned i=0; i<atat_bank->data_size/4; i++) {
      //printf("ATAT[%d]: 0x%08x (%d)\n", i, p32[i], p32[i]);
      e->udpData.push_back(p32[i]);
   }

   if (e->udpData.size() < 9) {
      e->complete = false;
      e->error = true;
      return e;
   }

   e->complete = true;
   e->error = false;

   static uint32_t gFirstCounter = 0;

   if (gFirstCounter == 0)
      gFirstCounter = e->udpData[0];
   
   e->counter = e->udpData[0] + 1 - gFirstCounter; // udp packet counter counts from 0, we want our counter to count from 1

   double ts_freq = 62.5*1e6; // timestamp is 62.5 MHz

   static double gFirstTime = 0;
   static double gPrevTime = 0;
   static uint32_t gLastTs = 0;
   static int gEpoch = 0;

   uint32_t ts = e->udpData[2];

   if (ts < gLastTs)
      gEpoch++;
   gLastTs = ts;

   if (gFirstTime == 0) {
      gFirstTime = ts/ts_freq;
      gPrevTime = 0;
   }
   
   e->time = ts/ts_freq - gFirstTime + gEpoch*2.0*0x80000000/ts_freq;
   e->timeIncr = e->time - gPrevTime;
   gPrevTime = e->time;

   return e;
}

Alpha16Event* UnpackAlpha16Event(Alpha16Asm* adcasm, TMEvent* me)
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
         int module = adc_bank_name_to_module(b->name.c_str());

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
               e = adcasm->NewEvent();
            }
            adcasm->AddChannel(e, p, c);
         } else {
            printf("unknown packet type %d, version %d\n", packetType, packetVersion);
         }
      }
   }

   if (e) {
      adcasm->CheckEvent(e);
      //e->Print();
      //printf("\n");
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
         int module = pwb_bank_name_to_module(b->name.c_str());

         if (module < 0 || module > 80) {
            fprintf(stderr, "UnpackFeamEventNoEvb: bank name [%s] has invalid module number %d\n", b->name.c_str(), module);
            abort();
         }

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

            int column = (i/8);
            int ring = (i%8);
            
            evb->AddPacket(i, module, column, ring, f, p, data + p->off, p->buf_len);
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
            int module = pwb_bank_name_to_module(b->name.c_str());

            if (module < 0 || module > 80) {
               fprintf(stderr, "UnpackFeamEvent: bank name [%s] has invalid module number %d\n", b->name.c_str(), module);
               abort();
            }

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

               int column = (i/8);
               int ring = (i%8);

               evb->AddPacket(i, module, column, ring, f, p, data + p->off, p->buf_len);
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
