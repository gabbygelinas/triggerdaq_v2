//
// AgAsm.cxx
//
// ALPHA-g event assembler
// K.Olchanski
//

#include "AgAsm.h"

#include <stdio.h> // NULL, printf()
//#include <math.h> // fabs()
//#include <assert.h> // assert()

AgAsm::AgAsm()
{
   // empty
}

AgAsm::~AgAsm()
{
   if (fAdcAsm) {
      delete fAdcAsm;
      fAdcAsm = NULL;
   }

   if (fPwbMap) {
      delete fPwbMap;
      fPwbMap = NULL;
   }

   if (fPwbAsm) {
      delete fPwbAsm;
      fPwbAsm = NULL;
   }

   if (fFeamAsm) {
      delete fFeamAsm;
      fFeamAsm = NULL;
   }
}

AgEvent* AgAsm::UnpackEvent(TMEvent* me)
{
   AgEvent* e = new AgEvent();

   me->FindAllBanks();

   for (unsigned i=0; i<me->banks.size(); i++) {
      const TMBank* b = &me->banks[i];
      //printf("bank %s\n", b->name.c_str());
      if (0) {
      } else if (b->name == "ATAT") {
         if (!fTrgAsm) {
            fTrgAsm = new TrgAsm();
         }
         
         const char* bkptr = me->GetBankData(b);
         int bklen = b->data_size;

         e->trig = fTrgAsm->UnpackBank(bkptr, bklen);

         if (1) {
            printf("Unpacked TRG event: ");
            e->trig->Print();
            printf("\n");
         }
      } else if (b->name[0] == 'A') {
         // ADC UDP packet bank from feudp
      } else if (b->name[0] == 'B' && b->name[1] == 'B') {
         // obsolete bank "BBnn" from FEAM rev0 board
      } else if (b->name[0] == 'B' || b->name[0] == 'C') {
         // ADC bank from feevb
         int c1 = b->name[1]-'0';
         int c2 = b->name[2]-'0';
         int module = c1*10 + c2;
         if (module < 1 || module > ADC_MODULE_LAST) {
            fprintf(stderr, "AgAsm::UnpackEvent: bank name [%s] has invalid module number %d\n", b->name.c_str(), module);
            e->error = true;
            continue;
         }

         const char* bkptr = me->GetBankData(b);
         int bklen = b->data_size;

         if (!fAdcAsm) {
            fAdcAsm = new Alpha16Asm();
            fAdcAsm->fMap.Init(fAdcMap);
            fAdcAsm->fMap.Print();
         }

         if (!e->a16) {
            e->a16 = fAdcAsm->NewEvent();
         }

         fAdcAsm->AddBank(e->a16, module, b->name.c_str(), bkptr, bklen);
      } else if (b->name[0] == 'P' && ((b->name[1] == 'A') || (b->name[1] == 'B'))) {
         // PWB bank
         int c1 = b->name[2]-'0';
         int c2 = b->name[3]-'0';
         int imodule = c1*10 + c2;

         if (imodule < 0 || imodule > PWB_MODULE_LAST) {
            fprintf(stderr, "UnpackAlpha16Event: bank name [%s] has invalid module number %d\n", b->name.c_str(), imodule);
            e->error = true;
            continue;
         }

         const PwbModuleMapEntry* map = fPwbMap->FindPwb(imodule);
         
         //int ii = 0; // FIXME
         //
         //int column = (ii/8);
         //int ring = (ii%8);
         //
         //bool short_tpc = false; // FIXME (runinfo->fRunNo < 1450);
         //
         //if (short_tpc) {
         //   column = ii;
         //   ring = 0;
         //}
            
         const char* p8 = me->GetBankData(b);
         const uint32_t *p32 = (const uint32_t*)p8;
         const int n32 = b->data_size/4;
         
         if (0) {
            unsigned nprint = b->data_size/4;
            nprint=10;
            for (unsigned i=0; i<nprint; i++) {
               printf("PB05[%d]: 0x%08x (%d)\n", i, p32[i], p32[i]);
            }
         }

         int f = 0;
         if (b->name[0] == 'P' && b->name[1] == 'A')
            f = 1;
         else if (b->name[0] == 'P' && b->name[1] == 'B')
            f = 2;
         else {
            fprintf(stderr, "AgAsm::UnpackEvent: invalid PWB bank name [%s]\n", b->name.c_str());
            e->error = true;
            continue;
         }
         
         if (f == 1 || f == 2) { // old UDP data format
            if (!fFeamAsm) {
               fFeamAsm = new FeamAsm();
            }

            if (b->data_size < 26) {
               fprintf(stderr, "AgAsm::UnpackEvent: bank name [%s] has invalid FEAM packet length %d\n", b->name.c_str(), b->data_size);
               e->error = true;
               continue;
            }
            
            FeamPacket* p = new FeamPacket();
            
            p->Unpack(p8, b->data_size);
            
            if (p->error) {
               fprintf(stderr, "AgAsm::UnpackEvent: cannot unpack FeamPacket from bank [%s]\n", b->name.c_str());
               delete p;
               p = NULL;
               e->error = true;
               continue;
            }
            
            fFeamAsm->AddPacket(imodule, map->fColumn, map->fRing, f, p, p8 + p->off, p->buf_len);
         } else { // new UDP data format
            if (!fPwbAsm) {
               fPwbAsm = new PwbAsm();
            }
            fPwbAsm->AddPacket(imodule, map->fColumn, map->fRing, p8, b->data_size);
         }
      } else {
         fprintf(stderr, "AgAsm::UnpackEvent: unknown bank [%s]\n", b->name.c_str());
      }
   }

   if (e->a16) {
      fAdcAsm->CheckEvent(e->a16);
   }

   if (fPwbAsm) {
      //if (fPwbAsm->CheckComplete()) {
      //printf("PwbAsm ---> complete !!!\n");
      if (!e->feam) {
         e->feam = new FeamEvent();
      }
      fPwbAsm->BuildEvent(e->feam);
      printf("PwbAsm built an event:\n");
      e->feam->Print();
      printf("\n");
      //PrintFeamChannels(e->feam->hits);
   }

   if (fFeamAsm) {
      if (!e->feam) {
         e->feam = new FeamEvent();
      }
      fFeamAsm->BuildEvent(e->feam);
      printf("FeamAsm built an event:\n");
      e->feam->Print();
      printf("\n");
      //PrintFeamChannels(e->feam->hits);
   }

   if (1) {
      printf("AgAsm::UnpackEvent: returning event: ");
      e->Print();
      printf("\n");
   }
  
   return e;
}

/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */


