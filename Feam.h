//
// Unpacking FEAM data
// K.Olchanski
//

#ifndef Feam_H
#define Feam_H

#include <stdint.h>
#include <vector>

class FeamPacket
{
public:
   uint32_t cnt;
   uint16_t n;
   uint16_t x511;
   uint16_t buf_len;
   uint32_t ts_start;
   uint32_t ts_trig;
   int off;
   bool error;

public:
   FeamPacket(); // ctor
   ~FeamPacket(); // dtor
   void Unpack(const char* data, int size);
   void Print() const;
};

class FeamModuleData
{
public:
   int module;

   uint32_t cnt;
   uint32_t ts_start;
   uint32_t ts_trig;

   bool error;

   uint32_t next_n;

   int fSize;
   char* fPtr;

   double   fTime;
   double   fTimeIncr;

public:
   FeamModuleData(const FeamPacket* p, int xmodule);
   ~FeamModuleData(); // dtor
   void AddData(const FeamPacket*p, int module, const char* ptr, int size);
   void Finalize();
   void Print() const;
};

#define MAX_FEAM   8
#define MAX_FEAM_SCA    4
#define MAX_FEAM_CHAN  80
#define MAX_FEAM_BINS 900

struct FeamAdcData
{
   int nsca;
   int nchan;
   int nbins;

   int adc[MAX_FEAM_SCA][MAX_FEAM_CHAN][MAX_FEAM_BINS];
};

struct FeamEvent
{
   bool   complete;
   bool   error;
   int    counter;
   double time;

   std::vector<FeamModuleData*> modules;
   std::vector<FeamAdcData*> adcs;

   FeamEvent(); // ctor
   ~FeamEvent(); // dtor
   void Print() const;
};

extern void Unpack(FeamAdcData* a, FeamModuleData* m);

#endif

/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */


