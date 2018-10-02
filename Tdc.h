//
// ALPHA-g TPC
//
// TRB3 TDC functions
//
// Class definitions
//

#ifndef TDC_H
#define TDC_H

#include <stdint.h>
#include <vector>
#include <string>

class TdcHit
{
 public:
   // NB: be careful with alignement of these short data types!
   uint32_t epoch = 0;       // 24 bits
   uint8_t  fpga = 0;        // 0..3
   uint8_t  chan = 0;        // 7 bits
   uint16_t fine_time = 0;   // 10 bits
   uint16_t rising_edge = 0; // 1 bit
   uint16_t coarse_time = 0; // 11 bits
   void Print() const;
};

class TdcEvent
{
 public:
   bool complete = false; // event is complete
   bool error = false;    // event has an error
   std::string error_message; // error message
   int  counter = 0;      // event sequential counter
   double time = 0;       // event time, sec
   double timeIncr = 0;   // time from previous event, sec

   std::vector<TdcHit*> hits;

 public:
   TdcEvent(); // ctor
   ~TdcEvent(); // dtor

 public:
   void Print(int level=0) const;
};

class TdcAsm
{
 public: // configuration
   //double fTsFreq = 125000000; // 125 MHz timestamp clock
   //Alpha16Map fMap;

 public:
   TdcAsm(); // ctor
   ~TdcAsm(); // dtor

 public: // member functions
   void Print() const;
   TdcEvent* UnpackBank(const void* bkptr, int bklen8);

 public: // internal state
   int fEventCount = 0; // event counter

   double   fFirstEventTime = 0;
   uint32_t fLastEventTs = 0;
   double   fLastEventTime = 0;
   int      fTsEpoch = 0;
};

#endif

/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */
