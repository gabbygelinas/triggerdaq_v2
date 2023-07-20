//
// ALPHA-g TPC
//
// Simulation unpacking functions
//
// Class definitions
//

#ifndef SIM_HA
#define SIM_HA

#include <stdint.h>
#include <vector>
#include <string>
#include<iostream>
#include<array>
#include<map>
#include "TPCSimEvent.h"

class TPCSimAsm
{
 public: // configuration

 public:
   TPCSimAsm(); // ctor
   ~TPCSimAsm(); // dtor

 public: // member functions
   void Print() const;
   void Reset();
   TPCSimEvent* UnpackBank(const void* bkptr);

 public: // internal state
   int fEventCount = 0; // event counter
};

#endif

/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */
