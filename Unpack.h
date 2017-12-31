//
// ALPHA-g experiment
//
// Unpack ALPHA-g data from MIDAS event banks into C++ structures
//
// K.Olchanski
//

#ifndef UNPACK_H
#define UNPACK_H

#include "Alpha16.h"
#include "Feam.h"
#include "FeamEVB.h"
#include "midasio.h"

Alpha16Event* UnpackAlpha16Event(Alpha16EVB* evb, TMEvent* me);
FeamEvent* UnpackFeamEvent(FeamEVB* evb, TMEvent* me, const std::vector<std::string>& banks);

#endif

/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */


