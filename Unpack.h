//
// ALPHA-g experiment
//
// Unpacking of ALPHA-g data from MIDAS events
//
//

#ifndef UNPACK_H
#define UNPACK_H

#include "Alpha16.h"
#include "TMidasEvent.h"

Alpha16Event* UnpackAlpha16Event(Alpha16EVB* evb, const TMidasEvent* me);

#endif

/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */


