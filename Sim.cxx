//
// ALPHA-g TPC
//
// Simulation unpacking functions
//
// Class functions for Sim.h
//

#undef NDEBUG // this program requires working assert()

#include "Sim.h"

#include <stdio.h>
//#include <string.h>
#include <assert.h> // assert()


static std::string toString(int v)
{
   char buf[256];
   sprintf(buf, "%d", v);
   return buf;
}

SimEvent::SimEvent() // ctor
{
}

SimEvent::~SimEvent() // dtor
{
}

void SimEvent::Print(int level) const
{

}

SimAsm::SimAsm() // ctor
{
}

SimAsm::~SimAsm() // dtor
{
   printf("SimAsm: %d events\n", fEventCount);
}

void SimAsm::Print() const
{
   printf("SimAsm::Print!\n");
}

void SimAsm::Reset()
{
   fEventCount = 0;
}

static uint32_t getUint32(const void* ptr, int offset)
{
  uint8_t *ptr8 = (uint8_t*)(((char*)ptr)+offset);
  return (ptr8[0]<<24) | (ptr8[1]<<16) | (ptr8[2]<<8) | ptr8[3];
}

// Take a vector of bytes and convert it into a double assuming it is little
// endian
// This should be platform independent as long as the `int` and `float`
// endianness is the same (which it is on all modern platforms)
double from_le_bytes(std::vector<uint8_t> bytes) {
	double value = 0;
	uint64_t* ptr = reinterpret_cast<uint64_t*>(&value);
	for (size_t i = 0; i < sizeof(double); i++) {
		*ptr <<= 8;
		*ptr |= bytes[sizeof(double) - i - 1];
	}

	return value;
}

// Take a `const char*` which is a pointer to some bytes. Convert those bytes
// into 3 consecutive doubles and return them in an array
// The elements of the array represent (in order) the x, y, and z coordinates
// of the MC vertex.
std::array<double, 3> read_mcvx(const char* bytes) {
	std::array<double, 3> values;
	for (size_t i = 0; i < 3; i++) {
		std::vector<uint8_t> bytes_vec(bytes + i * sizeof(double), bytes + (i + 1) * sizeof(double));
		values[i] = from_le_bytes(bytes_vec);
	}

	return values;
}

SimEvent* SimAsm::UnpackBank(const void* bkptr, int bklen)
{
   SimEvent* e = new SimEvent();

   e->counter = ++fEventCount;

   std::array<double, 3> vertex = read_mcvx((const char *)bkptr);
   e->SetVertex(vertex);
   e->complete = true;
   
   return e;
}

/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */
