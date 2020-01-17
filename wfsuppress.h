// wfsuppress.h - waveform suppression

#include <stdint.h>

class WfSuppress
{
 private:
  uint16_t a0;
  uint16_t a1;
  uint16_t a2;
  uint16_t a3;
  uint16_t a4;
  uint16_t a5;
  uint16_t a6;
  uint16_t a7;
  uint16_t abase;
  int16_t  amp;
  int16_t  ampMin;
  int16_t  ampMax;
  uint16_t amin;
  uint16_t amax;
  int16_t  threshold;
  bool     above_threshold;
  bool     below_threshold;
  bool     keep;

 public:
  WfSuppress(); // ctor
  ~WfSuppress(); // dtor

 public:
  void Init(uint16_t a, int16_t t);
  bool Add(uint16_t a);
  uint16_t GetBase() const;
  int16_t GetAmp() const;
  void Print() const;
};

/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */
