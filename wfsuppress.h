// wfsuppress.h - waveform suppression

#include <stdint.h>

class WfSuppress
{
 private:
  int16_t a0;
  int16_t a1;
  int16_t a2;
  int16_t a3;
  int16_t a4;
  int16_t a5;
  int16_t a6;
  int16_t a7;
  int16_t abase;
  int16_t amp;
  int16_t ampMin;
  int16_t ampMax;
  int16_t amin;
  int16_t amax;
  int16_t threshold;
  bool    clipped;
  bool    above_threshold;
  bool    below_threshold;
  bool    keep;

 public:
  WfSuppress(); // ctor
  ~WfSuppress(); // dtor

 public:
  void Init(int16_t a, int16_t t);
  bool Add(int16_t a);
  bool    GetClipped() const { return clipped; };
  int16_t GetBase()    const { return abase; };
  int16_t GetAdcMin()  const { return amin; };
  int16_t GetAdcMax()  const { return amax; };
  int16_t GetAmp()     const { return amp; };
  int16_t GetAmpMin()  const { return ampMin; };
  int16_t GetAmpMax()  const { return ampMax; };
  void Print() const;
};

/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */
