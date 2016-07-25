#ifndef WAVEFORMH
#define WAVEFORMH

struct Waveform
{
  bool overflow;
  int nsamples;
  int offset;
  double *samples;

  Waveform(int num_samples) {
    overflow = false;
    nsamples = num_samples;
    offset = 0;
    samples = new double[nsamples];
  }

  ~Waveform() {
    nsamples = 0;
    delete samples;
    samples = NULL;
  }
};

#endif
/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */
