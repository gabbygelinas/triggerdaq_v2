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

#include "Avg.h"

static void smooth(Waveform* w)
{
  for (int i=1; i<w->nsamples-1; i++) {
    double s = 0.25*w->samples[i-1] + 0.5*w->samples[i] + 0.25*w->samples[i+1];
    w->samples[i] = s;
  }
}

static void add(Waveform* w, const Waveform* w1, const Waveform* w2, double weight1, double weight2)
{
  assert(w1->nsamples == w2->nsamples);
  assert(w->nsamples == w1->nsamples);
  for (int i=0; i<w->nsamples; i++) {
    w->samples[i] = w1->samples[i]*weight1 + w2->samples[i]*weight2;
  }
}

static double max(const Waveform* w)
{
  double max = w->samples[0];
  for (int i=1; i<w->nsamples; i++)
    if (w->samples[i] > max)
      max = w->samples[i];
  return max;
}

static double min(const Waveform* w)
{
  double min = w->samples[0];
  for (int i=1; i<w->nsamples; i++)
    if (w->samples[i] < min)
      min = w->samples[i];
  return min;
}

static int maxbin(const Waveform* w)
{
  int m = 0;
  for (int i=1; i<w->nsamples; i++)
    if (w->samples[i] > w->samples[m])
      m = i;
  return m;
}

static int minbin(const Waveform* w)
{
  int m = 0;
  for (int i=1; i<w->nsamples; i++)
    if (w->samples[i] < w->samples[m])
      m = i;
  return m;
}

static double avgbins(const Waveform* w, int bin1, int bin2)
{
  if (bin1 < 0)
    bin1 = 0;
  if (bin2 >= w->nsamples)
    bin2 = w->nsamples-1;

  Avg a;
  for (int i=bin1; i<=bin2; i++)
    a.Add(w->samples[i]);
  return a.Mean();
}

static double avgmin(const Waveform* w)
{
  int b = minbin(w);
  return avgbins(w, b-2, b+2);
}

static double avgmax(const Waveform* w)
{
  int b = maxbin(w);
  return avgbins(w, b-2, b+2);
}

static double baseline(const Waveform* w, int n1, int n2, double *varp = NULL, double *rmsp = NULL)
{
   if (n1 < 0)
      n1 = 0;
   if (n2 <= n1)
      n2 = w->nsamples;
   
   if (n2 > w->nsamples)
      n2 = w->nsamples;
   
   double min = w->samples[0];
   double max = w->samples[0];
   
   double sum0 = 0;
   double sum1 = 0;
   double sum2 = 0;
   
   for (int i=n1; i<n2; i++) {
      double s = w->samples[i];
      sum0 += 1.0;
      sum1 += s;
      sum2 += s*s;
      if (s < min)
         min = s;
      if (s > max)
         max = s;
      //printf("%d %d %f %f %f %f %f\n", i, w->samples[i], sum0, sum1, sum2, min, max);
   }
   
   double mean = 0;
   double var = 0;
   double rms = 0;
   
   if (sum0 > 0) {
      mean = sum1/sum0;
      var = sum2/sum0 - mean*mean;
      if (var>0)
         rms = sqrt(var);
   }
   
   //printf("sum: %f %f %f, mean %f, var %f, rms %f\n", sum0, sum1, sum2, mean, var, sqrt(var));
   
   if (varp)
      *varp = var;
   if (rmsp)
      *rmsp = rms;
   return mean;
}

static void printbins(const Waveform*w, int bin1, int bin2)
{
  if (bin1 < 0)
    bin1 = 0;
  if (bin2 >= w->nsamples)
    bin2 = w->nsamples-1;

  printf("bins %d..%d: ", bin1, bin2);
  for (int i=bin1; i<=bin2; i++)
    printf(" %.1f", w->samples[i]);
}

static void checkmonotony(const Waveform*w, int bin1, int bin2, double dir)
{
  if (bin1 < 0)
    bin1 = 0;
  if (bin2 >= w->nsamples)
    bin2 = w->nsamples-1;

  printf(" M");
  for (int i=bin1; i<bin2; i++) {
    if ((w->samples[i+1] - w->samples[i])*dir < 0)
      printf(" %d", i);
  }
}

static int led(const Waveform* w, double baseline, double gain, double threshold)
{
   for (int i=0; i<w->nsamples; i++)
      if ((w->samples[i]-baseline)*gain > threshold) {
         //printf("thr %.1f, bin %d: ", threshold, i);
         //printbins(w, i-10, i+10);
         //checkmonotony(w, i-10, i+10, gain);
         //printf("\n");
         return i;
      }
   
   return 0;
}

static void fit(const Waveform* w, double baseline, double gain, double threshold, int le_bin, int nbins, double *posp, double *slopep, double *ampp)
{
  double sum0 = 0;
  double sum1x = 0;
  double sum1y = 0;
  double sum2x = 0;
  double sum2y = 0;

  for (int i=le_bin-nbins; i<=le_bin+nbins; i++) {
    double x = i;
    double y = (w->samples[i]-baseline)*gain;
    sum0 += 1;
    sum1x += x;
    sum2x += x*x;
    sum1y += y;
    sum2y += y*y;
  }

  double meanx = sum1x/sum0;
  double meany = sum1y/sum0;
  double varx = sum2x/sum0 - meanx*meanx;
  double vary = sum2y/sum0 - meany*meany;
  double rmsx = sqrt(varx);
  double rmsy = sqrt(vary);
  double slope = rmsy/rmsx;

  double pos = meanx - (threshold - meany)/slope;

  printf("bin %d, thr %f, mean %f %f, rms %f %f, slope %f, pos %f\n", le_bin, threshold, meanx, meany, rmsx, rmsy, slope, pos);

  if (posp)
    *posp = pos;
  if (slopep)
    *slopep = slope;
  if (ampp)
    *ampp = meany;
}

static void zoom(Waveform* wz, const Waveform* w, int start, double baseline, double gain)
{
  wz->offset = w->offset + start;
  for (int i=0, j=start; i<wz->nsamples && j<w->nsamples; i++, j++) {
    double s = 0;
    if (j>=0 && j<w->nsamples)
      s = (w->samples[j] - baseline)*gain;
    wz->samples[i] = s;
  }
}

static void zoom0(Waveform* wz, const Waveform* w)
{
  int s = 0;

  for (int i=0; i<w->nsamples; i++)
    if (w->samples[i] < 3600) {
      s = i - 20;
      break;
    }

  if (s<0)
    s = 0;

  double sum0 = 0;
  double sum1 = 0;

  for (int i=0; i<s; i++) {
    sum0 += 1.0;
    sum1 += w->samples[i];
  }
    
  int offset = 0;
  if (sum0 > 0)
    offset = sum1/sum0;

  for (int i=0, j=s; i<wz->nsamples && j<w->nsamples; i++, j++)
    wz->samples[i] = offset - w->samples[j];
}

static void diff(Waveform* wd, const Waveform* w)
{
  wd->samples[0] = 0;
  for (int i=1; i<wd->nsamples-1; i++)
    wd->samples[i] = (w->samples[i+1] - w->samples[i-1]);
  wd->samples[wd->nsamples-1] = 0;
}

#endif
/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */
