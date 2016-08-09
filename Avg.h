// Avg.h

struct Avg
{
  double sum0;
  double sum1;
  double sum2;

  Avg() // ctor
  {
    Reset();
  }

  void Reset()
  {
    sum0 = 0;
    sum1 = 0;
    sum2 = 0;
  }

  void Add(double v)
  {
    sum0 += 1;
    sum1 += v;
    sum2 += v*v;
  }

  double N()
  {
    return sum0;
  }

  double Mean()
  {
    return sum1/sum0;
  }

  double Variance()
  {
    double mean = Mean();
    return sum2/sum0 - mean*mean;
  }

  double RMS()
  {
    double var = Variance();
    if (var > 0)
      return sqrt(var);
    else
      return 0;
  }
};

/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */
