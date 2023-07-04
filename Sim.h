//
// ALPHA-g TPC
//
// Simulation unpacking functions
//
// Class definitions
//

#ifndef SIM_H
#define SIM_H

#include <stdint.h>
#include <vector>
#include <string>
#include<iostream>
#include<array>
#include<map>


class SimEvent
{
 public:
   bool fComplete = false; // event is complete
   bool fError = false;    // event has an error
   std::string fError_message; // error message
   int  fCounter = 0;      // event sequential counter
   std::array<double, 3> fVertex;  // simulated verex position
   std::map<std::string, double> fSimDoubles; // other data from simulation, left open for the future
   std::map<std::string, double> fSimInts; // other data from simulation, left open for the future

 public:
   SimEvent(); // ctor
   ~SimEvent(); // dtor

 public:
   void Print(int level=0) const;

public:
    void SetDouble(const std::string& key, const double value) {fSimDoubles[key] = value;}
    double GetDouble(const std::string& key) {
      if (fSimDoubles.count(key)>0) return fSimDoubles[key];
      else return 0;
    }
    void SetInt(const std::string& key, const int value) {fSimInts[key] = value;}
    double GetInt(const std::string& key) {
      if (fSimInts.count(key)>0) return fSimInts[key];
      else return 0;
    }
    void SetVertex(std::array<double, 3> a_vertex) { fVertex = a_vertex; }
    void SetVertexX(double x) { fVertex[0] = x;}
    void SetVertexY(double y) { fVertex[1] = y;}
    void SetVertexZ(double z) { fVertex[2] = z;}
    std::array<double, 3> GetVertex() const { return fVertex;}
    double GetVertexX() const { return fVertex[0]; }
    double GetVertexY() const { return fVertex[1]; }
    double GetVertexZ() const { return fVertex[2]; }

};

class SimAsm
{
 public: // configuration

 public:
   SimAsm(); // ctor
   ~SimAsm(); // dtor

 public: // member functions
   void Print() const;
   void Reset();
   SimEvent* UnpackBank(const void* bkptr, int bklen8);

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
