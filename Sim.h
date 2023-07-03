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
   bool complete = false; // event is complete
   bool error = false;    // event has an error
   std::string error_message; // error message
   int  counter = 0;      // event sequential counter
   std::array<double, 3> vertex;  // simulated verex position
   std::map<std::string, double> sim_doubles; // other data from simulation, left open for the future
   std::map<std::string, double> sim_ints; // other data from simulation, left open for the future

 public:
   SimEvent(); // ctor
   ~SimEvent(); // dtor

 public:
   void Print(int level=0) const;

public:
    void SetDouble(std::string key, double value) {sim_doubles[key] = value;}
    double GetDouble(std::string key) {
      if (sim_doubles.count(key)>0) return sim_doubles[key];
      else return 0;
    }
    void SetInt(std::string key, int value) {sim_ints[key] = value;}
    double GetInt(std::string key) {
      if (sim_ints.count(key)>0) return sim_ints[key];
      else return 0;
    }
    void SetVertex(std::array<double, 3> a_vertex) { vertex = a_vertex; }
    void SetVertexX(double x) { vertex[0] = x;}
    void SetVertexY(double y) { vertex[1] = y;}
    void SetVertexZ(double z) { vertex[2] = z;}
    std::array<double, 3> GetVertex() { return vertex;}
    double GetVertexX() { return vertex[0]; }
    double GetVertexY() { return vertex[1]; }
    double GetVertexZ() { return vertex[2]; }

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
