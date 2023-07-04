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


class TPCSimEvent
{
 public:
   bool fComplete = false; // event is complete
   bool fError = false;    // event has an error
   std::string fError_message; // error message
   int  fCounter = 0;      // event sequential counter
   std::array<double, 3> fVertex;  // simulated verex position

 public:
   TPCSimEvent(); // ctor
   ~TPCSimEvent(); // dtor

 public:
   void Print(int level=0) const;

public:
    void SetVertex(std::array<double, 3> a_vertex) { fVertex = a_vertex; }
    void SetVertexX(double x) { fVertex[0] = x;}
    void SetVertexY(double y) { fVertex[1] = y;}
    void SetVertexZ(double z) { fVertex[2] = z;}
    std::array<double, 3> GetVertex() const { return fVertex;}
    double GetVertexX() const { return fVertex[0]; }
    double GetVertexY() const { return fVertex[1]; }
    double GetVertexZ() const { return fVertex[2]; }

};

#endif

/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */
