//
// ncfm.h - CFM configuration database
//
// Konstantin Olchanski
//

#ifndef NCFM_H
#define NCFM_H

#include <string>
#include <vector>
#include <map>

class NcfmParser
{
 public:
   std::map<std::string,std::string> fMap;

 public:
   NcfmParser(const std::vector<std::string>& file);

 public:
   int GetInt(const char* varname, int defaultValue = 0);
   double GetDouble(const char* varname, double defaultValue = 0);
   std::string GetString(const char* varname, const char* defaultValue = NULL);
};

class NcfmData;

class Ncfm
{
public:
   std::string fRoot;
   std::vector<NcfmData*> fData;

public:
   Ncfm(const char* root_dir); // ctor
   ~Ncfm(); // dtor
   
   void Print() const; // print the database
   
   int GetRev(const char* system, const char* subsystem, int runno);
   std::string GetFilename(const char* system, const char* subsystem, int runno);
   std::vector<std::string> ReadFile(const char* system, const char* subsystem, int runno);
   NcfmParser* ParseFile(const char* system, const char* subsystem, int runno);

public:
   std::string MakeFilename(const char* system, const char* subsystem, int rev) const;
   std::vector<std::string> ReadFile(const char* filename) const;

public:
   NcfmData* LoadIndexFile(const char* system, const char* subsystem) const;
   NcfmData* GetIndex(const char* system, const char* subsystem);
};

#endif
/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */
