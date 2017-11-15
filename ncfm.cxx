//
// ncfm.cxx - CFM configuration database
//
// Konstantin Olchanski
//

#include "ncfm.h"

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <assert.h> // assert()

#include <string>
#include <vector>

struct Index
{
  int index1;
  int index2;
  int rev;
};

class NcfmData
{
public:
   std::string sys;
   std::string subsys;
   std::vector<Index> fIndex;

public:
   void Print() const;
   int FindRev(int runno);
};

NcfmData* Ncfm::LoadIndexFile(const char* system, const char* subsystem) const
{
   std::string f = fRoot + "/" + system + "_" + subsystem + ".txt";
   
   FILE* fp = fopen(f.c_str(),"r");
   if (!fp) {
      fprintf(stderr,"Cannot read CFM index file \'%s\', errno %d (%s)\n", f.c_str(), errno, strerror(errno));
      return NULL;
   }

   NcfmData* c = new NcfmData();
   c->sys = system;
   c->subsys = subsystem;
  
   while (1) {
      char buf[128];
      if (!fgets(buf,sizeof(buf)-1,fp))
         break;

      //printf("read \'%s\'\n",buf);
      
      if (buf[0] == '#')
         continue;
      
      char* s = buf;
      int index1 = strtol(s,&s,0);
      int index2 = strtol(s,&s,0);
      int rev    = strtol(s,&s,0);

      Index ind;
      ind.index1 = index1;
      ind.index2 = index2;
      ind.rev = rev;

      c->fIndex.push_back(ind);
    }

   printf("CFM: Loaded %d entries from \'%s\'\n", (int)c->fIndex.size(), f.c_str());

   return c;
}

void NcfmData::Print() const
{
   printf("Index for system \'%s\' subsystem \'%s\': %d entries\n", sys.c_str(), subsys.c_str(), (int)fIndex.size());
   for (unsigned i=0; i<fIndex.size(); i++) {
      printf("entry %d: index: %d %d, revision: %d\n", i, fIndex[i].index1, fIndex[i].index2, fIndex[i].rev);
   }
}

NcfmData* Ncfm::GetIndex(const char* system, const char* subsystem)
{
   // find already loaded index file

   for (unsigned i=0; i<fData.size(); i++) {
      if (fData[i]->sys == system && fData[i]->subsys == subsystem)
         return fData[i];
   }
   
   // load index file

   NcfmData *c = LoadIndexFile(system, subsystem);
   if (c == NULL)
      return NULL;

   c->Print();

   fData.push_back(c);

   return c;
}

int NcfmData::FindRev(int runno)
{
   for (int i=fIndex.size()-1; i>=0; i--) {
      //printf("try entry %d\n",i);
      
      // check for open-end entries: ttt1 0 rev
      if ((runno >= fIndex[i].index1) && (fIndex[i].index2 == 0))
         return fIndex[i].rev;
      
      // check for interval entries: ttt1 ttt2 rev
      if ((runno >= fIndex[i].index1) && (runno <= fIndex[i].index2))
         return fIndex[i].rev;
   }
   
   // found no matching entries...
   return -1;
}

Ncfm::Ncfm(const char* root_dir) // ctor
{
   if (!root_dir)
      root_dir = ".";
   fRoot = root_dir;
}

Ncfm::~Ncfm() // dtor
{
   for (unsigned i=0; i<fData.size(); i++) {
      if (fData[i]) {
         delete fData[i];
         fData[i] = NULL;
      }
   }
}

int Ncfm::GetRev(const char* system, const char* subsystem, int runno)
{
   NcfmData *ctrl = GetIndex(system, subsystem);
   if (ctrl == NULL)
      return -1;
   return ctrl->FindRev(runno);
}

static std::string RevToString(int rev)
{
  char buf[128];
  sprintf(buf,"%06d",rev);
  return buf;
}

std::string Ncfm::MakeFilename(const char* system, const char* subsystem, int rev) const
{
   return fRoot + "/" + system + "_" + subsystem + "_" + RevToString(rev) + ".txt";
}

std::string Ncfm::GetFilename(const char* system, const char* subsystem, int runno)
{
   int rev = GetRev(system, subsystem, runno);
   if (rev < 0)
      return "";
   return MakeFilename(system, subsystem, rev);
}

std::vector<std::string> Ncfm::ReadFile(const char* system, const char* subsystem, int runno)
{
   std::vector<std::string> v;

   std::string f = GetFilename(system, subsystem, runno);
   if (f.length() < 1)
      return v;

   return ReadFile(f.c_str());
}

std::vector<std::string> Ncfm::ReadFile(const char* filename) const
{
   std::vector<std::string> v;
  
   FILE* fp = fopen(filename, "r");

   if (!fp) {
      fprintf(stderr,"Cannot read CFM file \'%s\', errno %d (%s)\n", filename, errno, strerror(errno));
      return v;
   }

   while (1) {
      char buf[1024];
      if (!fgets(buf,sizeof(buf)-1,fp))
         break;

      //printf("read \'%s\'\n",buf);
      
      if (buf[0] == '#')
         continue; // skip comments

      // strip trailing \n and \r

      int len = strlen(buf);
      if (len > 0) {
         if (buf[len-1] == '\n')
            buf[len-1] = 0;
      }

      len = strlen(buf);

      if (len > 0) {
         if (buf[len-1] == '\r')
            buf[len-1] = 0;
      }

      len = strlen(buf);

      if (len > 0) {
         if (buf[len-1] == '\n')
            buf[len-1] = 0;
      }

      len = strlen(buf);

      if (len < 1)
         continue; // skip empty lines

      v.push_back(buf);
   }

   fclose(fp);

   return v;
}

#if 0
int test(Ncfm* cfm, const char*sys, const char*subsys, int runno)
{
   int rev = cfm->GetRev(sys, subsys, runno);
   std::string f = cfm->GetFilename(sys, subsys, rev);
   printf("Test \'%s\' \'%s\' index %d: rev %d, file \'%s\'\n", sys, subsys, runno, rev, f.c_str());
   return 0;
}

int main(int argc,char*argv[])
{
   setbuf(stdout,0);
   setbuf(stderr,0);
   
   Ncfm* cfm = new Ncfm(NULL);
   
   test(cfm, "test", "sub", 10);
   for (int i=0; i<35; i++)
      test(cfm, "test", "sub", i);

   delete cfm;
   return 0;
}
#endif

/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */
