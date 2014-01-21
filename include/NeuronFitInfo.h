#ifndef NeuronFitInfo_h
#define NeuronFitInfo_h



#include "io_functions.h"
#include "ParameterDescription.h"
#include "NeuronGeometry.h"
#include "Trace.h"

#include <string>
#include <map>
#include <list>
#include <vector>
#include <set>
#include <sstream>



struct TraceInfo
{
  std::string name;
  std::string units;
  size_t numT;
  double dT;
  size_t traceNum;
};


struct OnePointInfo
{
  std::string name;
  std::string compartmentName;
};
// need to define < operator to insert into sets
inline bool operator<(const OnePointInfo& lhs, const OnePointInfo& rhs) {
  return (lhs.name + lhs.compartmentName) < (rhs.name + rhs.compartmentName);
}


struct TwoPointInfo
{
  std::string name;
  std::string preCompName;
  std::string postCompName;
};
// need to define < operator to insert into sets
inline bool operator<(const TwoPointInfo& lhs, const TwoPointInfo& rhs) {
  return (lhs.name + lhs.preCompName + lhs.postCompName) <
         (rhs.name + rhs.preCompName + rhs.postCompName);
}



class NeuronFitInfo
{
  public:
    virtual ~NeuronFitInfo() {}
    // read text of startup file
    void readStartupFile(const std::string & fileName);
    // output recorded traces to outFile
    void outputRecordedTraces(const std::string & fileName) const;
    // return a set describing where to insert OnePointInjectors
    const std::set<OnePointInfo> & getOnePointInjectorInfo(void) const {
      return onePointInfo;
    }
    // return a set describing where to insert TwoPointInjectors
    const std::set<TwoPointInfo> & getTwoPointInjectorInfo(void) const {
      return twoPointInfo;
    }
  
    std::list<Trace> electrophysData;
    ParameterDescriptionList parameterDescriptionList;
    NeuronGeometry::Geometry geometry;
    
    // directories containing channel shared objects
    std::set<std::string> injectorDirs;
    
    double tFinal;  // ms, final simulated/recorded time
    
  protected:
    // parse startup information and build member data
    void buildFitInfo(void);

    bool dataIsLoaded;
    std::vector<char> startupFileText;
    std::string startupFileName;
  
  private:
    // parse startup information to extract geometry info first
    void scanForGeometry(void);
    void addTime(const SplitLine & splitLine);
    void addGeometry(const SplitLine & splitLine);
    void addChannelDir(const SplitLine & splitLine);
    void checkValidity(const Trace & trace) const;
    void addRecordTrace(const SplitLine & splitLine);
    void addClampTrace(const SplitLine & splitLine);
    void addFitTrace(const SplitLine & splitLine);
    std::set<std::string>
      getCompartmentNamesByTag(const std::string & tag) const;
    void addChannels(std::stringstream & fileIn, size_t & lineNum,
                     SplitLine & splitLine);
    void addConnections(std::stringstream & fileIn, size_t & lineNum,
                     SplitLine & splitLine);
    void addChannel(const SplitLine & splitLine);
    void addConnection(const SplitLine & splitLine);
    void addParameters(std::stringstream & fileIn, size_t & lineNum,
                       SplitLine & splitLine);
    void addParameter(const SplitLine & splitLine);
    void scheduleLoadTrace(const std::string & fileName,
                           const size_t & traceNum);
    void loadTraces(void);
    std::vector<TraceInfo> getTraceHeader(std::ifstream & fileIn,
                                                size_t & lineNum) const;
    void getTraces(std::ifstream & fileIn, size_t & lineNum,
                   const std::vector<TraceInfo> & header,
                   std::vector<Trace*> & columns) const;

    std::map<std::string, std::vector<Trace*> > tracesToLoadMap;
    std::set<OnePointInfo> onePointInfo;
    std::set<TwoPointInfo> twoPointInfo;
};



#endif
