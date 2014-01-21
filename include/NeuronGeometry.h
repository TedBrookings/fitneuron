#ifndef _NEURON_GEOMETRY__H_
#define _NEURON_GEOMETRY__H_



#include "io_functions.h"


#include <vector>
#include <string>
#include <map>
#include <set>



// To do:
//  -create functions to write geometry file, given the lists of nodes and
//   segments.  This will be complicated to do *well* since it involves
//   analyzing tags and sorting the segments (indenting by distance from soma
//   would be nice, as would minimizing opening and closing of tags)
//  -create function to decrease the number of segments/nodes intelligently



namespace NeuronGeometry {
  struct Node;
  class Segment;
  
  class Geometry
  {
    public:
      virtual ~Geometry() {}
      // read geometry information from specified file
      void readGeometry(const std::string & fileName);
      // check if the geometry has been read in
      bool isEmpty(void) const;
      // return the geometryFileName
      const std::string & getGeometryFileName(void) const;
    
      std::vector<Node> nodes;
      std::vector<Segment> segments;
      size_t numCompartments;
      std::map<std::string, size_t> tags;
      std::set<std::string> compartmentNames;
      double surfaceArea;         // mm^2
      double volume;              // mm^3
    
    protected:
      // build geometry member data from information in geometry file
      void buildGeometry(void);
      // activate or inactivate a tag for subsequent segments, add to set tags
      void handleTags(std::string & line, std::set<std::string> & segTags);
      // define a node
      void defineNode(const SplitLine & splitLine);
      // define a segment with a circular cross-section
      void defineCircularSegment(const SplitLine & splitLine,
                                 const std::set<std::string> & segTags);
      // define a segment with an eliptical cross-section
      void defineElipticalSegment(const SplitLine & splitLine,
                                  const std::set<std::string> & segTags);
      // update summary member data with information from the last segment
      void updateGeometrySummary(void);
      
      std::vector<char> geometryFileContent;
      std::string geometryFileName;
  };
  
  
  
  struct Node
  {
    Node(void) { locSpecified = false; }
    Node(const double & x0, const double & y0, const double & z0) :
      x(x0), y(y0), z(z0) { locSpecified = true;}
    
    double x;   // um
    double y;   // um
    double z;   // um
    std::vector<size_t> segments;
    std::vector<size_t> compartments;
    
    bool locSpecified;
  };
  
  
  
  class Segment
  {
    public:
      Segment(void) {}
      Segment(const size_t & node0,
              const size_t & node1,
              const size_t & numCompartments,
              const double & length,
              const double & radius,
              const std::set<std::string> & tags,
              const std::map<std::string, size_t> geoTags);
      Segment(const size_t & node0,
              const size_t & node1,
              const size_t & numCompartments,
              const double & length,
              const double & semiMajor,
              const double & semiMinor,
              const double & elipseAngle,
              const std::set<std::string> & tags,
              const std::map<std::string, size_t> geoTags);
      
      size_t node0;
      size_t node1;
      size_t numCompartments;
      double length;          // um
      double semiMajor;       // um
      double semiMinor;       // um
      double elipseAngle;     // radians
      
      double crossSectionArea;// um^2
      double surfaceArea;     // mm^2
      double volume;          // mm^3
      // hold a set of tags that apply to every compartment in the Segment:
      std::set<std::string> tags;
      // hold a vector of valid names for each compartment in the Segment:
      std::vector<std::vector<std::string> > compartmentNames;
      
    protected:
      // assign names and aliases for each compartment
      void assignNames(const std::map<std::string, size_t> geoTags);
  };
};



 
 
#endif
