#include "../include/NeuronGeometry.h"
#include "../include/NeuronException.h"
#include "../include/constants.h"

#include <sstream>
#include <algorithm>



using namespace NeuronGeometry;
using namespace std;



//////////////////////////////  Helper functions  /////////////////////////////



template <class A, class B>
const B &
get(const map<A,B> & mapAB, const A & a)
{
  typename map<A,B>::const_iterator mapItr = mapAB.find(a);
  if(mapItr == mapAB.end())
    // this should never happen, so don't worry too much about detailed message
    throw NeuronGeometryException("Invalid map key");
  
  return mapItr->second;
}



//////////////////////  Segment public member functions  //////////////////////



Segment::Segment(const size_t & node0New,
                 const size_t & node1New,
                 const size_t & numCompartmentsNew,
                 const double & lengthNew,
                 const double & radiusNew,
                 const set<string> & tagsNew,
                 const std::map<std::string, size_t> geoTags) :
  node0(node0New), node1(node1New),
  numCompartments(numCompartmentsNew), length(lengthNew), semiMajor(radiusNew),
  semiMinor(radiusNew), elipseAngle(0), tags(tagsNew)
{
  crossSectionArea = Pi * semiMajor * semiMinor; // um^2
  surfaceArea = 2.0e-6 * Pi * radiusNew * length;   // mm^2
  volume = 1.0e-9 * crossSectionArea * length;   // mm^3
  
  // assign names and aliases for each compartment
  assignNames(geoTags);
}



Segment::Segment(const size_t & node0New,
                 const size_t & node1New,
                 const size_t & numCompartmentsNew,
                 const double & lengthNew,
                 const double & semiMajorNew,
                 const double & semiMinorNew,
                 const double & elipseAngleNew,
                 const set<string> & tagsNew,
                 const std::map<std::string, size_t> geoTags) :
  node0(node0New), node1(node1New),
  numCompartments(numCompartmentsNew), length(lengthNew),
  semiMajor(semiMajorNew), semiMinor(semiMinorNew),
  elipseAngle(elipseAngleNew), tags(tagsNew)
{
  // Approximate formula for perimeter of elipse, by David Cantrell
  // (Accurate to within 0.02%)
  const double s = 0.825056; // optimal power for arbitrary elipse eccentricity
  const double perimeter =
    4 * (semiMajor + semiMinor) - 2 * (4 - Pi) * semiMajor * semiMinor * pow(
      0.5 * (pow(semiMajor, s) + pow(semiMinor, s)), -1.0/s );
  crossSectionArea = Pi * semiMajor * semiMinor; // um^2
  surfaceArea = 1.0e-6 * perimeter * length;     // mm^2
  volume = 1.0e-9 * crossSectionArea * length;   // mm^3
  
  // assign names and aliases for each compartment
  assignNames(geoTags);
}



////////////////////  Segment protected member functions  /////////////////////



void
Segment::assignNames(const std::map<std::string, size_t> geoTags)
{
  // add names for each compartment
  
  // resize compartmentNames. It will hold a vector of valid names for
  //   each compartment in the Segment
  compartmentNames.resize(numCompartments);
  // each compartment will have:
  //    1 name (from being the nth compartment in the geometry)
  //    numTags additional names (being the nth compartment to have that tag)
  const size_t numNames = 1 + tags.size();
  
  vector<vector<string> >::iterator compItr = compartmentNames.begin();
  for(size_t compNum = 0; compNum < numCompartments; ++compItr, ++compNum) {
    // loop through compartments, adding names appropriately
    
    // resize the list of names for this compartment
    compItr->resize(numNames);
    
    // the first name is the true name, just the compartment number
    vector<string>::iterator nameItr = compItr->begin();
    *nameItr = numToString(compNum + get(geoTags, string("*")));
    
    // add on any aliases
    for(set<string>::const_iterator tagItr = tags.begin();
        tagItr != tags.end();
        ++tagItr) {
      ++nameItr;
      // create an alias for each tag, as [tag]_[#] where [#] is the number
      //  of previous compartments with that tag
      *nameItr = *tagItr + "_" + numToString(compNum + get(geoTags, *tagItr));
    }
  }
}



/////////////////////  Geometry public member functions  //////////////////////



void
Geometry::readGeometry(const string & fileName)
{
  ifstream fileIn (fileName.c_str());
  if(fileIn.fail())
    throw NeuronGeometryException("Error opening geometry data file: " +
                                  fileName);
                                  
  // store the geometry file name for later use
  geometryFileName = fileName;
  
  // figure out the length of the file
  fileIn.seekg(0, ios::end);
  const streampos endPos = fileIn.tellg();
  fileIn.seekg(0, ios::beg);
  const streampos beginPos = fileIn.tellg();
  const size_t fileLength = endPos - beginPos;
  
  // allocate space in geometryFileContent
  geometryFileContent.resize(fileLength);
  
  // read in the file data
  fileIn.read(&geometryFileContent[0], fileLength);
  
  // close the file
  fileIn.close();
  
  // build the geometry member data
  buildGeometry();
}



bool
Geometry::isEmpty(void) const
{
  // check if the geometry has been read in
  
  // the only way to specify geometryFileName is via readGeometry, so just
  // check that the file name has been specified
  return geometryFileName.empty();
}



const std::string &
Geometry::getGeometryFileName(void) const
{
  return geometryFileName;
}



////////////////////  Geometry protected member functions  ////////////////////



void
Geometry::buildGeometry(void)
{
  // initialize summary member data
  numCompartments = 0;
  tags["*"] = 0;
  surfaceArea = 0;
  volume = 0;

  // define some variables for parsing the geometry info
  string line;
  vector<string> splitLine;
  set<string> segmentTags;
  size_t lineNum = 0;

  try {
    // read from the vector<char> geometryFileContent just as though it was an
    //  ordinary file
    //stringstream fileIn (geometryFileContent);
    // try this (to avoid copying geometryFileContent)
    stringstream fileIn;
    fileIn.rdbuf()->pubsetbuf(&geometryFileContent[0],
                              geometryFileContent.size());

    while(fileIn.good()) {
      // get the next nontrivial line (removing comments that start with #)
      line = getNextLine(fileIn, lineNum, '#');
      
      // split line up into words
      splitLine = split(line);
      
      if(splitLine.empty())
        // empty line or just comment
        continue;
      
        switch(splitLine.size()) {
          case 1:  // start or end applying a tag
            handleTags(line, segmentTags);
            break;
            
          case 3:  // define a node
            defineNode(splitLine);
            break;
            
          case 5:  // define a segment with a circular cross-section
            defineCircularSegment(splitLine, segmentTags);
            break;
            
          case 7:  // define a segment with an eliptical cross-section
            defineElipticalSegment(splitLine, segmentTags);
            break;
            
          default: // syntax error
            throw NeuronGeometryException("Syntax error");
        }
    }
    
    // stringstream objects do not need to be closed
  }
  catch( exception & e ) {
    // caught an exception, add some detail and rethrow it
    const string errorStr =
      "Error reading geometry data from file: \"" + geometryFileName +
      "\" line " + numToString(lineNum) + ": " + e.what();
    throw NeuronGeometryException(errorStr);
  }
}



void
Geometry::handleTags(string & line, set<string> & segTags)
{
  // first remove whitespace
  removeWhiteSpace(line);

  if(line.length() < 3 || line[0] != '<' || *(line.rbegin()) != '>')
    // this isn't the beginning or end of a tag, it's an error
    throw NeuronGeometryException("Syntax error");

  if(line[1] == '/') {
    // this is the end of a tag, inactivate it
    
    // extract the actual tag
    string tag = string(line.begin() + 2, line.end() - 1);
    if(tag.empty())
      throw NeuronGeometryException("Empty tag");
    if(tag == "*")
      throw NeuronGeometryException("Can't tag with *");
  
    if(segTags.count(tag) == 0)
      // this tag not currently active, can't remove it
      throw NeuronGeometryException("Tag \"" + tag +"\" not currently active");
    
    // remove the tag, as requested
    segTags.erase(tag);
  }
  else {
    // this is the beginning of a tag, activate it
  
    // extract the actual tag
    string tag = string(line.begin() + 1, line.end() - 1);
    if(tag.empty())
      throw NeuronGeometryException("Empty tag");
    if(tag == "*")
      throw NeuronGeometryException("Can't tag with *");

    if(segTags.count(tag) != 0)
      // this tag is already active, can't add it
      throw NeuronGeometryException("Tag \"" + tag + "\" is already active");

    // add the tag, as requested
    segTags.insert(tag);
    
    // add this tag to the Geometry's list of all tags
    if(tags.count(tag) == 0)
      tags[tag] = 0;
  }
}



void
Geometry::defineNode(const SplitLine & splitLine)
{
  // define a node
  const double x = stringToDouble(splitLine[0]);
  const double y = stringToDouble(splitLine[1]);
  const double z = stringToDouble(splitLine[2]);
  
  for(vector<Node>::iterator nodeItr = nodes.begin();
      nodeItr != nodes.end();
      ++nodeItr)
    if(!nodeItr->locSpecified) {
      // specify the location of the first node that isn't specified
      nodeItr->x = x;
      nodeItr->y = y;
      nodeItr->z = z;
      nodeItr->locSpecified = true;
      return;
    }
    
  // no unspecified nodes, create a new node
  nodes.push_back( Node (x,y,z) );
}



void
Geometry::defineCircularSegment(const SplitLine & splitLine,
                                const set<string> & segmentTags)
{
  // define a segment with a circular cross-section
  const size_t node0 = stringToInt(splitLine[0]);
  const size_t node1 = stringToInt(splitLine[1]);
  const size_t segmentCompartments = stringToInt(splitLine[2]);
  const double length = stringToDouble(splitLine[3]);
  const double radius = stringToDouble(splitLine[4]);

  // define and add the segment
  segments.push_back(
    Segment (node0, node1, segmentCompartments, length, radius, segmentTags,
             tags) );
  
  // make sure there are enough nodes in the geometry
  const size_t numNeededNodes = max(node0, node1) + 1;
  if(numNeededNodes > nodes.size())
    nodes.resize(numNeededNodes);

  // add this segment to the list of segments in node0, node1
  const size_t segmentNum = segments.size();
  nodes[node0].segments.push_back(segmentNum);
  nodes[node0].compartments.push_back(numCompartments);
  numCompartments += segmentCompartments;
  nodes[node1].segments.push_back(segmentNum);
  nodes[node1].compartments.push_back(numCompartments - 1);
  // update summary member data (e.g. surface area, volume, etc)
  updateGeometrySummary();
}



void
Geometry::defineElipticalSegment(const SplitLine & splitLine,
                                 const set<string> & segmentTags)
{
  // define a segment with an eliptical cross-section
  const size_t node0 = stringToInt(splitLine[0]);
  const size_t node1 = stringToInt(splitLine[1]);
  const size_t segmentCompartments = stringToInt(splitLine[2]);
  const double length = stringToDouble(splitLine[3]);
  const double semiMajor = stringToDouble(splitLine[4]);
  const double semiMinor = stringToDouble(splitLine[5]);
  const double elipseAngle = stringToDouble(splitLine[6]);
  
  // define and add the segment
  segments.push_back(
    Segment (node0, node1, segmentCompartments, length, semiMajor,
             semiMinor, elipseAngle, segmentTags, tags) );
  
  // make sure there are enough nodes in the geometry
  const size_t numNeededNodes = max(node0, node1) + 1;
  if(numNeededNodes > nodes.size())
    nodes.resize(numNeededNodes);
                   
  // add this segment to the list of segments in node0, node1
  const size_t segmentNum = segments.size();
  nodes[node0].segments.push_back(segmentNum);
  nodes[node0].compartments.push_back(numCompartments);
  numCompartments += segmentCompartments;
  nodes[node1].segments.push_back(segmentNum);
  nodes[node1].compartments.push_back(numCompartments - 1);
  
  // update summary member data (e.g. surface area, volume, etc)
  updateGeometrySummary();
}



void
Geometry::updateGeometrySummary(void)
{
  // update summary member data with information from the last segment
  
  // make iterator to last segment
  const vector<Segment>::const_iterator segItr = segments.end() - 1;
  
  // update surface area
  surfaceArea += segItr->surfaceArea;
  // update volume
  volume += segItr->volume;
  // don't update numCompartments, the defineSegment function does that
  //numCompartments += segItr->numCompartments;
  
  // update number of compartments for each active tag
  set<string>::const_iterator tagItr = segItr->tags.begin();
  const set<string>::const_iterator end = segItr->tags.end();
  for(; tagItr != end; ++tagItr)
    tags[*tagItr] += segItr->numCompartments;
  // update the number of compartments for the "*" tag
  tags["*"] += segItr->numCompartments;
  
  // update the set of compartment names
  for(vector<vector<string> >::const_iterator compNameItr =
        segItr->compartmentNames.begin();
      compNameItr != segItr->compartmentNames.end();
      ++compNameItr) {
    // loop through compartments in segment
  
    // insert all the aliases for this compartment
    compartmentNames.insert(compNameItr->begin(), compNameItr->end());
  }
}
