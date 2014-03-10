#include "../include/NeuronFitInfo.h"
#include "../include/constants.h"
#include "../include/NeuronException.h"

#include <algorithm>
#include <limits>
#include <cmath>



using namespace std;



/////////////////////////////  Helper functions  //////////////////////////////



inline bool
endsSooner(const Trace & t1, const Trace & t2)
{
  const double tF1 = t1.deltaT * ((double)t1.trace.size() - 1.0);
  const double tF2 = t2.deltaT * ((double)t2.trace.size() - 1.0);
  return (tF1 < tF2);
}



//////////////////  NeuronFitInfo public member functions  ////////////////////



void
NeuronFitInfo::readStartupFile(const string & fileName)
{
  ifstream fileIn (fileName.c_str());
  if(fileIn.fail())
    throw FitInfoException("Error opening startup file: " + fileName);
  
  // store the startup file name for later use
  startupFileName = fileName;
  
  // figure out the length of the file
  fileIn.seekg(0, ios::end);
  const streampos endPos = fileIn.tellg();
  fileIn.seekg(0, ios::beg);
  const streampos beginPos = fileIn.tellg();
  const size_t fileLength = endPos - beginPos;
  
  // allocate space in geometryFileContent
  startupFileText.resize(fileLength);
  
  // read in the file data
  fileIn.read(&startupFileText[0], fileLength);
  
  // close the file
  fileIn.close();
  
  // electrophys data still needs to be loaded in
  dataIsLoaded = false;
  
  // parse the startup information and build member data
  buildFitInfo();
}



void
NeuronFitInfo::outputRecordedTraces(const std::string & fileName) const
{
  // output recorded traces to outFile
  
  // make a vector of indices to Record traces
  vector<size_t> recordInds;
  
  // first look to find the number of traces that are Record or Fit
  size_t numTraces = 0;
  for(const Trace & trace : electrophysData) {
    if(trace.traceAction == Trace::Record ||
       trace.traceAction == Trace::Fit) {
      ++numTraces;
    }
  }
  
  if(numTraces == 0)
    // quietly return
    return;
  
  // open fileName for writing
  ofstream out (fileName.c_str());
  if(!out.good())
    throw FitInfoException("Unable to open " + fileName + " for writing");
  
  // write the header
  out << "# number of simulated traces\n";
  out << numTraces << '\n';
  out << "# name units numT deltaT\n";
  for(const Trace & trace : electrophysData) {
    // loop through Record / Fit traces and output some header info

    // use the targetName if the descriptive name is empty
    const string traceName = trace.targetDescriptiveName.empty() ?
      trace.targetName : trace.targetDescriptiveName;

    if(trace.traceAction == Trace::Record) {
      // write the header info for this trace
      out << traceName << ' ' << trace.units << ' '
          << trace.trace.size() << ' ' << trace.deltaT << '\n';
    }
    else if(trace.traceAction == Trace::Fit) {
      // write the header info for this trace
      out << "error_" << traceName << ' ' << trace.units
          << ' ' << trace.errorTrace.size() << ' ' << trace.deltaT
          << '\n';
    }
  }
  
  // calculate a few widths to ensure all the significant info in the trace
  //  is output
  const streamsize doublePrecision = numeric_limits<double>::digits10 + 1;
  const streamsize expWidth = (streamsize)ceil(
    log10((double)numeric_limits<double>::max_exponent10));
  const streamsize doubleWidth = doublePrecision + expWidth + 4;
  
  // write the trace data
  for(const Trace & trace : electrophysData) {
    // loop through Record / Fit traces and output trace data / error data

    // use the targetName if the descriptive name is empty
    const string traceName = trace.targetDescriptiveName.empty() ?
      trace.targetName : trace.targetDescriptiveName;

    if(trace.traceAction == Trace::Record) {
      // mark the beginning of the trace
      out << left << "# " + traceName + "\n";
      // write the trace data out
      for(const double & traceValue : trace.trace) {
        out.precision(doublePrecision);
        out.width(doubleWidth);
        out << left << traceValue << '\n';
      }
    }
    else if(trace.traceAction == Trace::Fit) {
      // mark the beginning of the trace
      out << left << "# error_" << traceName + "\n";
      // write the error data out
      for(const double & traceValue : trace.errorTrace) {
        out.precision(doublePrecision);
        out.width(doubleWidth);
        out << left << traceValue << '\n';
      }
    }
    else {
      // don't output the trace
      continue;
    }
  }
  
  out.close();
}



/////////////////  NeuronFitInfo protected member functions  //////////////////



void
NeuronFitInfo::buildFitInfo(void)
{
  dataIsLoaded = electrophysData.size() > 0;
  
  if(geometry.isEmpty()) {
    // need to get the geometry before we can properly build fit info
    scanForGeometry();
  }

  if(!dataIsLoaded) {
    // set default value for tFinal
    tFinal = Inf;
  }
  
  // read from the vector<char> startupFileText just as though it was an
  //  ordinary file
  //stringstream fileIn (startupFileText);
  // try this (to avoid copying text content)
  stringstream fileIn;
  fileIn.rdbuf()->pubsetbuf(&startupFileText[0],
                            startupFileText.size());
  
  // define a few variables to help read in the file
  string line;
  string keyword;
  SplitLine splitLine;
  size_t lineNum = 0;
  while(fileIn.good()) {
    // get the next nontrivial line (removing comments that start with #)
    line = getNextLine(fileIn, lineNum, '#');

    // split line up into words
    splitLine = split(line);
    
    // convert the keyword (first word) to lower case
    keyword = splitLine[0];
    transform(keyword.begin(), keyword.end(), keyword.begin(), ::tolower);
    
    try {
      if( keyword == "time"   )
        // specify the fit/simulation time information
        addTime(splitLine);
      else if( keyword == "geometry" )
        // specify the name of the file containing geometry information
        addGeometry(splitLine);
      else if( keyword == "channeldir" )
        // specify directory of channel model shared object files
        addChannelDir(splitLine);
      else if( keyword == "record" )
        // record a simulated trace
        addRecordTrace(splitLine);
      else if( keyword == "clamp" )
        // clamp a value to a specified trace
        addClampTrace(splitLine);
      else if( keyword == "fit" )
        // produce an error estimate from deviations from a specified trace
        addFitTrace(splitLine);
      else if( keyword == "<channel>" )
        // add a series of ion channels to the model
        addChannels(fileIn, lineNum, splitLine);
      else if( keyword == "<connect>" )
        // add a series of synapses to the model
        addConnections(fileIn, lineNum, splitLine);
      else if( keyword == "channel" )
        // add an ion channel to the model
        addChannel(splitLine);
      else if( keyword == "connect" )
        // add a synapse to the model
        addConnection(splitLine);
      else if( keyword == "<parameter>" )
        // specify a series of parameters of the model
        addParameters(fileIn, lineNum, splitLine);
      else if( keyword == "parameter" )
        // specify a parameter of the model
        addParameter(splitLine);
      else
        // invalid keyword
        throw FitInfoException("Invalid keyword \"" + splitLine[0] + "\"");
    }
    catch( exception & e ) {
      // caught an exception, add some detail and rethrow it
      const string errorStr =
        "Error reading startup file: " + startupFileName + "\n" +
        "\tError in line " + numToString(lineNum) + ": " + e.what();
      throw FitInfoException(errorStr);
    }
  }

  // Done reading startup file
  //  (don't need to close stringstream)
  
  // wrap up by actually loading any needed data and allocating space for
  // simulated traces
  loadTraces();

  // sort the traces so that traces ending soonest are first
  electrophysData.sort(endsSooner);
}



void
NeuronFitInfo::scanForGeometry(void)
{
  // read from the vector<char> startupFileText just as though it was an
  //  ordinary file
  //stringstream fileIn (startupFileText);
  // try this (to avoid copying text content)
  stringstream fileIn;
  fileIn.rdbuf()->pubsetbuf(&startupFileText[0],
                            startupFileText.size());
  size_t lineNum = 0;
  
  while(fileIn.good()) {
     // get the next nontrivial line (removing comments that start with #)
    const string line = getNextLine(fileIn, lineNum, '#');
    
    // split line up into words
    const SplitLine splitLine = split(line);
    
    // convert the keyword (first word) to lower case
    string keyword = splitLine[0];
    transform(keyword.begin(), keyword.end(), keyword.begin(), ::tolower);
    
    if( keyword == "geometry" )
      try {
        // specify the name of the file containing geometry information
        addGeometry(splitLine);
      }
      catch( exception & e ) {
        // caught an exception, add some detail and rethrow it
        const string errorStr =
          "Error reading startup file: " + startupFileName + "\n" +
          "\tError in line " + numToString(lineNum) + ": " + e.what();
        throw FitInfoException(errorStr);
      }
  }
  
  // Done reading startup file
  //  (don't need to close stringstream)
}



void
NeuronFitInfo::addTime(const SplitLine & splitLine)
{
  // specify the fit/simulation time information

  if(splitLine.size() != 2)
    throw FitInfoException("Invalid use of \"time\" keyword");

  const double tF = stringToDouble(splitLine[1]);
  if(!(tF > 0 ))
    throw FitInfoException("Invalid tFinal");
  if(tF < tFinal)
    tFinal = tF;
}



void
NeuronFitInfo::addGeometry(const SplitLine & splitLine)
{
  // specify the name of the file containing geometry information and
  // read the geometry information in

  if(splitLine.size() != 2)
    throw FitInfoException("Invalid use of \"geometry\" keyword");
  
  if(geometry.isEmpty()) {
    // geometry hasn't yet been loaded, get it from the geometry file specified
    // in splitLine[1]
    geometry.readGeometry(userExpand(splitLine[1]));
  } else {
    // geometry is already loaded, just check to make sure that only one
    //  geometry file was supplied
    if(userExpand(splitLine[1]) != geometry.getGeometryFileName())
      throw FitInfoException("Multiple geometry files supplied.");
  }
}



void
NeuronFitInfo::addChannelDir(const SplitLine & splitLine)
{
  // specify directory of channel model shared object files
  
  if(splitLine.size() != 2)
    throw FitInfoException("Invalid use of \"channelDir\" keyword");
  
  // add this directory to the set
  injectorDirs.insert(userExpand(splitLine[1]));
}
        


inline void
NeuronFitInfo::checkValidity(const Trace & trace) const
{
  // check some trace values to make sure the trace is okay

  // ensure the target is valid
  const SplitLine splitLine = split(trace.targetName, "_");
  if(!splitLine.empty()) {
    // The trace targets a compartment-specific quantity
    const string & compartmentName = splitLine.back();
    if(!geometry.compartmentNames.count(compartmentName)) {
      throw FitInfoException("Trace target has invalid compartment: "
                             + trace.targetName);
    }
  }
  
  // ensure that fitTau and waitTime are sensible
  if(trace.traceAction == Trace::Fit) {
    if(trace.fitTau <= 0) {
      throw FitInfoException("Trace fitTau must be strictly > 0");
    }
    if(trace.waitTime < 0) {
      throw FitInfoException("Trace waitTime must be >= 0");
    }
  }
}



void
NeuronFitInfo::addRecordTrace(const SplitLine & splitLine)
{
  // record a simulated trace
  
  // line should be of form:
  //    record target dT
  // or:
  //    record target_compartmentName dT
  // (target may have underscores in name)
  
  if(splitLine.size() != 3)
    throw FitInfoException(
      "Incorrect number of arguments with \"record\" keyword");
  
  
  const string & target = splitLine[1];
  const double dT = stringToDouble(splitLine[2]);
  if(!(dT > 0))
    throw FitInfoException("Invalid trace dT for trace with target: " +target);

  // Note:
  // this is different from "clamp" or "fit" in that someone might for some
  // reason want multiple recordings of the same thing, so the presence of
  // duplicates isn't an error condition.


  if(dataIsLoaded) {
    // find an example of a matching trace
    for(const Trace & trace : electrophysData) {
      // loop through loaded data
      if(trace.matches(target, Trace::Record, dT)) {
        // if trace is a recording Trace that matches the target and dT,
        // we've found a match
        
        // check some to make sure the trace is okay
        checkValidity(trace);
        return;
      }
    }
    // no match, very bizarre error, should never happen
    throw FitInfoException("Data loaded but missing!");
  }
  else {
    // create a new trace and put it at the back of the list
    electrophysData.push_back(Trace());

    // make iterator to it
    Trace & trace = electrophysData.back();

    // set member data
    trace.targetName = target;
    trace.traceAction = Trace::Record;
    trace.deltaT = dT;
    
    // check some to make sure the trace is okay
    checkValidity(trace);
  }
}



void
NeuronFitInfo::addClampTrace(const SplitLine & splitLine)
{
  // clamp a value to a specified trace  
  
  // line should be of form:
  //    clamp target fileName traceNumber
  // or:
  //    clamp target_compartmentName fileName traceNumber
  // (target may have underscores in name)
  
  if(splitLine.size() != 4)
    throw FitInfoException(
      "Incorrect number of arguments with \"clamp\" keyword");
  
  const string & target = splitLine[1];
  const string & fileName = splitLine[2];
  const size_t traceNumber = stringToSize(splitLine[3]);
  
  // search for any traces that match the requested one (could be duplicates,
  // or could be that data is already loaded)
  size_t numMatches = 0;
  for(const Trace & trace : electrophysData) {
    // loop through traces
    if(trace.matches(target, Trace::Clamp)) {
      // existing trace matches the one we're trying to add
      numMatches++;  // inc the number of matches
      // check some to make sure the trace is okay
      checkValidity(trace);
    }
  }
  
  if(numMatches >= 2)
    // duplicate trace
    throw FitInfoException("Duplicate clamp trace");
  
  if(dataIsLoaded) {
    if(numMatches == 0)
      // no match, very bizarre error, should never happen
      throw FitInfoException("Data loaded but missing!");
  }
  else {
    if(numMatches > 0)
      // duplicate trace
      throw FitInfoException("Duplicate clamp trace");

    // schedule trace to be loaded
    scheduleLoadTrace(userExpand(fileName), traceNumber);
  
    // trace is put at end of list
    Trace & trace = electrophysData.back();

    // set member data
    trace.traceAction = Trace::Clamp;
    trace.targetName = target;

    // check some to make sure the trace is okay
    checkValidity(trace);
  }
}



void
NeuronFitInfo::addFitTrace(const SplitLine & splitLine)
{
  // produce an error estimate from deviations from a specified trace
  
  // line should be of form:
  //  fit target fileName traceNum fitTau waitTime
  // or:
  //  fit target_compartmentName fileName traceNum fitTau waitTime
  // (target may have underscores in name)
  if(splitLine.size() != 6)
    throw FitInfoException(
      "Incorrect number of arguments with \"fit\" keyword");
  
  const string & target = splitLine[1];
  const string & fileName = splitLine[2];
  const size_t traceNumber = stringToSize(splitLine[3]);
  const double fitTau = stringToDouble(splitLine[4]);
   const double waitTime = stringToDouble(splitLine[5]);
  
  // search for any traces that match the requested one (could be duplicates,
  // or could be that data is already loaded)
  size_t numMatches = 0;
  for(const Trace & trace : electrophysData) {
    // loop through traces
    if(trace.matches(target, Trace::Fit)) {
      // existing trace matches the one we're trying to add
      ++numMatches;  // inc the number of matches
      // check some to make sure the trace is okay
      checkValidity(trace);
    }
  }
  if(numMatches >= 2)
    // duplicate trace
    throw FitInfoException("Duplicate fit trace (data preloaded)");
  
  if(dataIsLoaded) {
    if(numMatches == 0)
      // no match, very bizarre error, should never happen
      throw FitInfoException("Data loaded but missing!");
  }
  else {
    if(numMatches > 0)
      // duplicate trace
      throw FitInfoException("Duplicate fit trace (data not preloaded)");

    // schedule trace to be loaded
    scheduleLoadTrace(userExpand(fileName), traceNumber);
  
    // trace is put at end of list
    Trace & trace = electrophysData.back();

    // set member data
    trace.traceAction = Trace::Fit;
    trace.targetName = target;
    trace.fitTau = fitTau;
    trace.waitTime = waitTime;
    // check some to make sure the trace is okay
    checkValidity(trace);
  }
}



void
NeuronFitInfo::addChannels(stringstream & fileIn, size_t & lineNum,
                           SplitLine & splitLine)
{
  // add a series of ion channels to the model
  
  // the <channel> keyword should be alone on the line
  if(splitLine.size() != 1)
    throw FitInfoException("Invalid use of <channel> keyword");

  // define line, keyword
  string line;
  string keyword = splitLine.front();
  
  while(fileIn.good() && keyword != "</channel>") {
    // get the next non-trivial line (removing comments indicated by #)
    line = getNextLine(fileIn, lineNum, '#');
    
    // split line up into words
    splitLine = split(line);
    
    // convert the keyword (first word) to lower case
    keyword = splitLine.front();
    transform(keyword.begin(), keyword.end(), keyword.begin(), ::tolower);

    if(keyword != "</channel>") {
      // insert the "channel" keyword
      splitLine.insert(splitLine.begin(), "channel");
      // add a channel with the resulting line
      addChannel(splitLine);
    }
  }
  
  if(keyword != "</channel>")
    throw FitInfoException("<channel> keyword never closed with </channel>");
  else if(splitLine.size() != 1)
    throw FitInfoException("Invalid use of <channel> keyword");
}



void
NeuronFitInfo::addConnections(stringstream & fileIn, size_t & lineNum,
                           SplitLine & splitLine)
{
  // add a series of ion channels to the model
  
  // the <connect> keyword should be alone on the line
  if(splitLine.size() != 1)
    throw FitInfoException("Invalid use of <connect> keyword");

  // define line, keyword
  string line;
  string keyword = splitLine.front();
  
  while(fileIn.good() && keyword != "</connect>") {
    // get the next non-trivial line (removing comments indicated by #)
    line = getNextLine(fileIn, lineNum, '#');
    
    // split line up into words
    splitLine = split(line);
    
    // convert the keyword (first word) to lower case
    keyword = splitLine.front();
    transform(keyword.begin(), keyword.end(), keyword.begin(), ::tolower);

    if(keyword != "</connect>") {
      // insert the "connect" keyword
      splitLine.insert(splitLine.begin(), "connect");
      // add a channel with the resulting line
      addConnection(splitLine);
    }
  }
  
  if(keyword != "</connect>")
    throw FitInfoException("<connect> keyword never closed with </connect>");
  else if(splitLine.size() != 1)
    throw FitInfoException("Invalid use of <connect> keyword");
}



set<string>
NeuronFitInfo::getCompartmentNamesByTag(const string & tag) const
{
  // return a set of compartment names that tag applies to (either as a tag or
  // as a compartment name)
  
  set<string> taggedCompartments;
  if(tag == "*") {
    // every compartment matches
    for(const NeuronGeometry::Segment & segment : geometry.segments) {
      for(const vector<string> & compartmentNames : segment.compartmentNames) {
        taggedCompartments.insert(compartmentNames.front());
      }
    }
    return taggedCompartments;
  }
  
  for(const NeuronGeometry::Segment & segment : geometry.segments) {
    if(segment.tags.count(tag)) {
      // all the compartments have this tag
      for(const vector<string> & compartmentNames : segment.compartmentNames) {
        taggedCompartments.insert(compartmentNames.front());
      }
    }
    else {
      // check compartment by compartment to see if any name matches
      for(const vector<string> & compartmentNames : segment.compartmentNames) {
        for(const string & alias : compartmentNames) {
          if(alias == tag) {
            // this compartment matches
            taggedCompartments.insert(compartmentNames.front());
            break;
          }
        }
      }
    }
  }
  return taggedCompartments;
}


void
NeuronFitInfo::addChannel(const SplitLine & splitLine)
{
  // add an ion channel to the model
  
  // line should be of form:
  //    channel channelName tagName
  if(splitLine.size() != 3)
    throw FitInfoException("Invalid use of \"channel\" keyword");
  
  const string & channelName = splitLine[1];
  const string & tagName = splitLine[2];
  
  // get a list of compartment names corresponding to compartments that have
  // tagName (either as a tag or a compartment name)
  const auto taggedCompartmentNames = getCompartmentNamesByTag(tagName);
  if(taggedCompartmentNames.empty()) {
    throw FitInfoException("Invalid tag or compartment name \"" + tagName
                           + "\"");
  }
  
  // add appropriate channel to onePointMap for all the specified compartments
  for(const string & compartmentName : taggedCompartmentNames) {
    onePointInfo.insert( {channelName, compartmentName} );
  }  
}



void
NeuronFitInfo::addConnection(const SplitLine & splitLine)
{
  // add a TwoPointInjector to the model
  
  // line should be of form:
  //    connect synapseName preTag postTag
  if(splitLine.size() != 4)
    throw FitInfoException("Invalid use of \"connect\" keyword");
  
  const string & synapseName = splitLine[1];
  const string & preTag = splitLine[2];
  const string & postTag = splitLine[3];
  
  // get a list of compartment names corresponding to compartments that have
  // preTag (either as a tag or a compartment name)
  const auto preCompartmentNames = getCompartmentNamesByTag(preTag);
  if(preCompartmentNames.empty()) {
    throw FitInfoException("Invalid tag or compartment name \"" + preTag
                           + "\"");
  }
  // get a list of compartment names corresponding to compartments that have
  // preTag (either as a tag or a compartment name)
  const auto postCompartmentNames = getCompartmentNamesByTag(postTag);
  if(postCompartmentNames.empty()) {
    throw FitInfoException("Invalid tag or compartment name \"" + postTag
                           + "\"");
  }
  
  // add appropriate channel to onePointMap for all the specified compartments
  for(const string & preCompName : preCompartmentNames) {
    for(const string & postCompName : postCompartmentNames) {
      twoPointInfo.insert( {synapseName, preCompName, postCompName} );
    }
  }  
}


void
NeuronFitInfo::addParameters(stringstream & fileIn, size_t & lineNum,
                           SplitLine & splitLine)
{
  // // specify a series of parameters of the model
  
  // the <parameter> keyword should be alone on the line
  if(splitLine.size() != 1)
    throw FitInfoException("Invalid use of <channel> keyword");

  // define line, keyword
  string line;
  string keyword = splitLine.front();
  
  while(fileIn.good() && keyword != "</parameter>") {
    // get the next non-trivial line (removing comments indicated by #)
    line = getNextLine(fileIn, lineNum, '#');
    
    // split line up into words
    splitLine = split(line);
        
    // convert the keyword (first word) to lower case
    keyword = splitLine.front();
    transform(keyword.begin(), keyword.end(), keyword.begin(), ::tolower);

    if(keyword != "</parameter>") {
      // insert the "parameter" keyword
      splitLine.insert(splitLine.begin(), "parameter");
      try {
        // add a channel with the resulting line
        addParameter(splitLine);
      }
      catch(exception & e) {
        throw FitInfoException("Invalid parameter");
      }
    }
  }
  
  if(keyword != "</parameter>")
    throw FitInfoException("<channel> keyword never closed with </parameter>");
  else if(splitLine.size() != 1)
    throw FitInfoException("Invalid use of <parameter> keyword");
}



void
NeuronFitInfo::addParameter(const SplitLine & splitLine)
{
  // specify a parameter of the model
  
  string paramName;
  switch(splitLine.size()) {
    case 3: {  // a constant value
      const double value = stringToDouble(splitLine[2]);
      ParameterDescription param (splitLine[1], value);
      parameterDescriptionList.add(param);
      break;
    }
    case 4: {  // parameter ranges between minVal, maxVal
      const double minVal = stringToDouble(splitLine[2]);
      const double maxVal = stringToDouble(splitLine[3]);
      ParameterDescription param (splitLine[1], minVal, maxVal);
      parameterDescriptionList.add(param);
    
      break;
    }
    case 6: {  // parameter ranges between minVal, maxVal; initial population
               // ranges between startMin and startMax
      const double minVal = stringToDouble(splitLine[2]);
      const double maxVal = stringToDouble(splitLine[3]);
      const double startMin = stringToDouble(splitLine[4]);
      const double startMax = stringToDouble(splitLine[5]);
      ParameterDescription param (splitLine[1], minVal, maxVal,
                                  startMin, startMax);
      parameterDescriptionList.add(param);
      break;
    }

    default:
      throw FitInfoException("Invalid use of \"parameter\" keyword");
  }
}



void
NeuronFitInfo::scheduleLoadTrace(const string & fileName,
                                 const size_t & traceNum)
{
  // get deltaT numT and numTraces from file
  ifstream fileIn (fileName.c_str());
  if(!fileIn)
    throw FitInfoException("Can't open file \"" + fileName + "\"");
  
  size_t numTraces, lineNum;
  TraceInfo traceInfo;
  try {
    vector<TraceInfo> infoVec = getTraceHeader(fileIn, lineNum);
    fileIn.close();
    
    numTraces = infoVec.size();
    // make sure traceNum is in range [1, numTraces]
    if(traceNum < 1 || traceNum > numTraces)
      throw FitInfoException("Invalid trace number: " + numToString(traceNum));
    
    traceInfo = infoVec[traceNum - 1];
  }
  catch (exception & e) {
    fileIn.close();
    const string errorStr =
      "While reading electrophys header from file \"" + fileName + "\", " +
      e.what();
    throw FitInfoException(errorStr);
  }
  
  // make sure dT is sensible
  if(!(traceInfo.dT > 0))
    throw FitInfoException("Invalid trace deltaT");
  // make sure numT is sensible
  if(traceInfo.numT < 2)
    throw FitInfoException("Invalid trace numT");
    
  // if fileName isn't a key in tracesToLoadMap, create an entry that is a
  //   vector of NULL
  if(tracesToLoadMap.count(fileName) == 0)
    tracesToLoadMap[fileName].resize(numTraces, NULL);
  
  // Create a new trace and set some member data
  Trace trace;
  //    -set deltaT
  trace.deltaT = traceInfo.dT;
  //    -set units
  trace.units = traceInfo.units;
  //    -allocate space in tracePtr->trace
  trace.trace.resize(traceInfo.numT);
  // add to the back of electrophysData
  electrophysData.push_back(std::move(trace));
  
  // add the new trace to the map of traces to load (later) from this file
  tracesToLoadMap[fileName][traceNum - 1] = &electrophysData.back();
  
  // shorten global tFinal if tF < tFinal
  const double tF = traceInfo.dT * (double)(traceInfo.numT - 1);
  if(tF < tFinal)
    tFinal = tF;
}



void
NeuronFitInfo::loadTraces(void)
{
  // first allocate memory for all recorded traces, fit error traces
  for(Trace & trace : electrophysData) {
    // loop through traces
    if(trace.traceAction == Trace::Record) {
      // this trace is recorded, so allocate space for data
      const size_t len = 1 + (size_t)(tFinal / trace.deltaT);
      
      trace.trace.resize(len);
    }
    else if(trace.traceAction == Trace::Fit) {
      // the error trace is recorded, so allocate space for data
      const size_t len = 1 + (size_t)(tFinal / trace.deltaT);
      
      trace.errorTrace.resize(len);
    }
  }
  
  // next load all the real traces in tracesToLoadMap
  size_t lineNum;
  for(auto & mapObj : tracesToLoadMap) {
    // loop through files that contain data
    try {
      lineNum = 0;
      
      // open the data file
      char const* const dataFileName = mapObj.first.c_str();
      ifstream fileIn (dataFileName);

      // read header info
      vector<TraceInfo> header = getTraceHeader(fileIn, lineNum);
      // read in the traces
      getTraces(fileIn, lineNum, header, mapObj.second);
      // close the file
      fileIn.close();
    }
    catch (exception & e) {
      const string errorStr =
        "While reading electrophys data from file \"" + mapObj.first +
        "\" line " + numToString(lineNum) + ": " + e.what();
      throw FitInfoException(errorStr);
    }
  }
}



vector<TraceInfo>
NeuronFitInfo::getTraceHeader(ifstream & fileIn, size_t & lineNum) const
{
  lineNum = 0;
  
  try {
    // get the next non-trivial line
    string line = getNextLine(fileIn, lineNum, '#');
    // that line should contain the number of traces
    const size_t numTraces = stringToSize(line);
    
    vector<TraceInfo> traceInfoVec (numTraces, TraceInfo ());
    
    // get all the header info
    size_t traceNum = 0;
    for(TraceInfo & info : traceInfoVec) {
      // for each trace, get the header info
      
      // get the next line
      line = getNextLine(fileIn, lineNum, '#');
      
      // split by whitespace into four columns (starting at end of line)
      SplitLine splitLine = rsplit(line, whitespace, 3);
      if(splitLine.size() != 4)
        throw FitInfoException("Incorrect number of columns in header");

      info.name = splitLine[0];
      info.units = splitLine[1];
      info.numT = stringToSize(splitLine[2]);
      info.dT = stringToDouble(splitLine[3]);
      info.traceNum = traceNum + 1;
    }
    
    return traceInfoVec;
  }
  catch(exception & e) {
    // close the file
    fileIn.close();
    // throw an exception
    const string errorStr =
      "Error reading line " + numToString(lineNum) + ": " + e.what();
    throw FitInfoException(errorStr);
  }
}



void
NeuronFitInfo::getTraces(ifstream & fileIn, size_t & lineNum,
                         const vector<TraceInfo> & header,
                         vector<Trace*> & columns) const
{
  // read data into the requested traces
  
  for(size_t traceNum = 0; traceNum < header.size(); ++traceNum) {
    // loop through the traces, loading in the requested ones
    
    Trace* tracePtr = columns[traceNum];
    const size_t numT = header[traceNum].numT;
    if(tracePtr == NULL) {
      // this trace isn't wanted, read in dummy values just to advance
      // through the file
      for(size_t numRead = 0; numRead < numT; ++numRead) {
        // read numT non-trivial lines in from file
        getNextLine(fileIn, lineNum, '#');
      }
    }
    else {
      // this trace is wanted, read data into it

      if(tracePtr->trace.size() != numT) {
        throw FitInfoException("Trace " + tracePtr->targetName
                             + " was not sized correctly for reading in data");
      }

      for(double & value : tracePtr->trace) {
        // read numT values into trace, one on each line
        value = stringToDouble( getNextLine(fileIn, lineNum, '#') );
      }
    }
  }
}


