function Neuron_CompareTraces(traceFile1, traceFile2, varargin)
% Neuron_CompareTraces(targetTraceFile, fitTraceFile, varargin)
%  Compare a fit to a target trace
% INPUTS:
%  -targetTraceFile: file name of target trace
%  -fitTraceFile: file name of trace from fit model
%  OPTIONS: [defaults]
%   -plotSubject: [false] if set to true or a string, make plots. If set to a
%                   string, alter the plot titles to encorporate string
%   -minT:        [0.0] Consider data from trace at or after this time
%                   (in seconds)
%   -maxT:        [Inf] Consider data from trace at or before this time
%                   (in seconds)
%   -matchGoal:   [2.0] When computing performance gamma, consider spikes a
%                   match if there error is <= this number (in ms)

% set the default options
defaultOptions = { ...
  'plotSubject', false, ...
  'minT', 0.0, ...
  'maxT', Inf ...
  'matchGoal', 2.0 ...
};
% get the options overrides from varargin
[options, modified] = GetOptions(defaultOptions, varargin); %#ok<NASGU>

traces1 = Neuron_LoadTraceFile(traceFile1);
traces2 = Neuron_LoadTraceFile(traceFile2);

v1 = getVoltageTrace(traces1, options, traceFile1);
v2 = getVoltageTrace(traces2, options, traceFile2);

if options.plotSubject
  spikes1 = GetSpikes(v1.dT, v1.trace, 'plotSubject', v1.name);
  fprintf('Target trace had %d spikes\n', length(spikes1.times))
  a1 = gca();
  if strcmp(v1.name, v2.name)
    name2 = [v1.name, 'B'];
  else
    name2 = v2.name;
  end
  spikes2 = GetSpikes(v2.dT, v2.trace, 'plotSubject', name2);
  a2 = gca();
  % link time and voltage axis together
  aHandles = [a1, a2];
  linkaxes(aHandles, 'xy');
  %linkaxes(aHandles, 'x');
else
  spikes1 = GetSpikes(v1.dT, v1.trace);
  fprintf('Target trace had %d spikes\n', length(spikes1.times))
  spikes2 = GetSpikes(v2.dT, v2.trace);
end

duration = v1.dT * (v1.numT - 1);
compareSpikes(spikes1, spikes2, duration, options);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vTrace = getVoltageTrace(traces, options, traceFileName)

for tracename = fieldnames(traces)'
  if any(strfind(tracename{1}, 'error'))
    continue
  end
  vTrace = traces.(tracename{1});
  if strcmp(vTrace.units, 'mV')
    if options.minT > 0 || options.maxT < Inf
      minInd = max(1, ...
                   ceil(1000 * options.minT / vTrace.dT) + 1);
      maxInd = min(floor(1000 * options.maxT / vTrace.dT) + 1, ...
                   length(vTrace.trace));
      vTrace.trace = vTrace.trace(minInd:maxInd);
    end
    return
  end
end

error('Could not find a voltage trace in %s', traceFileName)
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function compareOK = compareSpikes(spikes1, spikes2, duration, options)

[missingSpikes, multipleA, matched1, matched2] = ...
  checkSpikeOverlap(spikes1, spikes2);
numWrong = 0;
if any(missingSpikes)
  missingTimes = spikes1.times(missingSpikes) / 1000;
  numWrong = numWrong + length(missingTimes);
  fprintf('%d missing spike(s) at time: ', length(missingTimes))
  fprintf(' %g', missingTimes)
  fprintf('\n')
end
if any(multipleA)
  multipleTimes = spikes1.times(multipleA) / 1000;
  fprintf('%d instances of multiple spikes found when there should be only one at time: ', ...
          length(multipleTimes))
  numWrong = numWrong + length(multipleTimes);
  fprintf(' %g', multipleTimes)
  fprintf('\n')
end
[extraSpikes, multipleB, ~, ~] = checkSpikeOverlap(spikes2, spikes1);
if any(extraSpikes)
  extraTimes = spikes2.times(extraSpikes) / 1000;
  numWrong = numWrong + length(extraTimes);
  fprintf('%d extra spike(s) at time: ', length(extraTimes))
  fprintf(' %g', extraTimes)
  fprintf('\n')
end
if any(multipleB)
  multipleTimes = spikes2.times(multipleB) / 1000;
  numWrong = numWrong + length(multipleTimes);
  fprintf('%d instances of one spike found when there should be many at time: ', ...
       length(multipleTimes))
  fprintf(' %g', multipleTimes)
  fprintf('\n')
end

numMatched = length(matched1);
fprintf('%d spikes matched.\n', numMatched)
dT = spikes1.times(matched1) - spikes2.times(matched2);
matchRate = numMatched / (numWrong + numMatched);
fprintf('Spike match rate is %f%%\n', 100 * matchRate)
fprintf('Mean time difference = %g +- %g ms\n', mean(dT), std(dT))
meanAbsDT = mean(abs(dT));
stdAbsDT = std(abs(dT));
fprintf('Mean abs(time difference) = %g +- %g ms\n', meanAbsDT, stdAbsDT)
dT = sort(dT);
lowDT = dT(round(0.05 * length(dT)));
highDT = dT(round(0.95 * length(dT)));
fprintf('95%% of matched spikes have dT between %g and %g ms\n', ...
        lowDT, highDT)
percent2ms = 100 * sum(dT >= -2.0 & dT <= 2.0) / length(dT);
fprintf('%.1f%% of matched spikes have dT between +- 2 ms\n', percent2ms)
      
widthErr = spikes2.width(matched2) / spikes1.width(matched1) - 1.0;
fprintf('Mean spike width = %g +- %g mV, fractional error = %g%%\n', ...
        mean(spikes1.width(matched1)), std(spikes1.width(matched1)), ...
        100 * widthErr)
heightErr = spikes2.height(matched2) / spikes1.height(matched1) - 1.0;
fprintf('Mean spike height = %g +- %g mV, fractional error = %g%%\n', ...
        mean(spikes1.height(matched1)), std(spikes1.height(matched1)), ...
        100 * heightErr)

gamma = computePerformance(spikes1, spikes2, duration, options);
fprintf('Performance gamma = %g\n', gamma)
      
if options.plotSubject
  h = NamedFigure('Spike Time Differences');
  set(h, 'WindowStyle', 'docked')
  hist(dT, max(round(sqrt(length(dT))), 50))
  title('Histogram of spike time differences')
  xlabel('Time difference (ms)')
  ylabel('Count')
end

compareOK = ~(any(missingSpikes) ||...
              any(multipleA) || ...
              any(extraSpikes) || ...
              any(multipleB));
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [missing, multiple, matchedA, matchedB] = ...
  checkSpikeOverlap(spikesA, spikesB)

n1A = spikesA.preMinV.ind;
n2A = spikesA.postMinV.ind;
tA = spikesA.times;
n1B = spikesB.preMinV.ind;
n2B = spikesB.postMinV.ind;
tB = spikesB.times;
missing = [];
multiple = [];
matchedA = [];
matchedB = [];
for n = 1:length(n1A)
  % find spikes that overlap with the nth A-spike
  m = find(n1B < n2A(n) & n2B > n1A(n));
  % exclude spikes that are closer to (n-1)th or (n+1)th A-spike
  exclude = [];
  for k = 1:length(m)
    m_k = m(k);
    if (n > 1 && tB(m_k) - tA(n-1) < tA(n) - tB(m_k)) || ...
       (n < length(n1A) && tA(n+1) - tB(m_k) < tB(m_k) - tA(n))
      exclude = [exclude, m_k]; %#ok<AGROW>
    end
  end
  if any(exclude)
    m = setdiff(m, exclude);
  end

  % exclude spikes that are closer to (n+1)th spike
  
  
  if isempty(m)
    missing = [missing, n]; %#ok<AGROW>
  elseif length(m) == 1
    matchedA = [matchedA, n]; %#ok<AGROW>
    matchedB = [matchedB, m]; %#ok<AGROW>
  else
    % record that there were multiple matches for n
    multiple = [multiple, n]; %#ok<AGROW>
    % find the smallest abs(dT)
    [~, k_min] = min(abs(tA(n) - tB(m)));
    m_min = m(k_min);
    % record the corresponding spikes as the match
    matchedA = [matchedA, n]; %#ok<AGROW>
    matchedB = [matchedB, m_min]; %#ok<AGROW>
  end
end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function gamma = computePerformance(spikesA, spikesB, duration, options)
% compute gamma, a measure of spike prediction performance, as listed in
%   Jolivet et al 2008 "A benchmark test for a quantitative assessment of
%   simple neuron models"
%   equation (2)
% Assumes that spikesA are from data, and spikesB a model

% store the number of A and spikes
numASpikes = length(spikesA.times);
numBSpikes = length(spikesB.times);

% compute the rate of A
rateA = numASpikes / duration;

% estimate the number of spikes that would be "predicted" by a Poisson
% process with rate of A
% Note: deviation from Jolivet et al, I use 1.0 - exp(-f * delta)
%  instead of f * delta
probSpikePredicted = 2 * (1.0 - exp(-rateA * options.matchGoal));
numPoisson = numASpikes * probSpikePredicted;

% compute the number of "predicted" spikes
numCoinc = 0;
for tA = spikesA.times
  tMin = tA - options.matchGoal;
  tMax = tA + options.matchGoal;
  if any(spikesB.times >= tMin & spikesB.times <= tMax)
    numCoinc = numCoinc + 1;
  end
end

gamma = (numCoinc - numPoisson) * 2.0 / (numASpikes + numBSpikes) ...
        / (1.0 - probSpikePredicted);
end
