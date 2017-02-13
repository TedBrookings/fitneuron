This repository contains the essential source code and data for the paper
Brookings T, Goeritz ML, Marder E (2014), "Automatic parameter estimation of multicompartmental neuron models via minimization of trace error with control adjustment" currently in revisions

All code is Copyright (c) 2014 Ted Brookings under the MIT license (see
LICENSE.TXT).

As is hopefully apparent from the accompanying paper, the code here is more proof-of-principle than practical tool. It would be wonderful if one day I had to the time to develop and improve this software until it is capable of reliably producing useful neuron models. No promises though. Some initial progress has been made on improving the minimization routine (which is the main source of difficulties). I would also like to re-do some of the design to make use of lambda functions (which have subsequently become available in C++).

All of this code has been developed on Ubuntu linux. It has been subsequently tested on RedHat linux. Although some effort has been made to avoid requiring linux, that hasn't been a major goal, and none of the code has been tested on Windows or Mac.

A non-exhaustive list of the software that is required to run this code:
g++ 4.8.1
OpenMpi 1.6.2
gnu scientific library (gsl) 1.16
matplotlib 1.2.1
scipy 0.12.0

Supplemental Contents:
scripts - Contains python and matlab code for plotting and analyzing recorded/simulated traces
models - Contains subdirectories for each neuron model from the paper
         Each subdirectory contains necessary data and scripts to fit or simulate models, as well as a log
         of the outcome of the fits we performed. Keep in mind that the fits are stochastic, so exact duplication of results is
         not guaranteed.
figure 1 - Contains a python script which generates the panel data for figure 1
figure 6 - Contains a script and the basic data necessary to generate the compensation plots from figure 6
