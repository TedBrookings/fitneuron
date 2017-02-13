Minor issues:
1. The original fit for this paper included h_KCa_Soma and h_KCa_PN as fixed parameters. Since they are not elements of the model, this is useless (but harmless). I have included them in startFit.txt and startSimResults.txt to preserve compatibility with the existing fit results.

2. The original fit read in a file with 10 seconds of data and only used 2 seconds for fitting. The paper plots with 4 seconds for the sake of visibility. Thus if you run the simulate.sh script here, you will see a simulated model with 4 seconds of output overlayed on a target model with 10 seconds of data.
