Fitness landscape (surface plots of model error) are made for the viscinity of the model fit to synthetic data  with mismatched conductance models. That is, perturbed models are created by making small changes to the best-fit model, and then their error is calculated.

Visualization of fitness landscape can be made by running
./perturb.py [param1] [param2]
where [param1] and [param2] are replaced by parameters from the model. Figure 6 was made with
./perturb.py gBar_Na_Axon gBar_Kd_Axon
./perturb.py gBar_KCA_Passive gBar_CaT_Passive