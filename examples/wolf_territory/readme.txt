Dean Koch
Feb 2023

INTRODUCTION 

This runs a simulation experiment using the NetLogo model from the 2022 Sells et al paper "Economical Defense of Resources Structures Territorial Space Use in a Cooperative Carnivore" in Proceedings of the Royal Society B.

We run simulations from this model to get a set of artificial wolf territory maps. Then we use a (biased) grid based occupancy sampling method to estimate the area occupied, and measure the resulting bias. 

To repeat the workflow, first download the data supplement from the Sells et al paper and unzip to a local directory. This will produce a file "Model.zip" containing the NetLogo model configuration files. There is no need to unzip "Model.zip".

Note that "Model.zip" contains an installer for NetLogo v6.1.1. You can install this if you want, but for this simulation study you will need an earlier version, NetLogo v6.0.4. The installer for v6.0.4 can be found online (or ask me and I'll send a copy).


CONTENTS

Our simulation study uses 4 R files

* "0_netlogo_helpers.R" : helper function definitions
* "1_run_netlogo.R" : runs a NetLogo simulation, saves results to .tif
* "2_process_netlogo_results.R" : runs analysis on all output from (2), saves results to .rds
* "3_display_netlogo_analysis.R" : creates three .png plots using analysis results

There is a subfolder "writeup" with a (.Rmd) markdown document summarizing everything, but this isn't required for the workflow.


INSTRUCTIONS

1) Install NetLogo v6.0.4 (not the newer version included in the Sells et al data supplement!)

2) Open the R files numbered 1-3 and in each one change the path `proj_dir` to something appropriate for your computer (this is where model files and all results are stored). Check also that `supplement_zip` and `helper_path` point to the right place at the beginning of "run_netlogo.R".

3) Run the script "run_netlogo.R" multiple times. A new output .tif file is written to `proj_dir` each time.

4) Run the script "process_netlogo_results.R" to run the analysis on all .tif files found in `proj_dir`

5) Run the script "display_netlogo_analysis.R" to summarize the analysis in two figures


NOTES

* NetLogo simulations are slow. Each repetition of "1_run_netlogo.R" takes about 2-3 hours
* NetLogo v6.0.4 is (currently) required for compatibility with the RNetLogo package
* NetLogo input/output is found in the "Model" subfolder of `proj_dir` after running a simulation