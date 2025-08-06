# Welcome!
This is the GitHub Repository for the study "Are Resting-State EEG Spatiotemporal Complexity or Microstate Features Associated With Fluid Intelligence?". This project has two preregistrations. You can find the preregistration for the complexity analysis [here](https://osf.io/bt5v8) and the one for the microstate analysis [here](https://osf.io/8qtwj).

To reproduce the results from the paper, simply clone the repository and request the data from the authors.

## Analysis Pipeline
Make sure to stay in the directory of the script you run (e.g. Matlab script $\rightarrow$ Matlab directory), otherwise the script will have problems finding the other files.

1. Start within the `/Matlab` folder and run the script `Main.m`. Adapt the folder paths so it matches your directory. This script will automatically preprocess and epoch the data, and calculate multivariate multiscale sample entropy (mMSE). Depending on your machine an the number of cores, this will take quite a while. So grab a coffee, or five, in the meantime.
   
   1.1 Optional: Check out the `Quality_Checks.m` script for several plots.
3. Change to the `/Python` folder and start `get_microstates.py`. We used VS Code for the analysis, but feel free to use any interpreter. You can create a virtual environment with the `requirements.txt` in the folder. Licenses are within the `/Licences` folder.
4. When the microstate analysis is finished, run `plot_microstates.py`. This will also combine the behavioral and demographic data with the microstate data so it is conceptually similar to the Matlab output.
5. Now switch to the `/R` folder. Run `PLSR.R` for partial least square regression analysis on both mMSE and microstate data. For plotting and calculation of retest-correlations, run `Correlations_and_Plots.R`.

And that's it.
