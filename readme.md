Related publication: Angeloni et al. (2023) Dynamics of cortical contrast adaptation predict perception of signals in noise.
https://doi.org/10.1101/2021.08.11.455845

This code was tested using MATLAB 2019a on Mac OSX. The script `run_analysis.m` generates all of the main figures and statistics used in the paper. To use this script, download the data from [DRYAD](**dryad link**), unzip into a folder called **_data** and place that folder in the same directory as `run_analysis.m`. 

For more detailed instructions go to the [analysis](#running-the-analysis) section.

From scratch, the entire analysis takes ~2-3 days to run (this approximate runtime is on an iMac with 4GHz Core i7 CPU, 16GB RAM). A large proportion of this time is taken fitting LN models to each neuron in the dataset (Figure 6), which takes approx 10-40 s for each of the ~5000 recorded cells.

## Data

The **_data** folder contains data from several experiments:

1. The `_behavior` folder contains code and data from the behavioral experiments. Each mouse has a folder in `_behavior/_data/` which contains a .mat file and a .txt log file for each session.

2. The `_glm` folder contains acute recording data used to fit the Poisson GLM in Figure 2. Each neuron has a .mat file in this folder. Because the fitting procedure was run on a computing cluster, we also provided the fit results in the `_res` folder.

3. The `_spectrograms` folder contains stimulus information for each unique session during the microdrive recordings. Each .mat file contains the full spectrogram generated from the concatenated trials, and other general stimulus information.

4. `spikeData.mat` contains the spiking activity from all neurons recorded during the behavioral task.

    `cellinfo` is a cell array of useful information about each neuron. There are 11 columns, corresponding to:

    >[mouse, session/date, cell_number, cellID, single/multi-unit, 1 or 2 for single or multi-unit, cellID, average spks/s, behavioral data file location, spike data file location, task type]

    `spikes` is a cell array containing spike times for each neuron.

    `waveform` is a struct containing information about the spike waveform for each neuron.

5. `sessionData.mat` is a struct containing information about each session.

    `behavior` contains information about task performance on each session. The most important fields are:

    - `trialType`, in which column 1 is the volume index of the target (`0` = noise only, `1+` are different target volumes), column 2 is the target time offset index, and column 3 is the background scene index (columns 1 and 2 index into `session.SNR` and `session.offsets`, respectively)
    - `abort`: which tells whether each trial was aborted early or not
    - `response`: `0` when no licks in the response window, `1` when a lick occurred in the response window
    - `goodTrials`: `1` when the mouse is actively responding, `0` when the mouse stopped responding

    <p></p>
    
    `events` contains information about events recorded from the electrophysiology system. `ts` are the raw event timestamps, `state` contains the channel up (+) or down (-) for each event, `stim` is all stimulus events, `lick` is all lick events, `reward` is all reward events, `trialOn` is the onset of each trial, `targOn` is the onset of the response window in each trial.

    `session` contains information about each session, including the mouse, session date and time, recording times, task type, number of target levels and offsets, their corresponding values for that session.  
    >NOTE: In some sessions the `offsets` (ie. time of the response window) came 25 ms early. This is accounted for in `session.offsets` but is not accounted for in other places where the target timing is defined.


6. `muscimolCortexData.mat` contains spiking activity acutely recorded from mice with topically applied saline or muscimol (Extended Data Figure 4). It contains a 2(mice)x1 struct array with fields `d`, which contains spike times, cell information, events, and stimulus information and `ops`, which contains analysis information.

7. `psych_sig_cells.mat` contains frozen bootstrap results defining which neurons were significantly responsive to targets in the psychometric task.

## Running the analysis

1) Make a new folder and change to that directory:

    ```
    mkdir contrast_study
    cd contrast_study
    ```

2) Clone this repo (behavioral analyis):

    ```
    git clone https://github.com/chris-angeloni/contrast_behavior
    ```

3) Clone the GLM repo:

    ```
    git clone https://github.com/chris-angeloni/contrast_glm
    ```

4) Download the data from [DRYAD](**dryad link**) and unzip the file. Place the `_data` folder in `./contrast_behavior` (ie. in the same directory as `run_analysis.m`).

5) If you wish to plot the GLM results, make sure line 178 of `run_analysis.m` is: `cd ../contrast_glm/`.

6) **From the contrast_behavior repo**, run the run_analysis script:

    ```
    cd ./contrast_behavior
    run_analysis
    ```

`run_analysis.m` will generate all of the main figures, most of the extended data figures, and the statistics. On the first run, this script will take several days, most of the time being occupied by fitting the cross-validated LN models to each of the ~5000 neurons in the dataset (it takes anywhere from 20-60 seconds per neuron). After first the analysis, results files will have been saved and subsequent runs will be much faster.

For your convenience, we have provided precomputed results files used for the paper. They are all in the `_data` folder with the prefix `_res`. If you would like to run the analysis from scratch, move or delete the files with `_res` at the beginning and the script will regenerate them.

## Notes

Code for the GLM in Figure 2 can be found here: https://github.com/chris-angeloni/contrast_glm

Results for the LN model may vary by run (due to crossvalidation splits), however, we have found the effects to largely remain stable over multiple runs.

Some supplementary figures were generated separately from this script, largely due to processing constraints.
- Supplementary Figure 2: GC-GLM simulations were run on a computing cluster to massively parallelize and speed up the analysis. The folder `_glm\_simulation` contains the necessary code to do this, but is specific to the cluster setup used by the authors.
- Supplementary Figure 6: This is histology, not code generated.
- Supplementary Figure 7: This analysis examines whether STRFs recorded during the task are affected by contrast. If you would like to run this analysis, you can do so using `run_strf.m` in the main folder. This has not been extensively tested for stable results.
