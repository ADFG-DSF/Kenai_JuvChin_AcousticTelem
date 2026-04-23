# Kenai River Chinook Salmon Smolt Acoustic Telemetry Study, 2026-2027

This study will be the final two years of a three-year collaboration between the 
Alaska Department of Fish and Game (ADF&G) and University of Alaska Fairbanks 
(UAF) using juvenile salmon acoustic telemetry system (JSATS) technology for the 
first time in the State of Alaska to investigate the smolt life stage of Kenai 
River Chinook salmon. The primary purpose of this project is to improve 
understanding of Chinook salmon smolt behavior and survival as they migrate from 
the Kenai River into the marine nearshore waters of Cook Inlet. Minnow traps, 
beach seines, and other devices will be deployed to capture Chinook salmon smolt 
in the Kenai River drainage to surgically implant with an acoustic transmitter 
and release to continue their downstream migration. Movements, behavior, and 
survival of tagged fish will be monitored with strategically placed hydrophone 
receiver arrays in the Kenai River and nearshore waters of Cook Inlet near the 
Kenai River mouth.

A previous iteration of this study was conducted in 2025, which represented a 
pilot effort.



# Modeling overview

Estimation of survival and detection probabilities will be possible from detection
histories for each instrumented fish, which may be formatted like the following 
matrix, in which each row corresponds to a unique fish and each column corresponds
to a given receiver array:

         [,1] [,2] [,3] [,4] [,5] [,6]
    [1,]    1    1    1    1    0    1
    [2,]    1    0    1    1    1    0
    [3,]    NA   1    1    1    0    1
    [4,]    NA   0    1    0    0    0
    ...

This may be conceptualized as an unobserved state process in which survival of 
fish $i$ at time $j$ is conditionally distributed with natural survival probability 
$\phi_j$ and handling survival probability $\lambda$, that is:

$X_{i,t_i} \sim Binom \left((\phi_{t_i} \times \lambda), 1 \right)$ at time of
entry $t_i$, and

$X_{i,j} \sim Binom(\phi_j, X_{i,j-1})$ thereafter.

This gives rise to a data process, in which fish $i$ is detected at time $j$ 
with conditional probability $p_j$ depending on survival, that is:

$Y_{i,j} \sim Binom(p_j, X_j)$

Bayesian models were constructed to reflect this probability model, with three 
candidate formulations:

* A fairly simple Hidden Markov model consisting solely of the state & data processes
outlined above, with weakly informative priors on the survival and detection 
probabilities

* The same Hidden Markov model, but with the state process for survival imputed
when survival can be logically inferred from the observation history

* A Multinomial model with probability parameters of all possible detection histories
as derived from the relative survival and detection probabilities.

To date, all three model formulations give fully equivalent inferences; the only
difference being run-time.


# Data processing overview

Due to the high degree of ambient acoustic noise in a riverine or nearshore 
environment, it can be expected that the vast majority of acoustic signals 
recorded by each listening array will be acoustic artifacts rather than 
instrumented fish.  It will therefore be necessary to employ a suite of data 
filters in order to ensure data quality and prevent compromising inferences.  
Within a set of receiver files spanning much of the 2025 field season, an 
overall filtering acceptance rate of substantially less than 1% was observed.

The following filters have been developed, and will likely be employed in this order:

- Prefilter: Allow entries for which some minimum number of entries (default of 2) exist for each tag number.

- Tag Code filter: Allow entries for which tag number exists in the library of tags that have been deployed.

- Interval filter: Allow entries for which the time interval between entries of a given tag is consistent with the beep rate (default of 2.5-3.5 seconds).

- Event filter: Allow entries for which at least some number of entries (default of 3) have been recorded for a given tag within some amount of time (default of 30 seconds).

- Multipath filter: Exclude entries that represent a reflected signal, that is, within a threshold value (default of 0.3 seconds) of another entry.


# Folder Contents

## /data_processing

### /data_processing/R

**receiver_data.processing.R** is the script **intended for all data processing.**

This script is outlined as follows:

  * Data directories to the appropriate locations
  * Reading and filtration algorithms are defined as functions
  * An interactive Shiny app is defined to allow for fine-tuning of filtration steps and parameters
  * Finally, a processing loop is defined which may be run with or without saving data.  It is recommended to run this first without saving data as a final test of filtration steps.

**inriver_passage.R** is a trial script for summarizing inriver passage from 
all *filtered* receiver data to date.

#### /data_processing/experimentation

This folder contains some early experiments with data processing, which were 
superceded by receiver_data.processing.R

### /data_processing/test_data

This folder contains a single receiver file intended for testing

### /data_processing/TESTING

This folder contains all raw files and subfolders associated with the full processing
script.

## /spatial

This folder contains an R script used to create a rivernetwork object for the 
Kenai river, and the associated .Rdata object that can be loaded directly. 

**Analysis_OP_draft.Rmd** This R Markdown file was created to produce text and
simulation output to be pasted into the draft Operational Plan.

## /OP_2025

Files in this folder pertain directly to the Operational Plan outlining methods 
for the 2025 iteration of this study.

### /OP_2025/R
* **entry_model.R** This script simulates a sequence of survival and detection
probability vectors according to anticipated values, simulates detection histories 
from each set of probabilities, and then runs the most current Hidden Markov model 
using the simulated data as input to estimate parameters.  The inferential precision
of the model is then evaluated by comparing the sequences of probability vectors
to the values estimated by the candidate models.  The principal difference of the 
**entry model** versus the previous set of candidate models is that it allows 
the entrance of tagged individuals at multiple locations relative to listening 
arrays, as well as incorporating handling survival (1-handling mortality) as 
an estimable parameter.

**This is the script that was used to estimate inferential precision as reported
in the 2025 Operational Plan.**

#### /OP_2025/R/obsolete 

This folder contains a record of work under a previous state of the fieldwork 
design, and does not reflect
the current state of the project.

* **Kenai_telem_rp.R** This script simulates a sequence of survival and detection
probability vectors according to anticipated values, simulates detection histories 
from each set of probabilities, and then runs each of the three candidate models 
using the simulated data as input to estimate parameters.  The inferential precision
of each model is then evaluated by comparing the sequences of probability vectors
to the values estimated by the candidate models.

#### /OP_2025/R/experimenation 

This folder contains a record of early, experimental work, and does not reflect
the current state of the project.

* **multinomial_model.R** In this script, a multinomial model was derived and 
implemented, in which the probabilities represent the probability of each possible
detection history, as a function of the vectors of survival probability and 
detection probability.  

* **Kenai_telem_sim.R** The purpose of this script was originally to create a 
Hidden Markov model to estimate the respective survival and detection probabilities.
The model file has quite a lot of unnecessary components and its current state 
represents several avenues of experimentation.  After developing the HM model, 
this script was then repurposed for an early comparison between three candidate 
models, by means of meta-simulation.  


### /OP_2025/data

This folder contains a sequence of .Rdata files representing specific simulations,
that may be read as needed.  One is read in **Analysis_OP_draft.Rmd** in order to
produce a table of anticpated inferential precision for each parameter of interest.

### /OP_2025/testing

This folder contains a few raw receiver files for initial testing.



## /OP_2026

This folder was created for files associated with the Operational Plan describing
2026-2027 research efforts, which is currently empty.