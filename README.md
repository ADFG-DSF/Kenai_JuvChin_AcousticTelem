# Kenai River Chinook Salmon Smolt Acoustic Telemetry Study, 2025

This pilot study will be the first year of a 3-year collaboration with the 
Alaska Department of Fish and Game (ADF&G) and University of Alaska Fairbanks 
(UAF) that will use juvenile salmon acoustic telemetry system (JSATS) technology 
to investigate the juvenile life stage of Kenai River Chinook salmon. A 
combination of minnow traps, beach seines, and other devices will be deployed to 
capture Chinook salmon smolt in the Kenai River drainage. A sample of up to 500 
Chinook salmon smolt of taggable size will be surgically implanted with acoustic 
transmitters and allowed to continue their downstream migration. Downstream 
movements and inriver behavior of tagged fish will be tracked by ADF&G using 
strategically placed hydrophone receiver arrays spanning the Kenai River. Early 
marine behavior and survival of tagged smolt will be studied by UAF and ADF&G 
using hydrophone receivers in the nearshore waters of Cook Inlet near the Kenai 
River mouth. Hydrophone receivers will be downloaded periodically to determine 
fish passage, timing, location, migratory behavior, and survivorship. Mobile 
boat tracking may also be conducted to track tagged smolt to further understand 
migratory behavior and survival.  In addition, a feasibility study will be 
conducted in late fall 2025 to determine if rearing Chinook salmon juveniles 
can be acoustic tagged and tracked to overwintering locations. This first year 
of the 3-year collaborative study will be focused on capacity building and 
determining successful techniques to acoustically tag and track fish. The 
study design will be refined for the second and third study years based on 
insights learned from this initial study in 2025. 

# /OP_2025

To date, all simulation and analysis work has been directly pertaining to the 
Operational Plan written in 2025, specifically:

* Construction of a modeling framework and validating model appropriateness

* Estimation of expected inferential precision, by means of simulation.

## Modeling overview

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

$X_{i,t_i} \sim Binom \right((\phi_{t_i} \times \lambda), 1 \left)$ at time of
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

## Folder Contents

### /R

* **Kenai_telem_rp.R** This script simulates a sequence of survival and detection
probability vectors according to anticipated values, simulates detection histories 
from each set of probabilities, and then runs each of the three candidate models 
using the simulated data as input to estimate parameters.  The inferential precision
of each model is then evaluated by comparing the sequences of probability vectors
to the values estimated by the candidate models.

#### /R/experimenation 

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