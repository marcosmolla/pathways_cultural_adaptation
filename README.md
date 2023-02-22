# Code and Material for Pathways to Cultural Adaptation

Code for [Pathways to Cultural Adaptation](https://biorxiv.org/cgi/content/short/2023.02.21.529416v1) manuscript

## File overview

The following files are found in the `model/` folder

* `run.jl` initates julia simulations after loading methods from `methods.jl` and the job array specific parameters from `main.jl`
* `summariseFiles_rout.R` opens all `*.Rdata` files from the output directory, summarises values, returns a summary file into the working directory (with the current setting it will create three files one with a general summary, one withj a record of all the individual repertoires, and one with all the time series data)

The `figures/` folder contains `R` code to generate the figures of the main text. 

The `material/` folder contains supplemental material.

## Working with the simulation model
### Requirements to run the model

Running the simulation code requires:

* `julia` (dependencies: `Distributions`, `StatsBase`, `Random`, `RCall`); tested and verified for Julia version 1.8.5
* `r` (dependencies: igraph)

### Simulations reported in our article

In our manuscript, [Pathways to Cultural Adaptation](https://biorxiv.org/cgi/content/short/2023.02.21.529416v1), we report results for the following versions of the model:

1. Homogeneous environments, evolving $p_n$, $p_r$
2. Homogeneous environments, evolving $p_n$, $p_r$, variable innovation and social learning success rate
3. Heterogeneous environments, evolving $p_n$, $p_r$
4. Heterogeneous environments, evolving $p_n$, $p_r$, variable population size

Furthermore, in the ESM we report results for simulations with

1. Neutral selection
2. Fixed linking parameters
3. Random graphs
4. Simple Contagion
5. Fertility selection

To run the individual simulations you need to adjust the parameters in the `grid` array that is defined in `main.jl`.

For **homogeneous environments**, set parameters 11 ($\sigma$) and 12 ($\tau$) to 0. They control the shape of the lognormal distribution and Lognormal(0,0)$ = 1$. For **heterogeneous environments** change the values of $\sigma$ and $\tau$ accordingly.

To let $p_n$ and $p_r$ evolve, set parameter 13 to `true`. In this case, you might want to initialise the population with random values for $p_n$ and $p_r$, which you can do by setting parameters 8 and 9 to `1000`. To keep $p_n$ and $p_r$ fixed, set 8 and 9 to the preferred values and change 13 to `false`.

Learning success rates $\alpha$ and $\beta$ never evolve but can be set to different values by changing parameters 6 and 7.

Population size can be adjusted by changing parameter 1.

How selection works is determined by parameter 17, whereby `1` indicates mortality selection, `0` indicates neutral selection, and `-1` indicates fertility selection.

The mode of social learning contagion is controlled by parameter 15, whereby `true` indicates complex contagion, and `false` indicates simple contagion.

To simulate learning on random graphs, set parameter 8 to `-1` (which sets up a fully connected network), and adjust parameter 14, which controls the number of randomly selected neighbours an individual can learn from, $k$ (where $k<N$).

The minimum fitness can be adjusted with parameter 16. And saturating payoffs (Michaelis-Menten model) can be turned on or off with parameter 18.

## Additional material

In the main text of our manuscript [Pathways to Cultural Adaptation](https://biorxiv.org/cgi/content/short/2023.02.21.529416v1), we report on populations cycling between the high connectivity state (with high payoffs) and sparse networks (with low payoffs). 

Here is an example of this pattern:

![Track](/material/220927_1_varEnv_noEnv_anmatin_1_KS92EEM.gif)

When looking at several simulations across time, we observe how there both pathways present (low payoff and high payoff) throughout the simulaitons:

![Track](/material/220927_1_varEnv_noEnv_animatin_YQ56UZQ.gif)