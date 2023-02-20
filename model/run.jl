#!/usr/bin/env julia

### Load libraries
using Distributions
using StatsBase
using Random
using RCall
R"library(igraph)";

# Set working directory 
global paths="PATH/TO/DIRECTORY";
# Load methods from file
include(string(paths, "/methods.jl"));
# Run all iterations of the simulation 
for q in 1:2
  global queue = q;
  # Initialise parameters for the given iteration
  include(string(paths, "/main.jl"));
# Run Simulation
@time runit();
print(".")
end
