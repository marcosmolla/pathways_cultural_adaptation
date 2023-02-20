### 1- Initialize parameters ...
grid = collect(Base.product(
100,              # 1, Number of individual nodes
5e2,              # 2, Number of generations
2,                # 3, Repetitions (local)
500,              # 4, Number of traits;
100,              # 5, Number of learning turns (each IL and SL)
1,                # 6, Social learning success; #[.01 .25 .5 .75 .95]
1e-2,             # 7, Individual learning probability
1000,             # 8, p_n; Set to 1000 for random initialisation or to -1 to create a fully connected network
1000,             # 9, p_r; Set to 1000 for random initialisation
.05,              # 10, Mutation rate
[.2 1],           # 11, SD of LogNormal payoff distribution use "0" for uniform
[.0001 1],        # 12, Environmental turnover (average number of patches changing each round) # CHANGE THIS TO proportion of utilties that are redrawn each generation (at the end of N time steps)
true,             # 13, Evolve linking parameters
5,                # 14, K, number of neighbours to observe when focal is in a complete graph
true,             # 15, Complex contagion (if false, ise simple contagion for social learning)
1,                # 16, Minimum fitness
1,                # 17, For 1: Mortality selection, -1: fertility selection, 0: neutral selection
true,             # 18, Michaelis-Menthen payoff saturation

collect(1:1)))[queue]; # additional repetitions (global; to spread jobs over more cores)

global report = true; # whether to store results or not

# Record parameters for export
s = Dict(
  # Setup for simulation
  "nRound"            => grid[1]*grid[2], # The number here determines the number of simulation rounds, i.e. one generation is equivalent to replacing all individuals (on average) once
  "repe"              => grid[3], # number of repetitions per execution of the code

  # Setup for graph and connection dynamics
  "nod"               => grid[1], # number of nodes
  "deg"               => 4,       # average degree of initialisation network
  "pn"                => grid[8], # probaility to connect with neighbours of parents; set to 1000 
  "pr"                => grid[9], # probaility to connect with strangers
  "k"                 => grid[14],# number of neighbours to copy from when in complete graph

  # Setup for culture
  "nTraits"           => grid[4], # total number of potential starting point for trait trees, no more than this can be innovated
  "nLearnTurns"       => grid[5], # number of each IL and SL learning turns
  "slSuccess"         => grid[6],
  "ilSuccess"         => grid[7], 
  "complexContagion"  => grid[15], # is social learning complex or simple contagion

  # Setup for evolution
  "evolvePN"          => grid[13], # Do values for linking propensity for social inheritance evolve?
  "evolvePR"          => grid[13], # Do values for linking propensity for random connections evolve?
  "mutRate"           => grid[10], # Mutation rate for inherited parameters

  "sigma"             => grid[11], # standard deviation for LogNormal distribution
  "envTurnover"       => grid[12], # environmental turnover as proportion of utilities to change per generation
  "minFit"            => grid[16], # minimum fitness
  "mortalitySelect"   => grid[17], # minimum fitness
  "MM"                => grid[18], # minimum fitness
  );