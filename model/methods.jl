### Declare functions
# A function to bound the upper and lower limits of a vector (or single value)
function bound(x,lower=0,upper=1)
    return min(upper, max(lower, x))
  end

# A function to set initial values for ind arrays
  function setValues(evolve::Bool,times,mean,sd)
    if evolve
      return bound.(rand(Normal(mean,sd), times))
    else
      return bound.(fill(mean, times))
    end
  end

# Dynamic network updating, takes adj.matrix (ADJM), probabilities to connect to parent(s), neighbours, and randoms (PB, PN, PR), the index of the parent (if 0, default, one individual is randomly choosen as parent)
  function pbpnprDyn(ADJM, PN, PR, NOD, DB=0)
      # one iteration of the basic model
      if DB==0
        DB = sample(1:NOD,2, replace=false); #sample a death and a birth, without replacement; 1 newborn, 2 parent
      end
      inherit = ADJM[:,DB[2]].*(1 .-(rand(NOD).-PN .> 0.)); # socially inherited connections
      randconn = (1 .- ADJM[:,DB[2]]).*(1 .-(rand(NOD).-PR .> 0.)); # random connections
      newconn = inherit+randconn #total connections
      newconn[DB[2]]=1; #connect to the parent
      ADJM[:,DB[1]]=newconn; #replace the dead individual with the newborn
      ADJM[DB[1],:]=newconn; #same
      ADJM[DB[1],DB[1]]=0; #set self-link equal to zero
      return ADJM
  end

# Select lower triangle of a matrix (excluding the diagonal)
  function lowerTri(x)
    res = zeros(Int,size(x)[1],size(x)[1]);
    for i in 2:size(x)[1]
      res[i,1:(i-1)] .= 1;
    end
    return Bool.(res)
  end

# Random Graph Generation, takes number of nodes (NOD) and average degree (DEG), runs pbpnprDyn 10 times N number of times as burn in for the network
  function initGraph(NOD, DEG, PN, PR)
    P = (DEG*NOD)/(NOD^2); # connection probability
    ADJM = zeros(Int, NOD, NOD); # Set adjacency matrix
    diag = lowerTri(ADJM); # Select one half of the matrix
    npairs = sum(diag); # Number of pairs without diagonal (using only one matrix traingle)
    ADJM[diag] = rand(npairs).<= P;
    ADJM = ADJM + transpose(ADJM); # Populate other half of the matrix
    for i in 1:(10*NOD) # Let PBPNPR Dynamics run for burning in period
      ADJM = pbpnprDyn(ADJM, PN, PR, NOD);
    end
    return ADJM # return adjacency matrix
  end

  # Function to set up adjaceny matrix for a fully connected graph without self-links
    function completeGraph(NOD)
      ADJM = ones(Int, NOD, NOD);
      for i in 1:NOD 
        ADJM[i,i] = Int(0); #set self-link equal to zero
      end
      return ADJM
    end
  
# A function to inherit parental values, or reuse the standard value if this trait does not evovle
function inherit(;evolve=evolve, parents=parents, standard=standard, sd=0.01, rmutation=0.01)
  if evolve
    if rand()<=rmutation
      return bound.(rand(Normal(parents, sd))) 
    else
      return parents 
    end
  else
    return standard
  end
end

# Michaelis-Menthen Model (if MM is set to false, it returns the payoff unaltered)
function mmkinectics(VMAX,K,PAY)
  return (VMAX * PAY) / (K + PAY)
end

# A function to calculate the gini index of a distribution
function gini(x)
  n = length(x);
  mu = mean(x);
  ox = sort(x[:]);
  return 	(1/((n^2) * mu)) * sum(ox.*((2 .* collect(1:n)) .- n .- 1))
end

# A function to set the correct linking parameters, if set to 1000 the function returns two random values from a uniform distribution, whereby the first is between 0 and 1, and the second between 0 and 0.1 (i.e., values for p_n, and p_r)
function setLinking(P)
  if P[1]==1000
    return [rand(),rand()/10]
  else
    return P
  end
end

# A function to set utility values based on a lognormal distribution, all = 1 if sigma = 0
function setUtilities(SIGMA,N)
    return rand(LogNormal(0,SIGMA), N)
end

# A function for mortality/fertility selection
function selectionMF(S, INDPAY)
  # MORTALITY
  if S["mortalitySelect"] == 1 
    # choose newborn inversely realted to their payoff
    if all(INDPAY .== 1) 
      local nb = sample(1:S["nod"])
    else 
      local nb = wsample(1:S["nod"], 1 .- (INDPAY[:] ./ maximum(INDPAY)));
    end
    # choose parent relative to highest payoff
    local p = sample((1:S["nod"])[1:end .!= nb]);
    return (nb,p)
  # FERTILITY 
  elseif S["mortalitySelect"] == -1 
    # choose a random newborn
    local nb = sample(1:S["nod"]);
    # choose a parent relative to their payoff
    local p = wsample((1:S["nod"])[1:end .!= nb], INDPAY[1:end .!= nb]);
    return (nb,p)
  # NEUTRAL
  elseif S["mortalitySelect"] == 0 
    # choose a random newborn and parent
    return sample(1:S["nod"],2,replace=false)
  end
end

# Simple innovation function that does not take own and others' knowledge into account
function innovate_simple(S,INDREPNB)
  # Number of successful innovations
  local innovateTraits = sample(1:S["nTraits"], # sample from all traits
                                sum(rand(S["nLearnTurns"]) .<= S["ilSuccess"]), # n times 
                                replace=true); # with replacement
  # Update repertoire by adding 
  map(x -> INDREPNB[x] +=1, innovateTraits);
  return(INDREPNB)
end

# Simple copying function that takes neighbourhood repertoire size into account
function copy_simple(S, ADJM, NB, INDREP)
  # In the case of a complex graph, I need to select all actual neighbours
  global neibs = findall(x->x!=0, ADJM[:, NB]); # neighbour ID
  # In the case of a complete graph, I need to select $k$ random neighbours
  if S["pn"] == -1
    global neibs = sample(neibs, S["k"], replace=false);
  end
  # ID of observable neighbours
  global neighbourID = neibs;
  global neighbourRep= INDREP[neighbourID,:]; # neighbour repertoire
  if any(neighbourRep.>0)
    global weight = (sum(neighbourRep.!=0,dims=1) / sum(neighbourRep.!=0)); # weights
    socialTraits = wsample(1:S["nTraits"], 
                          Weights(weight[:]),
                          sum(rand(S["nLearnTurns"]) .<= S["slSuccess"]), # n times SL success
                          replace=true); # sample a trait relative to observation probability
    
    if S["complexContagion"]
      slSuccess = rand(length(socialTraits)) .<= (weight[socialTraits])[:];
      socialTraits = socialTraits[ slSuccess ];
    end

    # Update repertoire by adding 
    map(x -> INDREP[NB,x] +=1, socialTraits);
    # Cap at neighbour maximum (cannot be better than demonstrator(s) through social learning)
    for u in unique(socialTraits)
      INDREP[NB,u] = ifelse( INDREP[NB,u] > maximum(neighbourRep[:,u]), maximum(neighbourRep[:,u]), INDREP[NB,u]);
    end
  end
  return(INDREP[NB,:])
end

  
### Function to run the simulation model
function runit(;adjm=nothing)
  println("Starting simulation 2022_Updated .")
  if report
    # Initialise data structure for reporting summarised data for each complete repetition
    mmaxLevel, mpay, mntraits, mbetadiv, mpn, mpr, mMedNTraits, mMedMaxTraitLevel, mdeg, mpath , mclust , mclustLocal, mclustWeightedAvg, mclustLargest = [zeros(s["repe"]) for _ = 1:14];
    
    # Initialise data structure for summarising results in R
    R"recPNl <- recPRl <- recDegl <- recClustl <- recPathl <- recTraitNl <- recTraitLl <- recNTraitsl <- recBetaDivl <- recPayl <- adjml <- indRepl <- recClustWeightedAvgl <- recClustLargestl <- recRepl <- list()"
  end

  # For loop for repetitions with the same set of parameters
  for reps in 1:s["repe"]
    # Initialise linking parameters (random if set to "rand", otherwise using fixed values)
    global (pn,pr) = setLinking([s["pn"],s["pr"]]);

    # Initialise adjacency matrix (make fully connected if pn = -1)
    if pn == -1
      global adjm = completeGraph(s["nod"]);
    else
      global adjm = initGraph(s["nod"], s["deg"], pn, pr);
    end
    
    # Initialise skill utilities (all=1 if set sigma set to "uniform") # NOTE: patPay > utility
    global utility = bound.( setUtilities(s["sigma"], s["nTraits"]), 0, 10);

    # Initialise individuals
    global indPay          = ones(s["nod"]);
    global indPN           = setValues(s["evolvePN"],s["nod"],pn,.05);
    global indPR           = setValues(s["evolvePR"],s["nod"],pr,.005);
    global indRep      = zeros(Int, s["nod"], s["nTraits"]);
    for i in 1:s["nod"]
      indRep[i, sample(1:s["nTraits"])] = 1;
    end

    # Initialise data structure for recording data, once every 10 generations
    global recgen = 1; # counter for how far into the generation we are
    global recid  = 1; # index for recording data
    if report
      recLength = Int(s["nRound"]/(s["nod"]*10));
      global recPay = zeros(recLength);
      global recNTraits = zeros(recLength);
      global recBetaDiv = zeros(recLength);
      global recPN = zeros(recLength);
      global recPR = zeros(recLength);
      global recDegree = zeros(recLength);
      global recPath = zeros(recLength);
      global recClustLocal = zeros(recLength);
      global recClustWeightedAvg = zeros(recLength);
      global recClustLargest = zeros(recLength);
      global recMedNTraits = zeros(recLength);
      global recMedMaxTraitLevel = zeros(recLength);
      global recRep = zeros(recLength,s["nTraits"]);
      global recPatPay = zeros(recLength,s["nTraits"]);
    end

    # For loop simulating all individual timesteps of a single repetition
    for times in 1:s["nRound"]
      # Moran replacement with mortality or fertility selection:
      global (newborn,parents) = selectionMF(s,indPay);

      # For all simulations with complex networks (excluding fully connected):
      if pn != -1
        # Inherit linking parameters from single parents
        indPN[newborn] = inherit(evolve=s["evolvePN"], parents=indPN[parents], standard=pn, sd=0.1, rmutation=s["mutRate"]);
        indPR[newborn] = inherit(evolve=s["evolvePR"], parents=indPR[parents], standard=pr, sd=0.01, rmutation=s["mutRate"]);

        # Disconnect and reconnect newborn from social network
        adjm[newborn,:] .= 0;
        adjm[:,newborn] .= 0;
        adjm = pbpnprDyn(adjm, indPN[newborn], indPR[newborn], s["nod"], [newborn, parents]); 
      end
      
      # Empty memory of newborn
      indRep[newborn,:] .= 0;

      # Social learning
      indRep[newborn,:] = copy_simple(s, adjm, newborn, indRep);
      
      # Individual learning
      indRep[newborn,:] = innovate_simple(s, indRep[newborn,:]);
      
      # Calculate payoff (only for newborn)
      if s["MM"]
        indPay[newborn] = s["minFit"] + mmkinectics(50,50,sum(indRep[newborn,:] .* utility));
      end
      
      # Every 10 generations ...
      if recgen == s["nod"]*10
              if report
                # ... record results
                repBoo                = indRep.≠0; # turn repertoire into boolean

                recPay[recid]         = mean(indPay);
                recNTraits[recid]     = sum(any(repBoo, dims=1)); # number of all known traits in the population
                recBetaDiv[recid]     = gini(sum(repBoo,dims=1));
                recPN[recid]          = mean(indPN);
                recPR[recid]          = mean(indPR);
                recDegree[recid]      = sum(adjm)/s["nod"];
                recMedNTraits[recid]  = mean(sum(repBoo,dims=2)); # median number of known traits per individual
                recMedMaxTraitLevel[recid] = mean(maximum(indRep, dims=2)); # median of highest trait level for each individual

                # Calculate network metrics (for complex networks only)
                if s["pn"] != -1
                  @rput adjm s;
                  R"""
                  nod <- s$nod
                  net <- graph_from_adjacency_matrix(adjm, mode="undirected")
                  path <- mean(average.path.length(net))
                  # calculate average weighted cluster/component size
                  cs <- clusters(graph=net)$csize
                  tcs <- table(cs) # count for each cluster size
                  csize <- as.numeric(names(tcs)) # cluster size
                  avg_cs <- sum(((tcs*csize) / nod) * csize)
                  """;
                  @rget path avg_cs;
                  recPath[recid] = path;
                  recClustWeightedAvg[recid] = avg_cs;
                end
              end

        recid += 1; # keep increasing
        recgen = 1; # reset 
      else
        recgen += 1;
      end

      # Update fitness for everyone 
      if s["MM"]
        indPay = s["minFit"] .+ mmkinectics.(50,50,indRep * utility); # multiplication of a matrix and vector creates a vector of product sums!
      end

      ## Update environment
      change = rand(s["nTraits"]) .< (s["envTurnover"]/s["nod"]);
      if any(change)
        utility[change] = bound.( setUtilities(s["sigma"], sum(change)), 0, 10);
      end
      
    end
    #END of single simulation 

    # Summarising results of the single simulation run
    if report
      # At the end of a simulation run summarise values as means across the last 20% of recorded generations
      subs=round(Int, (recLength-recLength*.2)):Int(recLength); # keep final 20% of rounds
      repCols = sum(indRep.>0,dims=1); # ,1] is colsums ,2] is rowsums
      mpay[reps]=mean(recPay[subs]);
      mntraits[reps]=mean(recNTraits[subs]);
      mbetadiv[reps]=mean(recBetaDiv[subs]);
      mpn[reps]=mean(recPN[subs]);
      mpr[reps]=mean(recPR[subs]);
      mMedNTraits[reps]=mean(recMedNTraits[subs]);
      mMedMaxTraitLevel[reps]=mean(recMedMaxTraitLevel[subs]);
      mdeg[reps]=mean(recDegree[subs]);
      mpath[reps]=mean(recPath[subs]);
      mclustWeightedAvg[reps] = mean(recClustWeightedAvg[subs]);

      # Store data in an R object
      @rput recPN recPR recRep recDegree recClustLocal recClustWeightedAvg recClustLargest recPath recMedNTraits recMedMaxTraitLevel recNTraits recBetaDiv recPay reps s adjm indRep indPN indPR; 
      R"""
      recPNl[[reps]] <- recPN
      recPRl[[reps]] <- recPR
      recDegl[[reps]] <- recDegree
      recClustWeightedAvgl[[reps]] <- recClustWeightedAvg
      recPathl[[reps]] <- recPath
      recTraitNl[[reps]] <- recMedNTraits
      recTraitLl[[reps]] <- recMedMaxTraitLevel
      recNTraitsl[[reps]] <- recNTraits
      recBetaDivl[[reps]] <- recBetaDiv
      recPayl[[reps]] <- recPay
      adjml[[reps]] <- adjm
      indRepl[[reps]] <- indRep
      """
    end


  end
  #END of single repetition
          if report
            # Summarising all data from all repetitions and storing them in an *.Rdata file
            pat=string(paths, "/output/");
            if !isdir(pat)
              mkdir(pat)
            end

            outpath = string(pat,"envturnover_rout_",queue);

            @rput outpath mmaxLevel mpay mntraits mbetadiv mpn mpr mMedNTraits mMedMaxTraitLevel mclustWeightedAvg mclustLargest mdeg mpath mclustLocal indRep indPN indPR adjm; 

            R"""
            setup <- s

            setup$recPay <- mpay
            setup$recNTraits <- mntraits
            setup$recBetaDiv <- mbetadiv
            setup$recPN <- mpn
            setup$recPR <- mpr
            setup$recMedNTraits <- mMedNTraits
            setup$recMedMaxTraitLevel <- mMedMaxTraitLevel
            setup$recDegree <- mdeg
            setup$path <- mpath
            setup$clustWeightedAvg <- mclustWeightedAvg

            setup <- do.call(cbind,setup)
            save(setup, indRepl, adjml, indPN, indPR, recPNl, recPRl, recRepl, recDegl, recClustl, recClustWeightedAvgl, recClustLargestl, recPathl, recTraitNl, recTraitLl, recNTraitsl, recBetaDivl, recPayl, file=outpath)
            """
          end

  print(".. done with ", queue,". ")
  return true
end