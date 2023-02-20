suppressMessages(library(dplyr))
library(reshape2)

path <- getwd()

if(file.exists(paste(path,"output/",sep="/"))){
  ## Define functions
  loadData <- function(PATH, PNPR=FALSE, JULIA=FALSE){
    summaryFilePath <- paste(path, paste(tail(strsplit(path, split="/")[[1]],1),"summary",sep="_"), sep="/") # Set path to save results file
    pnpr2_data <- do.call(rbind, lapply(list.files(PATH), function(f){ # for all simulations do
      load(paste(PATH, f, sep="")) # Load file
      return(setup) # Return the setup file (contains summarised results for all files)
    }))

    pnprpay <- do.call(rbind, lapply(list.files(PATH), function(f){ # for all simulations do
      load(paste(PATH, f, sep="")) # Load file
      do.call("rbind", 
        lapply(1:length(recPNl), FUN = function(i){
          pnl <- recPNl[[i]]
          prl <- recPRl[[i]]
          pal <- recPayl[[i]]

          subs <- 1: length(pnl)
          dats <- cbind(
            sigma = setup[1,"sigma"],
            envTurnover = setup[1,"envTurnover"],
            nod = setup[1,"nod"],
            seed = round(runif(1,0.1,.9)*1e6),
            time = subs,
            pn = pnl[subs],
            pr = prl[subs], 
            pay = pal[subs], row.names=NULL)
            return(dats)
          })
        )
    }))

    save(pnpr2_data, pnprpay, file=summaryFilePath) # Save new file or append to previous
  }

  loadRepertoireData <- function(PATH){
    summaryFilePath3 <- paste(path, paste(tail(strsplit(path, split="/")[[1]],1),"traitHistogram",sep="_"), sep="/") # Set path to save results file
    traitHistogram_c <- do.call(rbind, lapply(list.files(PATH), function(f){ # for all simulations do
      load(paste(PATH, f, sep="")) # Load file
      
      # Results
      lapply(indRepl, function(ind){
        tabletable <- as.data.frame(apply(ind!=0, 1, which) |> unlist() |> table() |> table())
        colnames(tabletable) <- c("rank","freq")
        cbind(
          tabletable,
          data.frame(
            # State values
            nod=setup[1,"nod"],
            sigma=setup[1,"sigma"],
            envTurnover=setup[1,"envTurnover"],
            avgMaxProf = mean(apply(ind, 1, max)),
            avgRepSize = mean(rowSums(ind!=0)),
            # ID values
            seed=sample(x=LETTERS, size=10, replace=T) %>% paste(collapse="")),
          row.names = NULL
        )
      }) %>% bind_rows() -> reps
      return(reps) # Return the setup file (containts summarised results for all files)
    }))

    save(traitHistogram_c, 
    file=summaryFilePath3) # Save new file or append to previous
  }

  loadTimeData <- function(PATH){
    summaryFilePath2 <- paste(path, paste(tail(strsplit(path, split="/")[[1]],1),"timeData",sep="_"), sep="/") # Set path to save results file
    timeData <- do.call(rbind, lapply(list.files(PATH), function(f){ # for all simulations do
      load(paste(PATH, f, sep="")) # Load file
      return(
        data.frame(
        row.names=NULL,
        # State values
          nod=setup[1,"nod"],
          sigma=setup[1,"sigma"],
          envTurnover=setup[1,"envTurnover"],
          ilSuccess=setup[1,"ilSuccess"],
          slSuccess=setup[1,"slSuccess"],
        # ID values
          time=1:length(recPNl[[1]]),
          seed=lapply(1:setup[1,"repe"], function(y) sample(x=LETTERS, size=10, replace=T)) %>% lapply(paste,collapse="") %>% unlist %>% rep(each=length(recPNl[[1]])),
        # Results
          recPN=melt(recPNl)[,"value"],
          recPR=melt(recPRl)[,"value"],
          recDeg=melt(recDegl)[,"value"],
          recClustWeightedAvg=melt(recClustWeightedAvgl)[,"value"],
          recPath=melt(recPathl)[,"value"],
          recTraitN=melt(recTraitNl)[,"value"],
          recTraitL=melt(recTraitLl)[,"value"],
          recNTraits=melt(recNTraitsl)[,"value"],
          recBetaDiv=melt(recBetaDivl)[,"value"],
          recPay=melt(recPayl)[,"value"]
        )
      ) # Return the setup file (contains summarised results for all files)
    }))
    save(timeData, file=summaryFilePath2) # Save new file or append to previous
  }


  ## Run summarising function
  loadData(PATH=paste(path,"output/",sep="/"), PNPR=TRUE, JULIA=TRUE)
  loadRepertoireData(PATH=paste(path,"output/",sep="/"))
  loadTimeData(PATH=paste(path,"output/",sep="/"))
  print("Done!")
} else {
  print("There is no output directory to summarise (probably already done that!) :D ")
}
