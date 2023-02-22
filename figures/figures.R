library(ggplot2)
library(cowplot); theme_set(theme_cowplot())
library(dplyr)
library(tidyr)
library(viridis)
library(ggsci)
library(reshape2)

WD <- "PATH/TO/WORKING_DIRECTORY"
folder <- "NAME_OF_FOLDER_WITH_SUMMARY_FILE"
setwd(paste(WD,folder,sep=""))

### Figure 1 (excluding networks from inset in b) ----
# This figure requires the summary file that is created after executing summariseFiles_rout.R, which summarises the results files from the output folder, which is created by the run.jl file. 
# Set parameters for this simulations as follows: $\alpha=0.01$, $\beta=1$, $N=100$, $M=500$, $5{,}000$ generations, $\tau=0$ and $\sigma=0$, and 200 repetitions.
load(paste(folder,"summary",sep="_"))
  
data <- as.data.frame(pnpr2_data) 
q <- ggplot(data=data) +
  theme(axis.line=element_blank(), 
        legend.position="",
        strip.text=element_text(face="bold"),
        strip.background=element_blank()
  ) + 
  panel_border(size=1, colour="black") + 
  labs(color="Proficiency") +
  scale_color_gradient(low="#3B4992",high="#EE0000",guide="colourbar")
## Plots A-D
pp1 <- q + geom_point(aes(x=recPR, y=recPN, color=recMedMaxTraitLevel), size=1.5) + 
  xlab(expression("Random linking,"~p[r])) + 
  ylab(expression("Socal inheritance,"~p[n])) +
  scale_x_log10() + annotation_logticks(sides = "b")
  ylab("Average payoff") + 
  xlab("Avg. weighted component size") +
  theme(legend.position = c(.15,.7))
pp2 <- q + geom_point(aes(x=clustWeightedAvg, y=recDegree, color=recMedMaxTraitLevel), size=1.5) + 
  xlab("Avg. weighted component size") + 
  ylab("Degree")
pp3 <- q + geom_point(aes(x=recMedNTraits, y=recMedMaxTraitLevel, color=recMedMaxTraitLevel), size=1.5) +
  xlab("Average repertoire size") + 
  ylab("Mean highest proficiency")
pp4 <-  q + geom_point(aes(y=recPay, x=recNTraits, color=recMedMaxTraitLevel), size=1.5) + 
  ylab("Average payoff") + 
  xlab("Population trait count") +
  theme(legend.position = c(.58,.75))
### Scatter plot for turnover = 0 and sigma = 0, i.e. for the case where there is no environmental difference
save_plot(plot_grid(pp1,pp2,pp3,pp4, labels="auto", nrow=1), filename=paste(folder,"Fig1.pdf",sep="_"), base_height=4, base_width=14)


### Figure 2 ----
# This figure also requires the summary file which is the output from summariseFiles_rout.R 
load(paste(folder,"summary",sep="_"))
pnprpay <- as.data.frame(pnprpay)

## A, Overview plot showing where populations have been in the alst 20% of generations:
dat_00_end <- filter(pnprpay, sigma==0 & envTurnover==0 & time>400)
ggplot(dat_00_end) + 
  geom_hline(linetype=2, col="black", yintercept = c(.2,.5,.8)) +
  geom_point(aes(x=pr, y=pn, col=pay), alpha=.5) + 
  geom_point(data=expand.grid(x=c(.001,.01,.02,.04,.06), y=c(.2,.5,.8)), aes(x=x, y=y), col="red", shape=15) +
  xlab(expression("Random linking,"~p[r])) +
  ylab(expression("Socal inheritance,"~p[n])) +
  scale_color_viridis_c("Payoff", limits=c(7,NA)) +
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1), limits=c(NA, .5)) +
  annotation_logticks(sides = "b") +
  background_grid() +
  theme(axis.line = element_blank(),
        legend.position = c(.05,.7),
        panel.border = element_rect(colour = "black",linewidth = 1.2)) -> pnpr_pay_all

## B, Plots from horizontal cuts through A
dat_00_cut <- dat_00_end %>% mutate(pn_r = round(pn, 1)) %>% filter(pn_r %in% c(.2,.5,.8)) 
dat_00_cut$pn_r <- factor(dat_00_cut$pn_r, levels=c(.8,.5,.2))
ggplot(dat_00_cut) + geom_point(aes(x=pr, y=pay, col=pay), alpha=.4) +
  scale_color_viridis_c("Payoff", limits=c(7,NA)) +
  ylab("Average payoff") +
  xlab(expression("Random linking,"~p[r])) +
  facet_grid("pn_r", labeller = label_both) +
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1), limits=c(NA, .5)) +
  annotation_logticks(sides = "b") +
  background_grid() +
  theme(axis.line = element_blank(),
        legend.position = "",
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black",size = 1.2)) -> pr_cutpn_pay
save_plot(
  plot_grid(pnpr_pay_all, pr_cutpn_pay, rel_widths = c(2/3,1/3), labels="auto"),
  filename=paste(folder,"Fig2.pdf",sep="_"), base_height = 6, base_width = 12, bg="white")


### Figure 4 ----
# This figure also requires the summary file which is the output from summariseFiles_rout.R 
# Set parameters for this simulations as follows: $N=100$, $M=500$, $5{,}000$ generations, $\beta = [0 .25 .5 .75 1]$, and $\alpha = [.1 .05 .01 .005 .001]$. 
dats <- as.data.frame(pnpr2_data, stringsAsFactors = F) 
dats %>% group_by(ilSuccess, slSuccess) %>% summarise(meanPay=mean(recPay), lci=quantile(recPay, .05), hci=quantile(recPay, .95)) -> pay_dat
q2 <- ggplot(data=dats) +
  theme(axis.line=element_blank(), 
        strip.text=element_text(face="bold"),
        strip.background=element_blank(),
        text = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12)
  ) + 
  panel_border(size=1, colour="black") + 
  scale_color_manual(values = c("#003f5c","#58508d","#bc5090","#ff6361","#ffa600"), name="SL success") + 
  scale_fill_manual(values = c("#003f5c","#58508d","#bc5090","#ff6361","#ffa600"))+#, name="IL success") + 
  guides(colour = guide_legend(override.aes = list(alpha = 1)))

## A, Payoff over alpha
q2 + 
  geom_point(data=pay_dat, aes(x=ilSuccess, y=meanPay, col=factor(slSuccess))) + 
  geom_line(data=pay_dat, aes(x=ilSuccess, y=meanPay, col=factor(slSuccess)), linewidth=.8) + 
  geom_errorbar(data=pay_dat, aes(ymin=lci, ymax=hci, x=ilSuccess, col=factor(slSuccess)), width=.002) +
  theme(legend.position=c(.75, .25)) +
  ylab("Average payoff") + xlab(expression("Innovation success rate"~alpha)) -> ilsl_pay_plot

## B, Weighted component size over alpha
q2 + 
  geom_point(aes(x=factor(ilSuccess), y=clustWeightedAvg, color=factor(slSuccess)), position = position_dodge(width = .5), alpha=.25) +
  ylab("Avg. weighted component size") + xlab(expression("Innovation success rate"~alpha))  + 
  theme(legend.position = "") +
  guides(col=guide_legend("SL success")) -> ilsl_clust_plot

ggdraw(ilsl_clust_plot) + 
  draw_label("<- Connected            Unconnected ->", x=.96, y=.58, size = 11, angle = 270) -> ilsl_clust_plot_labels

save_plot(plot_grid(ilsl_pay_plot, ilsl_clust_plot_labels, labels="auto"), 
          filename = paste(folder,"Fig4.pdf",sep="_"), base_width = 10, base_height = 4)


### Figure 5 ----
# This figure also requires the summary file which is the output from summariseFiles_rout.R 
# Set parameters for this simulations as follows: $N=100$, $M=500$, $5{,}000$ generations, $\sigma = [.2 1]$, $\tau = [.001 1]$. 
data <- as.data.frame(pnpr2_data)

## A, Raster plots
data %>% group_by(envTurnover,sigma) %>% summarise(
  ClustWeightedAvg=mean(clustWeightedAvg),
  TotalTraits=mean(recNTraits),
  Repertoire=mean(recMedNTraits),
  Level=mean(recMedMaxTraitLevel),
) -> pnpr2_s

p <- ggplot(data=pnpr2_s, aes(y=factor(sigma), x=factor(envTurnover))) +
  ylab(expression(sigma)) +
  xlab(expression("Turnover "~tau)) +
  coord_equal() +
  theme(axis.line=element_blank()) +
  panel_border(size=1, colour="black") +
  scale_x_discrete(expand=c(0,0), breaks=10^seq(-4,0,1), labels=expression(10^-4,10^-3,10^-2,10^-1,10^0)) + scale_y_discrete(expand=c(0,0), breaks=seq(0,1,.2))

q1 <- p + geom_raster(aes(fill=Level)) + scale_fill_viridis(limits=c(1,max(pnpr2_s[, "Level"])), name="") + labs(title="Proficiency")
q2 <- p + geom_raster(aes(fill=Repertoire)) + scale_fill_viridis(limits=c(1,max(pnpr2_s[, "Repertoire"])),name="") + labs(title="Repertoire size")
q3 <- p + geom_raster(aes(fill=ClustWeightedAvg)) + scale_fill_viridis(name="") + labs(title="Avg. weigh. component size")
q4 <- p + geom_raster(aes(fill=TotalTraits)) + scale_fill_viridis(name="") + labs(title="Total skill count")#breaks=c(70,75,80),

fig5a <- plot_grid(q1,q2,q3,q4, nrow=2, labels="auto")


## B, Scatter plot with false colour scale 
data$sigma <- factor(data$sigma, levels=sort(unique(data$sigma), decreasing=T), labels=c("high variance", "low variance"))
data$envTurnover <- factor(data$envTurnover, levels=sort(unique(data$envTurnover), decreasing=F), labels=c("slow turnover", "fast turnover"))
## Calculate combined measure
data$combMeasure <- (-data[,"recMedNTraits"]/max(data[,"recMedNTraits"])) + (data[,"recMedMaxTraitLevel"]/max(data[,"recMedMaxTraitLevel"]))
q <- ggplot() + 
  facet_grid(sigma~envTurnover) +
  scale_colour_gradientn(colours=c("#3B4992","grey","#EE0000"), 
                         name="", 
                         breaks=c(min(data$combMeasure), max(data$combMeasure)),
                         label=c("large\nrepertoire","high\nproficiency")) +
  theme(axis.text=element_text(size=17),
        axis.title=element_text(size=20),
        axis.line=element_blank(),
        legend.text=element_text(size=8),
        strip.background=element_blank(),
        strip.text=element_text(face="bold", size=17)) +
  panel_border(size=1, colour="black")

q + 
  geom_point(data=data, aes(x=recPR, y=recPN, col=combMeasure)) + 
  theme(legend.position=c(.05,.8), 
        legend.direction = "vertical", 
        legend.key.width=unit(3, "mm"),
        legend.key.height=unit(6, "mm")) +
  xlab(expression("Random linking,"~p[r])) + 
  ylab(expression("Socal inheritance,"~p[n])) +
  scale_x_log10(breaks=c(0.01,0.1)) + annotation_logticks(sides = "b") -> fig5b

save_plot(plot_grid(fig5a, fig5b, labels=c("","e")), filename = paste(folder,"Fig5.png",sep="_"), base_height=8*.8, base_width=14, bg="white")