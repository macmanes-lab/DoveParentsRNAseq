# custom themes for plots

# hex colors for last full day of parents life

colorlastday <-  c("control" = "#77787A", 
                   "nest building" = "#DA733E", 
                   "eggs lay" = "#559FD9",
                   "eggs early" = "#9AD0C3",
                   "eggs middle" = "#68C2A8", 
                   "eggs later" = "#37A589",  
                   "eggs delay" = "#19977D",
                   "chicks hatch" = "#BAB5D2",   
                   "chicks early" = "#9993C3", 
                   "chicks later" = "#7874B0",
                   "empty nest" = "#754E27")

# hex colors for last 2nd to last day of parents life
colorpenultimate <- c("control" = "#77787A", 
                      "nest building" = "#DA733E", 
                      "eggs lay" = "#559FD9",
                      "eggs early" = "#9AD0C3",
                      "eggs middle" = "#68C2A8", 
                      "eggs later" = "#37A589",  
                      "eggs delay" = "#19977D",
                      "chicks hatch" = "#BAB5D2",   
                      "chicks early" = "#9993C3", 
                      "chicks later" = "#7874B0")


myannotationscolors <- list(lastday = c("control" = "#77787A", 
                                        "nest building" = "#DA733E", 
                                        "eggs lay" = "#559FD9",
                                        "eggs early" = "#9AD0C3",
                                        "eggs middle" = "#68C2A8", 
                                        "eggs later" = "#37A589",  
                                        "eggs delay" = "#19977D",
                                        "chicks hatch" = "#BAB5D2",   
                                        "chicks early" = "#9993C3", 
                                        "chicks later" = "#7874B0",
                                        "empty nest" = "#754E27"),
                       penultimate = c("control" = "#77787A", 
                                       "nest building" = "#DA733E", 
                                       "eggs lay" = "#559FD9",
                                       "eggs early" = "#9AD0C3",
                                       "eggs middle" = "#68C2A8", 
                                       "eggs later" = "#37A589",  
                                       "eggs delay" = "#19977D",
                                       "chicks hatch" = "#BAB5D2",   
                                       "chicks early" = "#9993C3", 
                                       "chicks later" = "#7874B0"))


# themes for whole plots

theme_noaxislabels <- theme(legend.position = "none", 
                            axis.title=element_blank(), axis.ticks = element_blank(), axis.text=element_blank(),
                            panel.background = element_blank()) 

theme_siggenes <- theme(legend.position = c(0.2,0.3), 
                            #axis.title.x=element_blank(), 
                        axis.ticks = element_blank(), axis.text=element_blank(),
                        panel.background = element_blank())

theme_nolegend <- theme(legend.position = "none",  axis.text.x=element_blank(), axis.ticks = element_blank()) 

theme_rmh <- function(){ 
  theme_bw(base_size=14) +
    theme(
      panel.grid.minor.x  = element_blank(),
      panel.grid.minor.y  = element_blank(),
      strip.background = element_rect(colour="white", fill="white"),
      legend.position = "right",           
      legend.margin=margin(t=-0.1, r=0, b=-0.1, l=-0.1, unit="cm"),
      legend.key.size = unit(0.5, "cm"))
}

theme_set(theme_classic())