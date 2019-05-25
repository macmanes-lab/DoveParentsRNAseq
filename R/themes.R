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

theme_noaxislabels <- theme(legend.position = "none", axis.title.x=element_blank(), axis.ticks = element_blank(), axis.text.x=element_blank()) 
theme_nolegend <- theme(legend.position = "none",  axis.text.x=element_blank(), axis.ticks = element_blank()) 
