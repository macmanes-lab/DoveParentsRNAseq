# custom theme

theme_B3 <- function () { 
  theme_classic(base_size = 7) +
    theme(
      panel.grid.major  = element_blank(),  # remove major gridlines
      panel.grid.minor  = element_blank(),  # remove minor gridlines
      plot.title = element_text(hjust = 0.5, face = "bold") # center & bold 
    )
}


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

mytheme <- function(){
  theme_minimal(base_size = 12) + 
    theme(panel.grid = element_blank(),
          panel.background = element_rect(fill = "transparent", colour = NA),
          plot.background = element_rect(fill = "transparent", colour = NA),
          legend.title = element_blank())
}


# custom levels

charlevels <- c("control", "bldg", "lay", "inc.d3", "inc.d9", "inc.d17", "hatch", "n5", "n9")
maniplevels <- c( "remove.d03" ,  "remove.d09" ,  "remove.d17" , "remove.d20" , "extend" , "prolong", "early")

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


