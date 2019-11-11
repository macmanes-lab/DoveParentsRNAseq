# custom themes

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
maniplevels <- c( "remove.d03" ,  "remove.d09" ,  "remove.d17" , "remove.d20" , 
                  "extend" , "prolong", "early")
maniplevels1 <- c( "m.inc.d3" ,  "m.inc.d8" ,  "m.inc.d9" , "m.inc.d17" , "m.n2", 
                   "prolong" , "extend")

alllevels <- c("control", "bldg", "lay", "inc.d3", "m.inc.d3" ,  
               "m.inc.d8" , "inc.d9", "m.inc.d9" ,
               "inc.d17",  "m.inc.d17", "hatch", 
               "m.n2",  "prolong" , "extend",
               "n5", "n9")


# custom colors

colorsmanip <- c( "m.inc.d3" = "#CDCDCD", 
                         "m.inc.d8" = "#6DA9D3", 
                         "m.inc.d9" = "#959595", 
                         "m.inc.d17" = "#626262",
                         "m.n2" = "#262625", 
                         "prolong" = "#1B86B9" , 
                         "extend" = "#01588C")

colorschar <-  c("control" = "#F8766D", 
                         "bldg"= "#D39200", 
                         "lay"= "#93AA00", 
                         "inc.d3"= "#00BA38", 
                         "inc.d9"= "#00C19F", 
                         "inc.d17"= "#00B9E3", 
                         "hatch"= "#619CFF", 
                         "n5"= "#DB72FB", 
                         "n9"= "#FF61C3")


colorscharmaip <-  c("control" = "#F8766D", 
                 "bldg"= "#D39200", 
                 "lay"= "#93AA00", 
                 "inc.d3"= "#00BA38", 
                 "inc.d9"= "#00C19F", 
                 "inc.d17"= "#00B9E3", 
                 "hatch"= "#619CFF", 
                 "n5"= "#DB72FB", 
                 "n9"= "#FF61C3",
                 "m.inc.d3" = "#CDCDCD", 
                 "m.inc.d8" = "#6DA9D3", 
                 "m.inc.d9" = "#959595", 
                 "m.inc.d17" = "#626262",
                 "m.n2" = "#262625", 
                 "prolong" = "#1B86B9" , 
                 "extend" = "#01588C")