# custom themes

theme_B3 <- function () { 
  theme_classic(base_size = 11) +
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

maniplevels <- c( "remove.d03" ,  "remove.d09" ,  "remove.d17" , "remove.d20" ,  "extend" , "prolong", "early")

maniplevels1 <- c( "m.inc.d3" ,  "m.inc.d8" ,  "m.inc.d9" , "m.inc.d17" , "m.n2", "prolong" , "extend")

levelstiming <- c( "m.inc.d8" , "prolong" , "extend")
levelsremoval <- c( "m.inc.d3" ,    "m.inc.d9" , "m.inc.d17" , "m.n2")

controlstiming <- c( "inc.d9" , "hatch" , "n5")
controlsremoval <- c( "inc.d3" ,    "inc.d9" , "inc.d17" , "hatch")


comparisonlevels <- c("control.bldg", "bldg.lay" ,"lay.inc_d3", "inc_d3.inc_d9",
                      "inc_d9.inc_d17", "inc_d17.hatch" , "hatch.n5", "n5.n9"   )

alllevels <- c("control", "bldg", "lay", "inc.d3", "m.inc.d3" ,  
                "inc.d9", "m.inc.d8" ,"m.inc.d9" ,
               "inc.d17",  "m.inc.d17","prolong" ,
               "hatch",  "m.n2",   "extend",
               "n5", "n9")

alllevels2 <- c("control", "bldg", "lay", "inc.d3", "inc.d9", 
                "inc.d17", "hatch", "n5", "n9",
                "m.inc.d3" ,  "m.inc.d8" ,  "m.inc.d9" , 
                "m.inc.d17" , "m.n2", "prolong" , "extend")

alllevels3 <- c("control", "bldg", "lay", "inc.d3", "inc.d9", 
                "inc.d17", "hatch", "n5", "n9",
                 "m.inc.d8" , "prolong" , "extend",
                "m.inc.d3" ,  "m.inc.d9" , "m.inc.d17" , "m.n2")

tissuelevels <- c("hypothalamus", "pituitary", "gonads")
sexlevels <- c("female", "male")


charlevelsvocano <- c("control", "bldg", "lay", "inc.d3", "inc.d9", "inc.d17", "hatch", "n5", "n9", "NS")

alllevels3volcano <- c("control", "bldg", "lay", "inc.d3", "inc.d9", 
                      "inc.d17", "hatch", "n5", "n9",
                      "m.inc.d3" ,  "m.inc.d9" , "m.inc.d17" , "m.n2",
                      "m.inc.d8" , "prolong" , "extend",
                      "NS")

levelsnewgrouping <-  c("control" , "bldg",  "earlyinc", 
                        "nearhatch",  "chickcare", "manip")

## comparisons

serialtimepoints <- c("control.bldg" , "bldg.lay", "lay.inc.d3", "inc.d3.inc.d9", 
                      "inc.d9.inc.d17", "inc.d17.hatch", "hatch.n5", "n5.n9")

manipVchar <- c("inc.d3.m.inc.d3", "inc.d9.m.inc.d8", "inc.d9.m.inc.d9","inc.d17.m.inc.d17",
                "hatch.prolong", "hatch.extend", "hatch.m.n2")

offspringremoval <- c("lay.inc.d3", "inc.d3.m.inc.d3", 
                      "inc.d3.inc.d9", "inc.d9.m.inc.d9",
                      "inc.d9.inc.d17", "inc.d17.m.inc.d17", 
                      "hatch.n5", "hatch.m.n2")

prolongdelay <- c( "inc.d3.inc.d9", "inc.d9.inc.d17", 
                   "inc.d9.m.inc.d8",
                   "inc.d17.hatch", "hatch.n5",
                   "hatch.prolong", "hatch.extend")

## custom shapes

myshapes = c("hypothalamus" = 20,  "pituitary" = 17,  "gonads" = 15)

# custom colors

colorsmanip <- c( "m.inc.d3" = "#CDCDCD", 
                         "m.inc.d8" = "#fcc5c0", 
                         "m.inc.d9" = "#959595", 
                         "m.inc.d17" = "#626262",
                         "m.n2" = "#262625", 
                         "prolong" = "#f768a1" , 
                         "extend" = "#f768a1")

colorschar <-  c("control" = "#F8766D", 
                         "bldg"= "#D39200", 
                         "lay"= "#93AA00", 
                         "inc.d3"= "#00BA38", 
                         "inc.d9"= "#00C19F", 
                         "inc.d17"= "#00B9E3", 
                         "hatch"= "#619CFF", 
                         "n5"= "#DB72FB", 
                         "n9"= "#FF61C3")

colorscharnew <-  c("control" = "#cc4c02", 
                 "bldg"= "#fe9929", 
                 "lay"= "#fed98e", 
                 "inc.d3"= "#78c679", 
                 "inc.d9"= "#31a354", 
                 "inc.d17"= "#006837", 
                 "hatch"= "#08519c",
                 "n5"= "#3182bd", 
                 "n9"= "#6baed6")

colorsnewgrouping <-  c("control" = "#F8766D", 
                 "bldg"= "#D39200", 
                 "earlyinc"= "#00C19F", 
                 "nearhatch"= "#619CFF", 
                 "chickcare"= "#FF61C3",
                 "manip" = "#959595")

colorsmanip <- c("m.inc.d3" = "#CDCDCD", 
                   "m.inc.d9" = "#959595", 
                   "m.inc.d17" = "#626262",
                   "m.n2" = "#262625", 
                   "m.inc.d8" = "#cbc9e2", 
                   "prolong" = "#9e9ac8" , 
                   "extend" = "#6a51a3" )

colorscharmaip <-  c(colorscharnew, colorsmanip)

colorscharmaip2 <-  c(colorscharnew, colorsmanip)

colorstissue <- c("hypothalamus" = "#d95f02",
                  "pituitary" = "#1b9e77",
                  "gonads" =  "#7570b3")

sexcolors <- c("female" = "#969696", "male" = "#525252")

colorsvolcano <-  c(colorscharnew,
                    
                    "m.inc.d8" = "#cbc9e2", 
                    "prolong" = "#9e9ac8" , 
                     "extend" = "#6a51a3",
                    
                    "male" = "#525252",
                    "female" = "#969696",
                    
                    "NS" = "#bdbdbd")

colorsvolcanochar <-  c(colorscharnew,
                    "NS" = "#bdbdbd")

myannotationcolors <- list(sex = sexcolors,
                        tissue = colorstissue,
                        treatment = colorscharmaip2)
