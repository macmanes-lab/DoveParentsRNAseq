# custom themes

## hex colors chosen using https://colorbrewer2.org/ 
## hex colors blended using from https://meyerweb.com/eric/tools/color-blend/#F8766D:D39200::hex 

theme_B3 <- function () { 
  theme_classic(base_size = 7,
                base_family = 'Helvetica') +
    theme(
      strip.background  = element_blank(),
      panel.grid.major  = element_blank(),  # remove major gridlines
      panel.grid.minor  = element_blank()  # remove minor gridlines
    )
}



# custom levels

# characterization
charlevels <- c("control", "bldg", "lay", "inc.d3", "inc.d9", "inc.d17", "hatch", "n5", "n9")
charlevelsnocontrol <- c( "bldg", "lay", "inc.d3", "inc.d9", "inc.d17", "hatch", "n5", "n9")
maniplevels <- c("m.inc.d3" ,  "early" ,  "m.inc.d9" , "prolong" , "m.inc.d17" , "extend", "m.n2")

# manipulation
timinglevels <- c("control", "bldg", "inc.d9", "early" , "inc.d17", "prolong" , "hatch", "extend", "n5")
removallevels <- c("control", "bldg","inc.d3",  "m.inc.d3" , "inc.d9", "m.inc.d9" ,
             "inc.d17", "m.inc.d17" , "hatch", "m.n2")


# mix of characterization and manipulation, ordered by time
alllevels <- c("control", "bldg", "lay", "inc.d3", "m.inc.d3" ,  
               "inc.d9", "m.inc.d9" , "early" ,
               "inc.d17",  "m.inc.d17", "prolong" , 
               "hatch",  "m.n2", "extend",
               "n5",  
               "n9")

# characteriation first, then manip by time
alllevels2 <- c("control", "bldg", "lay", "inc.d3", "inc.d9", 
                "inc.d17", "hatch", "n5", "n9",
                "m.inc.d3" ,  
                "early" ,  "m.inc.d9" , 
                "prolong" , "m.inc.d17" , 
                "extend", "m.n2")

hypothesislevels <- c(  "nest",  "eggs", "chicks" ,  "lo" , "hi" , "early", "late")

# tissue
tissuelevels <- c("hypothalamus", "pituitary", "gonads")
tissuelevel <- c("hypothalamus", "pituitary", "gonad")
myshapes = c("hypothalamus" = 20,  "pituitary" = 17,  "gonads" = 15)

# sex
sexlevels <- c("female", "male")

levelsnewgrouping <-  c("control" , "bldg",  "earlyinc", 
                        "nearhatch",  "chickcare", "manip")

levelsextint <- c("nest", "eggs" , "chicks" )

# DESeq
comparisonlevelschar <- c("bldg_lay",  "lay_inc.d3",  
                          "inc.d3_inc.d9",   "inc.d9_inc.d17",  
                          "inc.d17_hatch", "hatch_n5", "n5_n9")

comparisonlevelscontrol <- c("control_bldg", "control_lay",  "control_inc.d3",  
                          "control_inc.d9",   "control_inc.d17",  
                          "control_hatch", "control_n5", "control_n9")

comparisonlevelsbldg <- c( "bldg_lay",  "bldg_inc.d3",  
                             "bldg_inc.d9",   "bldg_inc.d17",  
                             "bldg_hatch", "bldg_n5", "bldg_n9")

comparisonlabelschar <- c("bldg vs.\nlay" ,"lay vs.\ninc.d3", 
                           "inc.d3 vs.\ninc.d9", "inc.d9 vs.\ninc.d17", 
                           "inc.d17 vs.\nhatch" , "hatch vs.\nn5", "n5 vs.\nn9")

comparisonlabelscharnobldg <- c("lay vs.\ninc.d3", 
                          "inc.d3 vs.\ninc.d9", "inc.d9 vs.\ninc.d17", 
                          "inc.d17 vs.\nhatch" , "hatch vs.\nn5", "n5 vs.\nn9")

comparisonlabelscontrol <- c("control\nvs. bldg" ,"control\nvs. lay", "control\nvs. inc.d3", 
                                "control\nvs. inc.d9", "control\nvs. inc.d17", 
                                "control\nvs. hatch" , "control\nvs. n5", "control\nvs. n9")

comparisonlabelsbldg <- c("bldg vs.\nlay", "bldg vs.\ninc.d3", 
                                "bldg vs.\ninc.d9", "bldg vs.\ninc.d17", 
                                "bldg vs.\nhatch" , "bldg vs.\nn5", "bldg vs.\nn9")


comparisonlevelsremoval <- c("ind.d3_m.inc.d3", "ind.d9_m.inc.d9", "ind.d19_m.inc.d17", "hatch.m.n2")

# custom colors

colorschar <-  c("control" = "#cc4c02", 
                 "bldg"= "#fe9929", 
                 "lay"= "#fed98e", 
                 "inc.d3"= "#78c679", 
                 "inc.d9"= "#31a354", 
                 "inc.d17"= "#006837", 
                 "hatch"= "#6baed6",
                 "n5"= "#3182bd", 
                 "n9"= "#08519c" )

colorsmanip <- c("m.inc.d3" = "#CDCDCD", 
                 "m.inc.d9" = "#959595", 
                 "m.inc.d17" = "#626262",
                 "m.n2" = "#262625", 
                 "early" = "#cbc9e2", 
                 "prolong" = "#9e9ac8" , 
                 "extend" = "#6a51a3" )

colorscharmaip <-  c(colorschar, colorsmanip)

colorsnewgrouping <-  c("control" = "#F8766D", 
                 "bldg"= "#D39200", 
                 "earlyinc"= "#00C19F", 
                 "nearhatch"= "#619CFF", 
                 "chickcare"= "#FF61C3",
                 "manip" = "#959595")

colorhypothesis <- c("none" = "#cc5500",
                     "loss" = "#8B4513",
                     "eggs" = "#7FA163",
                     "chicks" = "#3A80B9",
                     "lo" = "#B5C4B2",
                     "hi" = "#19757A",
                     "early" = "#B5C4B2",
                     "late" = "#19757A",
                     "reference" = "#F8766D" )

levelhypothesis <- c("reference", "eggs", "chicks", "loss", "early", "late")

colorstissue <- c("hypothalamus" = "#d95f02",
                  "pituitary" = "#1b9e77",
                  "gonads" =  "#7570b3",
                  "gonad" =  "#7570b3")

sexcolors <- c("female" = "#969696", "male" = "#525252")

colorsvolcano <-  c(colorschar,
                    
                    "early" = "#cbc9e2", 
                    "prolong" = "#9e9ac8" , 
                     "extend" = "#6a51a3",
                    
                    "male" = "#525252",
                    "female" = "#969696",
                    
                    "NS" = "#bdbdbd")

 
allcolors <- c(colorscharmaip, sexcolors, 
               colorstissue, colorhypothesis)
 
colorsvolcanochar <-  c(colorschar,
                    "NS" = "#bdbdbd")

myannotationcolors <- list(sex = sexcolors,
                        tissue = colorstissue,
                        treatment = colorscharmaip)
