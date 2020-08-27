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
maniplevels <- c("m.inc.d3" ,  "early" ,  "m.inc.d9" , "m.inc.d17", "prolong" ,  "m.n2", "extend")

# manipulation
timinglevels <- c("inc.d9", "early" , "inc.d17", "prolong" , "hatch", "extend", "n5")
removallevels <- c("inc.d3",  "m.inc.d3" , "inc.d9", "m.inc.d9" ,
             "inc.d17", "m.inc.d17" , "hatch", "m.n2")

# mix of characterization and manipulation, ordered by time
alllevels <- c("control", "bldg", "lay", "inc.d3", "m.inc.d3" ,  
               "inc.d9", "m.inc.d9" , "early" ,
               "inc.d17",  "m.inc.d17", "prolong", 
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
                             "bldg_inc.d9", "bldg_inc.d17",  
                             "bldg_hatch", "bldg_n5", "bldg_n9")

comparisonlabelschar <- c("bldg\nlay" ,"lay\ninc.d3", 
                           "inc.d3\ninc.d9", "inc.d9\ninc.d17", 
                           "inc.d17\nhatch" , "hatch\nn5", "n5\nn9")

comparisonlevelscharnobldg <- c("lay_inc.d3",  
                          "inc.d3_inc.d9",   "inc.d9_inc.d17",  
                          "inc.d17_hatch", "hatch_n5", "n5_n9")

comparisonlabelscharnobldg <- c("lay\ninc3", 
                          "inc3\ninc9", "inc.d9\ninc17", 
                          "inc17\nhatch" , "hatch\nn5", "n5\nn9")

comparisonlabelscontrol <- c("bldg\ncontrol" ,"lay\ncontrol", "inc3\ncontrol", 
                                "inc9\ncontrol", "inc17\ncontrol", 
                                "hatch\ncontrol" , "n5\ncontrol", "n9\ncontrol")

comparisonlabelsbldg <- c("lay\nbldg", "inc3\nbldg", 
                                "inc9\nbldg", "inc17\nbldg", 
                                "hatch\nbldg" , "n5\nbldg", "n9\nbldg")

comparisonlevelsremoval <- c("inc.d3_m.inc.d3", "inc.d9_m.inc.d9", "inc.d17_m.inc.d17", "hatch_m.n2")

comparisonlevelsreplace <- c("hatch_early", "inc.d17_prolong",  "hatch_prolong", 
                             "extend_hatch", "n5_extend")

comparisonlablelssreplace <- c("hatch\nearly", "inc.d17\nprolong",  "hatch\nprolong", 
                             "extend\nhatch", "n5\nextend")



comparisonlevelsmanip <- c("inc.d3_m.inc.d3", "inc.d9_m.inc.d9",  "early_inc.d9" ,
                           "hatch_early",  "inc.d17_m.inc.d17", 
                           "hatch_prolong", "hatch_m.n2",
                            "extend_hatch", "n5_extend")

comparisonlabelssmanip <- c("inc.d3vs.\nm.inc.d3", "inc.d9vs.\nm.inc.d9", "earlyvs.\ninc.d9" ,
                            "hatchvs.\nearly", 
                            "inc.d17vs.\nm.inc.d17", "inc.d17vs.\nprolong",
                            "hatchvs.\nprolong", "hatchvs.\nm.n2",
                             "extendvs.\nhatch", "n5vs.\nextend")

comparisonlevelscontrolreplace <- c("control_bldg", "control_lay",  "control_inc.d3",  
                             "control_inc.d9",   "control_inc.d17",  
                             "control_hatch", "control_n5", "control_n9",
                             "control_early", "control_prolong", "control_extend")
comparisonlabelscontrolreaplce  <- c("bldg\ncontrol" ,"lay\ncontrol", "inc3\ncontrol", 
                             "inc9\ncontrol", "inc17\ncontrol", 
                             "hatch\ncontrol" , "n5\ncontrol", "n9\ncontrol",
                             "early\ncontrol" , "prolong\ncontrol", "extend\ncontrol")


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
