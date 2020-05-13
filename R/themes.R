# custom themes

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

# manipulation
maniplevels <- c( "rm inc.d3" ,  "rm inc.d9" ,  "rm inc.d17" , "rm inc.d20" ,  "extend" , "prolong", "early")
maniplevels1 <- c( "m.inc.d3" ,  "early" ,  "m.inc.d9" , "m.inc.d17" , "m.n2", "prolong" , "extend")

# mix
levelstiming <- c( "early" , "prolong" , "extend")
levelsremoval <- c( "m.inc.d3" ,    "m.inc.d9" , "m.inc.d17" , "m.n2")
controlstiming <- c( "inc.d9" , "hatch" , "inc.d17")
controlsremoval <- c( "inc.d3" ,    "inc.d9" , "inc.d17" , "hatch")
allmaniplevels <- c(maniplevels1, controlsremoval, "n5", "lay", "bldg")

# mix of characterization and manipulation, ordered by time
alllevels <- c("control", "bldg", "lay", "inc.d3", "m.inc.d3" ,  
               "inc.d9", "m.inc.d9" , "early" ,
               "inc.d17",  "m.inc.d17","prolong" ,
               "hatch",  "m.n2",   "extend",
               "n5", "n9")

# characteriation first, then manip by time
alllevels2 <- c("control", "bldg", "lay", "inc.d3", "inc.d9", 
                "inc.d17", "hatch", "n5", "n9",
                "m.inc.d3" ,  "early" ,  "m.inc.d9" , 
                "m.inc.d17" , "m.n2", "prolong" , "extend")

# tissue
tissuelevels <- c("hypothalamus", "pituitary", "gonads")
tissuelevel <- c("hypothalamus", "pituitary", "gonad")
myshapes = c("hypothalamus" = 20,  "pituitary" = 17,  "gonads" = 15)

# sex
sexlevels <- c("female", "male")

levelsnewgrouping <-  c("control" , "bldg",  "earlyinc", 
                        "nearhatch",  "chickcare", "manip")

levelsextint <- c("eggs" , "chicks" , "loss")

# DESeq
comparisonlevelschar <- c("control_bldg",  "bldg_lay",  "lay_inc.d3",  
                          "inc.d3_inc.d9",   "inc.d9_inc.d17",  
                          "inc.d17_hatch", "hatch_n5", "n5_n9")

comparisonlabelschar = c("control vs bldg", "bldg vs lay" ,"lay vs inc.d3", 
                         "inc.d3 vs inc.d9", "inc.d9 vs inc.d17", 
                         "inc.d17 vs hatch" , "hatch vs n5", "n5 vs n9")

# custom colors

colorschar <-  c("control" = "#cc4c02", 
                 "bldg"= "#fe9929", 
                 "lay"= "#fed98e", 
                 "inc.d3"= "#78c679", 
                 "inc.d9"= "#31a354", 
                 "inc.d17"= "#006837", 
                 "hatch"= "#08519c",
                 "n5"= "#3182bd", 
                 "n9"= "#6baed6")

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

colorhypothesis <- c("eggs" = "#7FA163",
                     "chicks" = "#3A80B9",
                     "lo" = "#B5C4B2",
                     "hi" = "#19757A",
                     "loss" = "#8B4513")

colorstissue <- c("hypothalamus" = "#d95f02",
                  "pituitary" = "#1b9e77",
                  "gonads" =  "#7570b3")

sexcolors <- c("female" = "#969696", "male" = "#525252")

colorsvolcano <-  c(colorschar,
                    
                    "early" = "#cbc9e2", 
                    "prolong" = "#9e9ac8" , 
                     "extend" = "#6a51a3",
                    
                    "male" = "#525252",
                    "female" = "#969696",
                    
                    "NS" = "#bdbdbd")

 
allcolors <- c(colorscharmaip, sexcolors, colorstissue, colorhypothesis,  "gonad" =  "#7570b3")
 
colorsvolcanochar <-  c(colorschar,
                    "NS" = "#bdbdbd")

myannotationcolors <- list(sex = sexcolors,
                        tissue = colorstissue,
                        treatment = colorscharmaip)
