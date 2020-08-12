bldg.lay <- createDEGdftreatment("lay", "bldg", i)
lay.inc.d3 <- createDEGdftreatment("inc.d3", "lay",  i) 
inc.d3.inc.d9 <- createDEGdftreatment("inc.d9", "inc.d3", i) 
inc.d9.inc.d17 <- createDEGdftreatment("inc.d17", "inc.d9", i)
inc.d17.hatch <- createDEGdftreatment( "hatch", "inc.d17", i) 
hatch.n5 <- createDEGdftreatment( "n5", "hatch",  i) 
n5.n9 <- createDEGdftreatment("n9", "n5",  i) 

control.lay <- createDEGdftreatment("lay", "control", i)
control.inc.d3 <- createDEGdftreatment("inc.d3", "control",  i) 
control.inc.d9 <- createDEGdftreatment("inc.d9", "control", i) 
control.inc.d17 <- createDEGdftreatment("inc.d17", "control", i)
control.hatch <- createDEGdftreatment( "hatch", "control", i) 
control.n5 <- createDEGdftreatment( "n5", "control",  i) 
control.n9 <- createDEGdftreatment("n9", "control",  i) 

control.extend <- createDEGdftreatment("extend", "control",  i) 
control.prolong <- createDEGdftreatment("prolong", "control",  i) 

control.m.inc.d3 <- createDEGdftreatment("m.inc.d3", "control",  i) 
control.m.inc.d9 <- createDEGdftreatment("m.inc.d9", "control",  i) 
control.m.inc.d17 <- createDEGdftreatment("m.inc.d17", "control",  i)
control.m.n2 <- createDEGdftreatment("m.n2", "control",  i) 

bldg.inc.d3 <- createDEGdftreatment("inc.d3", "bldg",  i) 
bldg.inc.d9 <- createDEGdftreatment("inc.d9", "bldg", i) 
bldg.inc.d17 <- createDEGdftreatment("inc.d17", "bldg", i)
bldg.hatch <- createDEGdftreatment( "hatch", "bldg", i) 
bldg.n5 <- createDEGdftreatment( "n5", "bldg",  i) 
bldg.n9 <- createDEGdftreatment("n9", "bldg",  i) 

bldg.extend <- createDEGdftreatment("extend", "bldg",  i) 
bldg.prolong <- createDEGdftreatment("prolong", "bldg",  i) 


bldg.m.inc.d3 <- createDEGdftreatment("m.inc.d3", "bldg",  i) 
bldg.m.inc.d9 <- createDEGdftreatment("m.inc.d9", "bldg",  i) 
bldg.m.inc.d17 <- createDEGdftreatment("m.inc.d17", "bldg",  i)
bldg.m.n2 <- createDEGdftreatment("m.n2", "bldg",  i) 

hatch.m.n2 <- createDEGdftreatment("m.n2", "hatch",  i) 
hatch.prolong <- createDEGdftreatment("prolong", "hatch",  i) 

inc.d17.m.n2 <- createDEGdftreatment("m.n2", "inc.d17",  i) 
inc.d17.prolong <- createDEGdftreatment("prolong", "inc.d17",  i) 

prolong.extend <- createDEGdftreatment("extend", "prolong",  i) 

m.inc.d3.m.inc.d9 <- createDEGdftreatment("m.inc.d9", "m.inc.d3",  i) 
m.inc.d9.m.inc.d17 <- createDEGdftreatment("m.inc.d17", "m.inc.d9",  i) 
m.inc.d3.m.inc.d17 <- createDEGdftreatment("m.inc.d17", "m.inc.d3",  i) 
m.inc.d17.m.n2 <- createDEGdftreatment("m.n2", "m.inc.d17",  i) 
m.inc.d9.m.n2 <- createDEGdftreatment("m.n2", "m.inc.d9",  i)
m.inc.d3.m.n2 <- createDEGdftreatment("m.n2", "m.inc.d3",  i)

inc.d3.m.inc.d3 <- createDEGdftreatment("m.inc.d3", "inc.d3",  i) 
inc.d9.m.inc.d9 <- createDEGdftreatment("m.inc.d9", "inc.d9",  i) 
inc.d17.m.inc.d17 <- createDEGdftreatment("m.inc.d17", "inc.d17",  i) 
hatch.m.n2 <- createDEGdftreatment("m.n2", "hatch",  i) 

early.inc.d9 <- createDEGdftreatment("inc.d9", "early",  i) 
hatch.early <- createDEGdftreatment("early", "hatch",  i) 

prolong.inc.d17 <- createDEGdftreatment("inc.d17", "prolong",  i) 
hatch.prolong <- createDEGdftreatment("prolong", "hatch",  i) 

extend.hatch <- createDEGdftreatment("hatch", "extend",  i) 
n5.extend <- createDEGdftreatment("extend", "n5",  i) 

