import ocrofst
from ocrofstre import *

fst1 = GEN(STR("hello"))
fst2 = GEN(STR("hello"))
assert sum(ocrofst.beam_search(fst1,fst2,10)[4])==0.0

fst1 = GEN(STR("hello",1.0))
fst2 = GEN(STR("hello"))
assert sum(ocrofst.beam_search(fst1,fst2,10)[4])==1.0

fst1 = GEN(STR("hello",1.0))
fst2 = GEN(STR("hello",3.0))
assert sum(ocrofst.beam_search(fst1,fst2,10)[4])==4.0

fst1 = GEN(STR("hello"))
fst2 = GEN(STR("hellx"))
assert sum(ocrofst.beam_search(fst1,fst2,10)[4])>99999.0

fst1 = GEN(STAR(STR("hello",1.0)))
fst2 = GEN(STR("hellohellohello"))
assert sum(ocrofst.beam_search(fst1,fst2,10)[4])==3.0

fst1 = GEN(STAR(STR("hello",1.0)))
fst2 = GEN(SEQ(STR("hello"),STR("hello")))
assert sum(ocrofst.beam_search(fst1,fst2,10)[4])==2.0

fst1 = GEN(STAR(STR("hello",1.0)))
fst2 = GEN(SEQ(STR("hello"),COST("hello",1.0)))
assert sum(ocrofst.beam_search(fst1,fst2,10)[4])==3.0

fst1 = GEN(ALT(STR("hello",1.0),STR("world",2.0)))
fst2 = GEN(SEQ(STR("hello")))
assert sum(ocrofst.beam_search(fst1,fst2,10)[4])==1.0

fst1 = GEN(ALT(STR("hello",1.0),STR("world",2.0)))
fst2 = GEN(SEQ(STR("world")))
assert sum(ocrofst.beam_search(fst1,fst2,10)[4])==2.0

anyfst = ocrofst.OcroFST()
state = anyfst.newState()
anyfst.setStart(state)
anyfst.setAccept(state)
for i in range(32,4096):
    anyfst.addTransition(state,state,i,0.0,i)
anyfst.save("any.fst")
fst2 = GEN(STR("hello"))
assert sum(ocrofst.beam_search(anyfst,fst2,10)[4])==0.0

countfst = ocrofst.OcroFST()
state = countfst.newState()
countfst.setStart(state)
countfst.setAccept(state)
for i in range(32,4096):
    countfst.addTransition(state,state,i,1.0,i)
fst2 = GEN(STR("hello"))
assert sum(ocrofst.beam_search(countfst,fst2,10)[4])==5.0
