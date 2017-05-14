import ocrofst

print "### bug1"
fst1 = ocrofst.OcroFST()
fst1.load("bug1_line.fst")
print fst1.bestpath()
fst2 = ocrofst.OcroFST()
fst2.load("bug1_dict.fst")
v1,v2,ins,outs,costs = ocrofst.beam_search(fst1,fst2,1000)
print sum(costs)
