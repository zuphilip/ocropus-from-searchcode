import iulib, ocrofstll
from ocrofstll import L_SIGMA,L_RHO,L_PHI,L_EPSILON

fst1=ocrofstll.make_OcroFST()
fst2=ocrofstll.make_OcroFST()

s0 = fst1.newState()
s1 = fst1.newState()
s2 = fst1.newState()
s3 = fst1.newState()
s4 = fst1.newState()

fst1.setAccept(s3)
#addTransition(int from,int to,int output,float cost,int input) 
fst1.addTransition(s0,s1,3,3.0,3)
fst1.addTransition(s1, s2, 1, 1.0, 1)
#fst1.addTransition(s2, s3, 1002, 10.0,1002)
fst1.addTransition(s2, s3, 14, 20.0, 14)

a0 = fst2.newState()
a1 = fst2.newState()
a2= fst2.newState()
a3 = fst2.newState()
a4 = fst2.newState()
a5 = fst2.newState()
a6 = fst2.newState()
#O=15, t=20  c=3 R=18
fst2.setAccept(a4);
fst2.setAccept(a5);
fst2.setAccept(a6);
fst2.addTransition(a0, a1, 3, 23.0, 3);#c
fst2.addTransition(a1, a3, 15, 1.0, 15);#O
fst2.addTransition(a1, a2, 1, 20.0, 1);#a
fst2.addTransition(a2, a4, 20, 40.0, 20);#T
fst2.addTransition(a2, a5, 18, 18.0, 18);#R
fst2.addTransition(a3, a6, 23, 13.0, 23);#W

s = iulib.ustrg()
v1 = iulib.intarray()
v2 = iulib.intarray()
ins = iulib.intarray()
outs = iulib.intarray()
costs = iulib.floatarray()
n=1000
ocrofstll.beam_search(v1,v2,ins,outs,costs,fst1,fst2,n)

