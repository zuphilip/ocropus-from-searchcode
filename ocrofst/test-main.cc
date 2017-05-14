#include <stdio.h>
#if 0
#include "openfst-pfst.h"
#endif
#include "ocr-pfst.h"

using namespace colib;
using namespace iulib;
using namespace ocrofst;

int main() {
    autodel<OcroFST> fst1(make_OcroFST());
    autodel<OcroFST> fst2(make_OcroFST());

    {
        // setup fst1
        int s0 = fst1->newState();
        int s1 = fst1->newState();
        int s2 = fst1->newState();
        int s3 = fst1->newState();
    
        fst1->setAccept(s3);
        fst1->addTransition(s0, s1, 1001, 1.0, 1001);
        fst1->addTransition(s1, s2, L_SIGMA, 2.0, 1002);
        fst1->addTransition(s2, s3, 1003, 3.0, 1003);
    }

    {
        // setup fst2
        int s0 = fst2->newState();
        int s1 = fst2->newState();
        int s2 = fst2->newState();
        int s3 = fst2->newState();
    
        fst2->setAccept(s3);
        fst2->addTransition(s0, s1, 1001, 10.0, 1001);
        //fst2->addTransition(s1, s2, 1005, 200.0, L_RHO);
        fst2->addTransition(s1, s2, L_RHO, 40.0, L_RHO);
        fst2->addTransition(s2, s3, 1003, 30.0, 1003);
    }

    intarray vertices1;
    intarray vertices2;
    intarray inputs;
    intarray outputs;
    floatarray costs;
    beam_search(vertices1, vertices2, inputs, outputs, costs, *fst1, *fst2);

    for (int i = 0; i < inputs.length(); i++){
        printf("%d : %d", inputs[i],outputs[i]);
        printf("\n");
    }
    
#if 0
    StdVectorFst *f1 = StdVectorFst::Read("010009.fst");
    StdVectorFst *f2 = StdVectorFst::Read("line9_compile.fst");

    vector<int> v1;
    vector<int> v2;
    vector<int> in;
    vector<int> out;
    vector<float> costs1;
    op_beamsearch(v1, v2, in, out, costs1, f1, f2,1000);
#endif
    
    return 0;
    
}
