// -*- C++ -*-

%{
#pragma GCC diagnostic ignored "-Wstrict-aliasing"
#pragma GCC diagnostic ignored "-Wuninitialized"
%}

%module ocrofstll
%include "typemaps.i"
%include "cstring.i"
%{
#include <malloc.h>
#include <string.h>
#include <colib/checks.h>
#include <colib/narray.h>
#include <colib/narray-ops.h>
#include <colib/narray-util.h>
#include <iulib/iulib.h>
#include <iulib/components.h>
/* #include <ocrofst/ocr-pfst.h> */
#include "ocr-pfst.h"
using namespace colib;
using namespace iulib;
using namespace narray_ops;
using namespace ocrofst;
%}
%import iulib.i

%feature("director") IGenericFst;
struct IGenericFst : virtual IComponent {
    virtual void clear() = 0;
    virtual int newState();
    virtual void addTransition(int from,int to,int output,float cost,int input);
    virtual void addTransition(int from,int to,int symbol,float cost);
    virtual void setStart(int node);
    virtual void setAccept(int node,float cost=0.0);
    virtual int special(const char *s);
    virtual void bestpath(ustrg &result);
    virtual ~IGenericFst();
    virtual void setString(ustrg &text,floatarray &costs,intarray &ids);
    virtual int nStates();
    virtual int getStart();
    virtual float getAcceptCost(int node);
    virtual bool isAccepting(int node);
    virtual void getTransitions(intarray &tos,intarray &symbols,floatarray &costs,intarray &inputs,int from);
    virtual void rescore(int from,int to,int output,float new_cost,int input);
    virtual void rescore(int from, int to, int symbol, float new_cost);
    virtual void load(const char *file);
    virtual void save(const char *file);
};


enum {
    L_SIGMA = -4,
    L_RHO = -3,
    L_PHI = -2,
    L_EPSILON = 0,
};

struct OcroFST : IGenericFst {
};

%newobject make_OcroFST;
OcroFST *make_OcroFST();

void scale_fst(OcroFST &,float);


%inline %{

    double beam_search(ustrg &s,IGenericFst &u,IGenericFst &v,int n) {
        OcroFST *uo = dynamic_cast<OcroFST*>(&u);
        CHECK(uo!=0);
        OcroFST *vo = dynamic_cast<OcroFST*>(&v);
        CHECK(vo!=0);
        return beam_search(s,*uo,*vo,n);
    }
    void beam_search(intarray &v1,intarray &v2,intarray &ins,intarray &outs,floatarray &costs,IGenericFst &u,IGenericFst &v,int n) {
        OcroFST *uo = dynamic_cast<OcroFST*>(&u);
        CHECK(uo!=0);
        OcroFST *vo = dynamic_cast<OcroFST*>(&v);
        CHECK(vo!=0);
        beam_search(v1,v2,ins,outs,costs,*uo,*vo,n);
    }    
    
        
%}

