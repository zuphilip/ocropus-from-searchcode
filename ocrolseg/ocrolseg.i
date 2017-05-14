// -*- C++ -*-

%{
#pragma GCC diagnostic ignored "-Wstrict-aliasing"
#pragma GCC diagnostic ignored "-Wuninitialized"
%}

%module(docstring="Layout analysis for OCRopus") ocrolseg;
%feature("autodoc",1);
%include "typemaps.i"
#ifdef SWIGPYTHON
%include "cstring.i"
#endif

%{
#include "ocrolseg.h"
%}

#ifdef SWIGPYTHON
%exception {
    try {
        $action
    }
    catch(const char *s) {
        PyErr_SetString(PyExc_IndexError,s);
        return NULL;
    }
    catch(...) {
        PyErr_SetString(PyExc_IndexError,"unknown exception in ocrolseg");
        return NULL;
    }
}
#endif

%import iulib.i

%feature("director") IComponent;
struct IComponent {
    IComponent();
    virtual void reinit();
    virtual const char *name() = 0; // FIXME -- making this non-abstract makes things fail
    virtual const char *description();
    virtual const char *interface();
    virtual void pset(const char *name,const char *value);
    virtual void pset(const char *name,double value);
    void pset(const char *name,int value);
    void pset(const char *name,bool value);
    const char *pget(const char *name);
    double pgetf(const char *name);
    void pprint(FILE *stream=stdout,int depth=0);
    virtual ~IComponent();
};

%feature("director") ISegmentLine;
struct ISegmentLine : IComponent {
    virtual void charseg(intarray &out,bytearray &in);
};

ISegmentLine *make_SegmentLineByCCS();
ISegmentLine *make_SegmentLineByGCCS();
ISegmentLine *make_ConnectedComponentSegmenter();
ISegmentLine *make_CurvedCutSegmenter1();
ISegmentLine *make_CurvedCutSegmenter();
ISegmentLine *make_CurvedCutWithCcSegmenter();


struct IDpSegmenter : ISegmentLine {
    float down_cost;
    float outside_diagonal_cost;
    float outside_diagonal_cost_r;
    float inside_diagonal_cost;
    float boundary_diagonal_cost;
    float inside_weight;
    float boundary_weight;
    float outside_weight;
    int min_range;
    float cost_smooth;
    float min_thresh;
    intarray dimage;
};

IDpSegmenter *make_DpSegmenter();

/* void make_line_segmentation_black(intarray &); */

%pythoncode %{
import iulib
import numpy

def checknp(a):
    """Checks whether the argument is a numpy array.  Raises an error if not."""
    if type(a) in [iulib.bytearray,iulib.intarray,iulib.floatarray,iulib.rectarray]:
        raise Exception("numpy array expected; an narray was passed")
    assert type(a)==numpy.ndarray
def checkna(a):
    """Checks whether the argument is an narray.  Raises an error if not."""
    if type(a) in [iulib.bytearray,iulib.intarray,iulib.floatarray,iulib.rectarray]:
        return
    if type(a)==numpy.array:
        raise Exception("narray expected; a numpy array was passed")
    raise Exception("expected an narray, got something different")

def ctype(a):
    """Return the numpy type character for an array."""
    if type(a)==str: return a
    if type(a)==iulib.floatarray: return 'f'
    if type(a)==iulib.intarray: return 'i'
    if type(a)==iulib.bytearray: return 'B'
    return a.dtype

def numpy2narray(page,type='B'):
    """Convert a numpy image to an narray. Flips from raster to
    mathematical coordinates.  When converting float to integer
    types, multiplies with 255.0, and when converting integer to
    float types, divides by 255.0."""
    checknp(page)
    if type is None: type = ctype(page)
    if isfp(page) and not isfp(type):
        page = numpy.array(255*page,dtype='B')
    elif not isfp(page) and isfp(type):
        page = page/255.0
    page = page.transpose([1,0]+range(2,page.ndim))[:,::-1,...]
    return iulib.narray(page,type=type)

def narray2numpy(na,type='B'):
    """Convert an narray image to a numpy image. Flips from mathematical
    coordinates to raster coordinates.  When converting integer to float
    types, multiplies with 255.0, and when converting integer to float
    types divides by 255.0"""
    checkna(na)
    if type is None: type = ctype(na)
    if isfp(na) and not isfp(type):
        page = iulib.numpy(na,'f')
        page = numpy.array(255.0*page,dtype=type)
    elif not isfp(na) and isfp(type):
        page = iulib.numpy(na,type=type)
        page /= 255.0
    else:
        page = iulib.numpy(na,type=type)
    return page.transpose([1,0]+range(2,page.ndim))[::-1,...]

def isfp(a):
    """Check whether the array is a floating point array."""
    if type(a)==str:
        if a in ['f','d']: return 1
        else: return 0
    if type(a)==iulib.floatarray: return 1
    try:
        if a.dtype in [dtype('f'),dtype('d')]: return 1
    except:
        pass
    return 0

def line2narray(line,type='B'):
    """Convert line images to narrays."""
    checknp(line)
    return numpy2narray(line,type=type)

def narray2line(line,type='B'):
    """Convert line narrays to line images."""
    checkna(line)
    return narray2numpy(line,type=type)

def narray2lseg(na):
    """Convert an narray to a line segmentation."""
    checkna(na)
    pseg = iulib.numpy(na,type='i')
    pseg = numpy.transpose(pseg,[1,0])
    pseg = pseg[::-1,...]
    return pseg

def lseg2narray(lseg):
    """Convert a line segmentation (rank 2, 'i') to an narray."""
    checknp(lseg)
    assert lseg.dtype=='i' and lseg.ndim==2,"wanted rank 2 'i' array, got %s"%lseg
    lseg = lseg[::-1,...].transpose()
    lseg = iulib.narray(lseg,type='i')
    return lseg

class SegmentLine:
    """Segment a page into columns and lines (layout analysis)."""
    def __init__(self):
        self.comp = self.c_class[0]()
    def info(self,depth=0):
        """Print information about this object."""
        return self.comp.info()
    def pexists(self,name):
        """Check whether parameter NAME exists."""
        return self.comp.pexists(name)
    def pset(self,name,value):
        """Set parameter NAME to VALUE."""
        return self.comp.pset(name,value)
    def pget(self,name):
        """Get the value of string parameter NAME."""
        return self.comp.pget(name)
    def pgetf(self,name):
        """Get the value of floating point parameter NAME."""
        return self.comp.pgetf(name)
    def charseg(self,line):
        """Segment a text line into potential character parts."""
        result = iulib.intarray()
        self.comp.charseg(result,line2narray(line,'B'))
        iulib.make_line_segmentation_black(result)
        iulib.renumber_labels(result,1)
        return narray2lseg(result)

class SegmentLineByCCS(SegmentLine):
    """Segment a page into columns and lines using the RAST algorithm."""
    c_class = [make_SegmentLineByCCS]

class SegmentLineByGCCS(SegmentLine):
    """Segment a page into columns and lines using the RAST algorithm."""
    c_class = [make_SegmentLineByGCCS]

class ConnectedComponentSegmenter(SegmentLine):
    """Segment a page into columns and lines using the RAST algorithm."""
    c_class = [make_ConnectedComponentSegmenter]

class CurvedCutSegmenter1(SegmentLine):
    """Segment a page into columns and lines using the RAST algorithm."""
    c_class = [make_CurvedCutSegmenter1]

class CurvedCutSegmenter(SegmentLine):
    """Segment a page into columns and lines using the RAST algorithm."""
    c_class = [make_CurvedCutSegmenter]

class CurvedCutWithCcSegmenter(SegmentLine):
    """Segment a page into columns and lines using the RAST algorithm."""
    c_class = [make_CurvedCutWithCcSegmenter]

class DpSegmenter(SegmentLine):
    """Segment a page into columns and lines using the RAST algorithm."""
    c_class = [make_DpSegmenter]

%}
