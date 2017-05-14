// -*- C++ -*-

%{
#pragma GCC diagnostic ignored "-Wstrict-aliasing"
#pragma GCC diagnostic ignored "-Wuninitialized"
%}

%module(docstring="Layout analysis for OCRopus") ocrorast;
%feature("autodoc",1);
%include "typemaps.i"
#ifdef SWIGPYTHON
%include "cstring.i"
#endif

%{
#include "ocrorast.h"
#include "wrappers.h"
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
        PyErr_SetString(PyExc_IndexError,"unknown exception in iulib");
        return NULL;
    }
}
#endif

struct TextLineRAST {
    virtual void setParams(double max_slope,double ymin,double ymax) = 0;
    virtual void setMaxLines(int i) = 0;
    virtual void addChar(int x0,int y0,int x1,int y1) = 0;
    virtual void addColumn(int x0,int y0,int x1,int y1) = 0;
    virtual void compute() = 0;
    virtual int nlines() = 0;
    virtual float getQuality(int i);
    virtual void getLine(int i,float &c,float &m,float &d) = 0;
    virtual float getLine_c(int i);
    virtual float getLine_m(int i);
    virtual float getLine_d(int i);
};
TextLineRAST *make_TextLineRAST();

struct TextLineRAST2 {
    virtual void setParams(double max_slope,double ymin,double ymax) = 0;
    virtual void setMaxLines(int i) = 0;
    virtual void addChar(int x0,int y0,int x1,int y1) = 0;
    virtual void addColumn(int x0,int y0,int x1,int y1) = 0;
    virtual void compute() = 0;
    virtual int nlines() = 0;
    virtual float getQuality(int i);
    virtual void getLine(int i,float &c,float &m,float &d) = 0;
    virtual float getLine_c(int i);
    virtual float getLine_m(int i);
    virtual float getLine_d(int i);
    virtual void getBbox(int i,int x0,int y0,int x1,int y1) = 0;
};
TextLineRAST2 *make_TextLineRAST2();

%import iulib.i

void align_segmentation(intarray &segmentation,rectarray &bboxes);
void make_line_segmentation_black(intarray &);
void make_line_segmentation_white(intarray &);
void make_page_segmentation_black(intarray &);
void make_page_segmentation_white(intarray &);
// void renumber_labels(intarray &,int);
void skeletal_features(bytearray &endpoints,
                       bytearray &junctions,
                       bytearray &image,
                       float presmooth=0.0,
                       float skelsmooth=0.0);
void skeletal_feature_counts(int &nendpoints,
                             int &njunctions,
                             bytearray &image,
                             float presmooth=0.0,
                             float skelsmooth=0.0);
int endpoints_counts(bytearray &image,float presmooth=0.0,float skelsmooth=0.0);
int junction_counts(bytearray &image,float presmooth=0.0,float skelsmooth=0.0);
int component_counts(bytearray &image,float presmooth=0.0);
int hole_counts(bytearray &image,float presmooth=0.0);


%inline %{
    int skeletal_feature_hack(bytearray &image,
                              float presmooth=0.0,
                              float skelsmooth=0.0) {
        int ne,nj;
        skeletal_feature_counts(ne,nj,image,presmooth,skelsmooth);
        return 1000*ne+nj;
    }
%}

/*
  The high-level interface to the layout analysis algorithms.
 */

%feature("director") IComponent;
struct IComponent {
    virtual void reinit();
    virtual const char *name() = 0;
    virtual const char *description() = 0;
    virtual const char *interface() = 0;
    virtual void pset(const char *name,const char *value);
    virtual void pset(const char *name,double value);
    void pset(const char *name,int value);
    void pset(const char *name,bool value);
    const char *pget(const char *name);
    double pgetf(const char *name);
    void pprint(FILE *stream=stdout,int depth=0);
};

%feature("director") ISegmentPage;
struct ISegmentPage : IComponent {
    virtual void segment(intarray &out,bytearray &input) = 0;
    virtual void segment(intarray &out,bytearray &input,rectarray &obstacles);
    virtual void ~ISegmentPage() = 0;
};

ISegmentPage *make_SegmentPageBy1CP();
ISegmentPage *make_SegmentPageByRAST();
ISegmentPage *make_SegmentPageByRAST1();
ISegmentPage *make_SegmentPageByMorphTrivial();
ISegmentPage *make_SegmentPageByWCUTS();
ISegmentPage *make_SegmentPageByXYCUTS();
ISegmentPage *make_SegmentWords();

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

def page2narray(page,type='B'):
    """Convert page images to narrays."""
    checknp(page)
    return numpy2narray(page,type=type)

def narray2page(page,type='B'):
    """Convert narrays to page images."""
    checkna(page)
    return narray2numpy(page,type=type)

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

def narray2pseg(na):
    """Convert an narray to a page segmentation (rank 3, RGB)."""
    checkna(na)
    pseg = iulib.numpy(na,type='i')
    pseg = numpy.array([pseg>>16,pseg>>8,pseg],'B')
    pseg = numpy.transpose(pseg,[2,1,0])
    pseg = pseg[::-1,...]
    return pseg

def pseg2narray(pseg):
    """Convert a page segmentation (rank 3, RGB) to an narray."""
    checknp(pseg)
    assert pseg.dtype=='B' and pseg.ndim==3
    r = numpy2narray(ascontiguousarray(pseg[:,:,0]))
    g = numpy2narray(ascontiguousarray(pseg[:,:,1]))
    b = numpy2narray(ascontiguousarray(pseg[:,:,2]))
    rgb = iulib.intarray()
    iulib.pack_rgb(rgb,r,g,b)
    return rgb

class SegmentPage:
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
    def segment(self,page,obstacles=None,black=0):
        page = page2narray(page,'B')
        # iulib.write_image_gray("_seg_in.png",page)
        result = iulib.intarray()
        if obstacles not in [None,[]]:
            raise Unimplemented()
        else:
            self.comp.segment(result,page)
        if black: ocropus.make_page_segmentation_black(result)
        # iulib.write_image_packed("_seg_out.png",result)
        return narray2pseg(result)

class SegmentPageByRAST(SegmentPage):
    """Segment a page into columns and lines using the RAST algorithm."""
    c_class = [make_SegmentPageByRAST]

class SegmentPageByRAST1(SegmentPage):
    """Segment a page into columns and lines using the RAST algorithm,
    assuming there is only a single column.  This is more robust for
    single column documents than RAST."""
    c_class = [make_SegmentPageByRAST1]

class SegmentPageBy1CP(SegmentPage):
    """A very simple page segmentation algorithm assuming a single column
    document and performing projection."""
    c_class = [make_SegmentPageBy1CP]

class SegmentPageByXYCUTS(SegmentPage):
    """An implementation of the XYCUT layout analysis algorithm.  Not
    recommended for production use."""
    c_class = [make_SegmentPageByXYCUTS]

%}
