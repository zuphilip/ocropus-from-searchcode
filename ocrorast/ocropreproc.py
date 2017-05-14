################################################################
# preprocessing
################################################################

import os,os.path
from pylab import *
import ocrorast
from scipy.ndimage import measurements,interpolation,filters,morphology
import math
ion(); gray()

def dread(fname):
    """Read an image, similar to imread.  However, imread flips JPEG images;
    this fixes that."""
    _,ext = os.path.splitext(fname)
    image = imread(fname)
    if ext.lower() in [".jpg",".jpeg"]:
        image = image[::-1,:,...]
    return image

################################################################
### Binarization
################################################################

def is_binary(image):
    """Check whether an input image is binary"""
    return sum(image==amin(image))+sum(image==amax(image)) > 0.99*image.size
    
def gsauvola(image,sigma=150.0,R=None,k=0.3,filter='uniform',scale=2.0):
    """Perform Sauvola-like binarization.  This uses linear filters to
    compute the local mean and variance at every pixel."""
    if image.dtype==dtype('uint8'): image = image / 256.0
    if len(image.shape)==3: image = mean(image,axis=2)
    if filter=="gaussian":
        filter = filters.gaussian_filter
    elif filter=="uniform":
        filter = filters.uniform_filter
    else:
        pass
    scaled = interpolation.zoom(image,1.0/scale,order=0,mode='nearest')
    s1 = filter(ones(scaled.shape),sigma)
    sx = filter(scaled,sigma)
    sxx = filter(scaled**2,sigma)
    avg_ = sx / s1
    stddev_ = maximum(sxx/s1 - avg_**2,0.0)**0.5
    s0,s1 = avg_.shape
    s0 = int(s0*scale)
    s1 = int(s1*scale)
    avg = zeros(image.shape)
    interpolation.zoom(avg_,scale,output=avg[:s0,:s1],order=0,mode='nearest')
    stddev = zeros(image.shape)
    interpolation.zoom(stddev_,scale,output=stddev[:s0,:s1],order=0,mode='nearest')
    if R is None: R = amax(stddev)
    thresh = avg * (1.0 + k * (stddev / R - 1.0))
    return array(255*(image>thresh),'uint8')

def inverse(image):
    return amax(image)-image

def autoinvert(image):
    """Automatically invert document images, so that the majority of pixels
    (background pixels) are black."""
    if median(image)>mean([amax(image),amin(image)]):
        image = amax(image)-image
    return image

################################################################
### Bounding-box operations.
################################################################

def bounding_boxes_math(image):
    """Compute the bounding boxes in the image; returns mathematical
    coordinates."""
    image = (image>mean([amax(image),amin(image)]))
    image,ncomponents = measurements.label(image)
    objects = measurements.find_objects(image)
    result = []
    h,w = image.shape
    for o in objects:
        y1 = h-o[0].start
        y0 = h-o[0].stop
        x0 = o[1].start
        x1 = o[1].stop
        c = (x0,y0,x1,y1)
        result.append(c)
    return result

def select_plausible_char_bboxes(bboxes,dpi=300.0):
    """Performs simple heuristic checks on character bounding boxes;
    removes boxes that are too small or too large, or have the wrong
    aspect ratio."""
    s = dpi/300.0
    result = []
    for b in bboxes:
        x0,y0,x1,y1 = b
        w = x1-x0
        if w<s*5: continue
        h = y1-y0
        if h<s*5: continue
        a = w*1.0/h
        if a>s or a<0.25: continue
        if w>s*100: continue
        if h>s*100: continue
        result.append(b)
    return result

def estimate_skew_angle(image):
    """Estimate the skew angle of a document image, by first finding
    character bounding boxes, then invoking the RAST text line finder
    (without constraints) in order to find the longest line."""
    assert is_binary(image)
    finder = ocrorast.make_TextLineRAST2()
    finder.setMaxLines(1)
    nobjects = 0
    bboxes = bounding_boxes_math(image)
    bboxes = select_plausible_char_bboxes(bboxes)
    assert bboxes!=[]
    for c in bboxes:
        finder.addChar(*c)
        nobjects += 1
    finder.compute()
    m = finder.getLine_m(0)
    del finder
    return math.atan(m)

def deskew(image):
    """Actually deskew an image by first estimating the skew angle, then
    performing the rotation."""
    a = estimate_skew_angle(image)
    return interpolation.rotate(image,-a*180/pi,mode='nearest',order=0)

def check_contains_halftones(image,dpi=300.0):
    """Heuristic method for determining whether we should apply a halftone removal
    algorithm."""
    bboxes = bounding_boxes_math(image)
    r = 4*dpi/300.0
    big = 0
    for b in bboxes:
        x0,y0,x1,y1 = b
        if x1-x0>r or y1-y0>r: big += 1
    return big<0.3*len(bboxes)

def remove_small_components(image,r=3):
    """Remove any connected components that are smaller in both dimension than r"""
    image,ncomponents = measurements.label(image)
    objects = measurements.find_objects(image)
    for i in range(len(objects)):
        o = objects[i]
        if o[0].stop-o[0].start>r: continue
        if o[1].stop-o[1].start>r: continue
        c = image[o]
        c[c==i+1] = 0
    return (image!=0)

def remove_big_components(image,r=100):
    """Remove any connected components that are smaller in any dimension than r"""
    image,ncomponents = measurements.label(image)
    objects = measurements.find_objects(image)
    for i in range(len(objects)):
        o = objects[i]
        if o[0].stop-o[0].start<r and o[1].stop-o[1].start<r: continue
        c = image[o]
        c[c==i+1] = 0
    return (image!=0)

def remove_small_any(image,r=3):
    """Remove both small connected components and small holes."""
    image = remove_small_components(image,r=r)
    image = amax(image)-image
    image = remove_small_components(image,r=r)
    image = amax(image)-image
    return image

def rectangular_cover(image,minsize=5):
    """Cover the set of regions with their bounding boxes.  This is
    an image-to-image transformation."""
    image,ncomponents = measurements.label(image)
    objects = measurements.find_objects(image)
    output = zeros(image.shape)
    for i in range(len(objects)):
        o = objects[i]
        if o[0].stop-o[0].start<minsize: continue
        if o[1].stop-o[1].start<minsize: continue
        output[o] = 1
    return output

def find_halftones(image,dpi=300.0,threshold=0.05,r=5,sigma=15.0,cover=1):
    """Find halftone regions in an image.  First, find small components and
    holes, then smooth their occurrences and threshold, finally compute
    a rectangular cover of the thresholded and smoothed image."""
    filtered = remove_small_any(image,r=r)
    diff = ((image!=0)!=(filtered!=0))
    density = filters.gaussian_filter(1.0*diff,sigma*dpi/300.0)
    if cover:
        return rectangular_cover(density>threshold)
    else:
        return maximum(diff,density>threshold)

def remove_halftones(image,dpi=300.0,threshold=0.05,r=5,sigma=15.0):
    """Perform halftone removal using find_halftones."""
    halftones = find_halftones(image,dpi=dpi,threshold=threshold,r=r,sigma=sigma)
    return maximum(image-amax(image)*halftones,0)

################################################################
### All preprocessing steps put together.
################################################################

def preprocess(raw,autoinv=1,clean=1,minsize=2,maxsize=100):
    assert amax(raw)>amin(raw),"empty image supplied"
    # binarize if the image isn't already binary
    if not is_binary(raw):
        bin = gsauvola(raw)
    else:
        bin = array(255*(raw!=amin(raw)),'B')
    assert amax(bin)>amin(bin),"something went wrong with binarization"
    if autoinv:
        bin = autoinvert(bin)
    # now clean up for skew estimation
    cleaned = bin
    if check_contains_halftones(bin):
        cleaned = remove_halftones(bin)
    cleaned = remove_small_components(cleaned,minsize)
    cleaned = remove_big_components(cleaned,maxsize)
    # perform skew estimation
    a = estimate_skew_angle(cleaned)
    # if the clean flag is given, output the cleaned image
    if clean:
        bin = cleaned
    # perform the deskewing
    deskewed = interpolation.rotate(bin,-a*180/pi,mode='nearest',order=0)
    deskewed = array(255*(deskewed==0),'B')
    return deskewed

################################################################
### Here is a class that can be dropped into the regular
### OCRopus processing pipeline.
################################################################

class CommonPreprocessing:
    def binarize(self,page):
        bin = preprocess(page)
        assert bin.dtype==dtype('B'),"wrong type: %s"%(bin.dtype,)
        return bin,bin
    def binarize_color(self,page):
        return self.binarize(mean(page,axis=2))
