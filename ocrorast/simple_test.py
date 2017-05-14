from pylab import *
from scipy.misc import imsave
import ocropreproc
import ocrorast

image = ocropreproc.dread("raw.png")
bin = ocropreproc.gsauvola(image)
imsave("_binarized.png",bin)
bin = ocropreproc.autoinvert(bin)
imsave("_inverted.png",bin)
bin = ocropreproc.deskew(bin)
imsave("_deskewed.png",bin)
rast = ocrorast.SegmentPageByRAST()
seg = rast.segment(bin)
imsave("_segmented.png",seg)
