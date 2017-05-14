from pylab import *
from scipy.misc import imsave
import ocrolseg
image = imread("line.png")
image = array(image,'B')
print image.shape
seg = ocrolseg.SegmentLineByCCS()
out = seg.charseg(image)
print out.shape,amin(out),amax(out)
imsave("out.png",out)
