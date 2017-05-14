#include "ocrorast.h"

#if 0
// A valid line segmentation may contain 0 or 0xffffff as the
// background, and otherwise numbers components starting at 1.
// The segmentation consists of segmented background pixels
// (0x80xxxx) and segmented foreground pixels (0x00xxxx).  The
// segmented foreground pixels should constitute a usable
// binarization of the original image.

void check_line_segmentation(intarray &cseg) {
    if(cseg.length1d()==0) return;
    CHECK_ARG(cseg.rank()==2);
    for(int i=0;i<cseg.length1d();i++) {
        int value = cseg.at1d(i);
        if(value==0) continue;
        if(value==0xffffff) continue;
        if(value&0x800000)
            CHECK_ARG((value&~0x800000)<100000);
        else
            CHECK_ARG(value<100000);
    }
}

// FIXME/mezhirov add comments --tmb

void make_line_segmentation_black(intarray &a) {
    check_line_segmentation(a);
    replace_values(a, 0xFFFFFF, 0);
    for(int i = 0; i < a.length1d(); i++)
        a.at1d(i) &= 0xFFF;
}

// FIXME/mezhirov add comments --tmb

void make_line_segmentation_white(intarray &a) {
    replace_values(a, 0, 0xFFFFFF);
    //for(int i = 0; i < a.length1d(); i++)
    //    a.at1d(i) = (a.at1d(i) & 0xFFF) | 0x1000;
    check_line_segmentation(a);
}
#endif
