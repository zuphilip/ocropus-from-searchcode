#include "ocrorast.h"

#define RAD_TO_DEG 57.3

param_string debug_deskew("debug_deskew", 0,
        "output deskewed document image as png");

double estimate_skew_by_rast(colib::bytearray &in){
    autodel<DeskewPageByRAST> deskewer(new DeskewPageByRAST());
    return deskewer->getSkewAngle(in);
}

double DeskewPageByRAST::getSkewAngle(bytearray &in) {
    bytearray binarized;
    binarize_simple(binarized,in);

    // Do connected component analysis
    intarray charimage;
    copy(charimage, binarized);
    make_page_binary_and_black(charimage);
    label_components(charimage, false);
    rectarray bboxes;
    bounding_boxes(bboxes, charimage);
    return getSkewAngle(bboxes);
}

double DeskewPageByRAST::getSkewAngle(rectarray &bboxes) {
    // Clean non-text and noisy boxes and get character statistics
    autodel<CharStats> charstats(make_CharStats());
    charstats->getCharBoxes(bboxes);
    charstats->calcCharStats();
    if(bboxes.length()>max_n)
        throw "too many character boxes for deskewing";

    // Extract textlines
    autodel<CTextlineRAST> ctextline(make_CTextlineRAST());
    narray<TextLine> textlines;
    ctextline->max_results=1;
    ctextline->min_gap = int(charstats->word_spacing*1.5);
    ctextline->setMaxSlope(0.5);
    ctextline->extract(textlines, charstats);
    if(textlines.length())
        return atan(textlines[0].m);
    else{
        fprintf(stderr,"Warning: no textlines found. ");
        fprintf(stderr,"Skipping deskewing ...\n");
        return 0;
    }
}

void DeskewPageByRAST::cleanup_gray(bytearray &image, bytearray &in) {
    cleanup(image,in);
}

void DeskewPageByRAST::cleanup(bytearray &image, bytearray &in) {
    makelike(image, in);
    float angle = (float) getSkewAngle(in);
    float cx = image.dim(0)/2.0;
    float cy = image.dim(1)/2.0;
    if(contains_only(in, byte(0), byte(255)))
        rotate_direct_sample(image, in, angle, cx, cy);
    else
        rotate_direct_interpolate(image, in, angle, cx, cy);
    if(debug_deskew) {
        fprintf(stderr, "Skew angle found = %.3f degrees\n", angle*RAD_TO_DEG);
        write_png(stdio(debug_deskew, "w"), image);
    }

}

ICleanupBinary *make_DeskewPageByRAST() {
    return new DeskewPageByRAST();
}
ICleanupGray *make_DeskewGrayPageByRAST() {
    return new DeskewPageByRAST();
}

