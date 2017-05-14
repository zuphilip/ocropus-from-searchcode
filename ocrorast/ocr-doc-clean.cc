#include "ocrorast.h"

param_string debug_cleanup("debug_cleanup",0,"output cleaned document image as png");

struct DocClean : ICleanupBinary {
    ~DocClean() {}

    const char *description() {
        return "Running black filter,white filter on the image  thus removing noise \n";
    }

    void init(const char **argv) {
        // nothing to be done
    }

    const char *name() {
        return "docclean";
    }

    void cleanup(bytearray &out,bytearray &in) {

        autodel<NoiseFilter> noisefilter(make_NoiseFilter());

        //Check if it is the complete image
        if(is_line_image(in)){
            copy(out,in);
            return ;
        }
        //run black filter first
        bytearray result_bf,result_cc,result_bf_inverted;
        makelike(result_cc,in);
        fill(result_cc,0xff);
        noisefilter->blackFilter(result_bf,in);

        // Do connected component analysis
        intarray charimage;
        copy(result_bf_inverted,result_bf);
        make_page_binary_and_black(result_bf_inverted);
        copy(charimage,result_bf_inverted);
        label_components(charimage,false);

        // Clean noisy boxes
        rectarray bboxes;
        bounding_boxes(bboxes,charimage);
        noisefilter->ccanalysis(result_cc,result_bf,bboxes);

        //run whitefilter on output of connected component analysis result
        noisefilter->whiteFilter(out,result_cc);
        if(debug_cleanup) {
            write_png(stdio(debug_cleanup,"w"), out);
        }
    }
};

ICleanupBinary *make_DocClean() {
    return new DocClean();
}
