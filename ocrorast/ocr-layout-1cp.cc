#include "ocrorast.h"

struct SegmentPageBy1CP : ISegmentPage {
    ~SegmentPageBy1CP() {}

    const char *description() {
        return "segment characters by horizontal projection (assumes single column)\n";
    }

    void init(const char **argv) {
        // nothing to be done
    }

    const char *name() {
        return "seg1cp";
    }

    void segment(intarray &image,bytearray &in) {
        param_int thigh("thigh",20,"projection threshold");
        param_int tlow("tlow",1,"projection threshold");

        optional_check_background_is_lighter(in);
        bytearray binarized;
        binarize_simple(binarized, in);
        copy(image,in);

        for(int i=0;i<image.length1d();i++)
            image.at1d(i) = !image.at1d(i);

        intarray projection(image.dim(1));
        for(int j=0;j<image.dim(1);j++) {
            projection(j) = 0;
            for(int i=0;i<image.dim(0);i++)
                projection(j) += image(i,j);
        }

        projection.reshape(projection.dim(0),1);
        floatarray temp;
        copy(temp,projection);
        hysteresis_thresholding(temp,tlow,thigh);
        copy(projection,temp);
        projection.reshape(projection.dim(0));

        int count = 0;
        for(int j=1;j<projection.length();j++) {
            if(projection(j-1)==0 && projection(j)==1)
                count++;
            if(projection(j))
                projection(j) = count;
        }
        for(int j=0;j<image.dim(1);j++) {
            for(int i=0;i<image.dim(0);i++) {
                if(projection(j)==0) {
                    image(i,j) = 0;
                } else if(image(i,j)) {
                    // the count-x+1 is there to reverse the order
                    int line = count-projection(j)+1;
                    image(i,j) = pseg_pixel(1,line);
                }
            }
        }
        make_page_segmentation_white(image);
        check_page_segmentation(image);
    }
};

ISegmentPage *make_SegmentPageBy1CP() {
    return new SegmentPageBy1CP();
}
