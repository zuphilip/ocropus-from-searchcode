#include "ocrorast.h"


struct DocCleanConComp : ICleanupBinary {
    ~DocCleanConComp() {}

    const char *description() {
        return "Remove large components from the image\n";
    }

    void init(const char **argv) {
        // nothing to be done
    }

    const char *name() {
        return "doccleancon";
    }


    void cleanup(bytearray &out,bytearray &in) {

        int page_width = in.dim(0);
        int page_height = in.dim(1);

        intarray charimage;
        makelike(charimage,in);
        for(int i=0,l=in.length1d(); i<l; i++)
            charimage.at1d(i) = !in.at1d(i);
        label_components(charimage,false);

        // Clean non-text and noisy boxes and get character statistics
        rectarray bboxes;
        bounding_boxes(bboxes,charimage);
        ASSERT(bboxes.length()!=0);

        // get char stats
        autodel<CharStats> charstats(make_CharStats());
        charstats->getCharBoxes(bboxes);

        int dotlength = charstats->dot_boxes.length();
        int charlength = charstats->char_boxes.length();

        makelike(out,in);
        fill(out,255);
        for(int i=0; i<dotlength; i++){
            rectangle b = charstats->dot_boxes[i];
            b.x1=min(b.x1,page_width-1);
            b.y1=min(b.y1,page_height-1);
            for(int x=b.x0; x<=b.x1; x++)
                for(int y=b.y0; y<=b.y1; y++)
                    out(x,y) = in(x,y);
        }

        for(int i=0; i<charlength; i++){
            rectangle b = charstats->char_boxes[i];
            b.x1=min(b.x1,page_width-1);
            b.y1=min(b.y1,page_height-1);
            for(int x=b.x0; x<=b.x1; x++)
                for(int y=b.y0; y<=b.y1; y++)
                    out(x,y) = in(x,y);
        }
    }
};

ICleanupBinary *make_DocCleanConComp() {
    return new DocCleanConComp();
}
