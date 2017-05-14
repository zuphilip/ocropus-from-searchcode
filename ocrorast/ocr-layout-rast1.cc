#include "ocrorast.h"

struct SegmentPageByRAST1 : ISegmentPage {
    SegmentPageByRAST1() {
        max_results.bind(this,"max_results",1000,"maximum number of results");
        debug_segm.bind(this,"debug_segm","","output segmentation file");
        debug_layout.bind(this,"debug_layout",0,"print the intermediate results to stdout");
        max_descender.bind(this,"max_descender",20,"maximum descender");
    }

    ~SegmentPageByRAST1() {}

    p_int max_results;
    p_string debug_segm;
    p_int debug_layout;
    p_float max_descender;

    const char *name() {
        return "segrast1";
    }

    const char *description() {
        return "Segment page by RAST, assuming single column.";
    }

    void segmentInternal(intarray &image,
            bytearray &in_not_inverted,
            rectarray &extra_obstacles) {

        CHECK(extra_obstacles.length()==0);

        const int zero = 0;
        const int yellow = 0x00ffff00;
        bytearray in;
        copy(in, in_not_inverted);
        make_page_binary_and_black(in);

        intarray charimage;
        copy(charimage,in);
        label_components(charimage,false);

        rectarray bboxes;
        bounding_boxes(bboxes,charimage);
        if(bboxes.length()==0){
            makelike(image,in);
            fill(image,0x00ffffff);
            return ;
        }

        autodel<CharStats> charstats(make_CharStats());
        charstats->getCharBoxes(bboxes);
        charstats->calcCharStats();
        if(int(debug_layout)>=2){
            charstats->print();
        }

        narray<TextLine> textlines;

        autodel<CTextlineRAST> ctextline(make_CTextlineRAST());
        ctextline->setDefaultParameters();
        ctextline->min_q = 2.0;
        ctextline->min_count = 2;
        ctextline->min_length= 30;
        ctextline->max_results = max_results;
        ctextline->min_gap = 20000;
        ctextline->use_whitespace = 0;
        ctextline->all_params[2] = interval(0,max_descender);
        rectarray textline_obstacles;
        ctextline->extract(textlines,textline_obstacles,charstats);

        autodel<ReadingOrderByTopologicalSort>
            reading_order(make_ReadingOrderByTopologicalSort());
        rectarray gutters,hor_rulings,vert_rulings;
        reading_order->sortTextlines(textlines,gutters,hor_rulings,vert_rulings,*charstats);

        rectarray textcolumns;
        rectarray paragraphs;
        rectarray textline_boxes;
        for(int i=0, l=textlines.length(); i<l; i++)
            textline_boxes.push(textlines[i].bbox);
        autodel<ColorEncodeLayout> color_encoding(make_ColorEncodeLayout());
        copy(color_encoding->inputImage,in);
        for(int i=0, l=textlines.length(); i<l; i++)
            color_encoding->textlines.push(textlines[i].bbox);
        for(int i=0, l=textcolumns.length(); i<l; i++)
            color_encoding->textcolumns.push(textcolumns[i]);

        color_encoding->all = true;
        color_encoding->encode();

        intarray &segmentation = color_encoding->outputImage;

        image.makelike(segmentation);
        image = 0xffffff;
        for(int i=1;i<bboxes.length();i++) {
            rectangle b = bboxes[i];
            int xc = b.xcenter();
            int yc = b.ycenter();
            if(xc<0 || xc>=image.dim(0) || yc>=image.dim(1) || yc<0) continue;
            int pixel = segmentation(xc,yc);
            for(int x=b.x0;x<b.x1;x++) {
                for(int y=b.y0;y<b.y1;y++) {
                    if(charimage(x,y)==i) image(x,y) = pixel;
                }
            }
        }
#if 0
        write_image_binary("_in.png",in);
        write_image_packed("_charimage.png",charimage);
        write_image_packed("_segmentation.png",segmentation);
#endif

        replace_values(image,zero,yellow);
    }

    void segment(intarray &result,
            bytearray &in_not_inverted,
            rectarray &obstacles) {
        intarray debug_image;
        segmentInternal(result, in_not_inverted, obstacles);
    }

    void segment(intarray &result,
            bytearray &in_not_inverted) {
        rectarray obstacles;
        segment(result,in_not_inverted,obstacles);
    }

};

ISegmentPage *make_SegmentPageByRAST1() {
    return new SegmentPageByRAST1();
}
