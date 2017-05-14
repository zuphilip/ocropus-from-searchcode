#include <stdlib.h>
#include "colib/colib.h"
#include "colib/iarith.h"
#include "iulib/iulib.h"
#include "iulib/ocrinterfaces.h"

using namespace colib;
using namespace iulib;
using namespace ocropus;

void invert(bytearray &a);
int binarize_simple(bytearray &result, bytearray &image);
int binarize_simple(bytearray &image);
void binarize_with_threshold(floatarray &out, floatarray &in, float threshold);
void binarize_with_threshold(bytearray &out, floatarray &in, float threshold);
void binarize_with_threshold(floatarray &in, float threshold);
void binarize_with_threshold(bytearray &out, bytearray &in, int threshold);
void binarize_with_threshold(bytearray &in, int threshold);
void optional_check_background_is_darker(bytearray &a);
int average_on_border(bytearray &a);
void optional_check_background_is_lighter(bytearray &a);
inline bool background_seems_black(bytearray &a) {
    return average_on_border(a) <= (min(a) + max(a) / 2);
}
inline bool background_seems_white(bytearray &a) {
    return average_on_border(a) >= (min(a) + max(a) / 2);
}

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
void extract_holes(bytearray &holes,bytearray &binarized);
void compute_troughs(floatarray &troughs,bytearray &binarized,float rsmooth=1.0);
void ridgemap(narray<floatarray> &maps,bytearray &binarized,
              float rsmooth=1.0,float asigma=0.7,float mpower=0.5,
              float rpsmooth=1.0);


#if 0
void check_line_segmentation(intarray &cseg);
void make_line_segmentation_black(intarray &a);
void make_line_segmentation_white(intarray &a);
void check_page_segmentation(intarray &cseg);
void make_page_segmentation_black(intarray &a);
void make_page_segmentation_white(intarray &a);
#endif

inline void write_line_segmentation(FILE *stream,intarray &a) {
    check_line_segmentation(a);
    make_line_segmentation_white(a);
    write_image_packed(stream,a,"png");
}

inline void read_line_segmentation(intarray &a,FILE *stream) {
    read_image_packed(a,stream,"png");
    check_line_segmentation(a);
    make_line_segmentation_black(a);
}

inline int pseg_pixel(int column,int paragraph,int line) {
    ASSERT((column > 0 && column < 32) || column == 254 || column == 255);
    ASSERT((paragraph >= 0 && paragraph < 64) || (paragraph >=251 && paragraph <= 255));
    ASSERT(line>=0 && line<256);
    return (column<<16) | (paragraph<<8) | line;
}

inline int pseg_pixel(int column,int line) {
    ASSERT(column>0 && column<32);
    ASSERT(line>=0 && line<64*256);
    return (column<<16) | line;
}

inline int pseg_column(int pixel) {
    return (pixel>>16)&0xff;
}

inline int pseg_paragraph(int pixel) {
    return (pixel>>8) & 0x3f;
}

inline int pseg_line(int pixel) {
    return pixel & 0xff;
}

inline int pseg_pline(int pixel) {
    return pixel & 0x3fff;
}

inline int cseg_pixel(int chr) {
    ASSERT(chr>0 && chr<4096);
    return (1<<12) | chr;
}

inline void pseg_columns(intarray &a) {
    for(int i=0;i<a.length1d();i++) {
        int value = a.at1d(i);
        if(value==0xffffff) value = 0;
        value = pseg_column(value);
        if(value>=32) value = 0;
        a.at1d(i) = value;
    }
}

inline void pseg_plines(intarray &a) {
    for(int i=0;i<a.length1d();i++) {
        int value = a.at1d(i);
        if(value==0xffffff) value = 0;
        if(pseg_column(value)>=32) value = 0;
        value = pseg_pline(value);
        a.at1d(i) = value;
    }
}

struct RegionExtractor {
    intarray segmentation;
    intarray ids;
    narray<rectangle> boxes;
    void idsAndBoxes() {
        intarray temp;
        temp = segmentation;
        renumber_labels(segmentation,1);
        int n = max(segmentation)+1;
        ids.resize(n);
        ids = -1;
        for(int i=0;i<temp.length();i++)
            ids[segmentation[i]] = temp[i];
        bounding_boxes(boxes,segmentation);
    }
    void setImage(intarray &image) {
        copy(segmentation,image);
        idsAndBoxes();
    }
    void setImageMasked(intarray &image,int mask,int lo,int hi) {
        makelike(segmentation,image);
        fill(segmentation,0);
        for(int i=0;i<image.length1d();i++) {
            int pixel = image.at1d(i);
            if(pixel<lo || pixel>hi) continue;
            segmentation.at1d(i) = (pixel & mask);
        }
        idsAndBoxes();
    }
    void setPageColumns(intarray &image) {
        makelike(segmentation,image);
        fill(segmentation,0);
        for(int i=0;i<image.length1d();i++) {
            int pixel = image.at1d(i);
            int col = pseg_column(pixel);
            if(col<1||col>=32) continue;
            int par = pseg_paragraph(pixel);
            if(par>=64) continue;
            segmentation.at1d(i) = col;
        }
        idsAndBoxes();
    }
    void setPageParagraphs(intarray &image) {
        makelike(segmentation,image);
        fill(segmentation,0);
        for(int i=0;i<image.length1d();i++) {
            int pixel = image.at1d(i);
            int col = pseg_column(pixel);
            if(col<1||col>=32) continue;
            int par = pseg_paragraph(pixel);
            if(par>=64) continue;
            segmentation.at1d(i) = (col<<8) | par;
        }
        idsAndBoxes();
    }
    void setPageLines(intarray &image) {
        makelike(segmentation,image);
        fill(segmentation,0);
        for(int i=0;i<image.length1d();i++) {
            int pixel = image.at1d(i);
            int col = pseg_column(pixel);
            if(col<1||col>=32) continue;
            int par = pseg_paragraph(pixel);
            if(par>=64) continue;
            segmentation.at1d(i) = pixel;
        }
        idsAndBoxes();
    }
    int length() {
        return boxes.length();
    }
    int id(int i) {
        return ids[i];
    }
    rectangle bbox(int i) {
        return boxes[i];
    }
    void bounds(int i,int *x0=0,int *y0=0,int *x1=0,int *y1=0) {
        *x0 = boxes[i].x0;
        *y0 = boxes[i].y0;
        *x1 = boxes[i].x1;
        *y1 = boxes[i].y1;
    }
    int x0(int i) {
        return boxes[i].x0;
    }
    int y0(int i) {
        return boxes[i].y0;
    }
    int x1(int i) {
        return boxes[i].x1;
    }
    int y1(int i) {
        return boxes[i].y1;
    }
    template <class S,class T>
        void extract(narray<S> &output,narray<T> &input,int index,int margin=0) {
        rectangle r = boxes[index].grow(margin);
        r.intersect(rectangle(0,0,input.dim(0),input.dim(1)));
        CHECK_CONDITION(!r.empty());
        extract_subimage(output,input,r.x0,r.y0,r.x1,r.y1);
    }
    template <class S>
    void mask(narray<S> &output,int index,int margin=0) {
        rectangle r = boxes[index].grow(margin);
        r.intersect(rectangle(0,0,segmentation.dim(0),segmentation.dim(1)));
        CHECK_CONDITION(!r.empty());
        output.resize(r.x1-r.x0,r.y1-r.y0);
        fill(output,0);
        for(int i=r.x0;i<r.x1;i++) for(int j=r.y0;j<r.y1;j++) {
                if(segmentation(i,j)==index)
                    output(i-r.x0,j-r.y0) = 255;
            }
    }
    void extract_masked(bytearray &output,bytearray &input,int index,int grow,int value,int margin=0) {
        extract(output,input,index,margin);
        if(grow<0) return;
        bytearray mask;
        this->mask(mask,index,margin);
        binary_dilate_circle(mask,grow);
        for(int i=0;i<mask.length();i++)
            output[i] = mask[i] ? output[i] : value;
    }
};
inline void write_page_segmentation(FILE *stream,intarray &a) {
    check_page_segmentation(a);
    write_image_packed(stream,a,"png");
}
inline void read_page_segmentation(intarray &a,FILE *stream) {
    read_image_packed(a,stream,"png");
    check_page_segmentation(a);
}

struct TextLineParam {
    float c,m,d; // c is y-intercept, m is slope, d is the line of descenders
    void print(FILE *stream=stdout){
        fprintf(stream,"%.3f %f %.2f\n",c,m,d);
    }
};

struct TextLine : TextLineParam{
    TextLine(){
    }
    TextLine(TextLineParam &tl){
        c = tl.c;
        m = tl.m;
        d = tl.d;
    }
    int xheight;
    rectangle bbox;
    void print(FILE *stream=stdout){
        fprintf(stream,"%d %d %d %d ",bbox.x0,bbox.y0,bbox.x1,bbox.y1);
        fprintf(stream,"%.3f %f %.2f %d\n",c,m,d,xheight);
    }
};

struct line {
    float c,m,d; // c is y-intercept, m is slope, d is the line of descenders
    float start,end,top,bottom; // start and end of line segment
    float istart,iend; //actual start and end of line segment in the image
    float xheight;

    line() {}
line(TextLine &tl):
    c(tl.c), m(tl.m), d(tl.d),
        start(tl.bbox.x0), end(tl.bbox.x1), top(tl.bbox.y0), bottom(tl.bbox.y1),
        istart(tl.bbox.x0), iend(tl.bbox.x1), xheight(tl.xheight) {}

    TextLine getTextLine(){
        TextLine tl;
        tl.c = c;
        tl.m = m;
        tl.d = d;
        tl.xheight = (int)xheight;
        //rectangle r((int)start, (int)top, (int)end, (int)bottom);
        rectangle r((int)istart, (int)top, (int)iend, (int)bottom);
        tl.bbox = r;
        return tl;
    }
};

struct TextLineParam4line {
    float c,m,d; // c is y-intercept, m is slope, d is the line of descenders
    float a,x; // a is ascender height, x is the x-height
    void print(FILE *stream=stdout){
        fprintf(stream,"%.3f %f %.2f\n",c,m,d);
    }
};

struct TextLineExtended : TextLineParam4line{
    TextLineExtended(){
    }
    TextLineExtended(TextLineParam4line &tl){
        c = tl.c;
        m = tl.m;
        d = tl.d;
        a = tl.a;
        x = tl.x;
    }
    rectangle bbox;
    void print(FILE *stream=stdout){
        fprintf(stream,"%d %d %d %d ",bbox.x0,bbox.y0,bbox.x1,bbox.y1);
        fprintf(stream,"%.3f %f %.2f %.2f %.2f\n",c,m,d,a,x);
    }
    TextLine getTextLine(){
        TextLine tl;
        tl.c = c;
        tl.m = m;
        tl.d = d;
        tl.xheight = (int)x;
        tl.bbox = bbox;
        return tl;
    }

};

using namespace colib;

void sort_boxes_by_x0(rectarray &boxes);
void sort_boxes_by_y0(rectarray &boxes);
int calc_xheight(rectarray &bboxes);

struct CharStats {
    // FIXME this violates coding conventions
    int img_height;
    int img_width;
    int xheight;
    int char_spacing;
    int word_spacing;
    int line_spacing;
    rectarray concomps;
    rectarray char_boxes;
    rectarray dot_boxes;
    rectarray large_boxes;

    CharStats();
    CharStats(CharStats &c);
    ~CharStats();
    void print();
    void getCharBoxes(rectarray &concomps);
    void setCharBoxes(rectarray &concomps) { getCharBoxes(concomps); }
    void calcCharStats();
    void calcCharStats(rectarray &cboxes);
    void calcCharStatsForOneLine();
    void calcCharStatsForOneLine(rectarray &cboxes);
};
CharStats *make_CharStats();
CharStats *make_CharStats(CharStats &c);

/////////////////////////////////////////////////////////////////////
///
/// \struct WhitespaceCover
/// Purpose: Whitespace Cover finding algorithm.
///
//////////////////////////////////////////////////////////////////////
enum qfunc {width, height, area};

class WhitespaceCover {
    // FIXME this violates coding conventions in several ways
private:
    int verbose;
    int max_results;
    float min_weight;
    float max_overlap;
    float min_aspect;
    float max_aspect;
    float min_width;
    float min_height;
    float logmin_aspect;
    bool greedy;
    rectangle bounds;
    qfunc quality_func;

    typedef shortarray Matches;
    typedef counted<Matches> CMatches;
    /////////////////////////////////////////////////////////////////////
    ///
    /// \struct WState
    /// Purpose: Current state of the Whitespace Cover finding algorithm.
    ///
    //////////////////////////////////////////////////////////////////////

    struct WState {
        int current_nrects;
        float weight;
        short top,left,bottom,right;
        rectangle bounds;
        CMatches matches;

        bool isDone(WhitespaceCover *env);
        void update(WhitespaceCover *env);
        int maxCentricity(WhitespaceCover *env);

    };

    typedef counted<WState> CState;

    rectarray rects;
    int initial_nrects;
    heap<CState> queue;
    narray<CState> results;
    void compute();
    bool goodDimensions(CState &result);
    void generateChildStates(CState &state, rectangle &pivot);
public:
    WhitespaceCover();
    WhitespaceCover(rectangle image_boundary);
    ~WhitespaceCover() {}
    void clear() {
        rects.clear();
        initial_nrects = 0;
        queue.clear();
        results.clear();
        bounds = rectangle(0,0,1,1);
    }
    void init();
    const char *description();
    void compute(rectarray &whitespaces, rectarray &obstacles);
    void addRect(rectangle r) {
        rects.push(r);
    }
    void setMaxResults(int value) {
        max_results = value;
    }
    void setMinWeight(float value) {
        min_weight = value;
    }
    void setMinWidth(float value) {
        min_width = value;
    }
    void setMinHeight(float value) {
        min_height = value;
    }
    void setBounds(int x0,int y0,int x1,int y1) {
        bounds = rectangle(x0,y0,x1,y1);
    }
    void setVerbose(int value) {
        verbose = value;
    }
    void setGreedy(bool value) {
        greedy = value;
    }
    void setMaxOverlap(float value) {
        max_overlap = value;
    }
    void setAspectRange(float min,float max) {
        min_aspect = min;
        max_aspect = max;
    }
    void setLogminAspect(float m) {
        logmin_aspect = m;
    }
    void setQfunc(qfunc t) {
        quality_func = t;
    }
    // Fit the bounds tightly to include all rectangles in the stack 'rects'.
    void snugBounds() {
        bounds = rectangle();
        for(int i=0;i<rects.length();i++) {
            bounds.include(rects[i]);
        }
    }
    int nSolutions() {
        return results.length();
    }
    void solution(int index,int &x0,int &y0,int &x1,int &y1) {
        rectangle &b = results[index]->bounds;
        x0 = b.x0;
        y0 = b.y0;
        x1 = b.x1;
        y1 = b.y1;
    }
};
WhitespaceCover *make_WhitespaceCover(rectangle &r);
WhitespaceCover *make_WhitespaceCover(int x0, int y0, int x1, int y1);

using namespace colib;
using namespace iulib;

/////////////////////////////////////////////////////////////////////
///
/// \struct ColSeparators
/// Purpose: Filtering columns from the whitespace cover
///
//////////////////////////////////////////////////////////////////////

struct ColSeparators {
    float min_aspect;
    int min_space;
    int max_space;
    float width;
    int max_boxes;
    int min_boxes;

    ColSeparators();
    void findGutters(rectarray &colboxes,
                     rectarray &whitespaceboxes,
                     CharStats &charstats);
    void filterOverlaps(rectarray &colboxes, rectarray &colcandidates);
};
ColSeparators *make_ColSeparators();

using namespace colib;

struct ExtractRulings {
    // Analyze user supplied obstacles to separate horizontal/vertical rulings
    // and image regions
    void analyzeObstacles(rectarray &hor_rulings,
                          rectarray &vert_rulings,
                          rectarray &images,
                          rectarray &obstacles,
                          int xheight);
};

ExtractRulings *make_ExtractRulings();

using namespace colib;


/////////////////////////////////////////////////////////////////////
///
/// \struct CTextlineRAST4line
/// Purpose: 4line implementation of the constrained textline finding
/// algorithm using RAST. Returns parameters of text-lines in
/// descending order of quality.
///
//////////////////////////////////////////////////////////////////////

static const int ntl4params = 5;
struct CTextlineRAST4line {
    CTextlineRAST4line();
    virtual ~CTextlineRAST4line(){
    }
    int generation;
    bool lsq;
    double epsilon;
    int maxsplits;
    double delta;
    double adelta;

    float min_length;
    int min_gap;
    double min_q;
    int min_count;
    int max_results;
    bool use_whitespace;

    // The parameters are:
    // all_params[0] = y-intercept of baseline
    // all_params[1] = slope
    // all_params[2] = descender distance from baseline
    // all_params[3] = xheight
    // all_params[4] = ascender distance from xheight
    typedef vecni<ntl4params> Parameters;
    double splitscale[ntl4params];
    Parameters all_params;
    Parameters empty_parameters;

    vec2i normalized(vec2i v) {
        interval a = atan2(v.y,v.x);
        return vec2i(cos(a),sin(a));
    }

    inline interval influence(bool lsq,interval d,double epsilon) {
        if(lsq) return sqinfluence(d,epsilon);
        else return rinfluence(d,epsilon);
    }

    typedef narray<int> Matches;
    rectarray cboxes;
    rectarray wboxes;
    narray<bool> used;

    bool final(interval q,const Parameters &p) {
        return p[0].width()<delta &&
            p[1].width()<adelta &&
            p[2].width()<delta &&
            p[3].width()<delta &&
            p[4].width()<delta;
    }

    struct TLState4line {
        short generation;
        short depth;
        short rank;
        signed char splits;
        bool splittable;
        interval quality;
        float priority;
        Parameters params;
        Matches matches;

        TLState4line();
        void set(CTextlineRAST4line &line,int depth,Parameters &params,
                 Matches &candidates,int splits);
        void reeval(CTextlineRAST4line &line);
        void update(CTextlineRAST4line &line, Matches &candidates);
        TextLineParam4line returnLineParam();
    };

    typedef counted<TLState4line> CState;
    heap<CState> queue;
    narray<CState> results;
    autodel<CharStats> linestats;
    Matches all_matches;

    void setDefaultParameters();
    void setMaxSlope(double max_slope);
    void setMaxYintercept(double ymin, double ymax);
    void prepare();
    void makeSubStates(narray<CState> &substates,CState &state);
    int wboxIntersection(CState &top);
    void search();
    virtual void pushResult(CState &result);
    virtual void extract(narray<TextLineParam4line> &textlines,
                         autodel<CharStats> &charstats);
    virtual void extract(narray<TextLineParam4line> &textlines,
                         rectarray &columns,
                         autodel<CharStats> &charstats);
};
CTextlineRAST4line *make_CTextlineRAST4line();

/////////////////////////////////////////////////////////////////////
///
/// \struct CTextlineRAST
/// Purpose: Constrained Textline finding using RAST
///
//////////////////////////////////////////////////////////////////////
struct CTextlineRASTExtended : CTextlineRAST4line{

    CTextlineRASTExtended();
    ~CTextlineRASTExtended(){ }
    // fraction of area covered by line bounding box
    // so that char_box is included in line_box
    float minoverlap;

    // rejection threshold for the height of a box = tr*xheight
    float min_box_height;

    // average distance between words
    int word_gap;

    int min_height;
    int assign_boxes;
    bool aggressive;
    int extend;
    int pagewidth;
    int pageheight;

    rectarray cboxes_all;
    narray<bool> used_all;
    narray<TextLineExtended> result_lines;

    void setDefaultParameters();
    void pushResult(CState &result);
    void extract(narray<TextLineExtended> &textlines,
                 autodel<CharStats> &charstats);
    void extract(narray<TextLineExtended> &textlines,
                 rectarray &columns,
                 autodel<CharStats> &charstats);
};
CTextlineRASTExtended *make_CTextlineRASTExtended();

/////////////////////////////////////////////////////////////////////
///
/// \struct CTextlineRASTBasic
/// Purpose: Basic implementation of the constrained textline finding
/// algorithm using RAST. Returns parameters of text-lines in
/// descending order of quality.
///
//////////////////////////////////////////////////////////////////////

static const int ntlparams = 3;
struct CTextlineRASTBasic {
    CTextlineRASTBasic();
    virtual ~CTextlineRASTBasic(){
    }
    int generation;
    bool lsq;
    double epsilon;
    int maxsplits;
    double delta;
    double adelta;

    float min_length;
    int min_gap;
    double min_q;
    int min_count;
    int max_results;
    bool use_whitespace;

    typedef vecni<ntlparams> Parameters;
    double splitscale[ntlparams];
    Parameters all_params;
    Parameters empty_parameters;

    vec2i normalized(vec2i v) {
        interval a = atan2(v.y,v.x);
        return vec2i(cos(a),sin(a));
    }

    inline interval influence(bool lsq,interval d,double epsilon) {
        if(lsq) return sqinfluence(d,epsilon);
        else return rinfluence(d,epsilon);
    }

    typedef narray<int> Matches;
    rectarray cboxes;
    rectarray wboxes;
    narray<bool> used;

    bool final(interval q,const Parameters &p) {
        return p[0].width()<delta &&
            p[1].width()<adelta &&
            p[2].width()<delta;
    }

    struct TLStateBasic {
        short generation;
        short depth;
        short rank;
        signed char splits;
        bool splittable;
        interval quality;
        float priority;
        Parameters params;
        Matches matches;

        TLStateBasic();
        void set(CTextlineRASTBasic &line,int depth,Parameters &params,
                 Matches &candidates,int splits);
        void reeval(CTextlineRASTBasic &line);
        void update(CTextlineRASTBasic &line, Matches &candidates);
        TextLineParam returnLineParam();
    };

    typedef counted<TLStateBasic> CState;
    heap<CState> queue;
    narray<CState> results;
    autodel<CharStats> linestats;
    Matches all_matches;

    void setDefaultParameters();
    void setMaxSlope(double max_slope);
    void setMaxYintercept(double ymin, double ymax);
    void prepare();
    void makeSubStates(narray<CState> &substates,CState &state);
    int wboxIntersection(CState &top);
    void search();
    virtual void pushResult(CState &result);
    virtual void extract(narray<TextLineParam> &textlines,
                         autodel<CharStats> &charstats);
    virtual void extract(narray<TextLineParam> &textlines,
                         rectarray &columns,
                         autodel<CharStats> &charstats);
    float getQuality(int i) {
        return results[i]->quality.center();
    }
};
CTextlineRASTBasic *make_CTextlineRASTBasic();

/////////////////////////////////////////////////////////////////////
///
/// \struct CTextlineRAST
/// Purpose: Constrained Textline finding using RAST
///
//////////////////////////////////////////////////////////////////////
struct CTextlineRAST : CTextlineRASTBasic{

    CTextlineRAST();
    ~CTextlineRAST(){ }
    // fraction of area covered by line bounding box
    // so that char_box is included in line_box
    float minoverlap;

    // rejection threshold for the height of a box = tr*xheight
    float min_box_height;

    // average distance between words
    int word_gap;

    int min_height;
    int assign_boxes;
    bool aggressive;
    int extend;
    int pagewidth;
    int pageheight;

    rectarray cboxes_all;
    narray<bool> used_all;
    narray<TextLine> result_lines;

    void setDefaultParameters();
    void pushResult(CState &result);
    void extract(narray<TextLine> &textlines,
                 autodel<CharStats> &charstats);
    void extract(narray<TextLine> &textlines,
                 rectarray &columns,
                 autodel<CharStats> &charstats);
};
CTextlineRAST *make_CTextlineRAST();

const int GUTTER_COLOR = 0x00ffff80;
const int IMAGE_COLOR = 0x0000ff00;
const int GRAPHICS_COLOR = 0x0000fe00;
const int RULING_COLOR = 0x0000fa00;

const int MULTI_COLUMN_ELEMENT_COLOR = 0x00fe0000;

struct ColorEncodeLayout{
    bool all;
    intarray outputImage;
    bytearray inputImage;
    rectarray textlines;
    rectarray paragraphs;
    rectarray textcolumns;
    rectarray gutters;
    rectarray rulings;
    rectarray sidebar;
    rectarray caption;
    rectarray table;
    rectarray graphics;
    rectarray images;
    rectarray pageNumber;
    rectarray header;
    rectarray footer;
    rectarray noise;

    ColorEncodeLayout() {
        all = false;
    }

    void encode();

private:
    void encode_textlines();
    void encode_gutters();

    void encode_zones(rectarray &zones, int zone_color);
    void encode_images();
    void encode_graphics();
    void encode_rulings();
};

ColorEncodeLayout *make_ColorEncodeLayout();

// Get text columns from an array of text-line bounding boxes and an array
// of whitespace column separators or vertical rulings
void get_text_columns(rectarray &textcolumns,
                      rectarray &textlines,
                      rectarray &gutters);

// Get text columns from an array of paragraphs or zones using their alignment
void get_text_columns(rectarray &columns,
                      rectarray &paragraphs);

void get_paragraphs(rectarray &paragraphs,
                    rectarray &textlines,
                    CharStats &charstats);

//Extend textlines to the nearest whitespace gutter or image start/end
void extend_lines(narray<line> &lines,
                  rectarray &wboxes,
                  int image_width);

struct ReadingOrderByTopologicalSort {
    ReadingOrderByTopologicalSort();
    ~ReadingOrderByTopologicalSort() {}

    // Extract reading order using whitespace gutters as separators
    void sortTextlines(narray<TextLine> &textlines,
                       rectarray &gutters,
                       CharStats &charstats);

    // Extract reading order using whitespace gutters and
    // vertical/horizontal rulings as separators
    void sortTextlines(narray<TextLine> &textlines,
                       rectarray &gutters,
                       rectarray &hor_rulings,
                       rectarray &vert_rulings,
                       CharStats &charstats);
private:
    int id;
    narray<int> val;
    narray<int> ro_index;

    void visit(int k, narray<bool> &lines_dag);
    void depthFirstSearch(narray<bool> &lines_dag);
};

ReadingOrderByTopologicalSort *make_ReadingOrderByTopologicalSort();

void visualize_layout(intarray &debug_image,
                      bytearray &in_not_inverted,
                      narray<TextLine> &textlines,
                      rectarray &gutters,
                      rectarray &extra_obstacles,
                      CharStats &charstats);

struct SegmentPageByRAST : ISegmentPage {
    SegmentPageByRAST();
    ~SegmentPageByRAST() {}

    const char *name() {
        return "segrast";
    }


    p_int max_results;
    p_int gap_factor;
    p_int use_four_line_model;
    p_int all_pixels;
    p_float max_descender;

    const char *description() {
        return "Segment page by RAST";
    }

    void init(const char **argv) {
    }

    void set(const char* var,double value){
#if 0
        if (strcmp(var,"max_results")==0){
            CHECK_ARG(value>=1.0 && value<=5000);
            max_results = int(value);
        }
        else if (strcmp(var,"gap_factor")==0){
            CHECK_ARG(value>=1.0 && value<=5000);
            gap_factor = int(value);
        }
        else if (strcmp(var,"use_four_line_model")==0)
            use_four_line_model = bool(value);
#else
        throw "unimplemented";
#endif
    }

    void segment(intarray &image,bytearray &in_not_inverted);
    void segment(intarray &image,bytearray &in_not_inverted,
                 rectarray &extra_obstacles);
    void visualize(intarray &result, bytearray &in_not_inverted,
                   rectarray &extra_obstacles);

private:
    void segmentInternal(intarray &visualization,
                         intarray &image,
                         bytearray &in_not_inverted,
                         bool need_visualization,
                         rectarray &extra_obstacles);


};


/// Get line information by character segmentation.
/// Having a true character segmentation is best, but CCs should also work.
///
/// @param intercept Y coordinate of intersection of the baseline with Oy axis.
/// @param slope Slope of the baseline (dy/dx, i.e. tangent)
/// @param xheight Height of letters.
/// @param descender_sink Distance between descender line and baseline (positive)
/// @param ascender_rise Distance between ascender line and baseline plus x-height
///
/// WARNING: this doesn't work reliably

bool get_extended_line_info(float &intercept, float &slope,
                            float &xheight, float &descender_sink,
                            float &ascender_rise, intarray &seg);

/// Get line information of a black-and-white line.
///
/// WARNING: this doesn't work reliably

bool get_extended_line_info_using_ccs(float &intercept, float &slope,
                                      float &xheight, float &descender_sink,
                                      float &ascender_rise, bytearray &img);

bool get_rast_info(float &intercept, float &slope,bytearray &img);

using namespace colib;

struct DeskewPageByRAST : virtual ICleanupBinary, virtual ICleanupGray {
    p_int max_n;
    DeskewPageByRAST() {
        max_n.bind(this,"max_n",10000,"maximum number of character boxes for deskewing");
    }
    ~DeskewPageByRAST() {
    }

    const char *interface() {
        return "ICleanupBinary";
    }

    const char *name() {
        return "deskewrast";
    }

    const char *description() {
        return "Deskew page image by RAST\n";
    }

    void init(const char **argv) {
        // nothing to be done
    }

    double getSkewAngle(bytearray &in);
    double getSkewAngle(rectarray &bboxes);
    void cleanup_gray(bytearray &image, bytearray &in);
    void cleanup(bytearray &image, bytearray &in);
};

const int white = 0xff;
const int black = 0x00;
bool is_line_image(bytearray &in);
struct NoiseFilter {
    int xstep;
    int ystep;
    int xwidth;
    int yheight;
    int max_width;
    int max_height;
    int border_margin;
    float bfthreshold;
    float wfthreshold;

    NoiseFilter();
    ~NoiseFilter();
    void blackFilter(bytearray &out,bytearray &in);
    void whiteFilter(bytearray &out,bytearray &in);
    void blackFilter(bytearray &out,bytearray &in,
                     float threshold,int xstep,int ystep,int xwidth,int yheight);
    void whiteFilter(bytearray &out,bytearray &original,
                     float threshold,int xstep,int ystep, int xwidth,int yheight);
    void remove(int x0,int y0,int x1,int y1,bytearray &image);
    void ccanalysis(bytearray &out,bytearray &in,
                    rectarray &bboxes);
    float blackRatio(int x0,int y0,int x1, int y1,bytearray &image);
    float whiteRatio(int x0,int y0, int x1, int y1, bytearray &image);
};
NoiseFilter *make_NoiseFilter();

ICleanupBinary *make_DocCleanConComp();

void remove_border_noise(bytearray &out,
                         bytearray &in,
                         rectangle &pageframe);

static const int npfparams = 4;
struct PageFrameRAST : ICleanupBinary {
    PageFrameRAST();
    virtual ~PageFrameRAST(){
    }
    const char *description() {
        return "remove marginal noise using page frame detection \n";
    }

    void init(const char **argv) {
        // nothing to be done
    }

    const char *name() {
        return "rastframe";
    }

    void cleanup(bytearray &out,bytearray &in);

    // FIXME/faisal do not expose instance variables in header files! --tmb

    int generation;
    bool lsq;
    double epsilon;
    int maxsplits;
    double delta;
    double adelta;

    float min_overlap;
    float min_length;
    int min_gap;
    double min_q;
    int min_count;
    int max_results;
    typedef vecni<npfparams> Parameters;
    double splitscale[npfparams];
    Parameters all_params;
    Parameters empty_parameters;

    vec2i normalized(vec2i v) {
        interval a = atan2(v.y,v.x);
        return vec2i(cos(a),sin(a));
    }

    inline interval influence(bool lsq,interval d,double epsilon) {
        if(lsq) return sqinfluence(d,epsilon);
        else return rinfluence(d,epsilon);
    }

    typedef narray<int> Matches;
    rectarray lineboxes;
    rectarray zoneboxes;
    narray<bool> used;

    bool final(interval q,const Parameters &p) {
        return p[0].width()<delta &&
            p[1].width()<delta &&
            p[2].width()<delta &&
            p[3].width()<delta;
    }

    struct PFStateBasic {
        short generation;
        short depth;
        short rank;
        signed char splits;
        bool splittable;
        interval quality;
        float priority;
        Parameters params;
        Matches matches;

        PFStateBasic();
        void set(PageFrameRAST &pf,int depth,Parameters &params,
                 Matches &candidates,int splits);
        void reeval(PageFrameRAST &pf);
        void update(PageFrameRAST &pf, Matches &candidates);
    };

    typedef counted<PFStateBasic> CState;
    heap<CState> queue;
    narray<CState> results;
    autodel<CharStats> linestats;
    Matches all_matches;

    void setDefaultParameters();
    void setSearchParameters(bytearray &in);
    void prepare();
    void makeSubStates(narray<CState> &substates,CState &state);
    void search();
    void pushResult(CState &result);
};

struct WhitespaceCuts : public ColSeparators{
    float max_aspect;

    WhitespaceCuts();

    void get_vborder(rectarray &columns,
                     CharStats &charstats);

    void find_hspaces(rectarray &hspaces,
                      rectarray &whitespaceboxes,
                      CharStats &charstats);
    void filter_hspaces(rectarray &hspaces,
                        CharStats &charstats);
    void filter_hanging_hspaces(rectarray &hspaces,
                                rectarray &columns,
                                CharStats &charstats);
    void get_whitespace_cuts(rectarray &wcuts,
                             CharStats &charstats);
};

WhitespaceCuts *make_WhitespaceCuts();

struct SegmentPageByWCUTS : ISegmentPage {
    ~SegmentPageByWCUTS() {}

    const char *name() {
        return "segwcuts";
    }
    const char *description() {
        return "segment characters by RAST\n";
    }

    void init(const char **argv) {
        // nothing to be done
    }

    void segment(intarray &image,bytearray &in_not_inverted);
};

ISegmentPage *make_SegmentPageByWCUTS();

enum {HORIZONTAL_CUT, VERTICAL_CUT};

struct SegmentPageByXYCUTS : ISegmentPage {
private:
    unsigned int tnx; // noise threshold on projection on x-axis
    unsigned int tny; // noise threshold on projection on y-axis
    unsigned int tcx; // min gap size on x-axis projection
    unsigned int tcy; // min gap size on y-axis projection
public:
    SegmentPageByXYCUTS();
    SegmentPageByXYCUTS(unsigned int itnx,
                        unsigned int itny,
                        unsigned int itcx,
                        unsigned int itcy);

    ~SegmentPageByXYCUTS() {}

    const char *name() {
        return "segxy";
    }

    const char *description() {
        return "segment page by XY-Cut algorithm\n"
            "Default parameters: \n"
            "\ttnx=78, tny=32, tcx=35, tcy=54\n"
            "\ttnx,tny: cleaning trhesholds\n"
            "\ttcx,tcy = min gap size hor. and ver.\n" ;
    }

    void setParameters(unsigned int itnx, unsigned int itny, unsigned int itcx, unsigned int itcy);

    void segment(intarray &image,bytearray &in);
};

ISegmentPage *make_SegmentPageByXYCUTS(unsigned int itnx,
                                       unsigned int itny,
                                       unsigned int itcx,
                                       unsigned int itcy);

ISegmentPage *make_SegmentWords();
bool segment_words_by_projection(intarray &seg,
                                 bytearray &in, int nwords);

/// Assert that all items of this array of arrays are not empty.
template<class T>
void assert_all_items_nonempty(narray<narray<T> > &a) {
    for(int i=0;i<a.length();i++)
        ASSERT(a[i].length1d() > 0);
}

/// Remove from the segmentation those pixels which are white in gray_image.
void binarize_in_segmentation(intarray &segmentation, /* const */ bytearray &gray_image);

/// Set line number for all foreground pixels in a character segmentation.
void set_line_number(intarray &a, int lnum);


/// Unpack page segmentation into separate line masks with bounding boxes.
void extract_lines(narray<bytearray> &lines,narray<rectangle> &rboxes,intarray &image);

/// If the line is too small or too large, rescale it (with the mask)
/// to a decent height (30-60 pixels).
void rescale_if_needed(bytearray &bin_line, bytearray &gray_line);

/// Make a binary image from a line segmentation.
void forget_segmentation(bytearray &image, intarray &segmentation);

/// Returns true if there's a mapping between s1's colors and s2's colors.
bool is_oversegmentation_of(intarray &s1, intarray &s2);

/// Return true if there are no zeros in the array.
bool has_no_black_pixels(intarray &);

void blit_segmentation_line(intarray &page,
                            rectangle bbox,
                            intarray &line,
                            int line_no);

/// Blit the segmentation of src onto dst shifted by (x,y) and shifted by
/// values by max(dst).
void concat_segmentation(intarray &dst, intarray &src,
                         int x, int y);

// Enlarge segmentation and AND it with line_mask.
// Don't pass binarized grayscale image as line_mask,
// otherwise you might get debris not from the line.
// (that means we cannot really call this from inside LineOCR)
void normalize_segmentation(intarray &segmentation, bytearray &line_mask);

int max_cnum(intarray &seg);

void get_recoloring_map(intarray &recolor, intarray &image);
void remove_gaps_by_recoloring(intarray &image);

void ocr_bboxes_to_charseg(intarray &cseg,narray<rectangle> &bboxes,intarray &segmentation);
void evaluate_segmentation(int &nover,int &nunder,int &nmis,
                           intarray &model_raw,intarray &image_raw,float tolerance);
void align_segmentation(intarray &segmentation,narray<rectangle> &bboxes);

extern int log_reg_class_num;
extern float log_reg_factor;
extern float log_reg_offset;
extern int log_reg_feature_len;
extern float log_reg_data[];

const int MAX_LEN = 128; // MAX_LEN = dimension of the histogram
const int COMP_START = 0;
const int COMP_LEN_START = 1;
const int COMP_LEN_INC_A = 2;
const int COMP_LEN_INC_B = 0;

enum zone_class {math=0, logo=1, text=2, table=3, drawing=4,
                 halftone=5, ruling=6, noise=7, undefined=-1};

struct ZoneFeatures{

    void extractFeatures(floatarray &features, bytearray &image);

    void horizontalRunLengths(floatarray &resulthist,
                              floatarray &resultstats,
                              const bytearray &image);
    void verticalRunLengths(floatarray &resulthist,
                            floatarray &resultstats,
                            const bytearray &image);
    void mainDiagRunLengths(floatarray &resulthist,
                            floatarray &resultstats,
                            const bytearray &image);
    void sideDiagRunLengths(floatarray &resulthist,
                            floatarray &resultstats,
                            const bytearray &image);

    void compressHist(intarray &histogram);
    void compress2DHist(intarray &histogram);

    void concompHist(floatarray &result,
                     rectarray &concomps);

    void concompNeighbors(floatarray &result,
                          rectarray &concomps);

};

ZoneFeatures *make_ZoneFeatures();

struct LogReg{
    int feature_len;
    int class_num;
    float factor;
    float offset;
    floatarray lambda;

    void loadData();
    zone_class classify(floatarray &feature);
    void getClassProbabilities(floatarray &prob,
                               floatarray &feature);
};

LogReg *make_LogReg();

const int math_color = 0x0001fa01;
const int logo_color = 0x0001fb01;
const int text_color = 0x00ff0101;
const int table_color = 0x0001fd01;
const int drawing_color = 0x0001fe01;
const int halftone_color = 0x0001ff01;
const int ruling_color = 0x0001fc01;
const int noise_color = 0x00ffff00;

struct TextImageSegByLogReg : ITextImageClassification {
    ~TextImageSegByLogReg() {}

    const char *description() {
        return "Get text/image probability map\n";
    }

    const char *name() {
        return "tiseglogreg";
    }

    void init(const char **argv) {
        // nothing to be done
    }

    // Get text-image map from a segmented image
    void textImageProbabilities(intarray &out, intarray &in);
    // Get text-image map from a binary image
    void textImageProbabilities(intarray &out, bytearray &in);

    void getProbabilityMap(floatarray &class_prob,
                           rectarray &boxes,
                           bytearray &image);

    int getColor(floatarray &prob_map, int index);

};

ITextImageClassification *make_TextImageSegByLogReg();

#ifdef HAVE_LEPTONICA
struct RemoveImageRegions : virtual ICleanupBinary, ICleanupGray{
    ~RemoveImageRegions() {}

    const char *name() {
        return "removeimageregions";
    }

    const char *description() {
        return "Remove text or non-text zones\n";
    }

    void init(const char **argv) {
        // nothing to be done
    }

    void cleanup(bytearray &out, bytearray &in);

};

ICleanupBinary *make_RemoveImageRegionsBinary();
ICleanupGray *make_RemoveImageRegionsGray();
#endif

ISegmentPage *make_SegmentPageBy1CP();
ISegmentPage *make_SegmentPageByRAST();
ISegmentPage *make_SegmentPageByRAST1();
ISegmentPage *make_SegmentPageByMorphTrivial();
ISegmentPage *make_SegmentPageByWCUTS();
ISegmentPage *make_SegmentPageByXYCUTS();
ISegmentPage *make_SegmentWords();

