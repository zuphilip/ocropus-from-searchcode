#include "colib/colib.h"
#include "iulib/iulib.h"
#include "iulib/ocrinterfaces.h"

using namespace colib;
using namespace iulib;
using namespace ocropus;

//void make_line_segmentation_black(intarray &);

ISegmentLine *make_SegmentLineByCCS();
ISegmentLine *make_SegmentLineByGCCS();
ISegmentLine *make_ConnectedComponentSegmenter();
ISegmentLine *make_CurvedCutSegmenter1();
ISegmentLine *make_CurvedCutSegmenter();
ISegmentLine *make_CurvedCutWithCcSegmenter();

struct IDpSegmenter : ISegmentLine {
    float down_cost;
    float outside_diagonal_cost;
    float outside_diagonal_cost_r;
    float inside_diagonal_cost;
    float boundary_diagonal_cost;
    float inside_weight;
    float boundary_weight;
    float outside_weight;
    int min_range;
    float cost_smooth;
    float min_thresh;
    intarray dimage;
};

IDpSegmenter *make_DpSegmenter();

