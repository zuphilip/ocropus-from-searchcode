#include "ocrolseg.h"
#include "queue.h"

bool bgcheck = 1;

void invert(bytearray &a) {
    int n = a.length1d();
    for (int i = 0; i < n; i++) {
        a.at1d(i) = 255 - a.at1d(i);
    }
}

int binarize_simple(bytearray &result, bytearray &image) {
    int threshold = (max(image)+min(image))/2;
    makelike(result,image);
    for(int i=0;i<image.length1d();i++)
        result.at1d(i) = image.at1d(i)<threshold ? 0 : 255;
    return threshold;
}

int binarize_simple(bytearray &image) {
    return binarize_simple(image, image);
}

template <class T>
void remove_small_components(narray<T> &bimage,int mw,int mh) {
    intarray image;
    copy(image,bimage);
    label_components(image);
    narray<rectangle> rects;
    bounding_boxes(rects,image);
    if(rects.length()==0) return;
    bytearray good(rects.length());
    for(int i=0;i<good.length();i++)
        good[i] = 1;
    for(int i=0;i<rects.length();i++)
        if(rects[i].width()<mw && rects[i].height()<mh)
            good[i] = 0;
    for(int i=0;i<image.length();i++)
        if(!good(image[i]))
            image[i] = 0;
    for(int i=0;i<image.length1d();i++)
        if(!image[i])
            bimage[i] = 0;
}

template void remove_small_components<byte>(narray<byte> &,int,int);
template void remove_small_components<int>(narray<int> &,int,int);

inline int cseg_pixel(int chr) {
    ASSERT(chr>0 && chr<4096);
    return (1<<12) | chr;
}

#if 0
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
#endif

int average_on_border(bytearray &a) {
    int sum = 0;
    int right = a.dim(0) - 1;
    int top = a.dim(1) - 1;
    for(int x = 0; x < a.dim(0); x++)
        sum += a(x, 0);
    for(int x = 0; x < a.dim(0); x++)
        sum += a(x, top);
    for(int y = 1; y < top; y++)
        sum += a(0, y);
    for(int y = 1; y < top; y++)
        sum += a(right, y);
    return sum / ((right + top) * 2);
}

void optional_check_background_is_darker(bytearray &a) {
    if(bgcheck) {
        if(!(average_on_border(a) <= (min(a) + max(a) / 2)))
            throw "line image probably not right in line segmenter";
    }
}

void optional_check_background_is_lighter(bytearray &a) {
    if(bgcheck) {
        if(!(average_on_border(a) >= (min(a) + max(a) / 2)))
            throw "line image probably not right in line segmenter";

    }
}

#if 0
void make_line_segmentation_black(intarray &a) {
    check_line_segmentation(a);
    replace_values(a, 0xFFFFFF, 0);
    for(int i = 0; i < a.length1d(); i++)
        a.at1d(i) &= 0xFFF;
}

void make_line_segmentation_white(intarray &a) {
    replace_values(a, 0, 0xFFFFFF);
    //for(int i = 0; i < a.length1d(); i++)
    //    a.at1d(i) = (a.at1d(i) & 0xFFF) | 0x1000;
    check_line_segmentation(a);
}
#endif

struct SegmentLineByProjection : ISegmentLine {
    ~SegmentLineByProjection() {}

    const char *name() {
        return "projseg";
    }

    const char *description() {
        return "segment characters by 1D projection\n";
    }

    void init(const char **argv) {
        // nothing to be done
    }

    void charseg(intarray &image,bytearray &in) {
        optional_check_background_is_lighter(in);
        param_int thigh("thigh",5,"projection threshold");
        param_int tlow("tlow",1,"projection threshold");

        copy(image,in);
        for(int i=0;i<image.length1d();i++)
            image.at1d(i) = !image.at1d(i);

        intarray projection(image.dim(0));
        for(int i=0;i<image.dim(0);i++) {
            projection(i) = 0;
            for(int j=0;j<image.dim(1);j++) {
                projection(i) += image(i,j);
            }
        }

        projection.reshape(projection.dim(0),1);
        floatarray temp;
        copy(temp,projection);
        hysteresis_thresholding(temp,tlow,thigh);
        copy(projection,temp);
        projection.reshape(projection.dim(0));

        int count = 0;
        for(int i=1;i<projection.length();i++) {
            if(projection(i-1)==0 && projection(i)==1)
                count++;
            if(projection(i))
                projection(i) = count;
        }
        for(int i=0;i<image.dim(0);i++) {
            for(int j=0;j<image.dim(1);j++) {
                if(!projection(i)) {
                    image(i,j) = 0;
                } else if(image(i,j)) {
                    image(i,j) = cseg_pixel(projection(i));
                }
            }
        }
        make_line_segmentation_white(image);
        check_line_segmentation(image);
    }
};

ISegmentLine *make_SegmentLineByProjection() {
    return new SegmentLineByProjection();
}

bool is_singular0(bytearray &image,int i,int j) {
    if(i<1||i>=image.dim(0)-1||j<1||j>image.dim(1)-1) return 0;
    bytearray surround(8);
    int k = 0;
    surround(k++) = image(i+1,j);
    surround(k++) = image(i+1,j+1);
    surround(k++) = image(i,j+1);
    surround(k++) = image(i-1,j+1);
    surround(k++) = image(i-1,j);
    surround(k++) = image(i-1,j-1);
    surround(k++) = image(i,j-1);
    surround(k++) = image(i+1,j-1);
    // count the number of transitions
    int count = 0;
    for(k=0;i<8;k++) if(surround((k+1)%8)&&!surround(k)) count++;
    return count>2;
}
inline int neighbors(bytearray &image,int i,int j) {
    if(i<1||i>=image.dim(0)-1||j<1||j>image.dim(1)-1) return 0;
    if(!image(i,j)) return 0;
    int count = -1;
    for(int k=-1;k<=1;k++)
        for(int l=-1;l<=1;l++)
            if(image(i+k,j+l)) count++;
    return count;
}
void count_neighbors(bytearray &result,bytearray &image) {
    makelike(result,image);
    fill(result,0);
    int d=1;
    for(int i=d;i<image.dim(0)-d-1;i++) {
        for(int j=d;j<image.dim(1)-d-1;j++) {
            int n = neighbors(image,i,j);
            result(i,j) = n;
        }
    }
}
inline bool is_singular(bytearray &image,int i,int j) {
    return neighbors(image,i,j)>2;
}
void find_endpoints(bytearray &result,bytearray &image) {
    int d = 1;
    makelike(result,image);
    fill(result,0);
    for(int i=d;i<image.dim(0)-d-1;i++) {
        for(int j=d;j<image.dim(1)-d-1;j++) {
            int n = neighbors(image,i,j);
            CHECK_ARG(n<5);
            if(n==1) result(i,j) = 255;
        }
    }
}
void find_junctions(bytearray &result,bytearray &image) {
    int d=1;
    makelike(result,image);
    fill(result,0);
    for(int i=d;i<image.dim(0)-d-1;i++) {
        for(int j=d;j<image.dim(1)-d-1;j++) {
            int n = neighbors(image,i,j);
            CHECK_ARG(n<5);
            if(n>2) result(i,j) = 255;
        }
    }
}
void remove_singular_points(bytearray &image,int d) {
    for(int i=d;i<image.dim(0)-d-1;i++) {
        for(int j=d;j<image.dim(1)-d-1;j++) {
            if(is_singular(image,i,j)) {
                for(int k=-d;k<=d;k++)
                    for(int l=-d;l<=d;l++)
                        image(i+k,j+l) = 0;
            }
        }
    }
}
class SkelSegmenter : public ISegmentLine {
    virtual const char *description() {
        return "skeleton segmenter";
    }
    const char *name() {
        return "skelseg";
    }
    virtual void init(const char **argv=0) {
    }

    virtual void charseg(intarray &segmentation,bytearray &image) {
        bytearray timage;
        copy(timage,image);
        for(int i=0;i<image.length();i++) image[i] = !image[i];
        thin(timage);
        //write_png(stdio("_thinned","w"),timage);
        remove_singular_points(timage,2);
        //write_png(stdio("_segmented","w"),timage);
        intarray tsegmentation;
        copy(tsegmentation,timage);
        label_components(tsegmentation);
        remove_small_components(tsegmentation,4,4);
        //write_png_rgb(stdio("_labeled","w"),tsegmentation);
        copy(segmentation,image);
        propagate_labels_to(segmentation,tsegmentation);
        //write_png_rgb(stdio("_propagated","w"),segmentation);
    }
};

ISegmentLine *make_SkelSegmenter() {
    return new SkelSegmenter();
}



struct SegmentLineByCCS : ISegmentLine {
    SegmentLineByCCS() {
        pdef("swidth",0,"smearing width");
        pdef("sheight",10,"smearing height");
    }

    ~SegmentLineByCCS() {}

    const char *description() {
        return "connected component segmenter using morphology for grouping\n";
    }

    const char *name() {
        return "segccs";
    }

    void init(const char **argv) {
        // nothing to be done
    }

    void charseg(intarray &image,bytearray &in) {
        using namespace narray_ops;
        int swidth = pgetf("swidth");
        int sheight = pgetf("sheight");
        bytearray temp;
        copy(image,in);
        for(int i=0;i<image.length();i++) image[i] = !image[i];
        copy(temp,image);
        if(swidth>0||sheight>0) binary_close_rect(temp,swidth,sheight);
        intarray labels;
        copy(labels,temp);
        label_components(labels);
        for(int i=0;i<image.length1d();i++)
            if(image.at1d(i))
                image.at1d(i) = cseg_pixel(labels.at1d(i));
        make_line_segmentation_white(image);
        check_line_segmentation(image);
    }
};

struct SegmentLineByGCCS : ISegmentLine {
    SegmentLineByGCCS() {
    }

    ~SegmentLineByGCCS() {}

    const char *description() {
        return "connected component segmenter using grouping by overlap\n";
    }

    const char *name() {
        return "seggccs";
    }

    void init(const char **argv) {
        // nothing to be done
    }

    void charseg(intarray &image,bytearray &in) {
        using namespace narray_ops;
        bytearray temp;
        copy(image,in);
        for(int i=0;i<image.length();i++)
            image[i] = !image[i];
        copy(temp,image);
        intarray labels;
        copy(labels,temp);
        label_components(labels);
        narray<rectangle> boxes;
        bounding_boxes(boxes,labels);
        intarray equiv(boxes.length());
        for(int i=0;i<boxes.length();i++)
            equiv[i] = i;
        for(int i=1;i<boxes.length();i++) {
            rectangle p = boxes[i];
            for(int j=1;j<boxes.length();j++) {
                if(i==j) continue;
                rectangle q = boxes[j];
                int x0 = max(p.x0,q.x0);
                int x1 = min(p.x1,q.x1);
                int iw = x1-x0;
                if(iw<=0) continue; // no overlap
                int ow = min(p.width(),q.width());
                float frac = iw/float(ow);
                if(frac<0.5) continue; // insufficient overlap
                // printf("%d %d : %d %d : %g\n",i,j,iw,ow,frac);
                equiv[max(i,j)] = min(i,j);
            }
        }
        for(int i=0;i<labels.length();i++)
            labels[i] = equiv[labels[i]];
        renumber_labels(labels,1);
        image = labels;
        make_line_segmentation_white(image);
        check_line_segmentation(image);
    }
};

ISegmentLine *make_SegmentLineByCCS() {
    return new SegmentLineByCCS();
}
ISegmentLine *make_SegmentLineByGCCS() {
    return new SegmentLineByGCCS();
}

class ConnectedComponentSegmenter : public ISegmentLine {
    virtual const char *description() {
        return "connected component segmenter";
    }
    const char *name() {
        return "segccs";
    }
    virtual void init(const char **argv=0) {
    }

    virtual void charseg(intarray &segmentation,bytearray &image) {
        bytearray temp_image;
        copy(temp_image, image);
        binary_autoinvert(temp_image);
        copy(segmentation,temp_image);
        label_components(segmentation);
    }
};

ISegmentLine *make_ConnectedComponentSegmenter() {
    return new ConnectedComponentSegmenter();
}

void combine_segmentations(intarray &dst, intarray &src) {
    CHECK_ARG(samedims(dst, src));
    int n = max(dst) + 1;
    for(int i = 0; i < dst.length1d(); i++)
        dst.at1d(i) += src.at1d(i) * n;
    renumber_labels(dst, 1);
}

void local_min(floatarray &result,floatarray &data,int r) {
    int n = data.length();
    result.resize(n);
    for(int i=0;i<n;i++) {
        float lmin = data(i);
        for(int j=-r;j<=r;j++) {
            int k = i+j;
            if(unsigned(k)>=unsigned(n)) continue;
            if(data(k)>=lmin) continue;
            lmin = data(k);
        }
        result(i) = lmin;
    }
}

void local_minima(intarray &result,floatarray &data,int r,float threshold) {
    int n = data.length();
    result.clear();
    floatarray lmin;
    local_min(lmin,data,r);
    for(int i=1;i<n-1;i++) {
        if(data(i)<=threshold && data(i)<=lmin(i) &&
                data(i)<=data(i-1) && data(i)<data(i+1)) {
            result.push(i);
        }
    }
}

void line_segmentation_sort_x(intarray &segmentation) {
    CHECK_ARG(max(segmentation)<100000);
    narray<rectangle> bboxes;
    bounding_boxes(bboxes,segmentation);
    floatarray x0s;
    x0s.push(-999999);
    for(int i=1;i<bboxes.length();i++) {
        if(bboxes(i).empty()) {
            x0s.push(999999);
        } else {
            x0s.push(bboxes(i).x0);
        }
    }
    // dprint(x0s,1000); printf("\n");
    intarray permutation,rpermutation;
    quicksort(permutation,x0s);
    rpermutation.resize(permutation.length());
    for(int i=0;i<permutation.length();i++)
        rpermutation[permutation[i]] = i;
    // dprint(rpermutation,1000); printf("\n");
    for(int i=0;i<segmentation.length1d();i++) {
        if(segmentation.at1d(i)==0) continue;
        segmentation.at1d(i) = rpermutation(segmentation.at1d(i));
    }
}

void remove_small_components(intarray &segmentation,int r=5) {
    CHECK_ARG(max(segmentation)<100000);
    narray<rectangle> bboxes;
    bounding_boxes(bboxes,segmentation);
    for(int i=1;i<bboxes.length();i++) {
        rectangle b = bboxes(i);
        if(b.width()<r && b.height()<r) {
            for(int x=b.x0;x<b.x1;x++)
                for(int y=b.y0;y<b.y1;y++)
                    if(segmentation(x,y)==i)
                        segmentation(x,y) = 0;
        }
    }
}

void extract_holes(bytearray &holes,bytearray &binarized) {
    using namespace narray_ops;
    intarray temp;
    temp.copy(binarized);
    sub(255,temp);
    label_components(temp);
    int background = -1;
    for(int i=0;i<temp.dim(0);i++) {
        if(temp(i,0)!=0) {
            background = temp(i,0);
            break;
        }
    }
    makelike(holes,temp);
    holes = 0;
    CHECK(background>0);
    for(int i=0;i<temp.dim(0);i++) {
        for(int j=0;j<temp.dim(1);j++) {
            if(temp(i,j)>0 && temp(i,j)!=background)
                holes(i,j) = 255;
        }
    }
    dsection("segholes");
    dshow(holes,"y");
}


param_int seg_cuts_merge("seg_cuts_merge",10,"merge components smaller than this in seg-cuts");

void line_segmentation_merge_small_components(intarray &segmentation,int r=10) {
    CHECK_ARG(max(segmentation)<100000);
    make_line_segmentation_black(segmentation);
    narray<rectangle> bboxes;
    bounding_boxes(bboxes,segmentation);
    bboxes(0) = rectangle();
    bool changed;
    do {
        changed = false;
        for(int i=1;i<bboxes.length();i++) {
            rectangle b = bboxes(i);
            if(b.empty()) continue;
            if(b.width()>=r || b.height()>=r) continue;
#if 0
            // merge any small component with nearby big components
            floatarray dists;
            for(int j=0;j<bboxes.length();j++) {
                if(i==j || bboxes(j).empty())
                    dists.push(9999999);
                else
                    dists.push(fabs(bboxes(i).xcenter()-bboxes(j).xcenter()));
            }
            int closest = argmin(dists);
#else
            // merge small components only with touching components
            int closest = 0;
            rectangle b1 = b.grow(1);
            b1.intersect(rectangle(0,0,segmentation.dim(0),segmentation.dim(1)));
            for(int x=b1.x0;x<b1.x1;x++) {
                for(int y=b1.y0;y<b1.y1;y++) {
                    int value = segmentation(x,y);
                    if(value==0) continue;
                    if(value==i) continue;
                    closest = value;
                    break;
                }
            }
            if(closest==0) continue;
#endif
            for(int x=b.x0;x<b.x1;x++)
                for(int y=b.y0;y<b.y1;y++)
                    if(segmentation(x,y)==i)
                        segmentation(x,y) = closest;
            bboxes(i) = rectangle();
            changed = true;
        }
    } while(changed);
}

struct ICurvedCutSegmenter {
    int down_cost;
    int outside_diagonal_cost;
    int inside_diagonal_cost;
    int boundary_diagonal_cost;
    int inside_weight;
    int boundary_weight;
    int outside_weight;
    int min_range;
    int fill_holes;
    float min_thresh;
    //virtual void params_for_chars() = 0;
    virtual void params_for_lines() = 0;
    virtual void findAllCuts() = 0;
    virtual void findBestCuts() = 0;
    // virtual void relabel_image(bytearray &image) = 0;
    // virtual void relabel_image(intarray &image) = 0;
    virtual void setImage(bytearray &image) = 0;
    virtual ~ICurvedCutSegmenter() {}
};

struct CurvedCutSegmenterImpl : ICurvedCutSegmenter {
    // input
    intarray wimage;
    int where;

    // output
    intarray costs;
    intarray sources;
    int direction;
    int limit;

    intarray bestcuts;

    strg debug;
    intarray dimage;

    narray< narray <point> > cuts;
    floatarray cutcosts;

    CurvedCutSegmenterImpl() {
        //params_for_chars();
        params_for_lines();
        //params_from_hwrec_c();
    }

    void params_for_lines() {
        down_cost = 0;
        outside_diagonal_cost = 4;
        inside_diagonal_cost = 4;
        boundary_diagonal_cost = 0;
        outside_weight = 0;
        boundary_weight = -1;
        inside_weight = 4;
        min_range = 3;
        //min_thresh = -2.0;
        min_thresh = 10.0;
        fill_holes = 1;
    }

    // this function calculates the actual costs
    void step(int x0,int x1,int y) {
        int w = wimage.dim(0),h = wimage.dim(1);
        Queue<point> queue(w*h);
        for(int i=x0;i<x1;i++) queue.enqueue(point(i,y));
        int low = 1;
        int high = wimage.dim(0)-1;

        while(!queue.empty()) {
            point p = queue.dequeue();
            int i = p.x, j = p.y;
            int cost = costs(i,j);
            int ncost = cost+wimage(i,j)+down_cost;
            if(costs(i,j+direction)>ncost) {
                costs(i,j+direction) = ncost;
                sources(i,j+direction) = i;
                if(j+direction!=limit) queue.enqueue(point(i,j+direction));
            }
            if(i>low) {
                if(wimage(i,j)==0)
                    ncost = cost+wimage(i,j)+outside_diagonal_cost;
                else if(wimage(i,j)>0)
                    ncost = cost+wimage(i,j)+inside_diagonal_cost;
                else if(wimage(i,j)<0)
                    ncost = cost+wimage(i,j)+boundary_diagonal_cost;
                if(costs(i-1,j+direction)>ncost) {
                    costs(i-1,j+direction) = ncost;
                    sources(i-1,j+direction) = i;
                    if(j+direction!=limit) queue.enqueue(point(i-1,j+direction));
                }
            }
            if(i<high) {
                if(wimage(i,j)==0)
                    ncost = cost+wimage(i,j)+outside_diagonal_cost;
                else if(wimage(i,j)>0)
                    ncost = cost+wimage(i,j)+inside_diagonal_cost;
                else if(wimage(i,j)<0)
                    ncost = cost+wimage(i,j)+boundary_diagonal_cost;
                if(costs(i+1,j+direction)>ncost) {
                    costs(i+1,j+direction) = ncost;
                    sources(i+1,j+direction) = i;
                    if(j+direction!=limit) queue.enqueue(point(i+1,j+direction));
                }
            }
        }
    }

    void findAllCuts() {
        int w = wimage.dim(0), h = wimage.dim(1);
        // initialize dimensions of cuts, costs etc
        cuts.resize(w);
        cutcosts.resize(w);
        costs.resize(w,h);
        sources.resize(w,h);

        fill(costs, 1000000000);
        for(int i=0;i<w;i++) costs(i,0) = 0;
        fill(sources, -1);
        limit = where;
        direction = 1;
        step(0,w,0);

        for(int x=0;x<w;x++) {
            cutcosts(x) = costs(x,where);
            cuts(x).clear();
            // bottom should probably be initialized with 2*where instead of
            // h, because where cannot be assumed to be h/2. In the most extreme
            // case, the cut could go through 2 pixels in each row
            narray<point> bottom;
            int i = x, j = where;
            while(j>=0) {
                bottom.push(point(i,j));
                i = sources(i,j);
                j--;
            }
            //cuts(x).resize(h);
            for(i=bottom.length()-1;i>=0;i--) cuts(x).push(bottom(i));
        }

        fill(costs, 1000000000);
        for(int i=0;i<w;i++) costs(i,h-1) = 0;
        fill(sources, -1);
        limit = where;
        direction = -1;
        step(0,w,h-1);

        for(int x=0;x<w;x++) {
            cutcosts(x) += costs(x,where);
            // top should probably be initialized with 2*(h-where) instead of
            // h, because where cannot be assumed to be h/2. In the most extreme
            // case, the cut could go through 2 pixels in each row
            narray<point> top;
            int i = x, j = where;
            while(j<h) {
                if(j>where) top.push(point(i,j));
                i = sources(i,j);
                j++;
            }
            for(i=0;i<top.length();i++) cuts(x).push(top(i));
        }

        // add costs for line "where"
        for(int x=0;x<w;x++) {
            cutcosts(x) += wimage(x,where);
        }

    }

    void findBestCuts() {
        for(int i=0;i<cutcosts.length();i++) ext(dimage,i,int(cutcosts(i)+10)) = 0xff0000;
        for(int i=0;i<cutcosts.length();i++) ext(dimage,i,int(min_thresh+10)) = 0x800000;
        floatarray temp;
        gauss1d(temp,cutcosts,3.0);
        cutcosts.move(temp);
        local_minima(bestcuts,cutcosts,min_range,min_thresh);
        for(int i=0;i<bestcuts.length();i++) {
            narray<point> &cut = cuts(bestcuts(i));
            for(int j=0;j<cut.length();j++) {
                point p = cut(j);
                ext(dimage,p.x,p.y) = 0x00ff00;
            }
        }
        if(!debug.empty()) write_image_packed(debug,dimage);
        // dshow1d(cutcosts,"Y");
        //dshow(dimage,"Y");
    }

    void setImage(bytearray &image_) {
        bytearray image;
        image = image_;
        if(fill_holes) {
            bytearray holes;
            extract_holes(holes,image);
            for(int i=0;i<image.length();i++)
                if(holes[i]) image[i] = 255;
        }
        copy(dimage,image);
        int w = image.dim(0), h = image.dim(1);
        wimage.resize(w,h);
        fill(wimage, 0);
        float s1 = 0.0, sy = 0.0;
        for(int i=1;i<w;i++) for(int j=0;j<h;j++) {
            if(image(i,j)) { s1++; sy += j; }
            if(!image(i-1,j) && image(i,j)) wimage(i,j) = boundary_weight;
            else if(image(i,j)) wimage(i,j) = inside_weight;
            else wimage(i,j) = outside_weight;
        }
        where = int(sy/s1);
        for(int i=0;i<dimage.dim(0);i++) dimage(i,where) = 0x008000;
    }
};

class CurvedCutSegmenter : public ISegmentLine {
    public:
        autoref<CurvedCutSegmenterImpl> segmenter;
        int small_merge_threshold;

        CurvedCutSegmenter() {
            small_merge_threshold = 1;
        }

        const char *name() {
            return "curvedcut";
        }

        virtual const char *description() {
            return "curved cut segmenter";
        }

        virtual void set(const char *key,const char *value) {
            if(!strcmp(key,"debug"))
                segmenter->debug = value;
            else
                throw "unknown key";
        }

        virtual void set(const char *key,double value) {
            if(!strcmp(key,"fill_holes"))
                segmenter->fill_holes = (int)value;
            else if(!strcmp(key,"down_cost"))
                segmenter->down_cost = (int)value;
            else if(!strcmp(key,"small_merge_threshold"))
                small_merge_threshold = (int)value;
            else if(!strcmp(key,"outside_diagonal_cost"))
                segmenter->outside_diagonal_cost = (int)value;
            else if(!strcmp(key,"inside_diagonal_cost"))
                segmenter->inside_diagonal_cost = (int)value;
            else if(!strcmp(key,"boundary_diagonal_cost"))
                segmenter->boundary_diagonal_cost = (int)value;
            else if(!strcmp(key,"outside_weight"))
                segmenter->outside_weight = (int)value;
            else if(!strcmp(key,"boundary_weight"))
                segmenter->boundary_weight = (int)value;
            else if(!strcmp(key,"inside_weight"))
                segmenter->inside_weight = (int)value;
            else if(!strcmp(key,"min_range"))
                segmenter->min_range = (int)value;
            else if(!strcmp(key,"min_thresh"))
                segmenter->min_thresh = value;
            else
                throw "unknown key";
        }

        virtual void charseg(intarray &segmentation,bytearray &raw) {
            enum {PADDING = 3};
            optional_check_background_is_lighter(raw);
            bytearray image;
            image.copy(raw);
            binarize_simple(image);
            invert(image);

            segmenter->setImage(image);
            segmenter->findAllCuts();
            segmenter->findBestCuts();

            intarray seg;
            seg.copy(image);

            for(int r=0;r<segmenter->bestcuts.length();r++) {
                int w = seg.dim(0);
                int c = segmenter->bestcuts(r);
                narray<point> &cut = segmenter->cuts(c);
                for(int y=0;y<image.dim(1);y++) {
                    for(int i=-1;i<=1;i++) {
                        int x = cut(y).x;
                        if(x<1||x>=w-1) continue;
                        seg(x+i,y) = 0;
                    }
                }
            }
            label_components(seg);
            // dshowr(seg,"YY"); dwait();
            segmentation.copy(image);
            propagate_labels_to(segmentation,seg);

            line_segmentation_merge_small_components(segmentation,small_merge_threshold);
            line_segmentation_sort_x(segmentation);

            make_line_segmentation_white(segmentation);
            // set_line_number(segmentation, 1);
        }
};

ISegmentLine *make_CurvedCutSegmenter1() {
    return new CurvedCutSegmenter();
}

// FIXME/faisal Messy implementation--get rid of this. --tmb

class CurvedCutSegmenterToISegmentLineAdapter : public ISegmentLine {
    public:
        autoref<CurvedCutSegmenterImpl> segmenter;
        int small_merge_threshold;

        CurvedCutSegmenterToISegmentLineAdapter() {
            small_merge_threshold = 0;
        }

        const char *name() {
            return "curvedcut";
        }

        virtual const char *description() {
            return "curved cut segmenter";
        }

        virtual void set(const char *key,const char *value) {
            if(!strcmp(key,"debug"))
                segmenter->debug = value;
            else
                throw "unknown key";
        }

        virtual void set(const char *key,double value) {
            if(!strcmp(key,"down_cost"))
                segmenter->down_cost = (int)value;
            else if(!strcmp(key,"small_merge_threshold"))
                small_merge_threshold = (int)value;
            else if(!strcmp(key,"outside_diagonal_cost"))
                segmenter->outside_diagonal_cost = (int)value;
            else if(!strcmp(key,"inside_diagonal_cost"))
                segmenter->inside_diagonal_cost = (int)value;
            else if(!strcmp(key,"boundary_diagonal_cost"))
                segmenter->boundary_diagonal_cost = (int)value;
            else if(!strcmp(key,"outside_weight"))
                segmenter->outside_weight = (int)value;
            else if(!strcmp(key,"boundary_weight"))
                segmenter->boundary_weight = (int)value;
            else if(!strcmp(key,"inside_weight"))
                segmenter->inside_weight = (int)value;
            else if(!strcmp(key,"min_range"))
                segmenter->min_range = (int)value;
            else if(!strcmp(key,"min_thresh"))
                segmenter->min_thresh = value;
            else
                throw "unknown key";
        }

        virtual void charseg(intarray &result_segmentation,bytearray &orig_image) {
            enum {PADDING = 3};
            bytearray image;
            copy(image, orig_image);
            optional_check_background_is_lighter(image);
            binarize_simple(image);
            invert(image);
            pad_by(image, PADDING, PADDING);
            intarray segmentation;
            // pass image to segmenter
            segmenter->setImage(image);
            // find all cuts in the image
            segmenter->findAllCuts();
            // choose the best of all cuts
            segmenter->findBestCuts();

            segmentation.resize(image.dim(0),image.dim(1));
            for(int i=0;i<image.dim(0);i++) for(int j=0;j<image.dim(1);j++)
                segmentation(i,j) = image(i,j)?1:0;
            for(int r=0;r<segmenter->bestcuts.length();r++) {
                int c = segmenter->bestcuts(r);
                narray<point> &cut = segmenter->cuts(c);
                for(int y=0;y<image.dim(1);y++) {
                    for(int x=cut(y).x;x<image.dim(0);x++)
                        if(segmentation(x,y)) segmentation(x,y)++;
                }
            }
            extract_subimage(result_segmentation,segmentation,PADDING,PADDING,
                    segmentation.dim(0)-PADDING,segmentation.dim(1)-PADDING);

            if(small_merge_threshold>0) {
                line_segmentation_merge_small_components(result_segmentation,small_merge_threshold);
                line_segmentation_sort_x(result_segmentation);
            }

            make_line_segmentation_white(result_segmentation);
            // set_line_number(result_segmentation, 1);
        }
};



ISegmentLine *make_CurvedCutSegmenter() {
    return new CurvedCutSegmenterToISegmentLineAdapter();
}

struct CurvedCutSegmenterToISegmentLineAdapterWithCc:
    CurvedCutSegmenterToISegmentLineAdapter {
        virtual void charseg(intarray &result_segmentation,bytearray &orig_image) {
            bytearray image;
            copy(image, orig_image);
            optional_check_background_is_lighter(image);
            binarize_simple(image);
            invert(image);

            intarray ccseg;
            copy(ccseg, image);
            label_components(ccseg);

            CurvedCutSegmenterToISegmentLineAdapter::charseg(result_segmentation, orig_image);
            combine_segmentations(result_segmentation, ccseg);
        }
    };

ISegmentLine *make_CurvedCutWithCcSegmenter() {
    return new CurvedCutSegmenterToISegmentLineAdapterWithCc();
}


// FIXME this should really work "word"-wise, centered on each word,
// otherwise it does the wrong thing for non-deskewed lines
// (it worked "word"-wise in the original version)

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <ctype.h>

void fix_diacritics(intarray &segmentation) {
    narray<rectangle> bboxes;
    bounding_boxes(bboxes,segmentation);
    if(bboxes.length()<1) return;
    intarray assignments(bboxes.length());
    for(int i=0;i<assignments.length();i++)
        assignments(i) = i;
    for(int j=0;j<bboxes.length();j++) {
        float dist = 1e38;
        int closest = -1;
        for(int i=0;i<bboxes.length();i++) {
            // j should overlap i in the x direction
            if(bboxes[j].x1<bboxes[i].x0) continue;
            if(bboxes[j].x0>bboxes[i].x1) continue;
            // j should be above i
            if(!(bboxes[j].y0>=bboxes[i].y1)) continue;
#if 0
            // j should be smaller than i
            if(!(bboxes[j].area()<bboxes[i].area())) continue;
#endif
            float d = fabs(bboxes[j].xcenter()-bboxes[i].xcenter());
            if(d>=dist) continue;
            dist = d;
            closest = i;
        }
        if(closest<0) continue;
        assignments(j) = closest;
    }
    for(int i=0;i<segmentation.length();i++)
        segmentation[i] = assignments(segmentation[i]);
    renumber_labels(segmentation,1);
}

struct DpSegmenter : IDpSegmenter {
    // input
    floatarray wimage;
    int where;

    // output
    intarray costs;
    intarray sources;
    int direction;
    int limit;

    intarray bestcuts;

    strbuf debug;

    narray< narray <point> > cuts;
    floatarray cutcosts;

    DpSegmenter() {
        pdef("down_cost",0,"cost of down step");
        pdef("outside_diagonal_cost",1,"cost of outside diagonal step to the left");
        pdef("outside_diagonal_cost_r",1,"cost of outside diagonal step to the right");
        pdef("inside_diagonal_cost",4,"cost of inside diagonal step");
        pdef("outside_weight",0,"cost of outside pixel");
        pdef("inside_weight",8,"cost of inside pixel");
        pdef("cost_smooth",2.0,"smoothing parameter for costs");
        pdef("min_range",1,"min range value");
        pdef("min_thresh",80.0,"min threshold value");
        pdef("component_segmentation",1,"also perform connected component segmentation");
        pdef("fix_diacritics",1,"group dots above characters back with those characters");
        pdef("fill_holes",1,"fill holes prior to dp segmentation (for cases like oo)");
        pdef("debug","none","debug output file");
    }
    const char *name() {
        return "dpseg";
    }
    void setParams() {
        down_cost = pgetf("down_cost");
        outside_weight = pgetf("outside_weight");
        inside_weight = pgetf("inside_weight");
        outside_diagonal_cost = pgetf("outside_diagonal_cost");
        outside_diagonal_cost_r = pgetf("outside_diagonal_cost_r");
        inside_diagonal_cost = pgetf("inside_diagonal_cost");
        cost_smooth = pgetf("cost_smooth");
        min_range = pgetf("min_range");
        min_thresh = pgetf("min_thresh");
        if(strcmp(pget("debug"),"none")) debug = pget("debug");
    }

    // this function calculates the actual costs
    void step(int x0,int x1,int y) {
        int w = wimage.dim(0),h = wimage.dim(1);
        Queue<point> queue(w*h);
        for(int i=x0;i<x1;i++) queue.enqueue(point(i,y));
        int low = 1;
        int high = wimage.dim(0)-1;

        while(!queue.empty()) {
            point p = queue.dequeue();
            int i = p.x, j = p.y;
            int cost = costs(i,j);
            int ncost = cost+wimage(i,j)+down_cost;
            if(costs(i,j+direction)>ncost) {
                costs(i,j+direction) = ncost;
                sources(i,j+direction) = i;
                if(j+direction!=limit) queue.enqueue(point(i,j+direction));
            }
            if(i>low) {
                if(wimage(i,j)==0)
                    ncost = cost+wimage(i,j)+outside_diagonal_cost;
                else
                    ncost = cost+wimage(i,j)+inside_diagonal_cost;
                if(costs(i-1,j+direction)>ncost) {
                    costs(i-1,j+direction) = ncost;
                    sources(i-1,j+direction) = i;
                    if(j+direction!=limit) queue.enqueue(point(i-1,j+direction));
                }
            }
            if(i<high) {
                if(wimage(i,j)==0)
                    ncost = cost+wimage(i,j)+outside_diagonal_cost_r;
                else
                    ncost = cost+wimage(i,j)+inside_diagonal_cost;
                if(costs(i+1,j+direction)>ncost) {
                    costs(i+1,j+direction) = ncost;
                    sources(i+1,j+direction) = i;
                    if(j+direction!=limit) queue.enqueue(point(i+1,j+direction));
                }
            }
        }
    }

    void findAllCuts() {
        int w = wimage.dim(0), h = wimage.dim(1);
        // initialize dimensions of cuts, costs etc
        cuts.resize(w);
        cutcosts.resize(w);
        costs.resize(w,h);
        sources.resize(w,h);

        fill(costs, 1000000000);
        for(int i=0;i<w;i++) costs(i,0) = 0;
        fill(sources, -1);
        limit = where;
        direction = 1;
        step(0,w,0);

        for(int x=0;x<w;x++) {
            cutcosts(x) = costs(x,where);
            cuts(x).clear();
            // bottom should probably be initialized with 2*where instead of
            // h, because where cannot be assumed to be h/2. In the most extreme
            // case, the cut could go through 2 pixels in each row
            narray<point> bottom;
            int i = x, j = where;
            while(j>=0) {
                bottom.push(point(i,j));
                i = sources(i,j);
                j--;
            }
            //cuts(x).resize(h);
            for(i=bottom.length()-1;i>=0;i--) cuts(x).push(bottom(i));
        }

        fill(costs, 1000000000);
        for(int i=0;i<w;i++) costs(i,h-1) = 0;
        fill(sources, -1);
        limit = where;
        direction = -1;
        step(0,w,h-1);

        for(int x=0;x<w;x++) {
            cutcosts(x) += costs(x,where);
            // top should probably be initialized with 2*(h-where) instead of
            // h, because where cannot be assumed to be h/2. In the most extreme
            // case, the cut could go through 2 pixels in each row
            narray<point> top;
            int i = x, j = where;
            while(j<h) {
                if(j>where) top.push(point(i,j));
                i = sources(i,j);
                j++;
            }
            for(i=0;i<top.length();i++) cuts(x).push(top(i));
        }

        // add costs for line "where"
        for(int x=0;x<w;x++) {
            cutcosts(x) += wimage(x,where);
        }

    }

    void findBestCuts() {
        for(int i=0;i<cutcosts.length();i++) ext(dimage,i,int(cutcosts(i)+10)) = 0xff0000;
        for(int i=0;i<cutcosts.length();i++) ext(dimage,i,int(min_thresh+10)) = 0x800000;
        floatarray temp;
        gauss1d(temp,cutcosts,cost_smooth);
        cutcosts.move(temp);
        local_minima(bestcuts,cutcosts,min_range,min_thresh);
        for(int i=0;i<bestcuts.length();i++) {
            narray<point> &cut = cuts(bestcuts(i));
            for(int j=0;j<cut.length();j++) {
                point p = cut(j);
                ext(dimage,p.x,p.y) = 0x00ff00;
            }
        }
        if(debug) write_image_packed(debug,dimage);
        //dshow1d(cutcosts,"Y");
        //dshow(dimage,"Y");
    }

    void setImage(bytearray &image_) {
        bytearray image;
        image = image_;
        copy(dimage,image);
        if(pgetf("fill_holes")) {
            bytearray holes;
            extract_holes(holes,image);
            for(int i=0;i<image.length();i++)
                if(holes[i]) image[i] = 255;
        }
        dsection("segholes");
        dshow(image,"y");
        int w = image.dim(0), h = image.dim(1);
        wimage.resize(w,h);
        fill(wimage, 0);
        float s1 = 0.0, sy = 0.0;
        for(int i=1;i<w;i++) {
            for(int j=0;j<h;j++) {
                if(image(i,j)) { s1++; sy += j; }
                if(image(i,j)) wimage(i,j) = inside_weight;
                else wimage(i,j) = outside_weight;
            }
        }
        if(s1==0) where = image.dim(1)/2;
        else where = int(sy/s1);
        for(int i=0;i<dimage.dim(0);i++) dimage(i,where) = 0x008000;
    }

    // ISegmentLine methods

    virtual const char *description() {
        return "curved cut segmenter";
    }

    virtual void charseg(intarray &segmentation,bytearray &raw) {
        setParams();
        enum {PADDING = 3};
        optional_check_background_is_lighter(raw);
        bytearray image;
        image.copy(raw);
        binarize_simple(image);
        invert(image);

        setImage(image);
        findAllCuts();
        findBestCuts();

        intarray seg;
        seg.copy(image);

#ifndef PROPAGATE
        seg = 255;
#endif

        for(int r=0;r<bestcuts.length();r++) {
            int w = seg.dim(0);
            int c = bestcuts(r);
            narray<point> &cut = cuts(c);
            for(int y=0;y<image.dim(1);y++) {
                for(int i=-1;i<=1;i++) {
                    int x = cut(y).x;
                    if(x<1||x>=w-1) continue;
                    seg(x+i,y) = 0;
                }
            }
        }
        label_components(seg);
        // dshowr(seg,"YY"); dwait();
        segmentation.copy(image);

        for(int i=0;i<seg.length();i++)
            if(!segmentation[i]) seg[i] = 0;

        propagate_labels_to(segmentation,seg);

        if(pgetf("component_segmentation")) {
            intarray ccseg;
            ccseg.copy(image);
            label_components(ccseg);
            combine_segmentations(segmentation,ccseg);
            if(pgetf("fix_diacritics")) {
                fix_diacritics(segmentation);
            }
        }

#if 0
        line_segmentation_merge_small_components(segmentation,small_merge_threshold);
        line_segmentation_sort_x(segmentation);
#endif

        make_line_segmentation_white(segmentation);
        // set_line_number(segmentation, 1);
    }
};

IDpSegmenter *make_DpSegmenter() {
    return new DpSegmenter();
}
