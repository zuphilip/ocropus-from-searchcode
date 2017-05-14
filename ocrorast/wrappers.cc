#include "colib/colib.h"
#include "colib/iarith.h"
#include "iulib/imggauss.h"
#include "iulib/imgmisc.h"
#include "ocrorast.h"

using namespace colib;
using namespace iulib;

#include "wrappers.h"

struct TextLineRAST_Impl : TextLineRAST {
    CTextlineRASTBasic comp;
    rectarray chars;
    rectarray columns;
    narray<TextLineParam> textlines;
    TextLineRAST_Impl() {
        comp.setDefaultParameters();
    }
    float getQuality(int i) {
        return comp.getQuality(i);
    }
    void setParams(double max_slope,double ymin,double ymax) {
        comp.setMaxSlope(max_slope);
        comp.setMaxYintercept(ymin,ymax);
    }
    void setMaxLines(int i) {
        comp.max_results = i;
    }
    void addChar(int x0,int y0,int x1,int y1) {
        chars.push(rectangle(x0,y0,x1,y1));
    }
    void addColumn(int x0,int y0,int x1,int y1) {
        columns.push(rectangle(x0,y0,x1,y1));
    }
    void compute() {
        autodel<CharStats> charstats;
        charstats = new CharStats();
        charstats->setCharBoxes(chars);
        textlines.clear();
        comp.extract(textlines,columns,charstats);
    }
    int nlines() {
        return textlines.length();
    }
    void getLine(int i,float &c,float &m,float &d) {
        c = textlines[i].c;
        m = textlines[i].m;
        d = textlines[i].d;
    }
};

TextLineRAST *make_TextLineRAST() {
    return new TextLineRAST_Impl();
}

struct TextLineRAST2_Impl : TextLineRAST2 {
    CTextlineRAST comp;
    rectarray chars;
    rectarray columns;
    narray<TextLine> textlines;
    TextLineRAST2_Impl() {
        comp.setDefaultParameters();
    }
    float getQuality(int i) {
        return comp.getQuality(i);
    }
    void setMaxLines(int i) {
        comp.max_results = i;
    }
    void setParams(double max_slope,double ymin,double ymax) {
        comp.setMaxSlope(max_slope);
        comp.setMaxYintercept(ymin,ymax);
    }
    void addChar(int x0,int y0,int x1,int y1) {
        chars.push(rectangle(x0,y0,x1,y1));
    }
    void addColumn(int x0,int y0,int x1,int y1) {
        columns.push(rectangle(x0,y0,x1,y1));
    }
    void compute() {
        autodel<CharStats> charstats;
        charstats = new CharStats();
        charstats->setCharBoxes(chars);
        textlines.clear();
        comp.extract(textlines,columns,charstats);
    }
    int nlines() {
        return textlines.length();
    }
    void getLine(int i,float &c,float &m,float &d) {
        c = textlines[i].c;
        m = textlines[i].m;
        d = textlines[i].d;
    }
    void getBbox(int i,int &x0,int &y0,int &x1,int &y1) {
        rectangle r = textlines[i].bbox;
        x0 = r.x0;
        y0 = r.y0;
        x1 = r.x1;
        y1 = r.y1;
    }
};

TextLineRAST2 *make_TextLineRAST2() {
    return new TextLineRAST2_Impl();
}
