struct TextLineRAST {
    TextLineRAST() {}
    virtual void setParams(double max_slope,double ymin,double ymax) = 0;
    virtual void setMaxLines(int i) = 0;
    virtual void addChar(int x0,int y0,int x1,int y1) = 0;
    virtual void addColumn(int x0,int y0,int x1,int y1) = 0;
    virtual void compute() = 0;
    virtual int nlines() = 0;
    virtual void getLine(int i,float &c,float &m,float &d) = 0;
    virtual float getQuality(int i) = 0;
    float getLine_c(int i) { float c,m,d; getLine(i,c,m,d); return c; }
    float getLine_m(int i) { float c,m,d; getLine(i,c,m,d); return m; }
    float getLine_d(int i) { float c,m,d; getLine(i,c,m,d); return d; }
};

TextLineRAST *make_TextLineRAST();

struct TextLineRAST2 {
    TextLineRAST2() {}
    virtual void setParams(double max_slope,double ymin,double ymax) = 0;
    virtual void setMaxLines(int i) = 0;
    virtual void addChar(int x0,int y0,int x1,int y1) = 0;
    virtual void addColumn(int x0,int y0,int x1,int y1) = 0;
    virtual void compute() = 0;
    virtual int nlines() = 0;
    virtual void getLine(int i,float &c,float &m,float &d) = 0;
    virtual void getBbox(int i,int &x0,int &y0,int &x1,int &y1) = 0;
    virtual float getQuality(int i) = 0;
    float getLine_c(int i) { float c,m,d; getLine(i,c,m,d); return c; }
    float getLine_m(int i) { float c,m,d; getLine(i,c,m,d); return m; }
    float getLine_d(int i) { float c,m,d; getLine(i,c,m,d); return d; }
};

TextLineRAST2 *make_TextLineRAST2();

