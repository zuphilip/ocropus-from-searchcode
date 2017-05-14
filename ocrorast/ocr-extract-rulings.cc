#include "ocrorast.h"

using namespace colib;

void ExtractRulings::analyzeObstacles(rectarray &hor_rulings,
        rectarray &vert_rulings,
        rectarray &images,
        rectarray &obstacles,
        int xheight){

    int hor_ruling_maxheight = 2*xheight;
    int vert_ruling_maxwidth = 2*xheight;
    for(int i=0, l=obstacles.length(); i<l; i++){
        int width = obstacles[i].width();
        int height = obstacles[i].height();
        if(width<vert_ruling_maxwidth && height<hor_ruling_maxheight)
            images.push(obstacles[i]);
        else if(width<vert_ruling_maxwidth && height>hor_ruling_maxheight)
            vert_rulings.push(obstacles[i]);
        else if(width>vert_ruling_maxwidth && height<hor_ruling_maxheight)
            hor_rulings.push(obstacles[i]);
        else
            images.push(obstacles[i]);
    }
}


ExtractRulings *make_ExtractRulings(){
    return new ExtractRulings();
}
