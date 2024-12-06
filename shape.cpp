#include "include/GShader.h"
#include "include/GBitmap.h"
#include "include/GMatrix.h"
#include "include/GPixel.h"
#include "include/GPoint.h"
#include "include/GColor.h"
#include <math.h>
#include <iostream>

class Voronoi : public GShader {
    public:
        Voronoi(const GPoint pts[], const GColor colors[], int count) {
            for(int i=0; i<count; i++){
                this->colors.push_back(colors[i]);
                this->points.push_back(pts[i]);
            }
            this->count = count;
        }
        bool setContext(const GMatrix& ctm) override {
        if (auto inverted = (ctm * mtx).invert()) {
            fInverse = *inverted;
            return true;
        }
        return false;
    }
        bool isOpaque() override {
        }
        void shadeRow(int x, int y, int count, GPixel out[]) override {
             GPoint dstPix;
        GPoint srcPix;
        int srcX;
        int srcY;
        int finX;
        int finY;

        for (int i = 0; i < count; i++) {

            dstPix.x = x + i + (0.5f);
            dstPix.y = y + (0.5f);

            srcPix = fInverse * dstPix;

            srcX = GFloorToInt(srcPix.x);
            srcY = GFloorToInt(srcPix.y);

            finX = std::max(0, std::min(srcX, fBitmap.width()-1));
            finY = std::max(0, std::min(srcY, fBitmap.height()-1));
            
            GColor c = sdc(srcPix.x, srcPix.y);
            out[i] = color2Pix(c);
        }
        }

        GColor sdc(float x, float y){
            float min = INFINITY;
            int min_i = 0;
            for(int i = 0; i<this->count; i++){
                float euclid = sqrt(pow((this->points[i].x - x), 2) + pow((this->points[i].y-y), 2));
                if(euclid < min){
                    min = euclid;
                    min_i = i;
                }
            }
            return this->colors[min_i];
        }

        GPixel color2Pix(const GColor& color) {

            unsigned int newA = GRoundToInt(color.a * 255);
            unsigned int newR = GRoundToInt(color.r * color.a * 255);
            unsigned int newG = GRoundToInt(color.g * color.a * 255);
            unsigned int newB = GRoundToInt(color.b * color.a * 255);

            return GPixel_PackARGB(newA, newR, newG, newB);
        }


    private:
    std::vector<GColor> colors;
        GBitmap fBitmap;
        GMatrix fInverse;
        GTileMode tileMode;
        GMatrix mtx;
    std::vector<GPoint> points;
        int count;
};