#include "include/GShader.h"
#include "include/GBitmap.h"
#include "include/GMatrix.h"
#include "include/GPixel.h"
#include "include/GPoint.h"
#include "include/GColor.h"
#include <iostream>
#include <typeinfo>
#include "stdlib.h"

class MyShader:public GShader {
public:
    MyShader(const GBitmap& bitmap, const GMatrix& localMatrix, GTileMode tMode) : fBitmap(bitmap), mtx(localMatrix), tileMode(tMode) {}

    bool isOpaque() override {
        // NOT USED
        
        for (int y=0; y<fBitmap.height(); y++) {
            for (int x=0; x<fBitmap.width(); x++) {
                GPixel pix = *fBitmap.getAddr(x, y);
                if (GPixel_GetA(pix) != 1) {
                    return false;
                }
            }
        }
        return true;
    }

    bool setContext(const GMatrix& ctm) override {
        if (auto inverted = (ctm * mtx).invert()) {
            fInverse = *inverted;
            return true;
        }
        return false;
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

            if (tileMode == GTileMode::kClamp) {
                finX = std::max(0, std::min(srcX, fBitmap.width()-1));
                finY = std::max(0, std::min(srcY, fBitmap.height()-1));
            }
            else if (tileMode == GTileMode::kRepeat) {
                finX = srcX;
                finY = srcY;
            
                while (finX > fBitmap.width()-1) {
                    finX -= fBitmap.width();
                }
                while (finX < 0) {
                    finX += fBitmap.width();
                }
                while (finY > fBitmap.height()-1) {
                    finY -= fBitmap.height();
                }
                while (finY < 0) {
                    finY += fBitmap.height();                    
                }
            }
            else if (tileMode == GTileMode::kMirror) {
                finX = srcX;
                finY = srcY;
                int xCounter = 0;
                int yCounter = 0;
            
                while (finX > fBitmap.width()-1) {
                    finX -= fBitmap.width();
                    xCounter++;
                }
                while (finX < 0) {
                    finX += fBitmap.width();
                    xCounter++;
                }
                while (finY > fBitmap.height()-1) {
                    finY -= fBitmap.height();
                    yCounter++;
                }
                while (finY < 0) {
                    finY += fBitmap.height();
                    yCounter++;
                }

                if (xCounter % 2 == 1) {
                    finX = fBitmap.width()-1 - finX;
                }
                if (yCounter % 2 == 1) {
                    finY = fBitmap.height()-1 - finY;
                }
            }

            out[i] = *fBitmap.getAddr(finX, finY);
        }
    }

private:
    GBitmap fBitmap;
    GMatrix fInverse;
    GMatrix mtx;
    GTileMode tileMode;
};

std::shared_ptr<GShader> GCreateBitmapShader(const GBitmap& bitmap, const GMatrix& localMatrix, GTileMode tMode) {
    return std::shared_ptr<GShader>(new MyShader(bitmap, localMatrix, tMode));
}