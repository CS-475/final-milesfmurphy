#include "include/GShader.h"
#include "include/GBitmap.h"
#include "include/GMatrix.h"
#include "include/GPixel.h"
#include "include/GPoint.h"
#include "include/GColor.h"
#include <iostream>

class MyGradient : public GShader {
    public:
        MyGradient(GPoint p0, GPoint p1, const GColor colors[], int count, GTileMode tMode) {

            for (int i=0; i<count; i++) {
                flatGradient.push_back(colors[i]);
            }

            gradMat = GMatrix(p1.x-p0.x, -(p1.y-p0.y), p0.x,
                              p1.y-p0.y, p1.x-p0.x, p0.y);

            tileMode = tMode;
        }

        bool setContext(const GMatrix& ctm) override {
            if (auto inverted = (ctm * gradMat).invert()) {
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
            int index;
            float t;
            GColor c;

            for (int i = 0; i < count; i++) {

                dstPix.x = x + i + (0.5f);
                dstPix.y = y + (0.5f);

                srcPix = fInverse * dstPix;

                if (tileMode == GTileMode::kClamp) {
                    srcPix.x = (srcPix.x > 1) ? 1 : srcPix.x;
                    srcPix.x = (srcPix.x < 0) ? 0 : srcPix.x;
                }
                else if (tileMode == GTileMode::kRepeat) {
                    while (srcPix.x > 1) {
                        srcPix.x -= 1;
                    }
                    while (srcPix.x < 0) {
                        srcPix.x += 1;
                    }
                }
                else if (tileMode == GTileMode::kMirror) {
                    int counter = 0;
                    while (srcPix.x > 1) {
                        srcPix.x -= 1;
                        counter++;
                    }
                    while (srcPix.x < 0) {
                        srcPix.x += 1;
                        counter++;
                    }
                    if (counter % 2 == 1) {
                        srcPix.x = 1 - srcPix.x;
                    }
                }

                index = floor(srcPix.x * (flatGradient.size()-1));
                t = (srcPix.x * (flatGradient.size()-1)) - index;
                c = (1 - t)*flatGradient[index] + t*flatGradient[index+1];
                out[i] = color2Pix(c);
            }
        }

        GPixel color2Pix(const GColor& color) {

            unsigned int newA = GRoundToInt(color.a * 255);
            unsigned int newR = GRoundToInt(color.r * color.a * 255);
            unsigned int newG = GRoundToInt(color.g * color.a * 255);
            unsigned int newB = GRoundToInt(color.b * color.a * 255);

            return GPixel_PackARGB(newA, newR, newG, newB);
        }


    private:
        std::vector<GColor> flatGradient;
        GMatrix gradMat;
        GMatrix fInverse;
        GTileMode tileMode;
};

std::shared_ptr<GShader> GCreateLinearGradient(GPoint p0, GPoint p1, const GColor colors[], int count, GTileMode tMode) {
    if (count < 1) {
        return nullptr;
    }

    return std::unique_ptr<GShader>(new MyGradient(p0, p1, colors, count, tMode));
}