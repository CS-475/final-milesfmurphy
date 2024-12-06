#include "include/GShader.h"
#include "include/GBitmap.h"
#include "include/GMatrix.h"
#include "include/GPixel.h"
#include "include/GPoint.h"
#include "include/GColor.h"
#include <iostream>

class TwoDGradient : public GShader {
    public:
        TwoDGradient(GPoint p0, GPoint p1, GPoint p2, GColor color00, GColor color1, GColor color2);
        bool setContext(const GMatrix& ctm) override;
        bool isOpaque() override;
        void shadeRow(int x, int y, int count, GPixel out[]) override;
        GPixel color2Pix(const GColor& color);

    private:
        GMatrix cart2Bary;
        GMatrix unit2UV;
        GColor color0, color1, color2;
        GColor DC1, DC2;
        GPoint U, V;
        GPoint e0, e1;
};

static inline std::shared_ptr<GShader> GCreateTwoDGradient(GPoint p0, GPoint p1, GPoint p2, GColor color00, GColor color1, GColor color2) {
    return std::unique_ptr<GShader>(new TwoDGradient(p0, p1, p2, color00, color1, color2));
}