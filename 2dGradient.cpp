#include "include/GShader.h"
#include "include/GBitmap.h"
#include "include/GMatrix.h"
#include "include/GPixel.h"
#include "include/GPoint.h"
#include "include/GColor.h"
#include "2DGradient.h"
#include <iostream>

TwoDGradient::TwoDGradient(GPoint p0, GPoint p1, GPoint p2, GColor color00, GColor color01, GColor color02) {
    U = p1-p0;
    V = p2-p0;

    color0 = color00;
    color1 = color01;
    color2 = color02;

    DC1 = color01 - color00;
    DC2 = color02 - color00;
    
    unit2UV = GMatrix(U.x, V.x, p0.x, U.y, V.y, p0.y);

    e0.x = 0;
    e0.y = 1;
    e1.x = 1;
    e1.y = 0;

    // std::cout << (unit2UV*e0 == U) << std::endl;
    // std::cout << (unit2UV*e0 == V) << std::endl;
    // std::cout << (unit2UV*e1 == U) << std::endl;
    // std::cout << (unit2UV*e1 == V) << std::endl;
}

bool TwoDGradient::setContext(const GMatrix& ctm) {
    if (auto inv = (ctm*unit2UV).invert()) {
        cart2Bary = *inv;
        return true;
    }
    // else if (auto inv = cart2Bary.invert()) {
    //     cart2Bary = *inv;
    // }
    return false;
}

bool TwoDGradient::isOpaque() {}

void TwoDGradient::shadeRow(int x, int y, int count, GPixel out[]) {
    for (int i = 0; i < count; ++i) {
        GPoint P;
        P.x = x + 0.5;
        P.y = y + 0.5;
        GPoint tP = cart2Bary*P;

        GColor C = (tP.x * color1) + (tP.y * color2) + ((1-tP.x-tP.y) * color0);
        if (C.a > 1) {
            C.a = 1;
        }
        if (C.r > 1) {
            C.r = 1;
        }
        if (C.g > 1) {
            C.g = 1;
        }
        if (C.b > 1) {
            C.b = 1;
        }
        if (C.a < 0) {
            C.a = 0;
        }
        if (C.r < 0) {
            C.r = 0;
        }
        if (C.g < 0) {
            C.g = 0;
        }
        if (C.b < 0) {
            C.b = 0;
        }
        out[i] = color2Pix(C);
        x += 1;
    }
}

GPixel TwoDGradient::color2Pix(const GColor& color) {

    unsigned int newA = GRoundToInt(color.a * 255);
    unsigned int newR = GRoundToInt(color.r * color.a * 255);
    unsigned int newG = GRoundToInt(color.g * color.a * 255);
    unsigned int newB = GRoundToInt(color.b * color.a * 255);

    return GPixel_PackARGB(newA, newR, newG, newB);
}

// std::shared_ptr<GShader> GCreateTwoDGradient(GPoint p0, GPoint p1, GPoint p2, GColor color00, GColor color1, GColor color2) {
//     return std::unique_ptr<GShader>(new TwoDGradient(p0, p1, p2, color00, color1, color2));
// }