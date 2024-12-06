/*
 *  Copyright 2024 <me>
 */

#ifndef _g_starter_canvas_h_
#define _g_starter_canvas_h_

#include "include/GCanvas.h"
#include "include/GRect.h"
#include "include/GColor.h"
#include "include/GBitmap.h"
#include "Edge.h"
#include <stack>

class MyCanvas : public GCanvas {
public:
    MyCanvas(const GBitmap& device) : fDevice(device) {
        fCTM.push_back(*new GMatrix());
    }
    GPixel color2Pix(const GColor& color);
    Edge makeEdge(const GPoint& point1, const GPoint& point2);
    void clear(const GColor& color) override;
    void drawRect(const GRect& rect, const GPaint& paint) override;
    void drawConvexPolygon(const GPoint* points, int count, const GPaint& paint) override;
    void save() override;
    void restore() override;
    void concat(const GMatrix& matrix) override;
    void drawPath(const GPath&, const GPaint&) override;
    void blit(int L, int R, int y, GPaint paint);
    int clippedX(int x);
    int clippedY(int y);
    void drawMesh(const GPoint verts[], const GColor colors[], const GPoint texs[],
                  int count, const int indices[], const GPaint&) override;
    void drawQuad(const GPoint verts[4], const GColor colors[4], const GPoint texs[4],
                  int level, const GPaint&) override;
    void quadDivider(const GPoint verts[4], const GColor colors[4], const GPoint texs[4],
                  int level, const GPaint&);
    GPoint getMidpoint(GPoint a, GPoint b);

    const GMatrix& ctm() const { return fCTM.back(); }
    static inline bool rect_stays_rect(const GMatrix& m) {
        return (m[0] > 0) && (m[3] > 0) && (m[1] == 0) && (m[2] == 0);
    }

private:
    // Note: we store a copy of the bitmap
    const GBitmap fDevice;

    // Add whatever other fields you need
    std::vector<GMatrix> fCTM;
};

#endif