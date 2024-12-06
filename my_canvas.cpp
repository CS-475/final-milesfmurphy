/*
 *  Copyright 2024 <me>
 */

#include "my_canvas.h"
#include <list>
#include <iostream>
#include "include/GBlendMode.h"
#include "include/GShader.h"
#include "include/GPaint.h"
#include "include/GPath.h"
#include "2DGradient.h"
#include "ProxyShader.h"


GPixel MyCanvas::color2Pix(const GColor& color) {

    unsigned int newA = GRoundToInt(color.a * 255);
    unsigned int newR = GRoundToInt(color.r * color.a * 255);
    unsigned int newG = GRoundToInt(color.g * color.a * 255);
    unsigned int newB = GRoundToInt(color.b * color.a * 255);

    return GPixel_PackARGB(newA, newR, newG, newB);
}

int MyCanvas::clippedX(int x) {
    if (x < 0) {
        return 0;
    }
    if (x > fDevice.width()) {
        return fDevice.width();
    }
    return x;
}

int MyCanvas::clippedY(int y) {
    if (y < 0) {
        return 0;
    }
    if (y > fDevice.height()) {
        return fDevice.height();
    }
    return y;
}

void MyCanvas::save() {
    fCTM.push_back(fCTM.back());
}
void MyCanvas::restore() {
    if (!fCTM.empty()) {
        fCTM.pop_back();
    }
}
void MyCanvas::concat(const GMatrix& matrix) {
    fCTM.back() = fCTM.back() * matrix;
}

GPixel kClear(GPixel src, GPixel dst){
    return GPixel_PackARGB(1, 0, 0, 0);
}

GPixel kSrc(GPixel src, GPixel dst){
    return src;
}

GPixel kDst(GPixel src, GPixel dst){
    return dst;
}

GPixel kSrcOver(GPixel src, GPixel dst){
    float coeff = (255.0f - GPixel_GetA(src))/255.0f;
    unsigned R = GRoundToInt(GPixel_GetR(src) + coeff * GPixel_GetR(dst));
    unsigned G = GRoundToInt(GPixel_GetG(src) + coeff * GPixel_GetG(dst));
    unsigned B = GRoundToInt(GPixel_GetB(src) + coeff * GPixel_GetB(dst));
    unsigned A = GRoundToInt(GPixel_GetA(src) + coeff * GPixel_GetA(dst));
    return GPixel_PackARGB(A, R, G, B);
}

GPixel kDstOver(GPixel src, GPixel dst){
    float coeff = (255.0f - GPixel_GetA(dst))/255.0f;
    unsigned R = GRoundToInt(GPixel_GetR(dst) + coeff * GPixel_GetR(src));
    unsigned G = GRoundToInt(GPixel_GetG(dst) + coeff * GPixel_GetG(src));
    unsigned B = GRoundToInt(GPixel_GetB(dst) + coeff * GPixel_GetB(src));
    unsigned A = GRoundToInt(GPixel_GetA(dst) + coeff * GPixel_GetA(src));
    return GPixel_PackARGB(A, R, G, B);
}

GPixel kSrcIn(GPixel src, GPixel dst){
    float coeff = GPixel_GetA(dst)/255.0f;
    unsigned R = GRoundToInt(coeff * GPixel_GetR(src));
    unsigned G = GRoundToInt(coeff * GPixel_GetG(src));
    unsigned B = GRoundToInt(coeff * GPixel_GetB(src));
    unsigned A = GRoundToInt(coeff * GPixel_GetA(src));
    return GPixel_PackARGB(A, R, G, B);
}

GPixel kDstIn(GPixel src, GPixel dst){
    float coeff = GPixel_GetA(src)/255.0f;
    unsigned R = GRoundToInt(coeff * GPixel_GetR(dst));
    unsigned G = GRoundToInt(coeff * GPixel_GetG(dst));
    unsigned B = GRoundToInt(coeff * GPixel_GetB(dst));
    unsigned A = GRoundToInt(coeff * GPixel_GetA(dst));
    return GPixel_PackARGB(A, R, G, B);
}

GPixel kSrcOut(GPixel src, GPixel dst){
    float coeff = (255.0f - GPixel_GetA(dst))/255.0f;
    unsigned R = GRoundToInt(coeff * GPixel_GetR(src));
    unsigned G = GRoundToInt(coeff * GPixel_GetG(src));
    unsigned B = GRoundToInt(coeff * GPixel_GetB(src));
    unsigned A = GRoundToInt(coeff * GPixel_GetA(src));
    return GPixel_PackARGB(A, R, G, B);
}

GPixel kDstOut(GPixel src, GPixel dst){
    float coeff = (255.0f - GPixel_GetA(src))/255.0f;
    unsigned R = GRoundToInt(coeff * GPixel_GetR(dst));
    unsigned G = GRoundToInt(coeff * GPixel_GetG(dst));
    unsigned B = GRoundToInt(coeff * GPixel_GetB(dst));
    unsigned A = GRoundToInt(coeff * GPixel_GetA(dst));
    return GPixel_PackARGB(A, R, G, B);
}

GPixel kSrcATop(GPixel src, GPixel dst){
    GPixel part1 = kSrcIn(src, dst);
    GPixel part2 = kDstOut(src, dst);
    GPixel r = GPixel_PackARGB(GPixel_GetA(part1) + GPixel_GetA(part2), 
        GPixel_GetR(part1) + GPixel_GetR(part2), GPixel_GetG(part1) + GPixel_GetG(part2),
        GPixel_GetB(part1) + GPixel_GetB(part2));
    return r;
}

GPixel kDstATop(GPixel src, GPixel dst){
    GPixel part1 = kDstIn(src, dst);
    GPixel part2 = kSrcOut(src, dst);
    GPixel r = GPixel_PackARGB(GPixel_GetA(part1) + GPixel_GetA(part2), 
        GPixel_GetR(part1) + GPixel_GetR(part2), GPixel_GetG(part1) + GPixel_GetG(part2),
        GPixel_GetB(part1) + GPixel_GetB(part2));
    return r;
}

GPixel kXor(GPixel src, GPixel dst){
    GPixel part1 = kDstOut(src, dst);
    GPixel part2 = kSrcOut(src, dst);
    GPixel r = GPixel_PackARGB(GPixel_GetA(part1) + GPixel_GetA(part2), 
        GPixel_GetR(part1) + GPixel_GetR(part2), GPixel_GetG(part1) + GPixel_GetG(part2),
        GPixel_GetB(part1) + GPixel_GetB(part2));
    return r;
}

static inline GPixel blend(GPixel src, GPixel dst, GBlendMode b){
    switch(b) {
        case GBlendMode::kClear:
            return GPixel_PackARGB(0, 0, 0, 0);
        case GBlendMode::kSrc:
            return src;
        case GBlendMode::kDst:
            return dst;
        case GBlendMode::kSrcOver:
            return kSrcOver(src, dst);
        case GBlendMode::kDstOver:
            return kDstOver(src, dst);
        case GBlendMode::kSrcIn:
            return kSrcIn(src, dst);
        case GBlendMode::kDstIn:
            return kDstIn(src, dst);
        case GBlendMode::kSrcOut:
            return kSrcOut(src, dst);
        case GBlendMode::kDstOut:
            return kDstOut(src, dst);
        case GBlendMode::kSrcATop:
            return kSrcATop(src, dst);
        case GBlendMode::kDstATop:
            return kDstATop(src, dst);
        case GBlendMode::kXor:
            return kXor(src, dst);
    }

    return 0;
}

double distanceToLine(const GPoint& P1, const GPoint& P2, const GPoint& P3) {
    double x1 = P1.x, y1 = P1.y;
    double x2 = P2.x, y2 = P2.y;
    double x3 = P3.x, y3 = P3.y;
    
    double numerator = std::abs((x3 - x2) * (y2 - y1) - (y3 - y2) * (x2 - x1));
    double denominator = std::sqrt((x3 - x2) * (x3 - x2) + (y3 - y2) * (y3 - y2));
    
    return numerator / denominator;
}

Edge MyCanvas::makeEdge(const GPoint& point1, const GPoint& point2) {

    Edge newEdge;

    newEdge.point1 = point1;
    newEdge.point2 = point2;
    newEdge.top = GRoundToInt(std::min(point1.y, point2.y));
    newEdge.bot = GRoundToInt(std::max(point1.y, point2.y));
    newEdge.left = std::min(point1.x, point2.x);
    newEdge.right = std::max(point1.x, point2.x);

    newEdge.slope = ((point2.x - point1.x)/(point2.y - point1.y));
    newEdge.yIntercept = point1.x - (newEdge.slope*point1.y);

    if (point1.y < point2.y) { newEdge.up = true; } else { newEdge.up = false; }

    return newEdge;
}

void MyCanvas::blit(int L, int R, int y, GPaint paint) {
        std::vector<GPixel> row((R-L)+1);
        int xDup = L;
        GPixel dstPix;
        if (paint.peekShader() != nullptr) {
            paint.peekShader()->setContext(fCTM.back());
            paint.peekShader()->shadeRow(L, y, (R-L), &row[0]);
            for (int x = L; x < R; x++) {
                dstPix = *fDevice.getAddr(x, y);
                *fDevice.getAddr(x, y) = blend(row[x-xDup], dstPix, paint.getBlendMode());
            }
        }
        else {

            GPixel srcPix = color2Pix(paint.getColor());

            for (int x = L; x < R; x++) {
                dstPix = *fDevice.getAddr(x, y);
                *fDevice.getAddr(x, y) = blend(srcPix, dstPix, paint.getBlendMode());
            }
        }
}

void MyCanvas::clear(const GColor& color) {
    // your code here
    // fill my canvas with one color

    GPixel p = color2Pix(color);

    int w = fDevice.width();
    int h = fDevice.height();

    for(int x=0; x < w; x++) {
        for (int y=0; y < h; y++) {
            *fDevice.getAddr(x,y) = p;
        }
    }
}

static inline bool rect_stays_rect(const GMatrix& m) {
    return m[0] > 0 && m[3] > 0 && m[1] == 0 && m[2] == 0;
}

static inline GPoint* to_quad(const GRect& r, GPoint pts[4]) {
    pts[0] = {r.left, r.top};
    pts[1] = {r.right, r.top};
    pts[2] = {r.right, r.bottom};
    pts[3] = {r.left, r.bottom};
    return pts;
}

void MyCanvas::drawRect(const GRect& rect, const GPaint& paint) {
    GPoint corners[4] = {{rect.left, rect.top}, {rect.left, rect.bottom}, {rect.right, rect.bottom}, {rect.right, rect.top}};
    drawConvexPolygon(corners, 4, paint);
    return;
}

void MyCanvas::drawConvexPolygon(const GPoint* points, int count, const GPaint& paint) {

    if (count < 3) {
        return;
    }

    GPoint* transformed = new GPoint[count];
    fCTM.back().mapPoints(transformed, points, count);

    int left = fDevice.width();
    int right = 0;
    int top = fDevice.height();
    int bottom = 0;

    GPoint first = transformed[0];
    GPoint last;

    std::list<Edge> edges;

    for (int i=0; i<count-1; i++) {
        GPoint point1 = transformed[i];
        GPoint point2 = transformed[i+1];

        if (point1.y < top) {
            top = GRoundToInt(point1.y);
        }
        if (point1.y > bottom) {
            bottom = GRoundToInt(point1.y);
        }
        if (point1.y < left) {
            left = GRoundToInt(point1.y);
        }
        if (point1.x > right) {
            right = GRoundToInt(point1.x);
        }
        if (point2.y < top) {
            top = GRoundToInt(point2.y);
        }
        if (point2.y > bottom) {
            bottom = GRoundToInt(point2.y);
        }
        if (point2.x < left) {
            left = GRoundToInt(point2.x);
        }
        if (point2.x > right) {
            right = GRoundToInt(point2.x);
        }

        if (point1.y != point2.y) {
            edges.push_back(makeEdge(point1, point2));
        }

        if (i == count-2) {
            last = transformed[i+1];
        }
    }

    if (first.y != last.y) {
        edges.push_back(makeEdge(first, last));
    }

    if (edges.size() < 2) {
        return;
    }

    if (top < 0) {
        top = 0;
    }
    if (bottom > fDevice.height()) {
        bottom = fDevice.height();
    }
    if (left < 0) {
        left = 0;
    }
    if (right > fDevice.width()) {
        right = fDevice.width();
    }

    edges.sort([](const Edge& a, const Edge& b) {
        return a.top < b.top;
    });

    Edge edge1 = edges.front();
    edges.pop_front();
    Edge edge2= edges.front();
    edges.pop_front();

    for (int y = top; y<bottom; y++) {
        if (y > edge1.bot) {
            if (!edges.empty()) {
                edge1 = edges.front();
                edges.pop_front();
            }
            else {
                return;
            }
        }
        if (y > edge2.bot) {
            if (!edges.empty()) {
                edge2 = edges.front();
                edges.pop_front();
            }
            else {
                return;
            }
        }

        if (edges.empty() && y > edge1.bot && y > edge2.bot) {
            return;
        }


        int bound1 = GRoundToInt(edge1.getX(y));
        int bound2 = GRoundToInt(edge2.getX(y));

        if (bound1 < 0) {
            bound1 = 0;
        }
        if (bound1 > fDevice.width()) {
            bound1 = fDevice.width();
        }
        
        if (bound2 < 0) {
            bound2 = 0;
        }
        if (bound2 > fDevice.width()) {
            bound2 = fDevice.width();
        }

        int xDup = GRoundToInt(std::min(bound1, bound2));

        int len = (GRoundToInt(std::max(bound1, bound2)) - xDup) + 1;
        std::vector<GPixel> row(len);

        blit(xDup, std::max(bound1, bound2), y, paint);
    }
}


void MyCanvas::drawPath(const GPath& path, const GPaint& paint) {
    std::vector<Edge> edges;

    GPoint pts[GPath::kMaxNextPoints];
    GPoint tpts[GPath::kMaxNextPoints];
    GPath::Edger edgy(path);
    int top = 999999, bottom = 0;

    while (auto v = edgy.next(pts)) {

        fCTM.back().mapPoints(tpts, pts, GPath::kMaxNextPoints);

        switch (v.value()) {
            case GPathVerb::kLine: {
                if (GRoundToInt(tpts[0].y) != GRoundToInt(tpts[1].y)) {
                    Edge e = makeEdge(tpts[0], tpts[1]);
                    if (e.top < top) {
                        top = e.top;
                    }
                    if (e.bot > bottom) {
                        bottom = e.bot;
                    }
                    edges.push_back(e);
                }
                break;
            }
            case GPathVerb::kQuad: {
                if (tpts[0].y < top) {
                    top = tpts[0].y;
                }
                if (tpts[0].y > bottom) {
                    bottom = tpts[0].y;
                }
                if (tpts[1].y < top) {
                    top = tpts[1].y;
                }
                if (tpts[1].y > bottom) {
                    bottom = tpts[1].y;
                }
                if (tpts[2].y < top) {
                    top = tpts[2].y;
                }
                if (tpts[2].y > bottom) {
                    bottom = tpts[2].y;
                }

                GPoint tempSrc[] = {tpts[0], tpts[1], tpts[2]};
                GPoint tempForError[5];

                GPath::ChopQuadAt(tempSrc, tempForError, 0.5);

                double error = distanceToLine(tempForError[2], tpts[0], tpts[1]);
                int num_segs = (int)ceil(sqrt(error/0.25));

                int numEdges = 2;

                std::vector<std::vector<GPoint>> points;
                std::vector<std::vector<GPoint>> newPoints;

                points.push_back({tpts[0], tpts[1], tpts[2]});

                while (numEdges < num_segs) {
                    for (std::vector<GPoint> group : points) {
                        GPoint currentPoints[3] = {group[0], group[1], group[2]};
                        GPoint choppedPoints[5];

                        GPath::ChopQuadAt(currentPoints, &choppedPoints[0], 0.5);

                        newPoints.push_back({choppedPoints[0], choppedPoints[1], choppedPoints[2]});
                        newPoints.push_back({choppedPoints[2], choppedPoints[3], choppedPoints[4]});

                        numEdges++;
                    }

                    points = newPoints;
                    newPoints.clear();
                }

                std::vector<GPoint> finalPoints;
                for (std::vector<GPoint> g : points) {
                    finalPoints.push_back(g[0]);
                    finalPoints.push_back(g[1]);
                }

                finalPoints.push_back(points.back()[2]);

                for (int i = 0; i < finalPoints.size()-1; i++) {
                    if (finalPoints[i].y != finalPoints[i+1].y){
                        edges.push_back(makeEdge(finalPoints[i], finalPoints[i+1]));
                    }
                }

                break;
            }
            case GPathVerb::kCubic: {
                if (tpts[0].y < top) {
                    top = tpts[0].y;
                }
                if (tpts[0].y > bottom) {
                    bottom = tpts[0].y;
                }
                if (tpts[1].y < top) {
                    top = tpts[1].y;
                }
                if (tpts[1].y > bottom) {
                    bottom = tpts[1].y;
                }
                if (tpts[2].y < top) {
                    top = tpts[2].y;
                }
                if (tpts[2].y > bottom) {
                    bottom = tpts[2].y;
                }
                if (tpts[3].y < top) {
                    top = tpts[3].y;
                }
                if (tpts[3].y > bottom) {
                    bottom = tpts[3].y;
                }

                GPoint E;
                GPoint e0 = tpts[0] - (2*tpts[1]) + tpts[2];
                GPoint e1 = tpts[1] - (2*tpts[2]) + tpts[3];

                E.x = std::max(abs(e0.x), abs(e1.x));
                E.y = std::max(abs(e0.y), abs(e1.y));

                int numCSegs = (int)ceil(sqrt(3*sqrt((E.x*E.x)+(E.y*E.y))/(4*0.25)));

                int numCEdges = 2;
                std::vector<std::vector<GPoint>> cPoints;
                std::vector<std::vector<GPoint>> newCPoints;

                cPoints.push_back({tpts[0], tpts[1], tpts[2], tpts[3]});

                while (numCEdges < numCSegs) {
                    for (std::vector<GPoint> group : cPoints) {
                        GPoint currentPoints[4] = {group[0], group[1], group[2], group[3]};
                        GPoint choppedPoints[7];

                        GPath::ChopCubicAt(currentPoints, &choppedPoints[0], 0.5);

                        newCPoints.push_back({choppedPoints[0], choppedPoints[1], choppedPoints[2], choppedPoints[3]});
                        newCPoints.push_back({choppedPoints[3], choppedPoints[4], choppedPoints[5], choppedPoints[6]});

                        numCEdges++;
                    }

                    cPoints = newCPoints;
                    newCPoints.clear();
                }

                std::vector<GPoint> finalCPoints;
                for (std::vector<GPoint> g : cPoints) {
                    finalCPoints.push_back(g[0]);
                    finalCPoints.push_back(g[1]);
                    finalCPoints.push_back(g[2]);
                }

                finalCPoints.push_back(cPoints.back()[3]);

                for (int i = 0; i < finalCPoints.size()-1; i++) {
                    if (finalCPoints[i].y != finalCPoints[i+1].y){
                        edges.push_back(makeEdge(finalCPoints[i], finalCPoints[i+1]));
                    }
                }
                break;
                }
        }
    }

    if (top < 0) {
        top = 0;
    }
    if (bottom > fDevice.height()) {
        bottom = fDevice.height();
    }

    std::sort(edges.begin(), edges.end(), [](const Edge &a, const Edge &b) {
        return a.top < b.top;
    });

    std::vector<Edge> currentEdges;
    int start = 0;
    int w = 0;

    for (int y = GRoundToInt(top); y < GRoundToInt(bottom); y++) {
        for (Edge e : edges) {
            if (e.existsAtY(y)) {
                currentEdges.push_back(e);
            }
        }

        std::sort(currentEdges.begin(), currentEdges.end(), [y](Edge &a, Edge &b) {
            return a.getX(y) < b.getX(y);
        });

        for (Edge e : currentEdges) {
            if (w == 0) {
                start = clippedX(GRoundToInt(e.getX(y)));
            }

            if (e.up) {w += 1;} else {w -= 1;}

            if (w == 0) {
                blit(start, clippedX(GRoundToInt(e.getX(y))), y, paint);
            }
        }
        // assert(w==0);
        currentEdges.clear();
        w = 0;
    }
}

void MyCanvas::drawMesh(const GPoint verts[], const GColor colors[], const GPoint texs[],
                        int count, const int indices[], const GPaint& paint) {

    int n = 0;
    GPoint point0, point1, point2;
    GColor color0, color1, color2;
    GPoint texs0, texs1, texs2;

    if (colors) {
        for (int i = 0; i < count; ++i) {
            point0 = verts[indices[n+0]];
            point1 = verts[indices[n+1]];
            point2 = verts[indices[n+2]];

            color0 = colors[indices[n+0]];
            color1 = colors[indices[n+1]];
            color2 = colors[indices[n+2]];

            GPaint gradPaint = GPaint(GCreateTwoDGradient(point0, point1, point2, color0, color1, color2));

            GPoint points[] = {point0, point1, point2};
            
            drawConvexPolygon(points, 3, gradPaint);

            n += 3;
        }
    }
    else if (texs) {
        for (int i = 0; i < count; ++i) {
            point0 = verts[indices[n+0]];
            point1 = verts[indices[n+1]];
            point2 = verts[indices[n+2]];

            texs0 = texs[indices[n+0]];
            texs1 = texs[indices[n+1]];
            texs2 = texs[indices[n+2]];


            GPoint U = texs1-texs0;
            GPoint V = texs2-texs0;

            GMatrix T = GMatrix(U.x, V.x, texs0.x, U.y, V.y, texs0.y);
            GMatrix tInv;
            if (auto inv = T.invert()) {
                tInv = *inv;
            }

            GPoint pU = point1-point0;
            GPoint pV = point2-point0;

            GMatrix P = GMatrix(pU.x, pV.x, point0.x, pU.y, pV.y, point0.y);
            
            GMatrix tex2Draw = P*tInv;

            GPaint texPaint = GPaint(GCreateProxyShader(paint.peekShader(), tex2Draw));

            GPoint points[] = {point0, point1, point2};

            drawConvexPolygon(points, 3, texPaint);
            
            n += 3;
        }
    }
    else {
        for (int i = 0; i < count; ++i) {
            point0 = verts[indices[n+0]];
            point1 = verts[indices[n+1]];
            point2 = verts[indices[n+2]];
            
            n += 3;
        }
    }
}

void MyCanvas::drawQuad(const GPoint verts[4], const GColor colors[4], const GPoint texs[4],
                        int level, const GPaint& paint) {

    GPoint A = verts[0];
    GPoint B = verts[1];
    GPoint C = verts[2];
    GPoint D = verts[3];

    GPoint E;
    double error;
    int numSubdivisions;

    E = (A-B+C-D)*0.5;
    // std::cout << std::sqrt((E.x*E.x) + (E.y*E.y)) << std::endl;
    error = std::sqrt((E.x*E.x) + (E.y*E.y));

    numSubdivisions = std::ceil(std::pow(error, 0.25));
    
    quadDivider(verts, colors, texs, numSubdivisions-1, paint);
}

void MyCanvas::quadDivider(const GPoint verts[4], const GColor colors[4], const GPoint texs[4], int level, const GPaint& paint) {
    GPoint A1 = verts[0];
    GPoint B1 = getMidpoint(verts[0], verts[1]);
    GPoint D1 = getMidpoint(verts[0], verts[3]);

    GPoint A2 = B1;
    GPoint B2 = verts[1];
    GPoint C2 = getMidpoint(verts[1], verts[2]);
    GPoint D2 = getMidpoint(D1, C2);

    GPoint C1 = D2;

    GPoint A3 = C1;
    GPoint B3 = C2;
    GPoint C3 = verts[2];
    GPoint D3 = getMidpoint(verts[3], verts[2]);

    GPoint A4 = D1;
    GPoint B4 = C1;
    GPoint C4 = D3;
    GPoint D4 = verts[3];

    if (texs) {
        GPoint tA1 = texs[0];
        GPoint tB1 = getMidpoint(texs[0], texs[1]);
        GPoint tD1 = getMidpoint(texs[0], texs[3]);

        GPoint tA2 = tB1;
        GPoint tB2 = texs[1];
        GPoint tC2 = getMidpoint(texs[1], texs[2]);
        GPoint tD2 = getMidpoint(tD1, tC2);

        GPoint tC1 = tD2;

        GPoint tA3 = tC1;
        GPoint tB3 = tC2;
        GPoint tC3 = texs[2];
        GPoint tD3 = getMidpoint(texs[3], texs[2]);

        GPoint tA4 = tD1;
        GPoint tB4 = tC1;
        GPoint tC4 = tD3;
        GPoint tD4 = texs[3];
    

        if (colors) {
            GColor cA1 = colors[0];
            GColor cB1 = (colors[0]*0.5) + (colors[1]*0.5);
            GColor cD1 = (colors[0]*0.5) + (colors[3]*0.5);

            GColor cA2 = cB1;
            GColor cB2 = colors[1];
            GColor cC2 = (colors[1]*0.5) + (colors[2]*0.5);
            GColor cD2 = (cD1*0.5) + (cC2*0.5);

            GColor cC1 = cD2;

            GColor cA3 = cC1;
            GColor cB3 = cC2;
            GColor cC3 = colors[2];
            GColor cD3 = (colors[3]*0.5) + (colors[2]*0.5);

            GColor cA4 = cD1;
            GColor cB4 = cA3;
            GColor cC4 = cD3;
            GColor cD4 = colors[3];
        

            if (level == 0) {
                GPoint tri1[] = {A1, B1, D1};
                GPoint ttri1[] = {tA1, tB1, tD1};
                GColor ctri1[] = {cA1, cB1, cD1};

                GPoint tri2[] = {C1, B1, D1};
                GPoint ttri2[] = {tC1, tB1, tD1};
                GColor ctri2[] = {cC1, cB1, cD1};

                GPoint tri3[] = {A2, B2, D2};
                GPoint ttri3[] = {tA2, tB2, tD2};
                GColor ctri3[] = {cA2, cB2, cD2};

                GPoint tri4[] = {C2, B2, D2};
                GPoint ttri4[] = {tC2, tB2, tD2};
                GColor ctri4[] = {cC2, cB2, cD2};

                GPoint tri5[] = {A3, B3, D3};
                GPoint ttri5[] = {tA3, tB3, tD3};
                GColor ctri5[] = {cA3, cB3, cD3};

                GPoint tri6[] = {C3, B3, D3};
                GPoint ttri6[] = {tC3, tB3, tD3};
                GColor ctri6[] = {cC3, cB3, cD3};

                GPoint tri7[] = {A4, B4, D4};
                GPoint ttri7[] = {tA4, tB4, tD4};
                GColor ctri7[] = {cA4, cB4, cD4};

                GPoint tri8[] = {C4, B4, D4};
                GPoint ttri8[] = {tC4, tB4, tD4};
                GColor ctri8[] = {cC4, cB4, cD4};

                int indices[] = {0, 1, 2};
                drawMesh(tri1, ctri1, ttri1, 1, indices, paint);
                drawMesh(tri2, ctri2, ttri2, 1, indices, paint);
                drawMesh(tri3, ctri3, ttri3, 1, indices, paint);
                drawMesh(tri4, ctri4, ttri4, 1, indices, paint);
                drawMesh(tri5, ctri5, ttri5, 1, indices, paint);
                drawMesh(tri6, ctri6, ttri6, 1, indices, paint);
                drawMesh(tri7, ctri7, ttri7, 1, indices, paint);
                drawMesh(tri8, ctri8, ttri8, 1, indices, paint);
            }
            else {
                GPoint verts1[] = {A1, B1, C1, D1};
                GPoint texs1[] = {tA1, tB1, tC1, tD1};
                GColor colors1[] = {cA1, cB1, cC1, cD1};
                quadDivider(verts1, colors1, texs1, level-1, paint);

                GPoint verts2[] = {A2, B2, C2, D2};
                GPoint texs2[] = {tA2, tB2, tC2, tD2};
                GColor colors2[] = {cA2, cB2, cC2, cD2};
                quadDivider(verts2, colors2, texs2, level-1, paint);

                GPoint verts3[] = {A3, B3, C3, D3};
                GPoint texs3[] = {tA3, tB3, tC3, tD3};
                GColor colors3[] = {cA3, cB3, cC3, cD3};
                quadDivider(verts3, colors3, texs3, level-1, paint);

                GPoint verts4[] = {A4, B4, C4, D4};
                GPoint texs4[] = {tA4, tB4, tC4, tD4};
                GColor colors4[] = {cA4, cB4, cC4, cD4};
                quadDivider(verts4, colors4, texs4, level-1, paint);
            }
        }
        else {
            if (level == 0) {
                GPoint tri1[] = {A1, B1, D1};
                GPoint ttri1[] = {tA1, tB1, tD1};

                GPoint tri2[] = {C1, B1, D1};
                GPoint ttri2[] = {tC1, tB1, tD1};

                GPoint tri3[] = {A2, B2, D2};
                GPoint ttri3[] = {tA2, tB2, tD2};

                GPoint tri4[] = {C2, B2, D2};
                GPoint ttri4[] = {tC2, tB2, tD2};

                GPoint tri5[] = {A3, B3, D3};
                GPoint ttri5[] = {tA3, tB3, tD3};

                GPoint tri6[] = {C3, B3, D3};
                GPoint ttri6[] = {tC3, tB3, tD3};

                GPoint tri7[] = {A4, B4, D4};
                GPoint ttri7[] = {tA4, tB4, tD4};

                GPoint tri8[] = {C4, B4, D4};
                GPoint ttri8[] = {tC4, tB4, tD4};

                int indices[] = {0, 1, 2};
                drawMesh(tri1, colors, ttri1, 1, indices, paint);
                drawMesh(tri2, colors, ttri2, 1, indices, paint);
                drawMesh(tri3, colors, ttri3, 1, indices, paint);
                drawMesh(tri4, colors, ttri4, 1, indices, paint);
                drawMesh(tri5, colors, ttri5, 1, indices, paint);
                drawMesh(tri6, colors, ttri6, 1, indices, paint);
                drawMesh(tri7, colors, ttri7, 1, indices, paint);
                drawMesh(tri8, colors, ttri8, 1, indices, paint);
            }
            else {
                GPoint verts1[] = {A1, B1, C1, D1};
                GPoint texs1[] = {tA1, tB1, tC1, tD1};
                quadDivider(verts1, colors, texs1, level-1, paint);

                GPoint verts2[] = {A2, B2, C2, D2};
                GPoint texs2[] = {tA2, tB2, tC2, tD2};
                quadDivider(verts2, colors, texs2, level-1, paint);

                GPoint verts3[] = {A3, B3, C3, D3};
                GPoint texs3[] = {tA3, tB3, tC3, tD3};
                quadDivider(verts3, colors, texs3, level-1, paint);

                GPoint verts4[] = {A4, B4, C4, D4};
                GPoint texs4[] = {tA4, tB4, tC4, tD4};
                quadDivider(verts4, colors, texs4, level-1, paint);
            }
        }
    }
    else if (colors) {
        GColor cA1 = colors[0];
        GColor cB1 = (colors[0]*0.5) + (colors[1]*0.5);
        GColor cD1 = (colors[0]*0.5) + (colors[3]*0.5);

        GColor cA2 = cB1;
        GColor cB2 = colors[1];
        GColor cC2 = (colors[1]*0.5) + (colors[2]*0.5);
        GColor cD2 = (cD1*0.5) + (cC2*0.5);

        GColor cC1 = cD2;

        GColor cA3 = cC1;
        GColor cB3 = cC2;
        GColor cC3 = colors[2];
        GColor cD3 = (colors[3]*0.5) + (colors[2]*0.5);

        GColor cA4 = cD1;
        GColor cB4 = cA3;
        GColor cC4 = cD3;
        GColor cD4 = colors[3];
    

        if (level == 0) {
            GPoint tri1[] = {A1, B1, D1};
            GColor ctri1[] = {cA1, cB1, cD1};

            GPoint tri2[] = {C1, B1, D1};
            GColor ctri2[] = {cC1, cB1, cD1};

            GPoint tri3[] = {A2, B2, D2};
            GColor ctri3[] = {cA2, cB2, cD2};

            GPoint tri4[] = {C2, B2, D2};
            GColor ctri4[] = {cC2, cB2, cD2};

            GPoint tri5[] = {A3, B3, D3};
            GColor ctri5[] = {cA3, cB3, cD3};

            GPoint tri6[] = {C3, B3, D3};
            GColor ctri6[] = {cC3, cB3, cD3};

            GPoint tri7[] = {A4, B4, D4};
            GColor ctri7[] = {cA4, cB4, cD4};

            GPoint tri8[] = {C4, B4, D4};
            GColor ctri8[] = {cC4, cB4, cD4};

            int indices[] = {0, 1, 2};
            drawMesh(tri1, ctri1, texs, 1, indices, paint);
            drawMesh(tri2, ctri2, texs, 1, indices, paint);
            drawMesh(tri3, ctri3, texs, 1, indices, paint);
            drawMesh(tri4, ctri4, texs, 1, indices, paint);
            drawMesh(tri5, ctri5, texs, 1, indices, paint);
            drawMesh(tri6, ctri6, texs, 1, indices, paint);
            drawMesh(tri7, ctri7, texs, 1, indices, paint);
            drawMesh(tri8, ctri8, texs, 1, indices, paint);
        }
        else {
            GPoint verts1[] = {A1, B1, C1, D1};
            GColor colors1[] = {cA1, cB1, cC1, cD1};
            quadDivider(verts1, colors1, texs, level-1, paint);

            GPoint verts2[] = {A2, B2, C2, D2};
            GColor colors2[] = {cA2, cB2, cC2, cD2};
            quadDivider(verts2, colors2, texs, level-1, paint);

            GPoint verts3[] = {A3, B3, C3, D3};
            GColor colors3[] = {cA3, cB3, cC3, cD3};
            quadDivider(verts3, colors3, texs, level-1, paint);

            GPoint verts4[] = {A4, B4, C4, D4};
            GColor colors4[] = {cA4, cB4, cC4, cD4};
            quadDivider(verts4, colors4, texs, level-1, paint);
        }
    }
}

GPoint MyCanvas::getMidpoint(GPoint a, GPoint b) {
    GPoint r;
    r.x = (a.x+b.x)*0.5;
    r.y = (a.y+b.y)*0.5;
    return r;
}

std::unique_ptr<GCanvas> GCreateCanvas(const GBitmap& device) {
    return std::unique_ptr<GCanvas>(new MyCanvas(device));
}

std::string GDrawSomething(GCanvas* canvas, GISize dim) {

    GPaint paintor = GPaint(GColor::RGBA(0.5, 0.5, 0.5, 1));

    GPoint p1, p2, p3;
    p1.x = 10;
    p1.y = 200;
    p2.x = 200;
    p2.y = 200;
    p3.x = 90;
    p3.y = 10;
    GPoint points[3] = {p1, p2, p3};

    int indices[3] = {0, 1, 2};

    GColor c1 = GColor::RGBA(1, 0.56, 0, 1);
    GColor c2 = GColor::RGBA(0, 0.1, 0.9, 1);

    GColor colors[] = {c1, c1, c2};

    canvas->drawMesh(points, colors, NULL, 1, indices, GPaint());

    return "Oblong"; //name of artwork
}