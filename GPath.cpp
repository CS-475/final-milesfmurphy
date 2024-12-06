#include "include/GPath.h"
#include <iostream>
#include "./include/GRect.h"


GPoint mid(const GPoint& p0, const GPoint& p1, float t) {
    return {p0.x + ((p1.x - p0.x) * t),
            p0.y + ((p1.y - p0.y) * t)};
}

std::vector<double> solveQuadratic(double A, double B, double C) {
    
    double discriminant = (B*B)-(4*A*C);
    std::vector<double> ret = {};

    if (A == 0) {
        if (B != 0) {
            double root1 = (-C)/B;

            if (root1 >= 0 && root1 <= 1) {
                ret.push_back(root1);
            }
        }
        return ret;
    }

    if (C == 0) {
        double root1 = (-B)/A;

        if (root1 >= 0 && root1 <= 1) {
            ret.push_back(root1);
        }
        return ret; 
    }

    if (discriminant > 0) {
        
        double root1 = ((-B) + sqrt(discriminant)) / (2*A);
        double root2 = ((-B) - sqrt(discriminant)) / (2*A);

        if (root1 >= 0 && root1 <= 1) {
            ret.push_back(root1);
        }
        if (root2 >= 0 && root2 <= 1) {
            ret.push_back(root2);
        }
    }

    else if (discriminant == 0) {
        double root1 = (-B)/(2*A);

        if (root1 >= 0 && root1 <= 1) {
            ret.push_back(root1);
        }
    }

    return ret;
}

float findQuadTX(GPoint A, GPoint B, GPoint C) {
    float tX = (A.x-B.x)/(A.x-(2*B.x)+C.x);

    if (tX <= 1.0 && tX >= 0.0) {
        return tX;
    }

    return 2.0;
}

float findQuadTY(GPoint A, GPoint B, GPoint C) {
    float tY = (A.y-B.y)/(A.y-(2*B.y)+C.y);

    if (tY <= 1.0 && tY >= 0.0) {
        return tY;
    }

    return 2.0;
}

GPoint getQuadExtrema(GPoint A, GPoint B, GPoint C, float t) {
    float x = A.x*((1-t)*(1-t)) + 2*B.x*t*(1-t) + C.x*t*t;
    float y = A.y*((1-t)*(1-t)) + 2*B.y*t*(1-t) + C.y*t*t;
    return {x, y};
}

GPoint getCubicExtrema(GPoint A, GPoint B, GPoint C, GPoint D, float t) {
    double x = (A.x*(1-t)*(1-t)*(1-t)) + (3*B.x*t*(1-t)*(1-t)) + (3*C.x*(1-t)*(t)*(t)) + D.x*(t)*(t)*(t);
    double y = (A.y*(1-t)*(1-t)*(1-t)) + (3*B.y*t*(1-t)*(1-t)) + (3*C.y*(1-t)*(t)*(t)) + D.y*(t)*(t)*(t);

    return {x, y};
}

GRect GPath::bounds() const {
    // int top = 999999, bot = 0, left = 999999, right = 0;

    // if (fPts.empty()) {
    //     return GRect::LTRB(0, 0, 0, 0);
    // }

    // int index = 0;

    // for (int i = 0; i < fVbs.size(); i++) {
    //     switch (fVbs[i]) {
    //         case GPathVerb::kMove:
    //             if (fVbs[i+1] != GPathVerb::kMove) {
    //                 if (fPts[index].x < left) {
    //                     left = GRoundToInt(fPts[index].x);
    //                 }
    //                 if (fPts[index].x > right) {
    //                     right = GRoundToInt(fPts[index].x);
    //                 }
    //                 if (fPts[index].y > bot) {
    //                     bot = GRoundToInt(fPts[index].y);
    //                 }
    //                 if (fPts[index].y < top) {
    //                     top = GRoundToInt(fPts[index].y);
    //                 }
    //             }

    //             index++;
    //             break;
    //         case GPathVerb::kLine:
    //             if (fPts[index].x < left) {
    //                 left = GRoundToInt(fPts[index].x);
    //             }
    //             if (fPts[index].x > right) {
    //                 right = GRoundToInt(fPts[index].x);
    //             }
    //             if (fPts[index].y > bot) {
    //                 bot = GRoundToInt(fPts[index].y);
    //             }
    //             if (fPts[index].y < top) {
    //                 top = GRoundToInt(fPts[index].y);
    //             }

    //             index++;
    //             break;
    //         case GPathVerb::kQuad:
    //             if (float t = findQuadTX(fPts[index], fPts[index+1], fPts[index+2]) <= 1.1) {
    //                 GPoint extrema = getQuadExtrema(fPts[index], fPts[index+1], fPts[index+2], t);

    //                 if (extrema.x < left) {
    //                     left = GRoundToInt(extrema.x);
    //                 }
    //                 if (extrema.x > right) {
    //                     right = GRoundToInt(extrema.x);
    //                 }
    //                 if (extrema.y > bot) {
    //                     bot = GRoundToInt(extrema.y);
    //                 }
    //                 if (extrema.y < top) {
    //                     top = GRoundToInt(extrema.y);
    //                 }
    //             }

    //             if (float t = findQuadTY(fPts[index], fPts[index+1], fPts[index+2]) <= 1.1) {
    //                 GPoint extrema = getQuadExtrema(fPts[index], fPts[index+1], fPts[index+2], t);

    //                 if (extrema.x < left) {
    //                     left = GRoundToInt(extrema.x);
    //                 }
    //                 if (extrema.x > right) {
    //                     right = GRoundToInt(extrema.x);
    //                 }
    //                 if (extrema.y > bot) {
    //                     bot = GRoundToInt(extrema.y);
    //                 }
    //                 if (extrema.y < top) {
    //                     top = GRoundToInt(extrema.y);
    //                 }
    //             }

    //             if (fPts[index+2].x < left) {
    //                 left = GRoundToInt(fPts[index+2].x);
    //             }
    //             if (fPts[index+2].x > right) {
    //                 right = GRoundToInt(fPts[index+2].x);
    //             }
    //             if (fPts[index+2].y > bot) {
    //                 bot = GRoundToInt(fPts[index+2].y);
    //             }
    //             if (fPts[index+2].y < top) {
    //                 top = GRoundToInt(fPts[index+2].y);
    //             }

    //             index++;
    //             index++;
    //             break;

    //         case GPathVerb::kCubic:
    //             GPoint A = fPts[index];
    //             GPoint B = fPts[index+1];
    //             GPoint C = fPts[index+2];
    //             GPoint D = fPts[index+3];

    //             std::vector<double> xRoots = solveQuadratic((-3*A.x)+(9*B.x)-(9*C.x)+(3*D.x), (6*A.x)-(12*B.x)+(6*C.x), (3*A.x)-(3*B.x));
    //             std::vector<double> yRoots = solveQuadratic((-3*A.y)+(9*B.y)-(9*C.y)+(3*D.y), (6*A.y)-(12*B.y)+(6*C.y), (3*A.y)-(3*B.y));

    //             for (double root : xRoots) {
    //                 GPoint extrema = getCubicExtrema(A, B, C, D, root);

    //                 if (extrema.x < left) {
    //                     left = GRoundToInt(extrema.x);
    //                 }
    //                 if (extrema.x > right) {
    //                     right = GRoundToInt(extrema.x);
    //                 }
    //                 if (extrema.y > bot) {
    //                     bot = GRoundToInt(extrema.y);
    //                 }
    //                 if (extrema.y < top) {
    //                     top = GRoundToInt(extrema.y);
    //                 }
    //             }
    //             for (double root : yRoots) {
    //                 GPoint extrema = getCubicExtrema(A, B, C, D, root);
                    
    //                 if (extrema.x < left) {
    //                     left = GRoundToInt(extrema.x);
    //                 }
    //                 if (extrema.x > right) {
    //                     right = GRoundToInt(extrema.x);
    //                 }
    //                 if (extrema.y > bot) {
    //                     bot = GRoundToInt(extrema.y);
    //                 }
    //                 if (extrema.y < top) {
    //                     top = GRoundToInt(extrema.y);
    //                 }
    //             }
    //             break;
    //     }
    // }
    // return GRect::LTRB(left, top, right, bot);

    
    if (fPts.empty()) {
        return GRect::LTRB(0, 0, 0, 0);
    }

    int top = GRoundToInt(fPts[0].y), bot = GRoundToInt(fPts[0].y), left = GRoundToInt(fPts[0].x), right = GRoundToInt(fPts[0].x);

    for (int i=0; i<fPts.size(); i++) {
        if (fPts[i].x < left) {
            left = GRoundToInt(fPts[i].x);
        }
        if (fPts[i].x > right) {
            right = GRoundToInt(fPts[i].x);
        }
        if (fPts[i].y > bot) {
            bot = GRoundToInt(fPts[i].y);
        }
        if (fPts[i].y < top) {
            top = GRoundToInt(fPts[i].y);
        }
    }

    return GRect::LTRB(left, top, right, bot);
}

void GPath::ChopQuadAt(const GPoint src[3], GPoint dst[5], float t) {

    GPoint AB = mid(src[0], src[1], t);
    GPoint BC = mid(src[1], src[2], t);

    GPoint midP = mid(AB, BC, t);

    dst[0] = src[0];
    dst[1] = AB;
    dst[2] = midP;
    dst[3] = BC;
    dst[4] = src[2];

    // std::cout << "(" << src[0].x << ", " << src[0].y << "), (" << AB.x << ", " << AB.y << "), (" << midP.x << ", " << midP.y << "), (" << BC.x << ", " << BC.y << "), (" << src[2].x << ", " << src[2].y << ")" << std::endl;
}

void GPath::ChopCubicAt(const GPoint src[4], GPoint dst[7], float t) {

    GPoint AB = mid(src[0], src[1], t);
    GPoint BC = mid(src[1], src[2], t);
    GPoint CD = mid(src[2], src[3], t);

    GPoint ABBC = mid(AB, BC, t);
    GPoint BCCD = mid(BC, CD, t);

    GPoint midP = mid(ABBC, BCCD, t);

    dst[0] = src[0];
    dst[1] = AB;
    dst[2] = ABBC;
    dst[3] = midP;
    dst[4] = BCCD;
    dst[5] = CD;
    dst[6] = src[3];
}