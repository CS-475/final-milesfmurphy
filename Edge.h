#include "include/GPoint.h"
#include <iostream>

struct Edge {
    float left, right;
    int top, bot;
    GPoint point1, point2;
    GPoint topPoint, botPoint, leftPoint, rightPoint;
    
    float slope;
    float yIntercept;

    bool up;

    bool existsAtY(int y) {
        if ((y<bot) && (y>=top)) {
            return true;
        }

        return false;
    }

    float getX(int y) {
        return float(((y+0.5)*slope) + yIntercept);
    }

    GPoint getTopPoint() {
        if (point1.y < point2.y) {
            return point1;
        }
        else {
            return point2;
        }
    }
    GPoint getBotPoint() {
        if (point1.y < point2.y) {
            return point2;
        }
        else {
            return point1;
        }
    }    
    GPoint getLeftPoint() {
        if (point1.x < point2.x) {
            return point1;
        }
        else {
            return point2;
        }
    }
    GPoint getRightPoint() {
        if (point1.x < point2.x) {
            return point2;
        }
        else {
            return point1;
        }
    }
};