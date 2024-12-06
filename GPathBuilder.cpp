#include "./include/GPathBuilder.h"
#include <cmath>
#include <numbers>
#include <iostream>

void GPathBuilder::addRect(const GRect& r, GPathDirection orientation) {
    switch (orientation) {
        case (GPathDirection::kCW):
            moveTo({r.left, r.top});
            lineTo({r.right, r.top});
            lineTo({r.right, r.bottom});
            lineTo({r.left, r.bottom});
            break;
        case (GPathDirection::kCCW):
            moveTo({r.left, r.top});
            lineTo({r.left, r.bottom});
            lineTo({r.right, r.bottom});
            lineTo({r.right, r.top});
    }
}

void GPathBuilder::addPolygon(const GPoint pts[], int count) {
    moveTo(pts[0]);
    for (int i = 1; i < count; i++) {
        lineTo(pts[i]);
    }
}

void GPathBuilder::addCircle(GPoint center, float radius, GPathDirection direction) {
    GMatrix mtx = GMatrix::Translate(center.x, center.y)*GMatrix::Scale(radius, radius);

    double pi = acos(-1.0);

    GPoint A = {1, 0};
    GPoint B = {1, tan(pi/8)};
    GPoint C = {sqrt(2)/2, sqrt(2)/2};
    GPoint D = {tan(pi/8), 1};
    GPoint E = {0, 1};

    GPoint F = {-tan(pi/8), 1};
    GPoint G = {-sqrt(2)/2, sqrt(2)/2};
    GPoint H = {-1, tan(pi/8)};
    GPoint I = {-1, 0};

    GPoint J = {-1, -tan(pi/8)};
    GPoint K = {-sqrt(2)/2, -sqrt(2)/2};
    GPoint L = {-tan(pi/8), -1};
    GPoint M = {0, -1};

    GPoint N = {tan(pi/8), -1};
    GPoint O = {sqrt(2)/2, -sqrt(2)/2};
    GPoint P = {1, -tan(pi/8)};
    
    A = mtx*A;
    B = mtx*B;
    C = mtx*C;
    D = mtx*D;
    E = mtx*E;
    F = mtx*F;
    G = mtx*G;
    H = mtx*H;
    I = mtx*I;
    J = mtx*J;
    K = mtx*K;
    L = mtx*L;
    M = mtx*M;
    N = mtx*N;
    O = mtx*O;
    P = mtx*P;

    if (direction == GPathDirection::kCW) {
        moveTo(A);
        quadTo(B, C);
        quadTo(D, E);
        quadTo(F, G);
        quadTo(H, I);
        quadTo(J, K);
        quadTo(L, M);
        quadTo(N, O);
        quadTo(P, A);
    }
    else {
        moveTo(A);
        quadTo(P, O);
        quadTo(N, M);
        quadTo(L, K);
        quadTo(J, I);
        quadTo(H, G);
        quadTo(F, E);
        quadTo(D, C);
        quadTo(B, A);
    }
}