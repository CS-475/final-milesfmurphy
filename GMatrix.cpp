#include "include/GMatrix.h"


GMatrix::GMatrix() : GMatrix(1, 0, 0,
                             0, 1, 0) {}

GMatrix GMatrix::Translate(float tx, float ty) {
    return GMatrix(1.0f, 0.0f, tx, 
                   0.0f, 1.0f, ty);
}

GMatrix GMatrix::Scale(float sx, float sy) {
    return GMatrix(sx, 0.0f, 0.0f,
                   0.0f, sy, 0.0f);
}

GMatrix GMatrix::Rotate(float radians) {
    float cos = cosf(radians);
    float sin = sinf(radians);
    return GMatrix(cos, -sin, 0.0f,
                   sin, cos, 0.0f);
}

GMatrix GMatrix::Concat(const GMatrix& a, const GMatrix& b) {
    float result[6] = {
        (a[0] * b[0]) + (a[2] * b[1]), (a[1] * b[0]) + (a[3] * b[1]), (a[0] * b[2]) + (a[2] * b[3]),
        (a[1] * b[2]) + (a[3] * b[3]), (a[0] * b[4]) + (a[2] * b[5]) + a[4], (a[1] * b[4]) + (a[3] * b[5]) + a[5]
    };
    return GMatrix(result[0], result[2], result[4], result[1], result[3], result[5]);
}

nonstd::optional<GMatrix> GMatrix::invert() const {
    double det = (fMat[0] * fMat[3]) - (fMat[1] * fMat[2]);
    if (det == 0) {
        return nonstd::nullopt;
    }

    det = (1.0f / det);

    double inverted[6] = {
        fMat[3]*det, -1.0f*fMat[1]*det, -1.0f*fMat[2]*det,
        fMat[0]*det, ((fMat[2]*fMat[5])-(fMat[3]*fMat[4]))*det, ((fMat[1]*fMat[4])-(fMat[0]*fMat[5]))*det
    };

    return GMatrix(inverted[0], inverted[2], inverted[4], inverted[1], inverted[3], inverted[5]);
}

void GMatrix::mapPoints(GPoint dst[], const GPoint src[], int count) const {
    for (int i=0; i<count; i++) {
        double x = src[i].x;
        double y = src[i].y;
        dst[i].x = (*this)[0]*x + (*this)[2]*y + (*this)[4];
        dst[i].y = (*this)[1]*x + (*this)[3]*y + (*this)[5];
    }
}