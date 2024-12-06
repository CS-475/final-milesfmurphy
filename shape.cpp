#include "include/GShader.h"
#include "include/GBitmap.h"
#include "include/GMatrix.h"
#include "include/GPixel.h"
#include "include/GPoint.h"
#include "include/GColor.h"
#include <math.h>
#include <iostream>
#include <vector>

class CustomVoronoi : public GShader {
public:
    CustomVoronoi(const GPoint vertices[], const GColor shades[], int numPoints) {
        for (int idx = 0; idx < numPoints; ++idx) {
            colorPalette.push_back(shades[idx]);
            coordSet.push_back(vertices[idx]);
        }
        totalPoints = numPoints;
    }

    bool setContext(const GMatrix& transformation) override {
        GMatrix combinedMatrix = transformation * localMatrix;
        if (auto maybeInvertible = combinedMatrix.invert()) {
            inverseMatrix = *maybeInvertible;
            return true;
        }
        return false;
    }

    bool isOpaque() override {
        return false; // Assuming default non-opaque behavior.
    }

    void shadeRow(int startX, int startY, int rowLength, GPixel pixelBuffer[]) override {
        GPoint targetPoint;
        GPoint mappedPoint;

        for (int idx = 0; idx < rowLength; ++idx) {
            targetPoint.x = startX + idx + 0.5f;
            targetPoint.y = startY + 0.5f;
            mappedPoint = inverseMatrix * targetPoint;

            float sourceX = mappedPoint.x;
            float sourceY = mappedPoint.y;

            GColor chosenColor = getClosestColor(sourceX, sourceY);
            pixelBuffer[idx] = convertColorToPixel(chosenColor);
        }
    }

private:
    GColor getClosestColor(float x, float y) {
        float nearestDistance = INFINITY;
        int closestIndex = 0;

        for (int idx = 0; idx < totalPoints; ++idx) {
            float dx = coordSet[idx].x - x;
            float dy = coordSet[idx].y - y;
            float distance = sqrt(dx * dx + dy * dy);

            if (distance < nearestDistance) {
                nearestDistance = distance;
                closestIndex = idx;
            }
        }
        return colorPalette[closestIndex];
    }

    GPixel convertColorToPixel(const GColor& color) {
        int alpha = GRoundToInt(color.a * 255);
        int red = GRoundToInt(color.r * color.a * 255);
        int green = GRoundToInt(color.g * color.a * 255);
        int blue = GRoundToInt(color.b * color.a * 255);

        return GPixel_PackARGB(alpha, red, green, blue);
    }

    std::vector<GColor> colorPalette;
    std::vector<GPoint> coordSet;
    GMatrix inverseMatrix;
    GMatrix localMatrix;
    int totalPoints;
};