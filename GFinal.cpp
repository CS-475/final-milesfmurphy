#include <memory>
#include <algorithm>
#include "include/GFinal.h"
#include "include/GBitmap.h"
#include "include/GMatrix.h"
#include "include/GPixel.h"
#include "include/GPoint.h"
#include "include/GColor.h"
#include "shape.cpp"
#include "ProxyShader.h"

class MatrixColorShader : public GShader {
public:
    MatrixColorShader(GShader* baseShader, const GColorMatrix& colorTransform)
        : shaderDelegate(baseShader), colorMatrix(colorTransform) {}

    bool isOpaque() override {
        return shaderDelegate->isOpaque();
    }

    bool setContext(const GMatrix& contextMatrix) override {
        return shaderDelegate->setContext(contextMatrix);
    }

    void shadeRow(int x, int y, int count, GPixel pixels[]) override {
        shaderDelegate->shadeRow(x, y, count, pixels);
        for (int idx = 0; idx < count; ++idx) {
            float alphaScale = GPixel_GetA(pixels[idx]) / 255.0f;
            float alpha = GPixel_GetA(pixels[idx]) / 255.0f;
            float red = GPixel_GetR(pixels[idx]) / (255.0f * alphaScale);
            float green = GPixel_GetG(pixels[idx]) / (255.0f * alphaScale);
            float blue = GPixel_GetB(pixels[idx]) / (255.0f * alphaScale);

            float newRed = (colorMatrix[0] * red + colorMatrix[4] * green + colorMatrix[8] * blue +
                            colorMatrix[12] * alpha + colorMatrix[16]);
            float newGreen = (colorMatrix[1] * red + colorMatrix[5] * green + colorMatrix[9] * blue +
                              colorMatrix[13] * alpha + colorMatrix[17]);
            float newBlue = (colorMatrix[2] * red + colorMatrix[6] * green + colorMatrix[10] * blue +
                             colorMatrix[14] * alpha + colorMatrix[18]);
            float newAlpha = (colorMatrix[3] * red + colorMatrix[7] * green + colorMatrix[11] * blue +
                              colorMatrix[15] * alpha + colorMatrix[19]);

            // Clamping to [0, 1]
            newAlpha = std::clamp(newAlpha, 0.0f, 1.0f);
            newRed = std::clamp(newRed, 0.0f, 1.0f);
            newGreen = std::clamp(newGreen, 0.0f, 1.0f);
            newBlue = std::clamp(newBlue, 0.0f, 1.0f);

            int pixelAlpha = GRoundToInt(newAlpha * 255.0f);
            int pixelRed = GRoundToInt(newRed * newAlpha * 255.0f);
            int pixelGreen = GRoundToInt(newGreen * newAlpha * 255.0f);
            int pixelBlue = GRoundToInt(newBlue * newAlpha * 255.0f);

            pixels[idx] = GPixel_PackARGB(pixelAlpha, pixelRed, pixelGreen, pixelBlue);
        }
    }

private:
    GShader* shaderDelegate;
    GColorMatrix colorMatrix;
};

class CustomFinal : public GFinal {
public:
    CustomFinal() {}

    std::shared_ptr<GShader> createVoronoiShader(const GPoint points[], const GColor colors[], int count) override {
        if (count < 1) {
            return nullptr;
        }
        return std::unique_ptr<GShader>(new CustomVoronoi(points, colors, count));
    }

    std::shared_ptr<GShader> createColorMatrixShader(const GColorMatrix& matrix, GShader* shader) override {
        return std::unique_ptr<GShader>(new MatrixColorShader(shader, matrix));
    }
};

std::unique_ptr<GFinal> GCreateFinal() {
    return std::unique_ptr<GFinal>(new CustomFinal());
}