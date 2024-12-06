#include <memory>
#include<algorithm>
#include "include/GFinal.h"
#include "include/GBitmap.h"
#include "include/GMatrix.h"
#include "include/GPixel.h"
#include "include/GPoint.h"
#include "include/GColor.h"
#include "shape.cpp"
#include "ProxyShader.h"

class ColorMatrixShader : public GShader {
    public:
        ColorMatrixShader(GShader* shader, const GColorMatrix& extraTransform)
        : fRealShader(shader), fMat(extraTransform) {}
        
        bool isOpaque() override { return fRealShader->isOpaque(); }

        bool setContext(const GMatrix& ctm) override {
            return fRealShader->setContext(ctm);
        }

        void shadeRow(int x, int y, int count, GPixel row[]) override {
            fRealShader->shadeRow(x, y, count, row);
            for(int i = 0; i<count; i++){
                GPixel p;
                float sa = GPixel_GetA(row[i]);
                float a =  GPixel_GetA(row[i]) / 255.0f;
                float r = GPixel_GetR(row[i])/(255.0*sa); 
                float g = GPixel_GetG(row[i])/(255.0*sa);
                float b = GPixel_GetB(row[i])/(255.0*sa);
                
                
                float r_new = (fMat[0] * r + fMat[4] * g + fMat[8] * b + fMat[12] * a + fMat[16]);
                float g_new = (fMat[1] * r + fMat[5] * g + fMat[9] * b + fMat[13] * a + fMat[17]);
                float b_new = (fMat[2] * r + fMat[6] * g + fMat[10] * b + fMat[14] * a + fMat[18]);
                float a_new = (fMat[3] * r + fMat[7] * g + fMat[11] * b + fMat[15] * a + fMat[19]);

                if(a_new < 0.0){
                    a_new = 0.0;
                }
                if(a_new > 1.0){
                    a_new = 1.0;
                }
                if(r_new < 0.0){
                    r_new = 0.0;
                }
                if(r_new > 1.0){
                    a_new = 1.0;
                }
                if(g_new < 0.0){
                    g_new = 0.0;
                }
                if(g_new > 1.0){
                    g_new = 1.0;
                }
                if(b_new < 0.0){
                    b_new = 0.0;
                }
                if(b_new > 1.0){
                    b_new = 1.0;
                }
                
                int newA = GRoundToInt(a_new * 255.0);
                int newR = GRoundToInt(r_new * a_new * 255.0);
                int newG = GRoundToInt(g_new * a_new * 255.0);
                int newB = GRoundToInt(b_new * a_new * 255.0);

                row[i] = GPixel_PackARGB(newA, newR, newG, newB);
            }
            fRealShader->shadeRow(x, y, count, row);
        }

    private:
        GShader* fRealShader;
        GColorMatrix  fMat;
};

class MyFinal : public GFinal{
    public:
    MyFinal(){};
    std::shared_ptr<GShader> createVoronoiShader(const GPoint points[], const GColor colors[], int count) {
        if (count < 1) {
        return nullptr;
    }
    return std::unique_ptr<GShader>(new CustomVoronoi(points, colors, count));
    }

     std::shared_ptr<GShader> createColorMatrixShader(const GColorMatrix& mat, GShader* realShader) {
        return std::unique_ptr<GShader>(new ColorMatrixShader(realShader, mat));
    }
};





std::unique_ptr<GFinal> GCreateFinal(){
    return std::unique_ptr<GFinal>(new MyFinal());
};