#ifndef _MY_GL_H_
#define _MY_GL_H_


#include "tgaimage.h"
#include "math.h"
#include <array>

extern Matrix Projection;
extern Matrix MV;
extern Matrix ViewPort;
extern Matrix MV;



struct MyShader {
    virtual ~MyShader(){};
    virtual void vertex(int i) = 0;
    virtual void fragment(float x, float y, Vec3f &c, std::vector<std::vector<double>>& zbuffer, TGAImage &image, TGAColor &color, Vec3f &intensity) = 0;
};



void line(Vec3f p0, Vec3f p1, TGAImage &image, TGAColor color);

Vec3f Centroid(std::array<Vec3f, 3> pts, float x, float y);

void triangle2(Vec2i t0, Vec2i t1, Vec2i t2, TGAImage &image, float intensity, TGAColor color);

void triangle( std::vector<std::array<Vec3f, 3>> node, MyShader &shader,std::vector<std::vector<double>> &zbuffer, TGAImage &image, Vec3f intensity);

Vec3f textpos(Vec3f v);

Vec3f m2v(Matrix m);

Matrix v2m(Vec3f v);

void viewport(int x, int y, int w, int h);

Matrix translation(Vec3f v);

Matrix zoom(float factor);

Matrix rotation_x(float cosangle, float sinangle);

Matrix rotation_y(float cosangle, float sinangle);

Matrix rotation_z(float cosangle, float sinangle);

void ModelView(Vec3f eye, Vec3f center, Vec3f up);

void projection(float v);



#endif


