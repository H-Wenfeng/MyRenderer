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

struct Material
{
    Vec3f diffusion_color;
    float diffusion_coefficient;
    float specular_coefficient;
    float reflect_coefficient;
    float refract_coefficient;
    float specular_exponent;
    float refract_exponent;
    Material(Vec3f diffusion_color, float diffusion_coefficient, float specular_coefficient, float reflect_coefficient, float refract_coefficient, float specular_eponent, float refract_exponent);
};



struct Sphere
{   
    Vec3f center;
    float r;  
    Material material;
    Sphere(Vec3f c, float radius, Material material);
    bool intersect(const Vec3f &light, const Vec3f &camera, Vec3f &intersectionPoint);
};

struct Light {
    Light(const Vec3f &p, const float &i) : position(p), intensity(i) {}
    Vec3f position;
    float intensity;
};




void Ray_Render(Vec3f &camera, std::vector<Sphere>  &Spheres, TGAImage &image, std::vector<std::vector<double>> &zbuffer, std::vector<Light>  &lights);

bool plane(Vec3f &dir, Vec3f &orig, Vec3f &intersectionPoint);

Vec3f cast_ray(Vec3f &dir, Vec3f &origin, size_t p, std::vector<Sphere> &Spheres, std::vector<Light> &lights, int &times);

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


