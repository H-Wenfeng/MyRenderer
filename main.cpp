#include <vector>
#include <cmath>
#include "tgaimage.h"
#include "model.h"
#include "math.h"
#include <iostream>
#include <array>
#include "mygl.h"
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red = TGAColor(255, 0, 0, 255);
const TGAColor green = TGAColor(0, 255, 0, 255);
const TGAColor blue = TGAColor(0, 0, 255, 255);
const TGAColor yellow = TGAColor(255, 255, 0, 255);

extern const int width = 800;
extern const int height = 800;

std::vector<std::vector<double>> shadowbuffer;
Model *model = NULL;

Vec3f eye(1, 0, 5);
Vec3f center(0, 0, 0);
Vec3f origin_light_dir(1, 1, 1);
int if_change_light = 1;
Matrix tobuffer = Matrix::identity(4);
Display *display = XOpenDisplay(NULL);
Window window = XCreateSimpleWindow(display, RootWindow(display, DefaultScreen(display)), 10, 10, 900, 900, 1, BlackPixel(display, DefaultScreen(display)), WhitePixel(display, DefaultScreen(display)));
XEvent event;
GC gc = XCreateGC(display, window, 0, nullptr);
bool using_tangent = false;
bool using_nm = true;
TGAColor nm;

void show_image(TGAImage &image)
{

    char image2[900][900][4];
    for (int i = 0; i <= 800; i++)
        for (int j = 0; j <= 800; j++)
        {
            image2[j][i][0] = image.get(i, j).b;
            image2[j][i][1] = image.get(i, j).g;
            image2[j][i][2] = image.get(i, j).r;
            image2[j][i][3] = 0;
        }

    XImage *xImage = XCreateImage(display, DefaultVisual(display, DefaultScreen(display)), 24,
                                  ZPixmap, 0, reinterpret_cast<char *>(image2), 900, 900, 32, 0);

    // 绘制图像
    XPutImage(display, window, gc, xImage, 0, 0, 0, 0, 900, 900);
    XFlush(display);
}
void camera_control(char keyBuffer[32])
{

    if (keyBuffer[0] == 'a' || keyBuffer[0] == 'A')
    {
        eye.x = eye.x - 1;
        event.type = Expose;
        XPutBackEvent(display, &event);
        XClearWindow(display, window);
        std::cout << "Simulated Expose event" << std::endl;
        if_change_light = 0;
    }
    if (keyBuffer[0] == 'd' || keyBuffer[0] == 'D')
    {
        eye.x = eye.x + 1;
        event.type = Expose;
        XPutBackEvent(display, &event);
        XClearWindow(display, window);
        std::cout << "Simulated Expose event" << std::endl;
        if_change_light = 0;
    }
    if (keyBuffer[0] == 'w' || keyBuffer[0] == 'W')
    {
        eye.y = eye.y + 1;
        event.type = Expose;
        XPutBackEvent(display, &event);
        XClearWindow(display, window);
        std::cout << "Simulated Expose event" << std::endl;
        if_change_light = 0;
    }
    if (keyBuffer[0] == 's' || keyBuffer[0] == 'S')
    {
        eye.y = eye.y - 1;
        event.type = Expose;
        XPutBackEvent(display, &event);
        XClearWindow(display, window);
        std::cout << "Simulated Expose event" << std::endl;
        if_change_light = 0;
    }
    if (keyBuffer[0] == 'i' || keyBuffer[0] == 'I')
    {
        eye.z = eye.z + 1;
        event.type = Expose;
        XPutBackEvent(display, &event);
        XClearWindow(display, window);
        std::cout << "Simulated Expose event" << std::endl;
        if_change_light = 0;
    }
    if (keyBuffer[0] == 'k' || keyBuffer[0] == 'K')
    {
        eye.z = eye.z - 1;
        event.type = Expose;
        XPutBackEvent(display, &event);
        XClearWindow(display, window);
        std::cout << "Simulated Expose event" << std::endl;
        if_change_light = 0;
    }
}

void light_control(char keyBuffer[32])
{

    if (keyBuffer[0] == '4')
    {
        origin_light_dir.x = origin_light_dir.x - 1;
        event.type = Expose;
        XPutBackEvent(display, &event);
        XClearWindow(display, window);
        std::cout << "Simulated Expose event" << std::endl;
        if_change_light = 1;
    }
    if (keyBuffer[0] == '6')
    {
        origin_light_dir.x = origin_light_dir.x + 1;
        event.type = Expose;
        XPutBackEvent(display, &event);
        XClearWindow(display, window);
        std::cout << "Simulated Expose event" << std::endl;
        if_change_light = 1;
    }
    if (keyBuffer[0] == '8')
    {
        origin_light_dir.y = origin_light_dir.y + 1;
        event.type = Expose;
        XPutBackEvent(display, &event);
        XClearWindow(display, window);
        std::cout << "Simulated Expose event" << std::endl;
        if_change_light = 1;
    }
    if (keyBuffer[0] == '2')
    {
        origin_light_dir.y = origin_light_dir.y - 1;
        event.type = Expose;
        XPutBackEvent(display, &event);
        XClearWindow(display, window);
        std::cout << "Simulated Expose event" << std::endl;
        if_change_light = 1;
    }
    if (keyBuffer[0] == '+')
    {
        origin_light_dir.z = origin_light_dir.z + 1;
        event.type = Expose;
        XPutBackEvent(display, &event);
        XClearWindow(display, window);
        std::cout << "Simulated Expose event" << std::endl;
        if_change_light = 1;
    }
    if (keyBuffer[0] == '-')
    {
        origin_light_dir.z = origin_light_dir.z - 1;
        event.type = Expose;
        XPutBackEvent(display, &event);
        XClearWindow(display, window);
        std::cout << "Simulated Expose event" << std::endl;
        if_change_light = 1;
    }
}

struct PhongShader : public MyShader
{

    std::vector<std::array<Vec3f, 3>> node;
    virtual void vertex(int i)
    {
        node.resize(3); // 为 node 向量分配大小为 3 的内存空间
        for (int j = 0; j < 3; j++)
        {
            node[0][j] = m2v(MV * v2m(model->vert(model->face(i)[j]))); // 顶点坐标
            node[1][j] = textpos(model->vts(model->ft(i)[j]));          // 顶点纹理坐标
            node[2][j] = m2v(MV * v2m(model->vns(model->fn(i)[j])));    // 顶点法向量
        }
    }

    void fragment(float x, float y, Vec3f &c, std::vector<std::vector<double>> &zbuffer, TGAImage &image, TGAColor &color, Vec3f &light_dir)
    {

        float drawx = c.x * node[1][0].x + c.y * node[1][1].x + c.z * node[1][2].x;
        float drawy = c.x * node[1][0].y + c.y * node[1][1].y + c.z * node[1][2].y;

        Vec3f forz(x, y, zbuffer[x][y]);

        // forz = m2v(ViewPort*Projection*MV*(ViewPort*Projection*MV).inverse()*v2m(forz));
        forz = m2v(tobuffer * (ViewPort * Projection * MV).inverse() * v2m(forz));
        float shadow = .3 + .7 * (shadowbuffer[forz.x][forz.y] < forz.z + 20);

        color = model->texture.get(drawx * model->texture.get_width(), drawy * model->texture.get_height());
        nm = model->normal.get(drawx * model->normal.get_width(), drawy * model->normal.get_height());
        float sp = model->spec.get(drawx * model->spec.get_width(), drawy * model->spec.get_height()).raw[0];

        Vec3f norm;
        // 使用法线贴图
        if (using_nm == true)
        {
            norm[0] = nm.r / 255.f * 2.f - 1.f;
            norm[1] = nm.g / 255.f * 2.f - 1.f;
            norm[2] = nm.b / 255.f * 2.f - 1.f;
            norm = m2v(MV * v2m(norm)); // 法向量的空间变换
            norm.normalize();
        }
        // 使用切线贴图
        else if (using_tangent == true)
        {
            // 如果不使用法线贴图，使用插值法向量
            norm[0] = c.x * node[2][0].x + c.y * node[2][1].x + c.z * node[2][2].x;
            norm[1] = c.x * node[2][0].y + c.y * node[2][1].y + c.z * node[2][2].y;
            norm[2] = c.x * node[2][0].z + c.y * node[2][1].z + c.z * node[2][2].z;
            // norm = (node[0][1]- node[0][2])^(node[0][1] - node[0][0]);//使用三角形法向量
            norm.normalize();

            Matrix A(3, 3);
            A[0][0] = node[0][1].x - node[0][0].x;
            A[0][1] = node[0][1].y - node[0][0].y;
            A[0][2] = node[0][1].z - node[0][0].z;
            A[1][0] = node[0][2].x - node[0][0].x;
            A[1][1] = node[0][2].y - node[0][0].y;
            A[1][2] = node[0][2].z - node[0][0].z;
            A[2][0] = norm[0];
            A[2][1] = norm[1];
            A[2][2] = norm[2];
            A = A.inverse();

            Matrix b(3, 1);
            b[0][0] = node[1][1].x - node[1][0].x;
            b[1][0] = node[1][2].x - node[1][0].x;
            b[2][0] = 0.f;

            Matrix b2(3, 1);
            b2[0][0] = (-node[1][1].y + node[1][0].y); // 变换回原始贴图坐标
            b2[1][0] = (-node[1][2].y + node[1][0].y);
            b2[2][0] = 0.f;

            Matrix i(3, 1);
            i = A * b;
            Vec3f T;
            T.x = i[0][0];
            T.y = i[1][0];
            T.z = i[2][0];
            T.normalize();

            Matrix j(3, 1);
            j = A * b2;
            Vec3f B;
            B.x = j[0][0];
            B.y = j[1][0];
            B.z = j[2][0];
            B.normalize();

            Matrix TBN(3, 3);
            TBN[0][0] = T.x;
            TBN[0][1] = T.y;
            TBN[0][2] = T.z;
            TBN[1][0] = B.y;
            TBN[1][1] = B.y;
            TBN[1][2] = B.z;
            TBN[2][0] = norm[0];
            TBN[2][1] = norm[1];
            TBN[2][2] = norm[2];
            Matrix nm_(3, 1);
            nm_[0][0] = double(nm.r) / 255.f * 2.f - 1.f;
            nm_[1][0] = double(nm.g) / 255.f * 2.f - 1.f;
            nm_[2][0] = double(nm.b) / 255.f * 2.f - 1.f;

            nm_ = TBN.transpose() * nm_;
            norm[0] = nm_[0][0];
            norm[1] = nm_[1][0];
            norm[2] = nm_[2][0];
            norm.normalize();
        }
        else
        {
            // 如果不使用法线贴图，使用插值法向量
            norm[0] = c.x * node[2][0].x + c.y * node[2][1].x + c.z * node[2][2].x;
            norm[1] = c.x * node[2][0].y + c.y * node[2][1].y + c.z * node[2][2].y;
            norm[2] = c.x * node[2][0].z + c.y * node[2][1].z + c.z * node[2][2].z;
            // norm = (node[0][1]- node[0][2])^(node[0][1] - node[0][0]);//使用三角形法向量
            norm.normalize();

        }
        // 镜面反射
        Vec3f r = (norm * (norm * light_dir * 2.f) - light_dir).normalize();
        sp = pow(std::max(r.z, 0.f), sp);
        float diff = std::max(0.f, norm * light_dir);
        float cr = std::min(255., 10. + shadow * color.r * (1.2 * diff + .8 * sp));
        float cg = std::min(255., 10. + shadow * color.g * (1.2 * diff + .8 * sp));
        float cb = std::min(255., 10. + shadow * color.b * (1.2 * diff + .8 * sp));
        color = TGAColor(cr, cg, cb, 255);
    }
};

struct DepthShader : public MyShader
{

    std::vector<std::array<Vec3f, 3>> node;
    virtual void vertex(int i)
    {
        node.resize(3); // 为 node 向量分配大小为 3 的内存空间
        for (int j = 0; j < 3; j++)
        {
            node[0][j] = m2v(MV * v2m(model->vert(model->face(i)[j]))); // 顶点坐标
            // node[1][j] = textpos(model->vts(model->ft(i)[j]));                             // 顶点纹理坐标
            // node[2][j] = m2v(MV.inverse().transpose() * v2m(model->vns(model->fn(i)[j]))); // 顶点法向量
        }
    }

    void fragment(float x, float y, Vec3f &c, std::vector<std::vector<double>> &zbuffer, TGAImage &image, TGAColor &color, Vec3f &light_dir)
    {
        float z = zbuffer[x][y];
        // std::cout<<z<<std::endl;
        float cr = z;
        float cg = z;
        float cb = z;
        color = TGAColor(cr, cg, cb, 255);
    }
};

void Shadow_rendering(Model *model, DepthShader &Shader, std::vector<std::vector<double>> &buffer, TGAImage &image, Vec3f light_dir)
{
    for (int i = 0; i < model->nfaces(); i++)
    {
        Shader.vertex(i);
        triangle(Shader.node, Shader, buffer, image, light_dir);
    }
}

void Color_rendering(Model *model, PhongShader &Shader, std::vector<std::vector<double>> &buffer, TGAImage &image, Vec3f light_dir)
{
    for (int i = 0; i < model->nfaces(); i++)
    {
        Shader.vertex(i);
        triangle(Shader.node, Shader, buffer, image, light_dir);
    }
}

// 坐标变换链:Viewport * Projection * View * Model * v
int main(int argc, char **argv)
{
    XSelectInput(display, window, ExposureMask | KeyPressMask);
    XMapWindow(display, window);
    XFlush(display);

    // glows.read_tga_file("./tga/african_head_glow.tga");
    int defaultValue = -1 * 0x3f3f3f3f;
    shadowbuffer.resize(width + 1);
    for (int i = 0; i <= width + 1; i++)
        shadowbuffer[i].resize(height + 1, defaultValue);
    std::vector<std::vector<double>> zbuffer;
    zbuffer.resize(width + 1);
    for (int i = 0; i <= width + 1; i++)
        zbuffer[i].resize(height + 1, defaultValue); // 只会给新添加的元素赋值
    TGAImage depth(width, height, TGAImage::RGB);    // 深度图
    TGAImage image(width, height, TGAImage::RGB);

    while (1)
    {

        XNextEvent(display, &event);
        /* 绘制窗口或者重新绘制 */
        if (event.type == Expose)
        {
            model = new Model("./tga/african_head.obj");
            // std::cout << eye.x << ' ' << eye.y << ' ' << eye.z << std::endl;
            // std::cout << origin_light_dir.x << ' ' << origin_light_dir.y << ' ' << origin_light_dir.z << std::endl;
            Vec3f light_dir = origin_light_dir;

            for (int i = 0; i <= width; i++)
                for (int j = 0; j <= height; j++)
                    zbuffer[i][j] = defaultValue;

            if (if_change_light == 1)
            {
                for (int i = 0; i <= width; i++)
                    for (int j = 0; j <= height; j++)
                        shadowbuffer[i][j] = defaultValue;

                light_dir.normalize();
                ModelView(light_dir, center, Vec3f(0, 1, 0));
                viewport(width / 8, height / 8, width * 3 / 4, height * 3 / 4);
                projection(0);

                DepthShader Dshader;
                std::cout << "Rendering shadow!" << std::endl;
                Shadow_rendering(model, Dshader, shadowbuffer, depth, light_dir);
                tobuffer = ViewPort * Projection * MV;
                if_change_light = 0;
            }

            ModelView(eye, center, Vec3f(0, 1, 0));
            viewport(width / 8, height / 8, width * 3 / 4, height * 3 / 4);
            projection(-1.f / (eye - center).norm());

            light_dir = m2v(MV * v2m(light_dir));
            light_dir.normalize();
            PhongShader shader;
            std::cout << "Rendering color!" << std::endl;
            Color_rendering(model, shader, zbuffer, image, light_dir);

            // model = new Model("./tga/african_head_eye_inner.obj");
            // std::cout << "Rendering color!" << std::endl;
            // Color_rendering(model, shader, zbuffer, image, light_dir);

            std::cout << "Rendering Finish!" << std::endl;

            image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
            show_image(image);
            image.clear();
            depth.clear();
            // image.write_tga_file("output.tga");
            // XDestroyImage(xImage);
            // 刷新窗口
        }

        if (event.type == KeyPress)
        {
            char keyBuffer[32];
            XLookupString(&event.xkey, keyBuffer, sizeof(keyBuffer), NULL, NULL);
            light_control(keyBuffer);
            camera_control(keyBuffer);
        }
    }

    // 释放图形上下文
    XFreeGC(display, gc);
    // 关闭显示连接
    XCloseDisplay(display);

    return 0;
}
