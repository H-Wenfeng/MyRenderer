#include <cmath>
#include <limits>
#include <cstdlib>
#include <iostream>
#include "geometry.h"
#include "mygl.h"

Matrix Projection = Matrix::identity(4);
Matrix MV = Matrix::identity(4);
Matrix ViewPort = Matrix::identity(4);

void line(Vec3f p0, Vec3f p1, TGAImage &image, TGAColor color)
{
	float x1 = p0.x;
	float x2 = p1.x;
	float y1 = p0.y;
	float y2 = p1.y;
	int trans = 0;
	if (fabs(x1 - x2) < fabs(y1 - y2))
	{
		std::swap(x1, y1);
		std::swap(x2, y2);
		trans = 1;
	}
	if (x1 > x2)
	{
		std::swap(y1, y2);
		std::swap(x1, x2);
	}
	if (((x2 - x1) == 0) && ((y1 - y2) == 0))
		image.set(x1, y1, color);
	else if ((y2 - y1) == 0)
		for (int i = x1; i <= x2; i += 1)
			if (trans)
				image.set(y1, i, color);
			else
				image.set(i, y1, color);
	else
	{

		float slope = (y2 - y1) / (x2 - x1);

		for (float i = x1; i <= x2; i += 1)
		{
			float drawx = i;
			float drawy = y1 + (i - x1) * slope;
			int tmp = drawy;
			if (drawy - tmp >= 0.5)
				tmp++;
			if (trans == 0)
				image.set(drawx, tmp, color);
			else
				image.set(tmp, drawx, color);
		}
	}
}

void triangle2(Vec2i t0, Vec2i t1, Vec2i t2, TGAImage &image, float intensity, TGAColor color)
{ // 上下两半给三角形着色
	if (t0.y > t1.y)
		std::swap(t0, t1);
	if (t0.y > t2.y)
		std::swap(t0, t2);
	if (t1.y > t2.y)
		std::swap(t1, t2);
	float halfx = t1.x;
	float halfy = t1.y;
	float a = t1.y - t0.y;
	float b = t2.y - t1.y;
	float px = (a * t2.x + b * t0.x) / (a + b);
	for (float y = t1.y; y <= t2.y; y += 1)
	{
		float lbound = px + (y - halfy) * ((t2.x - px) / (t2.y - halfy));
		float rbound = halfx + (y - halfy) * ((t2.x - halfx) / (t2.y - halfy));

		for (float drawx = std::min(lbound, rbound); drawx <= std::max(lbound, rbound); drawx += 1)
			image.set(drawx, y, color);
	}

	for (float y = t0.y; y <= t1.y; y += 1)
	{
		float lbound = px + (y - halfy) * ((t0.x - px) / (t0.y - halfy));
		float rbound = halfx + (y - halfy) * ((t0.x - halfx) / (t0.y - halfy));

		for (float drawx = std::min(lbound, rbound); drawx <= std::max(lbound, rbound); drawx += 1)
			image.set(drawx, y, color);
	}
}

Vec3f Centroid(std::array<Vec3f, 3> pts, float x, float y)
{
	// 	硬解求重心坐标
	// 	二元方程组：
	//	ua + vb + c =0
	// 	ud + ve + f = 0
	Vec3f c;
	c.x = ((pts[2].x - x) * (pts[1].y - pts[2].y) - (pts[1].x - pts[2].x) * (pts[2].y - y)) / ((pts[1].x - pts[2].x) * (pts[0].y - pts[2].y) - (pts[1].y - pts[2].y) * (pts[0].x - pts[2].x));
	c.y = -1 * (pts[2].x - x + (pts[0].x - pts[2].x) * c.x) / (pts[1].x - pts[2].x);
	c.z = 1. - c.x - c.y;

	return c;
}

Vec3f barycentric(Vec3f A, Vec3f B, Vec3f C, Vec3f P)
{
	// 用法向量求
	Vec3f s[2];
	for (int i = 2; i--;)
	{
		s[i][0] = C[i] - A[i];
		s[i][1] = B[i] - A[i];
		s[i][2] = A[i] - P[i];
	}
	Vec3f u = s[0] ^ s[1];
	if (std::abs(u[2]) > 1e-2) 
		return Vec3f(1.f - (u.x + u.y) / u.z, u.y / u.z, u.x / u.z);
	return Vec3f(-1, 1, 1);
}

// 光栅化,使用引用传参以修改zbuffer
void triangle(std::vector<std::array<Vec3f, 3>> node, MyShader &Shader, std::vector<std::vector<double>> &zbuffer, TGAImage &image, Vec3f intensity)
{
	std::array<Vec3f, 3> pts = node[0];
	pts[0] = m2v(ViewPort * Projection * v2m(pts[0]));
	pts[1] = m2v(ViewPort * Projection * v2m(pts[1]));
	pts[2] = m2v(ViewPort * Projection * v2m(pts[2]));

	float lbound = std::min(pts[0].x, std::min(pts[1].x, pts[2].x));
	float rbound = std::max(pts[0].x, std::max(pts[1].x, pts[2].x));
	float ubound = std::min(pts[0].y, std::min(pts[1].y, pts[2].y));
	float bbound = std::max(pts[0].y, std::max(pts[1].y, pts[2].y));
	for (int x = lbound; x <= rbound; x++)
		for (int y = ubound; y <= bbound; y++)
		{

			Vec3f c = Centroid(pts, x, y);
			float v = c.x;
			float u = c.y;
			float w = c.z;
			float z = pts[0].z * v + pts[1].z * u + pts[2].z * w; // 插值Z坐标

			if ((v >= 0.) && (u >= 0.) && (w >= 0) && zbuffer[x][y] < z) //该点重心在三角形内，并且该深度没有被渲染
			{
				zbuffer[x][y] = z;
				TGAColor color;
				Shader.fragment(x, y, c, zbuffer, image, color, intensity);				
				image.set(x, y, color);
				
			}
		}
}

Vec3f textpos(Vec3f v)
{
	return Vec3f(float(v.x), 1. - float(v.y), v.z); // 因为贴图是从左上角开始索引，所以我们转换成左下角
}

Vec3f m2v(Matrix m)
{

	return Vec3f(m[0][0] / m[3][0], m[1][0] / m[3][0], m[2][0] / m[3][0]); // 齐次坐标->笛卡尔坐标
}

Matrix v2m(Vec3f v)
{
	/*
	齐次坐标
	| 1  0  0  2 |
	| 0  1  0  3 |
	| 0  0  1  4 |
	| 0  0  0  1 |
	*/
	Matrix m = Matrix::identity(4);
	m[0][0] = v.x; // 笛卡尔坐标->齐次坐标
	m[1][0] = v.y;
	m[2][0] = v.z;
	m[3][0] = 1.f;
	return m;
}
// 视口变换
void viewport(int x, int y, int w, int h)
{
	Matrix m = Matrix::identity(4); // 创建一个单位矩阵

	// 视口变换的缩放和平移操作
	m[0][3] = x + w / 2.0f; // 原点平移
	m[1][3] = y + h / 2.0f;
	m[2][3] = 2000.f / 2.0f; // 缩放深度

	m[0][0] = w / 2.0f; // 缩放x轴
	m[1][1] = h / 2.0f;
	m[2][2] = 2000.f / 2.0f;
	ViewPort = m;
}

Matrix zoom(float factor)
{
	/*
	| 1     0      0     0 |
	| 0  0.7071 -0.7071  0 |
	| 0  0.7071  0.7071  0 |
	| 0     0      0     1 |
	*/
	Matrix Z = Matrix::identity(4);
	Z[0][0] = Z[1][1] = Z[2][2] = factor;
	return Z;
}

Matrix rotation_x(float cosangle, float sinangle)
{
	Matrix R = Matrix::identity(4);
	R[1][1] = R[2][2] = cosangle;
	R[1][2] = -sinangle;
	R[2][1] = sinangle;
	return R;
}

Matrix rotation_y(float cosangle, float sinangle)
{
	Matrix R = Matrix::identity(4);
	R[0][0] = R[2][2] = cosangle;
	R[0][2] = sinangle;
	R[2][0] = -sinangle;
	return R;
}

Matrix rotation_z(float cosangle, float sinangle)
{
	/*
	| 0.866 -0.5    0    0 |
	|  0.5   0.866  0    0 |
	|   0     0     1    0 |
	|   0     0     0    1 |
	*/
	Matrix R = Matrix::identity(4);
	R[0][0] = R[1][1] = cosangle;
	R[0][1] = -sinangle;
	R[1][0] = sinangle;
	return R;
}

void ModelView(Vec3f eye, Vec3f center, Vec3f up)
{

	MV = Matrix::identity(4);
	// 变换相机的视角，变换坐标系
	/*
	|x'|        |x|	  |	   |
	|y'| = M^-1(|y| - |Oxyz|）
	|z'|	    |z|   |    |
	*/
	// z是从视角指向新的原点
	Vec3f z = (eye - center).normalize(); // 将3个坐标轴归一化成单位向量
	Vec3f x = (up ^ z).normalize();		  // 通过叉乘求法向量,这里被重载成了^
	Vec3f y = (z ^ x).normalize();
	Matrix Minv = Matrix::identity(4);
	for (int i = 0; i < 3; i++)
	{
		Minv[0][i] = x[i];
		Minv[1][i] = y[i];
		Minv[2][i] = z[i];
		Minv[i][3] = -center[i];
	}
	MV = Minv;
	// 现在只需要求出M^-1
	/*
	如果我们的变换矩阵 M 是均匀缩放、旋转和平移（欧几里德空间的等距）的组合，
	则 M 等于它的逆转置，因为在这种情况下逆和转置相互抵消。但是由于我们的矩
	阵包括透视变形，通常这个技巧没有帮助。*/
}

void projection(float v)
{
	Projection = Matrix::identity(4);
	Projection[3][2] = v;
}