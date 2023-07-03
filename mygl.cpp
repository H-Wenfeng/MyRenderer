#include <cmath>
#include <limits>
#include <cstdlib>
#include <iostream>
#include "math.h"
#include "mygl.h"

Matrix Projection = Matrix::identity(4);
Matrix MV = Matrix::identity(4);
Matrix ViewPort = Matrix::identity(4);

Vec3f back_board(0.4, 0.4, 0.4);
TGAImage env(14576, 7100, TGAImage::RGB);

Material::Material(Vec3f diffusion_color, float diffusion_coefficient, float specular_coefficient, float reflect_coefficient, float refract_coefficient, float specular_exponent, float refract_exponent) : diffusion_color(diffusion_color), diffusion_coefficient(diffusion_coefficient),																																																			specular_coefficient(specular_coefficient), reflect_coefficient(reflect_coefficient), refract_coefficient(refract_coefficient), specular_exponent(specular_exponent), refract_exponent(refract_exponent) {}
Sphere::Sphere(Vec3f c, float radius, Material material) : center(c), r(radius), material(material) {}

bool Sphere::intersect(const Vec3f &light, const Vec3f &camera, Vec3f &intersectionPoint)
{
	Vec3f oc = camera - center;
	float a = light * light;
	float b = 2.0f * (oc * light);
	float c = oc * oc - r * r;

	float discriminant = b * b - 4 * a * c;
	if (discriminant < 0)
	{
		// 没有交点
		return false;
	}

	// 计算两个交点的参数化值
	float t1 = (-b + std::sqrt(discriminant)) / (2.0f * a);
	float t2 = (-b - std::sqrt(discriminant)) / (2.0f * a);
	// 选择更接近射线原点的交点
	float t;
	if (t1 <= 0.f)
		t = t2;
	else if (t2 <= 0.f)
		t = t1;
	else
		t = std::min(t1, t2);
	// 计算交点的位置向量
	if (t <= 0.f)
		return false;
	intersectionPoint = camera + light * t;

	return true;
}

void Ray_Render(Vec3f &camera, std::vector<Sphere> &Spheres, TGAImage &image, std::vector<std::vector<double>> &zbuffer, std::vector<Light> &lights)
{
	float FOV = M_PI/3.;
	float width = image.get_width();
	float height = image.get_height();
	env.read_tga_file("./tga/galaxy3.tga");
	for (size_t p = 0; p < Spheres.size(); p++)
	{
		Sphere sph = Spheres[p];
		for (int i = 0; i <= width; i++)
			for (int j = 0; j <= height; j++)
			{
				float x = (2 * (i + 0.5) / width - 1) * tan(FOV / 2.) * width / height; // x坐标：将像素的水平位置映射到范围[-1, 1]内，然后乘以tan(fov/2)得到相对于宽度的偏移量。
				float y = (2 * (j + 0.5) / height - 1) * tan(FOV / 2.);					// 将像素的水平位置(i)加上0.5是为了将光线的起始点定位在像素中心，而不是像素的边界。这是因为光线的方向是从摄像机通过像素中心的，如果起始点在像素边界，可能会导致采样不准确或产生锯齿状的图像。通过将像素位置加上0.5，可以将起始点移动到像素中心，从而更准确地采样光线的方向。这种技巧被称为"偏移0.5"，常用于图像渲染和采样过程中。
				Vec3f light(x, y, camera.z - 1);
				
				Vec3f interselectionP(0, 0, -0x3f3f3f3f);
				Vec3f PinterselectionP(0, 0, -0x3f3f3f3f);
				
				light.normalize();
				float theta = acos(light.y);
				float phi = atan2(light.z, light.x);
				float normalizedTheta = std::max(0., std::min(1., theta / M_PI)) * env.get_height();
				float normalizedPhi = std::max(0., std::min(1., (phi + M_PI) / (2 * M_PI))) * env.get_width();
				back_board = Vec3f(env.get(normalizedPhi, normalizedTheta).r / 255., env.get(normalizedPhi, normalizedTheta).g / 255., env.get(normalizedPhi, normalizedTheta).b / 255.);
				// std::cout<<normalizedTheta<<" "<<normalizedPhi<<" "<<back_board.z<<std::endl;

				if (sph.intersect(light, camera, interselectionP) == true || plane(light, camera, PinterselectionP) == true)
				{
					int z = std::max(interselectionP.z, PinterselectionP.z);
					if (zbuffer[i][j] < z)
					{
						zbuffer[i][j] = z;
						int times = 4;
						Vec3f color = cast_ray(light, camera, p, Spheres, lights, times);
						image.set(i, j, TGAColor(color.x * 255., color.y * 255., color.z * 255., 255));
					}	
						
				}
				else if (zbuffer[i][j] == -0x3f3f3f3f)
				{
					zbuffer[i][j] = -0x3f3f3f3f + 1;
					image.set(i, j, TGAColor(back_board.x * 255., back_board.y * 255., back_board.z * 255., 255));
				}
			}
	}
}

Vec3f reflect(Vec3f &norm, Vec3f &dir, Vec3f &interselectionP, size_t p, std::vector<Sphere> &Spheres, std::vector<Light> &lights, int &times)
{
	Vec3f normt = norm;
	float costt = normt * dir * -1.f;
	if (costt < 0)
	{
		normt = normt * -1;
		costt = normt * dir * -1.f;
	}
	Vec3f reflect_dir = dir + normt * 2.f * costt;
	reflect_dir.normalize();
	Vec3f reflect_orig = interselectionP;
	Vec3f reflect_color;
	if (times >= 0)
	{
		Vec3f tmp;
		int dis = 0x3f3f3f3f;
		int which = -1;
		for (size_t q = 0; q < Spheres.size(); q++)
		{
			if (Spheres[q].intersect(reflect_dir, reflect_orig, tmp) == true && p != q)
				if ((tmp - reflect_orig).norm() < dis)
				{
					which = q;
					dis = (tmp - reflect_orig).norm();
				}
		}		
		if (which != -1 || plane(reflect_dir, reflect_orig, tmp))
			return cast_ray(reflect_dir, reflect_orig, which, Spheres, lights, times);
		else
		{
				float theta = acos(reflect_dir.y);
				float phi = atan2(reflect_dir.z, reflect_dir.x);
				float normalizedTheta = std::max(0., std::min(1., theta / M_PI)) * env.get_height();
				float normalizedPhi = std::max(0., std::min(1., (phi + M_PI) / (2 * M_PI))) * env.get_width();
				back_board = Vec3f(env.get(normalizedPhi, normalizedTheta).r / 255., env.get(normalizedPhi, normalizedTheta).g / 255., env.get(normalizedPhi, normalizedTheta).b / 255.);
				return back_board;
		}
	}
	return back_board;
}
Vec3f refract(Vec3f &norm, Vec3f &dir, Vec3f &interselectionP, size_t p, std::vector<Sphere> &Spheres, std::vector<Light> &lights, int &times)
{

	Vec3f norm2 = norm;
	float sin = sqrt(1 - pow((dir * norm2) / (norm2.norm() * dir.norm()), 2));
	float theta_crit = 1 / Spheres[p].material.refract_exponent;
	Vec3f refract_origin;
	Vec3f refract_dir;

	float cost1 = norm2 * dir * -1.f;
	if (cost1 < 0 && sin > 1 / theta_crit)
		refract_dir = Vec3f(0, 0, 0);
	else
	{
		if (cost1 < 0)
		{
			norm2 = norm2 * -1;
			cost1 = norm2 * dir * -1.f;
			theta_crit = 1 / theta_crit;
		}
		float cost2 = sqrt(1 - theta_crit * theta_crit * (1 - cost1 * cost1));

		refract_dir = dir * theta_crit + norm2 * (theta_crit * cost1 - cost2);
	}
	refract_dir.normalize();
	Vec3f interselectionP2;
	refract_origin = refract_dir * norm2 < 0 ? interselectionP - norm2 * 1e-3 : interselectionP + norm2 * 1e-3;
	// refract_origin = interselectionP2;
	Spheres[p].intersect(refract_dir, refract_origin, interselectionP2);
	norm2 = (interselectionP2 - Spheres[p].center).normalize();
	sin = sqrt(1 - pow((refract_dir * norm2) / (norm2.norm() * refract_dir.norm()), 2));
	cost1 = norm2 * refract_dir * -1.f;
	if (cost1 < 0 && sin > 1 / theta_crit)
		refract_dir = Vec3f(0, 0, 0);
	else
	{
		if (cost1 < 0)
		{
			norm2 = norm2 * -1;
			cost1 = norm2 * refract_dir * -1.f;
			theta_crit = 1 / theta_crit;
		}
		float cost2 = sqrt(1 - theta_crit * theta_crit * (1 - cost1 * cost1));

		refract_dir = refract_dir * theta_crit + norm2 * (theta_crit * cost1 - cost2);
	}
	refract_origin = interselectionP2;
	refract_dir.normalize();

	if (times >= 0)
	{
		Vec3f tmp;
		int dis = 0x3f3f3f3f;
		int which = -1;
		for (size_t q = 0; q < Spheres.size(); q++)
		{
			if (Spheres[q].intersect(refract_dir, refract_origin, tmp) == true && p != q)
				if ((tmp - refract_origin).norm() < dis)
				{
					which = q;
					dis = (tmp - refract_origin).norm();
				}
		}
		if (which != -1 || plane(refract_dir, refract_origin, tmp))
		{
			return cast_ray(refract_dir, refract_origin, which, Spheres, lights, times);
		}
		else
		{
			float theta = acos(refract_dir.y);
			float phi = atan2(refract_dir.z, refract_dir.x);
			float normalizedTheta = std::max(0., std::min(1., theta / M_PI)) * env.get_height();
			float normalizedPhi = std::max(0., std::min(1., (phi + M_PI) / (2 * M_PI))) * env.get_width();
			back_board = Vec3f(env.get(normalizedPhi, normalizedTheta).r / 255., env.get(normalizedPhi, normalizedTheta).g / 255., env.get(normalizedPhi, normalizedTheta).b / 255.);
			return back_board;

		}	
	}
	return back_board;
}

bool plane(Vec3f &dir, Vec3f &orig, Vec3f &intersectionPoint)
{
	if(fabs(dir.y) > 1e-3)
	{
		float planed = -(orig.y+4)/dir.y;
		// (orig.y - (-4)) / dir.y;
		intersectionPoint = orig + dir*planed;
		if (planed>0 && fabs(intersectionPoint.x)<8 && intersectionPoint.z<-10 && intersectionPoint.z>-30 ) 		
			return true;
	}
	return false;
}

Vec3f cast_ray(Vec3f &dir, Vec3f &origin, size_t p, std::vector<Sphere> &Spheres, std::vector<Light> &lights, int &times)
{
	times--;
	Vec3f interselectionP(0, 0, 0x3f3f3f3f);
	Vec3f PinterselectionP(0, 0, 0x3f3f3f3f);
	plane(dir, origin, PinterselectionP);

	Sphere sph = Spheres[p];
	sph.intersect(dir, origin, interselectionP);

	float dsphere = (interselectionP - origin).norm();
	float dplane = (PinterselectionP - origin).norm();

	if(plane(dir, origin, PinterselectionP) !=true || dsphere < dplane)
	{

		Vec3f norm = interselectionP - sph.center;
		norm.normalize();
		Vec3f reflect_color = reflect(norm, dir, interselectionP, p, Spheres, lights, times);
		Vec3f refract_color = refract(norm, dir, interselectionP, p, Spheres, lights, times);
		float diff = 0;
		float spec = 0;
		for (size_t k = 0; k < lights.size(); k++)
		{
			Vec3f light_dir = lights[k].position - interselectionP;
			light_dir.normalize();
			int flag = 0;
			for (size_t g = 0; g < Spheres.size(); g++)
			{
				if (g != p)
				{
					Vec3f tmp;
					Vec3f ndir = lights[k].position - interselectionP;
					float dis = ndir.norm();
					if (Spheres[g].intersect(ndir.normalize(), interselectionP, tmp) == true && (dis >= (lights[k].position - tmp).norm()))
						flag = 1;
				}
			}
			if (flag == 1)
				continue;
			Vec3f r = (norm * (norm * light_dir) * 2.f - light_dir).normalize();
			float sp = sph.material.specular_exponent;
			Vec3f v = origin - interselectionP;
			v.normalize();
			spec += std::max(0., pow((r * v) / (r.norm() * v.norm()), sp)) * lights[k].intensity;
			diff += std::max(0.f, norm * light_dir) * lights[k].intensity;
		}
		float cr = std::min(1., sph.material.diffusion_color.x * diff * sph.material.diffusion_coefficient + spec * sph.material.specular_coefficient * 1. + reflect_color.raw[0] * sph.material.reflect_coefficient + refract_color.raw[0] * sph.material.refract_coefficient);
		float cg = std::min(1., sph.material.diffusion_color.y * diff * sph.material.diffusion_coefficient + spec * sph.material.specular_coefficient * 1. + reflect_color.raw[1] * sph.material.reflect_coefficient + refract_color.raw[1] * sph.material.refract_coefficient);
		float cb = std::min(1., sph.material.diffusion_color.z * diff * sph.material.diffusion_coefficient + spec * sph.material.specular_coefficient * 1. + reflect_color.raw[2] * sph.material.reflect_coefficient + refract_color.raw[2] * sph.material.refract_coefficient);
		return Vec3f(cr, cg, cb);
	}
	
	else if(plane(dir, origin, PinterselectionP) == true)
	{
		// std::cout<<PinterselectionP.x<<" "<<PinterselectionP.y<<" "<<PinterselectionP.z<<std::endl;
		Vec3f norm(0, 1, 0);			
		float diff = 0;
		for (size_t k = 0; k < lights.size(); k++)
			{
				Vec3f light_dir = lights[k].position - PinterselectionP;
				light_dir.normalize();
				int flag = 0;

				for (size_t g = 0; g < Spheres.size(); g++)
				{				
					{
						Vec3f tmp;
						Vec3f ndir = lights[k].position - PinterselectionP;
						float dis = ndir.norm();
						if (Spheres[g].intersect(ndir.normalize(), PinterselectionP, tmp) == true && (dis >= (lights[k].position - tmp).norm()))
							flag = 1;
					}
				}
				if (flag == 1)
					continue;
				diff += std::max(0.f, norm * light_dir) * lights[k].intensity;
			}
			Vec3f diffusion;
			diffusion = (int(.5*PinterselectionP.x+1000) + int(.5*PinterselectionP.z)) & 1 ? Vec3f(1.0, 1.0, 1.0) : Vec3f(.9, .4, .1);
			float cr = diffusion.x * diff *0.2;
			float cg = diffusion.y * diff *0.2;
			float cb = diffusion.z * diff *0.2;
			return Vec3f(cr, cg, cb);
	}
	return back_board;
}

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
	// 保护边界
	float lbound = std::max(0.f, std::min(pts[0].x, std::min(pts[1].x, pts[2].x)));
	float rbound = std::min(image.get_width() * 1.f, std::max(pts[0].x, std::max(pts[1].x, pts[2].x)));
	float ubound = std::max(0.f, std::min(pts[0].y, std::min(pts[1].y, pts[2].y)));
	float bbound = std::min(image.get_height() * 1.f, std::max(pts[0].y, std::max(pts[1].y, pts[2].y)));
	for (int x = lbound; x <= rbound; x++)
		for (int y = ubound; y <= bbound; y++)
		{

			Vec3f c = Centroid(pts, x, y);
			float v = c.x;
			float u = c.y;
			float w = c.z;
			float z = pts[0].z * v + pts[1].z * u + pts[2].z * w; // 插值Z坐标

			if ((v >= 0.) && (u >= 0.) && (w >= 0) && zbuffer[x][y] < z) // 该点重心在三角形内，并且该深度没有被渲染
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
	m[2][3] = 0.5f; // 缩放深度

	m[0][0] = w / 2.0f; // 缩放x轴
	m[1][1] = h / 2.0f; // 一般调整这两个缩放来测试新模型
	m[2][2] = 0.5f;
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
