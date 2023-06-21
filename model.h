#ifndef __MODEL_H__
#define __MODEL_H__

#include <vector>
#include "math.h"
#include "tgaimage.h"
#include <string>

class Model {
private:
	std::vector<Vec3f> verts_;
	std::vector<Vec3f> vts_;
	std::vector<Vec3f> vns_;
	std::vector<std::vector<int> > faces_;
	std::vector<std::vector<int> > ft_;
	std::vector<std::vector<int> > vn_;
public:
	Model(std::string filename);
	~Model();
	int nverts();
	int nvts();
	int nfaces();
	int nvn();
	Vec3f vert(int i);
	Vec3f vts(int i);
	Vec3f vns(int i);
	std::vector<int> face(int idx);
	std::vector<int> ft(int idx);
	std::vector<int> fn(int idx);
	TGAImage texture;
    TGAImage normal;
    TGAImage spec;

};

#endif //__MODEL_H__