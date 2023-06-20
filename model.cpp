#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include "model.h"

Model::Model(std::string filename) : verts_(), faces_()
{
    std::ifstream in;
    in.open(filename, std::ifstream::in);
    if (in.fail())
        return;
    std::string line;
    while (!in.eof())
    {
        std::getline(in, line);
        std::istringstream iss(line.c_str());
        char trash;
        if (!line.compare(0, 2, "v "))
        {
            iss >> trash;
            Vec3f v;
            for (int i = 0; i < 3; i++)
                iss >> v.raw[i];
            verts_.push_back(v);
        }
        else if (!line.compare(0, 2, "f "))
        {
            std::vector<int> f;
            std::vector<int> ft;
            std::vector<int> vn;
            int idx;
            int secondidx;
            int thirdidx;
            iss >> trash; // f 448/437/448 449/438/449 450/439/450
            while (iss >> idx >> trash >> secondidx >> trash >> thirdidx)
            {
                idx--; // in wavefront obj all indices start at 1, not zero
                f.push_back(idx);
                // std::cout<<idx<<' ';
                secondidx--;
                ft.push_back(secondidx);
                thirdidx--;
                vn.push_back(thirdidx);
                // std::cout<<secondidx<<' ';
            }
            // std::cout<<std::endl;
            ft_.push_back(ft);
            faces_.push_back(f);
            vn_.push_back(vn);
        }
        else if (!line.compare(0, 3, "vt "))
        {

            iss >> trash >> trash;
            Vec3f vt;
            // std::cout<<trash<<' ';
            for (int i = 0; i < 3; i++)
            {
                iss >> vt.raw[i];
                // std::cout<<vt.raw[i]<<' ';
            }
            // std::cout<<std::endl;
            vts_.push_back(vt);
        }
        else if (!line.compare(0, 3, "vn "))
        {

            iss >> trash >> trash;
            Vec3f norm;
            // std::cout<<trash<<' ';
            for (int i = 0; i < 3; i++)
            {
                iss >> norm.raw[i];
                // std::cout<<vt.raw[i]<<' ';
            }
            // std::cout<<std::endl;
            vns_.push_back(norm.normalize());
        }
    }
    std::cout << "# v# " << verts_.size() << " f# " << faces_.size() << " vt# " << vts_.size() << " vn# " << vns_.size() << std::endl;
    std::string filehead = filename.replace(filename.length() - 4, 4, "");
    texture.read_tga_file(filehead + "_diffuse.tga"); // normal.read_tga_file("./tga/african_head_nm.tga");
    normal.read_tga_file(filehead + "_nm.tga");
    spec.read_tga_file(filehead + "_spec.tga");
}

Model::~Model()
{
}

int Model::nverts()
{
    return (int)verts_.size();
}
int Model::nvts()
{
    return (int)vts_.size();
}
int Model::nfaces()
{
    return (int)faces_.size();
}

int Model::nvn()
{
    return (int)vn_.size();
}

std::vector<int> Model::face(int idx)
{
    return faces_[idx];
}

std::vector<int> Model::ft(int idx)
{
    return ft_[idx];
}

std::vector<int> Model::fn(int idx)
{
    return vn_[idx];
}

Vec3f Model::vert(int i)
{
    return verts_[i];
}

Vec3f Model::vts(int i)
{
    return vts_[i];
}

Vec3f Model::vns(int i)
{
    return vns_[i];
}
