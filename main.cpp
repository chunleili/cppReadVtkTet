#include <iostream>
#include <fstream>
#include <algorithm>
#include <list>
#include <tuple>
#include <vector>
#include <filesystem>
#include "MeshIO.hpp"

using vec3f = std::array<float,3>;
using vec4i = std::array<int,4>;
using vec3i = std::array<int,3>;

struct ReadVtkTet{
private:
    void extractSurf();
    std::vector<vec3f> pos;
    std::vector<vec4i> quads;
    std::vector<vec3i> surfs;
public:
    void apply() {
    std::string path = "E:\\codes\\read_vtk_tet\\cppReadVtkTet\\Armadillo13K.vtk";

    zs::Mesh<float, 3, int, 4> tet;
    read_tet_mesh_vtk(path, tet);
    const auto numVerts = tet.nodes.size();
    const auto numEles = tet.elems.size();
    pos.resize(numVerts);
    quads.resize(numEles);

    for (int i = 0; i < numVerts; i++)
    {
        pos[i] = tet.nodes[i];
    }
    for (int i = 0; i < numEles; i++)
    {
        quads[i] = tet.elems[i];
    }

    extractSurf();
  }
};



void ReadVtkTet::extractSurf()
{
    int numTets = quads.size();
    int numFaces = quads.size()*4;
    
    using vec4i = std::array<int,4>;

    //list_faces
    std::vector<vec4i> faces;
    for (int i = 0; i < quads.size(); i++)
    {
        std::array<int,4> tet=quads[i];

        std::sort(tet.begin(),tet.end());
        
        int t0 = tet[0];
        int t1 = tet[1];
        int t2 = tet[2];
        int t3 = tet[3];

        vec4i f0{t0, t1, t2, i};
        vec4i f1{t0, t1, t3, i};
        vec4i f2{t0, t2, t3, i};
        vec4i f3{t1, t2, t3, i};

        faces.push_back(f0);
        faces.push_back(f1);
        faces.push_back(f2);
        faces.push_back(f3);
    }

    auto myLess = [](auto a, auto b){
        return std::tie(a[0], a[1], a[2]) < std::tie(b[0], b[1], b[2]);
    };

    //sort faces
    std::sort(faces.begin(),faces.end(), myLess);


    /* -------------------------------------------------------------------------- */
    /*                       use list to remove shared faces                      */
    /* -------------------------------------------------------------------------- */
    auto myEqual = [](auto a, auto b){
        if((a[0]==b[0])&&(a[1]==b[1])&&(a[2]==b[2])) return true; 
        else return false;};
    //copy to a list
    std::list<vec4i> facelist;
    for(auto it:faces)
        facelist.push_back(it);
    
    //remove shared faces from the list
    auto f_prev = facelist.begin();
    auto f = facelist.begin();  f++;
    for(;f!=facelist.end();)
    {
        if(myEqual(*f, *f_prev))
        {
            f++; //move to the next
            if(f==facelist.end())//if f move to the end(), should break
            {
                facelist.erase(f_prev,f);
                break;
            }
            f_prev = facelist.erase(f_prev,f); //return the next
            f++; // move f to the next of next
        }
        else
        {
            f_prev = f;
            f++;
        }
    }

    //recontruct the surf with orders
    for(auto x:facelist)
    {
        int tetId = x[3];
        std::array<int, 4> vert = quads[tetId];
        std::array<bool, 4> hasVert={false, false, false, false};
        
        for (int j = 0; j < 4; j++)
            for (int k = 0; k < 3; k++)
                if(x[k] == vert[j])
                    hasVert[j] = true;
            
        if (hasVert[0] &&  hasVert[2] && hasVert[1])
            surfs.push_back(vec3i{vert[0],vert[2],vert[1]}); 

        if (hasVert[0] &&  hasVert[3] && hasVert[2])
            surfs.push_back(vec3i{vert[0],vert[3],vert[2]});

        if (hasVert[0] &&  hasVert[1] && hasVert[3])
            surfs.push_back(vec3i{vert[0],vert[1],vert[3]});

        if (hasVert[1] &&  hasVert[2] && hasVert[3])
            surfs.push_back(vec3i{vert[1],vert[2],vert[3]});
    }
    
    std::cout<<"before extractSurf numFaces: "<<numFaces<<std::endl;
    std::cout<<"after extractSurf numSurfs: "<<surfs.size()<<std::endl;
    
    std::ofstream f_surf, f_pos;
    f_surf.open("surf.txt");
    for(const auto& f:surfs)
    {
        for(const auto& x:f)
            f_surf<<x<<"\t";
        f_surf<<"\n";
    }

    f_pos.open("pos.txt");
    for(const auto& f:pos)
    {
        for(const auto& x:f)
            f_pos<<x<<"\t";
        f_pos<<"\n";
    }
    f_pos.close();
}


int main(){
    ReadVtkTet a;
    a.apply();
}