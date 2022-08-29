#include <iostream>
#include <fstream>
#include <algorithm>
#include <list>
#include <set>
#include <tuple>
#include <vector>
#include <filesystem>
#include "MeshIO.hpp"

using vec3f = std::array<float,3>;
using vec4i = std::array<int,4>;
using vec3i = std::array<int,3>;
using vec2i = std::array<int,2>;

struct ReadVtkTet{
private:
    void extractSurf();
    void extractEdge();
    std::vector<vec3f> pos;
    std::vector<vec4i> quads;
    std::vector<vec3i> surfs;
    std::vector<vec2i> edges;
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
    extractEdge();
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
    f_surf.close();

    f_pos.open("pos.txt");
    for(const auto& f:pos)
    {
        for(const auto& x:f)
            f_pos<<x<<"\t";
        f_pos<<"\n";
    }
    f_pos.close();
}


void ReadVtkTet::extractEdge()
{   
    //use the std::set to maintain the uniqueness of edges
    std::set<std::set<int>> edgeSet;
    for (int i = 0; i < quads.size(); i++)
    {
        std::array<int,4> tet=quads[i];
        
        int t0 = tet[0];
        int t1 = tet[1];
        int t2 = tet[2];
        int t3 = tet[3];

        std::set<int> l0{t0, t1}, l1{t0, t2}, l2{t0, t3},
                      l3{t1, t2}, l4{t1, t3}, l5{t2, t3};
        
        edgeSet.insert(l0);
        edgeSet.insert(l1);
        edgeSet.insert(l2);
        edgeSet.insert(l3);
        edgeSet.insert(l4);
        edgeSet.insert(l5);
    }

    //copy back from set
    size_t numEdges = edgeSet.size();
    edges.reserve(numEdges);
    for(auto &&line : edgeSet)
    {
        vec2i to;
        std::copy(line.begin(), line.end(), to.begin());
        edges.push_back(to);
    }
    

    std::cout<<"after extractEdge numEdges: "<<edges.size()<<std::endl;

    std::ofstream f_lines;
    f_lines.open("edges.txt");
    for(const auto& f:edges)
    {
        for(const auto& x:f)
            f_lines<<x<<"\t";
        f_lines<<"\n";
    }
    f_lines.close();
}

int main(){
    ReadVtkTet a;
    a.apply();
}