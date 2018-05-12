#ifndef TRIMESH_HPP
#define TRIMESH_HPP

#include <set>
#include <string>

#include "ScalarFunction.hpp"



namespace contourtree {

class TriMesh : public ScalarFunction
{
public:
    struct Vertex {
        std::set<uint32_t> adj;
    };

public:
    TriMesh();

    int getMaxDegree();
    int getVertexCount();
    int getStar(int64_t v, std::vector<int64_t> &star);
    bool lessThan(int64_t v1, int64_t v2);
    scalar_t getFunctionValue(int64_t v);

public:
    void loadData(std::string fileName);

public:
    int nv;
    std::vector<scalar_t> fnVals;
    std::vector<Vertex> vertices;
    int maxStar;
};

}

#endif // TRIMESH_HPP
