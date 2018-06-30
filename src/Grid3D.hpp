#ifndef GRID3D_H
#define GRID3D_H

#include <set>
#include <stdint.h>
#include <vector>
#include <string>


#include "ScalarFunction.hpp"

namespace contourtree {

struct Tet {
    int64_t v[4];
};

class Grid3D : public ScalarFunction
{
public:
    Grid3D (int resx, int resy, int resz);
	Grid3D (scalar_t *fnVals, int resx, int resy, int resz);

public:
    int getMaxDegree();
    int getVertexCount();
    int getStar(int64_t v, std::vector<int64_t> &star);
    bool lessThan(int64_t v1, int64_t v2);
    scalar_t getFunctionValue(int64_t v);

public:
	scalar_t *data(){return fnVals;}
    inline size_t dimX() const{return dimx;}
    inline size_t dimY() const{return dimy;}
    inline size_t dimZ() const{return dimz;}

public:
    void loadGrid(std::string fileName);

protected:
    void updateStars();

private:
    int dimx, dimy, dimz;
    int nv;
    std::vector<Tet> tets;
    int starin[14][3];
    int64_t star[14];
	scalar_t *fnVals = nullptr;
	std::vector<scalar_t> fnVals_;

public:
    inline int64_t index(int64_t x, int64_t y, int64_t z) {
        return (x + y * dimx + z * dimx * dimy);
    }
};

} // namespace 

#endif // GRID3D_H
