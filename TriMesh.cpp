#include "TriMesh.hpp"

#include <fstream>
//#include <windows.h>
#include <iostream>
#include <set>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>

namespace contourtree {

TriMesh::TriMesh()
{

}

int TriMesh::getMaxDegree() {
    return maxStar;
}

int TriMesh::getVertexCount() {
    return nv;
}

int TriMesh::getStar(int64_t v, std::vector<int64_t> &star) {
    int ct = 0;
   for(uint32_t vv :vertices[v].adj) {
        star[ct] = vv;
        ct ++;
    }
    return ct;
}

bool TriMesh::lessThan(int64_t v1, int64_t v2) {
    if(fnVals[v1] < fnVals[v2]) {
        return true;
    } else if(fnVals[v1] == fnVals[v2]) {
        return (v1 < v2);
    }
    return false;
}

unsigned char TriMesh::getFunctionValue(int64_t v) {
    return this->fnVals[v];
}

void TriMesh::loadData(std::string fileName)
{
    std::ifstream ip(fileName);
    if(!ip.is_open()) {
        std::cout << "could not read file" << fileName;
    }
	std::string line;

	std::getline(ip, line);

	std::stringstream ss(line);

	int nt = 0;

	ss >> nv;
	ss >> nt;

    vertices.resize(nv);
    fnVals.resize(nv);
	
	for(int i = 0; i < nv; i++) {
		std::getline(ip, line);
		int fn = line[3];
		fnVals[i] = fn;
	}

    //for(int i = 0;i < nv;i ++) {
        //line = text.readLine().split(" ");
       // int fn = std::string(line[3]).toInt();
        //fnVals[i] = fn;

    for(int i = 0; i < nt; i++) {
        std::getline(ip,line);
        int v1 = line[1];
        int v2 = line[2];
        int v3 = line[3];

        vertices[v1].adj.insert(v2);
        vertices[v1].adj.insert(v3);
        vertices[v2].adj.insert(v1);
        vertices[v2].adj.insert(v3);
        vertices[v3].adj.insert(v2);
        vertices[v3].adj.insert(v1);
    }

    maxStar = 0;
    for(int i = 0;i < nv;i ++) {
		maxStar = std::max<int>(maxStar, vertices[i].adj.size());
    }
    std::cout << maxStar;
}

}