#include "ContourTreeData.hpp"
//#include <QTextStream>
#include <fstream>
#include <cassert>
#include <iostream>
#include "constants.h"
#include <sstream>

#include <algorithm> 
#include <functional> 
#include <cctype>
#include <locale>


//#ifdef WIN32
#include <string>
#include <windows.h>

std::string workingdir()
{
	char buf[256];
	GetCurrentDirectoryA(256, buf);
	return std::string(buf) + '\\';
}
//#else
//std::string workingdir() { assert(false); }
//
//#endif

namespace contourtree {

ContourTreeData::ContourTreeData() {

}

// trim from start (in place)
inline void ltrim(std::string &s) {
	s.erase(s.begin(), std::find_if(s.begin(), s.end(),
		std::not1(std::ptr_fun<int, int>(std::isspace))));
}

// trim from end (in place)
inline void rtrim(std::string &s) {
	s.erase(std::find_if(s.rbegin(), s.rend(),
		std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
}

// trim from both ends (in place)
inline std::string trim(const std::string &s) {
	std::string sc = s;
	ltrim(sc);
	rtrim(sc);

	return sc;
}


void ContourTreeData::loadBinFile(std::string fileName) {
    // read meta data
    {
		std::cout << fileName << std::endl;
		std::cout << workingdir() << std::endl;

		std::ifstream ip(fileName + ".rg.dat");
		if (!ip.is_open()) {
			std::cout << "could not read file" << fileName + ".rg.dat";
			assert(false);
        }
       // QTextStream text(&ip);
		std::string line;

		std::getline(ip, line);
		std::stringstream ss(line);
        ss >> noNodes;

		std::getline(ip, line);
		ss = std::stringstream(line);
		ss >> noArcs;

        assert(noNodes == noArcs + 1);
        ip.close();
    }
    std::cout << noNodes << noArcs;

    std::vector<int64_t> nodeids(noNodes);
    std::vector<scalar_t> nodefns(noNodes);
    std::vector<char> nodeTypes(noNodes);
    std::vector<int64_t> arcs(noArcs * 2);

    // read the tree
    std::string rgFile = fileName + ".rg.bin";
    std::ifstream ip(rgFile, std::ios::binary);
    ip.read((char *)nodeids.data(),nodeids.size() * sizeof(int64_t));
    ip.read((char *)nodefns.data(),nodeids.size());
    ip.read((char *)nodeTypes.data(),nodeids.size());
    ip.read((char *)arcs.data(),arcs.size() * sizeof(int64_t));
    ip.close();

    std::cout << "finished reading data";
    this->loadData(nodeids,nodefns,nodeTypes,arcs);
}

void ContourTreeData::loadTxtFile(std::string fileName) {
	std::ifstream ip("fileName");
	if (!ip.is_open()) {
		std::cout << "could not read file" << "fileName";

    }
//std::stringstream text(&ip);
	std::string line;

	std::getline(ip, line);

	std::stringstream ss(line);

	int noArcs = 0;

	ss >> noNodes;
	ss >> noArcs;

    std::vector<int64_t> nodeids(noNodes);
	std::vector <scalar_t> nodefns(noNodes);
    std::vector<char> nodeTypes(noNodes);
    std::vector<int64_t> arcs(noArcs * 2);

    for(size_t i = 0;i < noNodes;i ++) {
        std::getline(ip,line);
		std::stringstream ss(line);


        int64_t v;
		float fn;
        char t;
		std::string type;

		ss >> v >> fn >> type;

		type = trim(type);
		
			
        if(type == "MINIMA") {
            t = MINIMUM;
        } else if(type == "MAXIMA") {
            t = MAXIMUM;
        } else if(type == "SADDLE") {
            t = SADDLE;
        } else {
            t = REGULAR;
        }
        nodeids[i] = v;
        nodefns[i] = (scalar_t)(fn);
        nodeTypes[i] = t;
    }
    for(size_t i = 0;i < noArcs;i ++) {
		std::getline(ip, line);
        int v1 = line[0];
        int v2 = line[1];
        arcs[i * 2 + 0] = v1;
        arcs[i * 2 + 1] = v2;
    }
    ip.close();
    std::cout << "finished reading data";
    this->loadData(nodeids,nodefns, nodeTypes,arcs);
}

void ContourTreeData::loadData(const std::vector<int64_t> &nodeids, const std::vector<scalar_t> &nodefns, const std::vector<char> &nodeTypes, const std::vector<int64_t> &iarcs) {
    nodes.resize(noNodes);
    nodeVerts.resize(noNodes);
    fnVals.resize(noNodes);
    type.resize(noNodes);
    arcs.resize(noArcs);

    for(uint32_t i = 0;i < noNodes;i ++) {
        nodeVerts[i] = nodeids[i];
        // TODO hard coding again.
        fnVals[i] = (float)(nodefns[i]) / 255.;
        type[i] = nodeTypes[i];
        nodeMap[nodeVerts[i]] = i;
    }

    for(uint32_t i = 0;i < noArcs;i ++) {
        arcs[i].from = nodeMap[iarcs[i * 2 + 0]];
        arcs[i].to = nodeMap[iarcs[i * 2 + 1]];
        arcs[i].id = i;
        nodes[arcs[i].from].next.push_back(i);
        nodes[arcs[i].to].prev.push_back(i);
    }
}

} // namespace
