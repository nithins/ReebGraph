#ifndef CONTOURTREE_HPP
#define CONTOURTREE_HPP

#include <string>
#include <vector>
#include "ScalarFunction.hpp"

namespace contourtree {

class MergeTree;

class ContourTree
{
public:
    struct Node {
        int64_t v;
        std::vector<int64_t> next;
        std::vector<int64_t> prev;
    };

public:
    ContourTree();

    void setup(const MergeTree * tree);
    void computeCT();
    void output(std::string fileName);

    void output(
        std::vector<int64_t>  &nodeids,
        std::vector<scalar_t> &nodefns,
        std::vector<char>     &nodeTypes,
        std::vector<int64_t>  &arcs,
		uint32_t *arcMap
     );

private:
    void remove(int64_t xi, std::vector<Node>& nodeArray);
    void removeAndAdd(std::vector<int64_t> &arr, int64_t rem, int64_t add);
    void remove(std::vector<int64_t> &arr, int64_t xi);
    void addArc(int64_t from, int64_t to);

public:
    const MergeTree * tree;
    std::vector<Node> nodesJoin;
    std::vector<Node> nodesSplit;
    std::vector<Node> ctNodes;

    int64_t nv;
};

}

#endif // CONTOURTREE_HPP
