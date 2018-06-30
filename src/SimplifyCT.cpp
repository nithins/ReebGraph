#include "SimplifyCT.hpp"
#include "utl.h"

#include <iostream>
#include <fstream>
#include <string>

namespace contourtree {

bool BranchCompare::operator() (uint32_t v1, uint32_t v2) {
     return sim->compare(v1,v2);
}

SimplifyCT::SimplifyCT() {
    queue = std::priority_queue<uint32_t,std::vector<uint32_t>,BranchCompare>(BranchCompare(this));
    order.clear();
}

void SimplifyCT::setInput(ContourTreeData *data) {
    this->data = data;
}

void SimplifyCT::addToQueue(uint32_t ano) {
    if(isCandidate(branches[ano])) {
        queue.push(ano);
        inq[ano] = true;
    }
}

bool SimplifyCT::isCandidate(const Branch &br) {
    uint32_t from = br.from;
    uint32_t to = br.to;
    if(nodes[from].prev.size() == 0) {
        // minimum
        if(nodes[to].prev.size() > 1) {
            return true;
        } else {
            return false;
        }
    }
    if(nodes[to].next.size() == 0) {
        // maximum
        if(nodes[from].next.size() > 1) {
            return true;
        } else {
            return false;
        }
    }
    return false;
}

void SimplifyCT::initSimplification(SimFunction* f) {
    branches.resize(data->noArcs);
    nodes.resize(data->noNodes);
    for(uint32_t i = 0;i < branches.size();i ++) {
        branches[i].from = data->arcs[i].from;
        branches[i].to = data->arcs[i].to;
        branches[i].parent = -1;
        branches[i].arcs.push_back(i);

        nodes[branches[i].from].next.push_back(i);
        nodes[branches[i].to].prev.push_back(i);
    }

    fn.resize(branches.size());
    removed.resize(branches.size(),false);
    invalid.resize(branches.size(),false);
    inq.resize(branches.size(), false);

    vArray.resize(nodes.size());

    simFn = f;
    if(f != NULL) {
        simFn->init(fn, branches);

        for(uint32_t i = 0;i < branches.size();i ++) {
            addToQueue(i);
        }
    }
}


bool SimplifyCT::compare(uint32_t b1, uint32_t b2) const {
    // If I want smallest weight on top, I need to return true if b1 > b2 (sort in descending order)
    if(fn[b1] > fn[b2]) {
        return true;
    }
    if(fn[b1] < fn[b2]) {
        return false;
    }
    float p1 = data->fnVals[branches[b1].to] - data->fnVals[branches[b1].from];
    float p2 = data->fnVals[branches[b2].to] - data->fnVals[branches[b2].from];
    if(p1 > p2) {
        return true;
    }
    if(p1 < p2) {
        return false;
    }
    int diff1 = branches[b1].to - branches[b1].from;
    int diff2 = branches[b2].to - branches[b2].from;
    if(diff1 > diff2) {
        return true;
    }
    if(diff1 < diff2) {
        return false;
    }
    return (branches[b2].from > branches[b1].from);
}

void SimplifyCT::removeArc(uint32_t ano) {
    Branch br = branches[ano];
    uint32_t from = br.from;
    uint32_t to = br.to;
    uint32_t mergedVertex = -1;
    if(nodes[from].prev.size() == 0) {
        // minimum
        mergedVertex = to;
    }
    if(nodes[to].next.size() == 0) {
        // maximum
        mergedVertex = from;
    }

    utl::deleteInstances(nodes[from].next, ano);
    utl::deleteInstances(nodes[to].prev, ano);
    removed[ano] = true;

    vArray[mergedVertex].push_back(ano);
    if(nodes[mergedVertex].prev.size() == 1 && nodes[mergedVertex].next.size() == 1) {
        mergeVertex(mergedVertex);
    }
    if(simFn != NULL)
        simFn->branchRemoved(branches, ano, invalid);
}

void SimplifyCT::mergeVertex(uint32_t v) {
    uint32_t prev = nodes[v].prev.at(0);
    uint32_t next = nodes[v].next.at(0);
    int a = -1;
    int rem = -1;
    if(inq[prev]) {
        invalid[prev] = true;
        removed[next] = true;
        branches[prev].to = branches[next].to;
        a = prev;
        rem = next;

        for(int i = 0;i < nodes[branches[prev].to].prev.size();i ++) {
            if(nodes[branches[prev].to].prev[i] == next) {
                nodes[branches[prev].to].prev[i] = prev;
            }
        }
    } else {
        invalid[next] = true;
        removed[prev] = true;
        branches[next].from = branches[prev].from;
        a = next;
        rem = prev;

        for(int i = 0;i < nodes[branches[next].from].next.size();i ++) {
            if(nodes[branches[next].from].next[i] == prev) {
                nodes[branches[next].from].next[i] = next;
            }
        }
        if(simFn != NULL && !inq[next]) {
            addToQueue(next);
        }
    }
    for(int i = 0;i < branches[rem].children.size();i ++) {
        int ch = branches[rem].children.at(i);
        branches[a].children.push_back(ch);
        ENSURES(branches[ch].parent == rem);
        branches[ch].parent = a;
    }
	for( const auto & arc: branches[rem].arcs)
		branches[a].arcs.push_back(arc);

    for(int i = 0;i < vArray[v].size();i ++) {
        uint32_t aa = vArray[v].at(i);
        branches[a].children.push_back(aa);
        branches[aa].parent = a;
    }
    branches[rem].parent = -2;
}

void SimplifyCT::simplify(contourtree::SimFunction *simFn) {
    std::cout << "init"<< std::endl;
    initSimplification(simFn);

    std::cout << "going over priority queue"<< std::endl;
    while(queue.size() > 0) {
        uint32_t ano = queue.top();
        queue.pop();
        inq[ano] = false;
        if(!removed[ano]) {
            if(invalid[ano]) {
                simFn->update(branches, ano);
                invalid[ano] = false;
                addToQueue(ano);
            } else {
                if(isCandidate(branches[ano])) {
                    removeArc(ano);
                    order.push_back(ano);
                }
            }
        }
    }
    std::cout << "pass over removed"<< std::endl;
    int root = 0;
    for(int i = 0;i < removed.size();i ++) {
        if(!removed[i]) {
            ENSURES(root == 0);
            order.push_back(i);
            root ++;
        }
    }
}


void SimplifyCT::simplify(const std::vector<uint32_t> &order, int topk, float th, const std::vector<float> &wts) {
    std::cout << "init"<< std::endl;
    initSimplification(NULL);

    std::cout << "going over order queue"<< std::endl;
    for(int i = 0;i < order.size();i ++) {
        inq[order.at(i)] = true;
    }
    if(topk > 0) {
        int ct = order.size() - topk;
		std::cout << order.size() << std::endl;
        for(int i = 0;i < ct;i ++) {
            uint32_t ano = order.at(i);
            ENSURES(isCandidate(branches[ano])) << "failing candidate test";
            inq[ano] = false;
            removeArc(ano);
        }
    } else if(th != 0) {
        for(int i = 0;i < order.size() - 1;i ++) {
            uint32_t ano = order.at(i);
            ENSURES(isCandidate(branches[ano])) << "failing candidate test";
            float fn = wts.at(i);
            if(fn > th) {
                break;
            }
            inq[ano] = false;
            removeArc(ano);
        }
    }
}

void SimplifyCT::outputOrder(std::vector<uint32_t> &order,std::vector<float> &wts)
{
    order = this->order;
    float pwt = 0;
    for(size_t i = 0;i < order.size();i ++) {
        uint32_t ano = order.at(i);
        float val = this->simFn->getBranchWeight(ano);
        wts.push_back(val);
        if(i > 0) {
            ENSURES(pwt <= val);
        }
        pwt = val;
    }


    if (wts.size() > 0) {
        // normalize weights
        float maxWt = wts.at(wts.size() - 1);
        if (maxWt == 0) maxWt = 1;
        for (int i = 0; i < wts.size(); i++) {
            wts[i] /= maxWt;
        }
    }

}

void SimplifyCT::outputOrder(std::string fileName) {
    std::cout << "Writing meta data"<< std::endl;
    {
		std::ofstream pr(fileName + ".order.dat");
		if (!pr.is_open()) {
			std::cout << "could not write file" << fileName + ".order.dat";

        }

		pr << order.size() << "\n";
		pr.close();
    }

    std::vector<float> wts;
    outputOrder(order,wts);

    std::cout << "writing tree output"<< std::endl;
    std::string binFile = fileName + ".order.bin";
    std::ofstream of(binFile,std::ios::binary);
    of.write((char *)order.data(),order.size() * sizeof(uint32_t));
    of.write((char *)wts.data(),wts.size() * sizeof(float));
//    of.write((char *)arcs.data(),arcs.size() * sizeof(uint32_t));
    of.close();
}

void SimplifyCT::computeFeatureHierarchy(const std::vector<uint32_t> &order,std::vector<uint32_t> &featureHierarchy){
    std::cout << "init" << std::endl;
    initSimplification(NULL);

    featureHierarchy.clear();

    std::cout << "going over order queue"<< std::endl;
    for(int i = 0;i < order.size();i ++) {
        inq[order.at(i)] = true;
    }

    for(int i = 0;i < order.size()-1;i ++) {
        uint32_t ano = order.at(i);
        ENSURES(isCandidate(branches[ano])) << "failing candidate test" << SVAR(ano) << SVAR(i);
        inq[ano] = false;


        Branch br = branches[ano];
        uint32_t from = br.from;
        uint32_t to = br.to;
        uint32_t mergedVertex = -1;
        uint32_t destroyedVertex = -1;

        uint32_t mergedVertexUp = -1;
        uint32_t mergedVertexDown = -1;

        if(nodes[from].prev.size() == 0) {
            // minimum
            mergedVertex = to;
            destroyedVertex = from;

            ENSURES(nodes[to].prev.size() == 2) << "Does not appaer to be a contour tree";

            uint32_t ano_0 = nodes[to].prev[0];
            uint32_t ano_1 = nodes[to].prev[1];
            if(ano_0 != ano) std::swap(ano_0,ano_1);

            ENSURES(nodes[to].next.size() == 1)<< "This arc could not be simplified";
            uint32_t ano_2 = nodes[to].next[0];

            ENSURES(ano_0 == ano && ano_1 != ano && ano_2 != ano);


            //std::cout <<SVAR(ano) << SVAR(ano_0) << SVAR(ano_1) << SVAR(ano_2) << std::endl;

            mergedVertexUp   = branches[ano_2].to;
            mergedVertexDown = branches[ano_1].from;

            //std::cout <<SVAR(mergedVertexUp) << SVAR(mergedVertexDown) << std::endl;


        }
        if(nodes[to].next.size() == 0) {
            // maximum
            mergedVertex = from;
            destroyedVertex = to;

            ENSURES(nodes[from].next.size() == 2)<< "Does not appaer to be a contour tree";

            uint32_t ano_0 = nodes[from].next[0];
            uint32_t ano_1 = nodes[from].next[1];
            if(ano_0 != ano) std::swap(ano_0,ano_1);

            ENSURES(nodes[from].prev.size() == 1);
            uint32_t ano_2 = nodes[from].prev[0];

            ENSURES(ano_0 == ano && ano_1 != ano && ano_2 != ano)"This arc could not be simplified";


            //std::cout <<SVAR(ano) << SVAR(ano_0) << SVAR(ano_1) << SVAR(ano_2) << std::endl;

            mergedVertexUp   = branches[ano_1].to;
            mergedVertexDown = branches[ano_2].from;

            //std::cout <<SVAR(mergedVertexUp) << SVAR(mergedVertexDown) << std::endl;

        }


//        std::cout <<"child1 = " << from<<", " << to << std::endl;
//        std::cout <<"child2 = " << mergedVertexDown<<", " << mergedVertex << std::endl;
//        std::cout <<"child3 = " << mergedVertex<<", " << mergedVertexUp << std::endl;
//        std::cout <<"parent = " << mergedVertexDown<<", " << mergedVertexUp << std::endl;

        featureHierarchy.push_back(destroyedVertex);
        featureHierarchy.push_back(mergedVertex);
        featureHierarchy.push_back(mergedVertexDown);
        featureHierarchy.push_back(mergedVertexUp);


        removeArc(ano);
    }


}

}
