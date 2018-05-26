#include "SimplifyCT2.hpp"
#include "constants.h"
#include "utl.h"

#include <algorithm>


using namespace  std;

typedef pair<int64_t,int64_t> arc_t;


std::vector<int64_t> contourtree::splitMonkeysAndNazis(
    std::vector<int64_t>   &nodeids,
    std::vector<scalar_t>  &nodeFuncs,
    std::vector<char>      &nodeType,
    std::vector<int64_t>   &arcs)
{

    int numNodes = nodeFuncs.size();
    int numArcs  = arcs.size()/2;


    std::vector<arc_t> oldArcNumToNodeIdPairs;
    for(int i = 0 ; i < arcs.size(); i+=2)
        oldArcNumToNodeIdPairs.push_back(std::make_pair(nodeids[arcs[i]],nodeids[arcs[i+1]]));


    vector<vector<int64_t>> nodeUp(numNodes);
    vector<vector<int64_t>> nodeDown(numNodes);

    auto makeDuplicateNode = [&](int64_t i)->int64_t {
        nodeFuncs.push_back(nodeFuncs[i]);
        nodeType.push_back(nodeType[i]);
        nodeids.push_back(nodeids[i]);

        nodeUp.resize(numNodes+1);
        nodeDown.resize(numNodes+1);
        return numNodes++;
    };

    for(int64_t i = 0 ; i < numArcs; ++i){
        auto a = arcs[i*2];
        auto b = arcs[i*2+1];

        nodeUp[a].push_back(b);
        nodeDown[b].push_back(a);
    }

    for(int64_t i = 0 ; i < numNodes; ++i) {

        for(int64_t j = 2 ; j < nodeDown[i].size(); ++j)
        {
//            std::cout << "Splitting Down Monkey " << SVAR(i) <<std::endl;

            auto ni = makeDuplicateNode(i);

            utl::replaceInstances(nodeUp[nodeDown[i][j]],i,ni);
            for(auto e : nodeUp[i])
                utl::replaceInstances(nodeDown[e],i,ni);

            nodeUp[ni] = nodeUp[i];
            nodeDown[ni] = {i,nodeDown[i][j]};
            nodeUp[i] = {ni};
        }

        if(nodeDown[i].size()>2)
            nodeDown[i] = {nodeDown[i][0],nodeDown[i][1]};

        for(int64_t j = 2 ; j < nodeUp[i].size(); ++j)
        {
//            std::cout << "Splitting up Monkey " << SVAR(i) <<std::endl;
            auto ni = makeDuplicateNode(i);

            utl::replaceInstances(nodeDown[nodeUp[i][j]],i,ni);
            for(auto e : nodeDown[i])
                utl::replaceInstances(nodeUp[e],i,ni);

            nodeDown[ni] = nodeDown[i];
            nodeUp[ni]   = {i,nodeUp[i][j]};
            nodeDown[i]  = {ni};
        }

        if(nodeUp[i].size()>2)
            nodeUp[i] = {nodeUp[i][0],nodeUp[i][1]};

        ENSURES(nodeUp[i].size()   <= 2);
        ENSURES(nodeDown[i].size() <= 2);

        if(nodeUp[i].size() == 2 && nodeDown[i].size() == 2) {
//            std::cout << "Splitting up Nazis " << SVAR(i) <<std::endl;
            auto ni = makeDuplicateNode(i);

            utl::replaceInstances(nodeDown[nodeUp[i][0]],i,ni);
            utl::replaceInstances(nodeDown[nodeUp[i][1]],i,ni);

            nodeUp[ni] = nodeUp[i];
            nodeDown[ni] = {i};
            nodeUp[i] = {ni};
        }

        ENSURES(nodeUp[i].size() + nodeDown[i].size() == 3 || nodeUp[i].size() + nodeDown[i].size() == 1);

    }

    for(int i = 0 ; i < numNodes; ++i) {
        ENSURES(nodeUp[i].size()   <= 2);
        ENSURES(nodeDown[i].size() <= 2);
        ENSURES(nodeUp[i].size() + nodeDown[i].size() == 3 || nodeUp[i].size() + nodeDown[i].size() == 1);
    }

    arcs.clear();
    for(int i = 0 ; i < numNodes; ++i) {
        for(int j = 0 ; j < nodeUp[i].size(); ++j)
        {
            arcs.push_back(i);
            arcs.push_back(nodeUp[i][j]);
        }
    }

    numArcs = arcs.size()/2;

    ENSURES(numArcs+1 == numNodes);

    std::map<arc_t,int64_t> nodeIdPairsToNewArcNum;
    for(int i = 0 ; i < arcs.size(); i+=2)
        nodeIdPairsToNewArcNum[std::make_pair(nodeids[arcs[i]],nodeids[arcs[i+1]])] = i/2;

    std::vector<int64_t> oldArcNumToNewArcNum;
    for(int i = 0 ; i < oldArcNumToNodeIdPairs.size(); ++i)
        oldArcNumToNewArcNum.push_back(nodeIdPairsToNewArcNum[oldArcNumToNodeIdPairs[i]]);

    return oldArcNumToNewArcNum;
}

void contourtree::simplifyPers(
        const std::vector<scalar_t>  &nodeFuncs,
        const std::vector<char>      &nodeType,
        const std::vector<int64_t>   &arcs,
        std::vector<uint32_t> &orderPairs,
        std::vector<float>    &wts,
        std::vector<uint32_t> &featureHierarchy,
        std::vector<int64_t>  &sarcs, // Surviving arcs
        int reqNumFeatures,
        float reqPers
        )
 {

//    // Debug Helper functions
//    auto nodeTypeStr= [&](int32_t a)->string {
//        if(nodeType[a] == MINIMUM) return "min";
//        if(nodeType[a] == SADDLE) return "sad";
//        if(nodeType[a] == MAXIMUM) return "max";
//        ENSURES(false);return "";
//    };


//    auto arcTypeStr= [&](arc_t a)->string {
//        return nodeTypeStr(a.first) +"-"+nodeTypeStr(a.first) +"-";
//    };

//    auto strArc = [&](arc_t a){
//        auto pers = getPersistence(a);
//        std::stringstream ss; ss <<a.first <<", " << a.second <<", " << SVAR(pers) << " " <<arcTypeStr(a) ;
//        return ss.str();
//    };

    size_t numNodes = nodeFuncs.size();
    size_t numArcs  = arcs.size()/2;

    ENSURES(nodeType.size() == numNodes && nodeType.size() == numNodes);
    ENSURES(numArcs +1 == numNodes);


    // Form a directed graph
    vector<vector<int64_t>> nodeUp(nodeType.size());
    vector<vector<int64_t>> nodeDown(nodeType.size());
    for(int i = 0 ; i < arcs.size(); i+= 2){
        auto a = arcs[i];
        auto b = arcs[i+1];

        nodeUp[a].push_back(b);
        nodeDown[b].push_back(a);
    }

    // Check there are no monkeys or nazis
    for(int i = 0 ; i < numNodes; ++i) {
        ENSURES(nodeUp[i].size()   <= 2);
        ENSURES(nodeDown[i].size() <= 2);
        ENSURES(nodeUp[i].size() + nodeDown[i].size() == 3 || nodeUp[i].size() + nodeDown[i].size() == 1);
    }

    // Persistence
    auto getPersistence = [&](arc_t a)->scalar_t {
        ENSURES(nodeFuncs[a.second] >= nodeFuncs[a.first]);
        return nodeFuncs[a.second] - nodeFuncs[a.first];
    };


    // Book keeping after a cancellation
    auto fmin=*min_element(nodeFuncs.begin(),nodeFuncs.end());
    auto fmax=*max_element(nodeFuncs.begin(),nodeFuncs.end());
    auto lpers=scalar_t(0);
    int curNumFeatures = numArcs;
    vector<bool> nodeIsCacelled(nodeType.size(),false);
    auto doneCancel = [&](arc_t arc,arc_t narc)
    {
        // Persistences NEED NOT be monotonic.
        // But they are monotonic in each dimension
        lpers = std::max(lpers,getPersistence(arc));

        orderPairs.push_back(arc.first);
        orderPairs.push_back(arc.second);

        wts.push_back(lpers/(fmax-fmin));

        bool isMaxSad = nodeType[arc.second] == MAXIMUM;

        featureHierarchy.push_back(isMaxSad);
        featureHierarchy.push_back((isMaxSad)?(arc.second):(arc.first)); // Extremum
        featureHierarchy.push_back((isMaxSad)?(arc.first):(arc.second)); // merged saddle
        featureHierarchy.push_back(narc.first); // merged saddle
        featureHierarchy.push_back(narc.second); // merged saddle

        curNumFeatures -=2;
    };

    // Check books if the procedure must quit before cancelling this arc
    auto shouldQuitBefore = [&](arc_t arc)->bool {
        return (curNumFeatures < reqNumFeatures) || (reqPers < getPersistence(arc) );
    };


    // Book keeping after all cacellations
    auto allCancelsDone = [&]()->bool {

        sarcs.clear();

        for(int i =0 ; i < nodeIsCacelled.size(); ++i)
            if(!nodeIsCacelled[i])
                for(auto j: nodeUp[i]){
                    ENSURES(!nodeIsCacelled[j]);
                    sarcs.push_back(i);
                    sarcs.push_back(j);
                }
    };



    // The actual cancellation procedure
    auto doCancel = [&](arc_t arc) {

        using namespace contourtree;

        auto e = (nodeType[arc.first] == MINIMUM)?(arc.first):(arc.second);
        auto s = (nodeType[arc.first] == MINIMUM)?(arc.second):(arc.first);


        //ENSURES(arc.first < arc.second);
        ENSURES(!nodeIsCacelled[e] && !nodeIsCacelled[s]);
        ENSURES((nodeType[e] == MINIMUM || nodeType[e] == MAXIMUM) && nodeType[s] == SADDLE);

        if(nodeType[e] == MINIMUM) ENSURES(nodeDown[e].size() == 0 && nodeUp[e].size() == 1);
        if(nodeType[e] == MAXIMUM) ENSURES(nodeDown[e].size() == 1 && nodeUp[e].size() == 0);
        if(nodeType[e] == MINIMUM) ENSURES(nodeDown[s].size() == 2 && nodeUp[s].size() == 1);
        if(nodeType[e] == MAXIMUM) ENSURES(nodeDown[s].size() == 1 && nodeUp[s].size() == 2);
        if(nodeType[e] == MINIMUM) ENSURES(utl::countInstances(nodeUp[e],s) == 1 && utl::countInstances(nodeDown[s],e) == 1 );
        if(nodeType[e] == MAXIMUM) ENSURES(utl::countInstances(nodeDown[e],s) == 1 && utl::countInstances(nodeUp[s],e) == 1 );


        nodeIsCacelled[e] = true;
        nodeIsCacelled[s] = true;

        if(nodeType[e] == MINIMUM) utl::deleteInstances(nodeDown[s],e);
        if(nodeType[e] == MAXIMUM) utl::deleteInstances(nodeUp[s],e);

        ENSURES(nodeDown[s].size() == 1 && nodeUp[s].size() == 1)
                << SVAR(nodeUp[s].size()) << SVAR(nodeDown[s].size());

        auto sd = nodeDown[s][0];
        auto su = nodeUp[s][0];


        ENSURES(utl::replaceInstances(nodeUp[sd],s,su) == 1);
        ENSURES(utl::replaceInstances(nodeDown[su],s,sd) == 1);

        nodeDown[s].clear();
        nodeDown[e].clear();
        nodeUp[s].clear();
        nodeUp[e].clear();

        doneCancel(arc,std::make_pair(sd,su));

        return std::make_pair(sd,su);
    };


    auto canCancel = [&](arc_t arc)->bool {

        using namespace contourtree;

        auto &a = arc.first;
        auto &b = arc.second;

 //       ENSURES(a < b);

        ENSURES((nodeType[a] == MINIMUM && nodeType[b] == SADDLE)||
                (nodeType[a] == MINIMUM && nodeType[b] == MAXIMUM)||
                (nodeType[a] == SADDLE && nodeType[b] == SADDLE)||
                (nodeType[a] == SADDLE && nodeType[b]  == MAXIMUM));

        if(nodeType[a] == MINIMUM && nodeDown[b].size() != 2)
            return false;

        if(nodeType[b] == MAXIMUM && nodeUp[a].size() != 2)
            return false;

        if(nodeIsCacelled[arc.first] || nodeIsCacelled[arc.second] )
            return false;

        if((nodeType[a] == MINIMUM && nodeType[b] == SADDLE)||
           (nodeType[a] == SADDLE  && nodeType[b] == MAXIMUM))
            return true;

        return false;
    };

    // Less than comparator.. It'll be properly inverted in the priority queue
    auto compareArcs = [&](arc_t a, arc_t b) {

//        //a < b is true
//        if(canCancel(a) && !canCancel(b))
//            return true;

//        //a < b is false
//        if(canCancel(b) && !canCancel(a))
//            return false;

        auto fda = getPersistence(a);
        auto fdb = getPersistence(b);

        if (fda != fdb)
            return fda < fdb;

        if(a.first != b.first)
            return a.first < b.first;

        if(a.second != b.second)
            return a.second < b.second;

        ENSURES(false) << "Dude. What the hell!!";

    };

    // Declare and prepare the pq
    priority_queue<arc_t,std::vector<arc_t>,std::function<bool(arc_t,arc_t)>>
            pq([&](arc_t a, arc_t b)->bool{return compareArcs(b,a);});

    for(int i = 0 ; i < arcs.size(); i += 2) {
        pq.push({arcs[i],arcs[i+1]});
    }

    // Go forth and cancel 'em all.
    while(pq.size() != 0) {

        auto arc = pq.top();pq.pop();

        if(canCancel(arc)) {

            if(shouldQuitBefore(arc))
                break;

            pq.push(doCancel(arc));
        }
    }

    allCancelsDone();
}

void contourtree::sqeezeCT(
        std::vector<int64_t>   &nodeIds,
        std::vector<scalar_t>  &nodeFuncs,
        std::vector<char>      &nodeTypes,
        std::vector<int64_t>   &sarcs
        )
{
    std::vector<int32_t> nodeRemap(nodeIds.size(),-1);
    int nsurv = 0;
    for(auto i : sarcs)
        if(nodeRemap[i] == -1)
            nodeRemap[i] = nsurv++;

    std::vector<int64_t> nodeids_(nsurv);
    std::vector<scalar_t> nodefns_(nsurv);
    std::vector<char> nodeTypes_(nsurv);

    for(int i = 0; i < nodeRemap.size(); ++i)
        if(nodeRemap[i] >= 0) {
            nodeids_[nodeRemap[i]] = nodeIds[i];
            nodefns_[nodeRemap[i]] = nodeFuncs[i];
            nodeTypes_[nodeRemap[i]] = nodeTypes[i];
        }

    for(auto &a : sarcs)
        a = nodeRemap[a];

    nodeIds   = nodeids_;
    nodeFuncs   = nodefns_;
    nodeTypes = nodeTypes_;
}


