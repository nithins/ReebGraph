#include "SimplifyCT2.hpp"
#include "constants.h"
#include "utl.h"

#include <algorithm>


using namespace  std;

std::vector<int64_t> contourtree::splitMonkeysAndNazis(
    std::vector<int64_t>   &nodeids,
    std::vector<scalar_t>  &nodeFuncs,
    std::vector<char>      &nodeType,
	std::vector<arc_t>     &arcs)
{

    int numNodes = nodeFuncs.size();
	int numArcs  = arcs.size();

	ENSURES(numNodes == numArcs + 1);

    std::vector<arc_t> oldArcNumToNodeIdPairs;
	for(int i = 0 ; i < arcs.size(); i++)
		oldArcNumToNodeIdPairs.push_back(std::make_pair(nodeids[arcs[i].first],nodeids[arcs[i].second]));


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
		auto a = arcs[i].first;
		auto b = arcs[i].second;

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
			arcs.push_back({i,nodeUp[i][j]});
        }
    }

	numArcs = arcs.size();

    ENSURES(numArcs+1 == numNodes);

    std::map<arc_t,int64_t> nodeIdPairsToNewArcNum;
	for(int i = 0 ; i < arcs.size(); i++)
		nodeIdPairsToNewArcNum[std::make_pair(nodeids[arcs[i].first],nodeids[arcs[i].second])] = i;

    std::vector<int64_t> oldArcNumToNewArcNum;
    for(int i = 0 ; i < oldArcNumToNodeIdPairs.size(); ++i)
        oldArcNumToNewArcNum.push_back(nodeIdPairsToNewArcNum[oldArcNumToNodeIdPairs[i]]);

    return oldArcNumToNewArcNum;
}

typedef contourtree::arc_t arc_t;


/// \brief the generic simplification kernel
struct contourTreeSimplificationKernel {

    /// \brief getFeature weight
    std::function<scalar_t(arc_t)>   getFeatureWeight;

    /// \brief called for book keepingafter every cancellation the.
    ///       the cancelled arc and the newly inserted arc are given as parameters.
    std::function<void(arc_t,arc_t)> doneCancellation;

    /// \brief should return true if cancellation should stop before the given weight
    std::function<bool(arc_t)>       stopCancellation;


    /// \brief the actual cancellation procedure
    /// \returns the set of arcs that survive after cancellation
	std::vector<arc_t>  operator()
	(const std::vector<char>  &nodeType,const std::vector<arc_t>   &arcs)
    {


        size_t numNodes = nodeType.size();
		size_t numArcs  = arcs.size();
        ENSURES(numArcs +1 == numNodes);


        // Form a directed graph
        vector<vector<int64_t>> nodeUp(nodeType.size());
        vector<vector<int64_t>> nodeDown(nodeType.size());
		for(auto arc: arcs){
			auto a = arc.first;
			auto b = arc.second;

            nodeUp[a].push_back(b);
            nodeDown[b].push_back(a);
        }

        // Check there are no monkeys or nazis
        for(int i = 0 ; i < numNodes; ++i) {
            ENSURES(nodeUp[i].size()   <= 2);
            ENSURES(nodeDown[i].size() <= 2);
            ENSURES(nodeUp[i].size() + nodeDown[i].size() == 3 || nodeUp[i].size() + nodeDown[i].size() == 1);
        }


        // Vector to indicate which nodes are cancelled
        vector<bool> nodeIsCacelled(nodeType.size(),false);

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

            auto fda = getFeatureWeight(a);
            auto fdb = getFeatureWeight(b);

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

		for(auto arc: arcs){
			pq.push(arc);
        }

        // Go forth and cancel 'em all.
        while(pq.size() != 0) {

            auto arc = pq.top();pq.pop();

            if(canCancel(arc)) {

                if(stopCancellation(arc))
                    break;

                auto narc = doCancel(arc);
                pq.push(narc);
                doneCancellation(arc,narc);
            }
        }

		std::vector<arc_t> sarcs;

        for(int i =0 ; i < nodeIsCacelled.size(); ++i)
            if(!nodeIsCacelled[i])
                for(auto j: nodeUp[i]){
                    ENSURES(!nodeIsCacelled[j]);
					sarcs.push_back({i,j});
                }

        return sarcs;
    }
};


struct PersistenceFunction {


    scalar_t operator()(arc_t a) const{
        ENSURES(nodeFuncs[a.second] >= nodeFuncs[a.first]);
        return (nodeFuncs[a.second] - nodeFuncs[a.first])/frange;
    }


    PersistenceFunction(const std::vector<scalar_t>  &nodeFuncs,bool normalize=false)
        :nodeFuncs(nodeFuncs)
    {
        if( normalize) {
            auto fmin=*min_element(nodeFuncs.begin(),nodeFuncs.end());
            auto fmax=*max_element(nodeFuncs.begin(),nodeFuncs.end());
            frange = fmax - fmin;
        }
    }

private:

    const std::vector<scalar_t>  &nodeFuncs;
    scalar_t frange = 1;
};

struct HierarchyRecorder {

    std::vector<uint32_t>        &featureHierarchy;
    const std::vector<char>      &nodeType;

    void operator()(arc_t arc,arc_t narc){

        bool isMaxSad = nodeType[arc.second] == contourtree::MAXIMUM;

        featureHierarchy.push_back(isMaxSad);
        featureHierarchy.push_back((isMaxSad)?(arc.second):(arc.first)); // Extremum
        featureHierarchy.push_back((isMaxSad)?(arc.first):(arc.second)); // merged saddle
        featureHierarchy.push_back(narc.first);
        featureHierarchy.push_back(narc.second);
    }
};

struct ArcParentRecorder {

    const std::vector<char> &nodeType;
	std::map<arc_t,arc_t>   &arcParent;

	ArcParentRecorder(const std::vector<char>  &nodeType,const std::vector<arc_t> & arcs)
		:nodeType(nodeType),arcParent(_arcParent){
		for(auto arc: arcs)
			arcParent[arc] = arc_t(-1,-1);
    }

    arc_t getParent(arc_t arc) const{
        while(arcParent.at(arc) != arc_t(-1,-1)) {arc = arcParent.at(arc);}
        return arc;
    }

    void operator()(arc_t arc,arc_t narc){
        auto sad = (nodeType[arc.second] == contourtree::MAXIMUM) ? (arc.first) :(arc.second);

        arc_t da1(narc.first,sad);
        arc_t da2(sad,narc.second);

        ENSURES(arcParent.count(arc) == 1 && arcParent[arc] == arc_t(-1,-1) &&
                arcParent.count(da1) == 1 && arcParent[da1] == arc_t(-1,-1) &&
                arcParent.count(da2) == 1 && arcParent[da2] == arc_t(-1,-1) &&
                arcParent.count(narc) == 0 );

        arcParent[da1]  = narc;
        arcParent[da2]  = narc;
        arcParent[arc]  = narc;
        arcParent[narc] = arc_t(-1,-1);

    }

private:
	std::map<arc_t,arc_t>   _arcParent;

};

struct OrderRecorder {

	std::vector<arc_t> &carcs;

    void operator()(arc_t arc){
		carcs.push_back(arc);
    }
};

struct WeightsRecorder {

    std::vector<float>                   &wts;
    scalar_t lWt = 0;

    void operator()(arc_t arc,scalar_t wt){
        lWt = max(lWt,wt);
        wts.push_back(lWt);
    }
};

void contourtree::simplifyPers(
        const std::vector<scalar_t>  &nodeFuncs,
        const std::vector<char>      &nodeType,
		const std::vector<arc_t>   &arcsInt64,
		std::vector<arc_t>    &carcs,
        std::vector<float>    &wts,
        std::vector<uint32_t> &featureHierarchy,
		std::vector<arc_t>    &sarcsInt64, // Surviving arcs
        int reqNumFeatures,
        float reqPers
        )
 {

    ENSURES(nodeType.size() == nodeFuncs.size());

    // Persistence
    auto getPersistence  = PersistenceFunction{nodeFuncs,true};
    auto recordHierarchy = HierarchyRecorder{featureHierarchy,nodeType};
	auto recordOrder     = OrderRecorder{carcs};
    auto recordWeight    = WeightsRecorder{wts};

    // Book keeping after a cancellation
	int curNumFeatures = arcsInt64.size()/2;
    auto doneCancel = [&](arc_t arc,arc_t narc)
    {
        recordHierarchy(arc,narc);
        recordOrder(arc);
        recordWeight(arc,getPersistence(arc));
        curNumFeatures -=2;
    };

    auto stopCancel = [&](arc_t arc)->bool {
        return (curNumFeatures < reqNumFeatures) || (reqPers < getPersistence(arc) );
    };

    auto simpKernel = contourTreeSimplificationKernel{getPersistence,doneCancel,stopCancel};

	sarcsInt64 = simpKernel(nodeType,arcsInt64);

}

void contourtree::sqeezeCT(
        std::vector<int64_t>   &nodeIds,
        std::vector<scalar_t>  &nodeFuncs,
        std::vector<char>      &nodeTypes,
		std::vector<arc_t>     &sarcs
        )
{
    std::vector<int32_t> nodeRemap(nodeIds.size(),-1);
    int nsurv = 0;
	for(auto sarc : sarcs)
		for(auto i: {sarc.first,sarc.second})
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
		a = {nodeRemap[a.first],nodeRemap[a.second]};

    nodeIds   = nodeids_;
    nodeFuncs   = nodefns_;
    nodeTypes = nodeTypes_;
}


std::vector<int64_t>  contourtree::preSimplifyPers(
        std::vector<int64_t>         &nodeIds,
        std::vector<scalar_t>        &nodeFuncs,
        std::vector<char>            &nodeType,
		std::vector<arc_t>           &arcs,
        float preSimpThresh
        )
{


    ENSURES(nodeType.size() == nodeFuncs.size() && nodeFuncs.size() == nodeType.size());

    auto getPersistence   = PersistenceFunction{nodeFuncs,true};
	auto recordArcParent  = ArcParentRecorder{nodeType,arcs};

	auto simpKernel = contourTreeSimplificationKernel{getPersistence,recordArcParent,
            [&](arc_t arc)->bool{return preSimpThresh < getPersistence(arc);}};


	std::vector<arc_t> oldArcs = arcs;

	arcs = simpKernel(nodeType,arcs);

    std::map<arc_t,int64_t> newArcToNum;
	for(int i = 0 ; i < arcs.size(); i++ )
		newArcToNum[arcs[i]] = i;


    std::vector<int64_t> oldArcNumToNewArcNum;
	for(size_t i =0; i < oldArcs.size(); ++i)
		oldArcNumToNewArcNum.push_back(newArcToNum.at(recordArcParent.getParent(oldArcs[i])));

	sqeezeCT(nodeIds,nodeFuncs,nodeType,arcs);

    return oldArcNumToNewArcNum;
}



