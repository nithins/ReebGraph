#include "SimplifyCT2.hpp"
#include "constants.h"
#include "utl.h"

#include <algorithm>


using namespace  std;
using namespace contourtree;

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


/// \brief the actual cancellation procedure
/// \returns the set of arcs that survive after cancellation
std::vector<arc_t>  contourtree::contourTreeSimplificationKernel::operator()
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

		if(canCancel(a) && canCancel(b)) {

			auto fda = getFeatureWeight(a);
			auto fdb = getFeatureWeight(b);

			if (fda != fdb)
				return fda < fdb;
		}

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
			doneCancellation(arc,narc);
			pq.push(narc);
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



scalar_t contourtree::PersistenceFunction::operator()(arc_t a) const{
	ENSURES(nodeFuncs[a.second] >= nodeFuncs[a.first]);
	return (nodeFuncs[a.second] - nodeFuncs[a.first])/frange;
}


contourtree::PersistenceFunction::PersistenceFunction
(const std::vector<scalar_t>  &nodeFuncs,bool normalize)
	:nodeFuncs(nodeFuncs)
{
	if( normalize) {
		auto fmin=*min_element(nodeFuncs.begin(),nodeFuncs.end());
		auto fmax=*max_element(nodeFuncs.begin(),nodeFuncs.end());
		frange = fmax - fmin;
	}
}

scalar_t  HyperVolumeFunction::fd(arc_t a) const{
	auto r = nodeFuncs[a.second] - nodeFuncs[a.first];
	ENSURES(r >= 0); return r;
}


HyperVolumeFunction::HyperVolumeFunction
(const std::vector<scalar_t>  &nodeFuncs,
 const std::vector<char>      &nodeType,
 const std::vector<arc_t>     &arcs,
 const std::vector<float>     &arcVols)
	:nodeFuncs(nodeFuncs)
	,nodeType(nodeType)
    ,arcHVol_(new arcHVol_t)
    ,arcHVol(*arcHVol_.get())
{
	int i=0;
	for(auto arc: arcs){
		auto v = arcVols[i++];
		arcHVol[arc] = {v,fd(arc)*v};
	}
}


scalar_t HyperVolumeFunction::operator()(arc_t a) const{
	return arcHVol.at(a)[1];
}

void HyperVolumeFunction::operator()(arc_t carc, arc_t narc){

	// Saddle Max carc
	arc_t larc = {carc.first,narc.second};
	arc_t marc = {narc.first,carc.first};

	// Saddle Min Carc
	if(nodeType[carc.first] == contourtree::MINIMUM){
		larc = {narc.first,carc.second};
		marc = {carc.second,narc.second};
	}

	ENSURES(arcHVol.count(carc) == 1 &&
			arcHVol.count(larc) == 1 &&
			arcHVol.count(marc) == 1 &&
			arcHVol.count(narc) == 0);

	float v = 0
			+ arcHVol[carc][0]
			+ arcHVol[larc][0]
			+ arcHVol[marc][0];

	float hv = 0
			+ arcHVol[carc][1] + arcHVol[carc][0]*fd(marc)
			+ arcHVol[larc][1] + arcHVol[larc][0]*fd(marc)
			+ arcHVol[marc][1];

	arcHVol[narc] = {v,hv};

	arcHVol.erase(carc);
	arcHVol.erase(larc);
	arcHVol.erase(marc);
}



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

/// \brief Record the parent of an arc after cancellation/deletion
///
/// \remarks The member arcParent and numToArc are references to private members. It is
///          done for the following reason.
///
///          When setting up the simplification kernel, the compiler will pass a copy of
///          this functor to the kernel. When this happens, all callbacks will be made
///          only on the copy and NOT the object created before setting up the kernel.
///          By keeping them as references, when the compiler creates a copy of this class,
///          the arcParent and numToArc members in the copy will reference the original
///          functor's members. Thus all updates made on the copy will reflect in the original.
///
struct ArcParentRecorder {

    const std::vector<char> &nodeType;
	std::map<arc_t,arc_t>   &arcParent;
	std::vector<arc_t>      &numToArc;

	ArcParentRecorder(const std::vector<char>  &nodeType,const std::vector<arc_t> & arcs)
		:nodeType(nodeType)
		,arcParent(_arcParent)
		,numToArc(_numToArc)
	{
		numToArc = arcs;
		for(auto arc: arcs)
			arcParent[arc] = arc_t(-1,-1);
    }

    arc_t getParent(arc_t arc) const{
        while(arcParent.at(arc) != arc_t(-1,-1)) {arc = arcParent.at(arc);}
        return arc;
    }

	std::vector<int64_t> getArcnumToSarcnum (const std::vector<arc_t> & sarcs) {

		std::map<arc_t,int64_t> sarcToNum;
		for(int i = 0 ; i < sarcs.size(); i++ )
			sarcToNum[sarcs[i]] = i;


		std::vector<int64_t> ArcnumToSarcnum;
		for(size_t i =0; i < numToArc.size(); ++i)
			ArcnumToSarcnum.push_back(sarcToNum.at(getParent(numToArc[i])));

		return ArcnumToSarcnum;

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
	std::vector<arc_t> _numToArc;

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
		const std::vector<arc_t>   &arcs,
		std::vector<arc_t>    &carcs,
        std::vector<float>    &wts,
        std::vector<uint32_t> &featureHierarchy,
		std::vector<arc_t>    &sarcs, // Surviving arcs
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
	int curNumFeatures = arcs.size();
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

	sarcs = simpKernel(nodeType,arcs);

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


std::vector<int64_t>  contourtree::preSimplify(
        std::vector<int64_t>         &nodeIds,
        std::vector<scalar_t>        &nodeFuncs,
        std::vector<char>            &nodeType,
        std::vector<arc_t>           &arcs,
        std::vector<float>           &arcVols,
        std::string smethod,
        float preSimpThresh,
        int reqNumFeatures)
{

    ENSURES(nodeType.size() == nodeFuncs.size() && nodeFuncs.size() == nodeType.size());
    ENSURES(smethod == "PersT" || smethod == "PersN" || smethod == "HvolN" ) << SVAR(smethod);

    auto arcParentRecord  = ArcParentRecorder{nodeType,arcs};
    int curNumFeatures = arcs.size();

    std::function<scalar_t(arc_t)> getFeatureWt
            = PersistenceFunction{nodeFuncs,true};

    std::function<bool(arc_t)>     stopCancel
            = [&](arc_t)->bool {return curNumFeatures <= reqNumFeatures;};

    std::function<void(arc_t arc,arc_t narc)>  doneCancel
            = [&](arc_t arc,arc_t narc){arcParentRecord(arc,narc);curNumFeatures -= 2;};

    if (smethod == "PersT")
        stopCancel = [&](arc_t arc)->bool {return preSimpThresh < getFeatureWt(arc);};

    if(smethod == "HvolN"){
        auto hvf = HyperVolumeFunction{nodeFuncs,nodeType,arcs,arcVols};
        doneCancel = [=](arc_t arc,arc_t narc) mutable {doneCancel(arc,narc);hvf(arc,narc);};
        getFeatureWt = hvf;
    }

    auto simpKernel = contourTreeSimplificationKernel{getFeatureWt,doneCancel,stopCancel};

	arcs = simpKernel(nodeType,arcs);

	auto rmap = arcParentRecord.getArcnumToSarcnum(arcs);

	sqeezeCT(nodeIds,nodeFuncs,nodeType,arcs);

	return rmap;
}






void contourtree::simplifyHyperVolume(
		const std::vector<scalar_t>  &nodeFuncs,
		const std::vector<char>      &nodeType,
		const std::vector<arc_t>   &arcs,
		const std::vector<float>   &arcVols,
		std::vector<arc_t>    &carcs,
		std::vector<float>    &wts,
		std::vector<uint32_t> &featureHierarchy,
		std::vector<arc_t>    &sarcs, // Surviving arcs
		int reqNumFeatures
		)
 {

	ENSURES(nodeType.size() == nodeFuncs.size());

	// Feature wts
	auto featureWtFunc   = HyperVolumeFunction{nodeFuncs,nodeType,arcs,arcVols};
	auto recordHierarchy = HierarchyRecorder{featureHierarchy,nodeType};
	auto recordOrder     = OrderRecorder{carcs};
	auto recordWeight    = WeightsRecorder{wts};

	// Book keeping after a cancellation
	int curNumFeatures = arcs.size();
	auto doneCancel = [&](arc_t arc,arc_t narc)
	{
		recordHierarchy(arc,narc);
		recordOrder(arc);
		recordWeight(arc,featureWtFunc(arc));
		featureWtFunc(arc,narc);
		curNumFeatures -=2;
	};

	auto stopCancel = [&](arc_t arc)->bool {
		return curNumFeatures <= reqNumFeatures ;
	};

	auto simpKernel = contourTreeSimplificationKernel{featureWtFunc,doneCancel,stopCancel};

	sarcs = simpKernel(nodeType,arcs);

    ENSURES(featureWtFunc.getArcHVol().size() == 1);

    for(auto &wt: wts)	wt/= featureWtFunc.getArcHVol().begin()->second[1];

}




