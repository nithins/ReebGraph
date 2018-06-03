/// \brief This file is meant to provide a simpler interface to simplification

#include "ContourTreeData.hpp"
#include "SimFunction.hpp"
#include <queue>
#include <vector>
#include <functional>
#include <array>

namespace contourtree {

typedef std::pair<int64_t,int64_t> arc_t;


/// \brief the generic simplification kernel
struct contourTreeSimplificationKernel {

	/// \brief getFeature weight
	std::function<scalar_t(arc_t)>   getFeatureWeight;

	/// \brief called for book keepingafter every cancellation the.
	///       the cancelled arc and the newly inserted arc are given as parameters.
	std::function<void(arc_t,arc_t)> doneCancellation;

	/// \brief should return true if cancellation should stop before the given arc
	std::function<bool(arc_t)>       stopCancellation;


	/// \brief the actual cancellation procedure
	/// \returns the set of arcs that survive after cancellation
	std::vector<arc_t>  operator()
	(const std::vector<char>  &nodeType,const std::vector<arc_t>   &arcs);
};


/// \brief functor object to compute the persistence
struct PersistenceFunction {

	/// \brief Ctor
	PersistenceFunction(const std::vector<scalar_t>  &nodeFuncs,bool normalize=false);

	/// \brief functor
	scalar_t operator()(arc_t a) const;

private:
	const std::vector<scalar_t>  &nodeFuncs;
	scalar_t frange = 1;
};

/// \brief functor object to compute HyperVolume function on arcs
struct HyperVolumeFunction {

	/// \brief Ctor
	HyperVolumeFunction(const std::vector<scalar_t>  &nodeFuncs,
						const std::vector<char>      &nodeType,
						const std::vector<arc_t> & arcs,
						const std::vector<float> & arcVols);

	/// \brief functor
	scalar_t operator()(arc_t a) const;

	/// \brief self update after cancellation
	void operator()(arc_t carc, arc_t narc);

	std::map<arc_t,std::array<float,2> > &arcHVol;

private:
	const std::vector<scalar_t>         &nodeFuncs;
	const std::vector<char>             &nodeType;
	std::map<arc_t,std::array<float,2> >  arcHVol_;
	scalar_t fd(arc_t a) const;
};



/// \brief Splitting Monkey saddles and Nazi configurations
/// \returns a mapping from old arc idxs to new arc idxs
std::vector<int64_t> splitMonkeysAndNazis(
        std::vector<int64_t>   &nodeids,
        std::vector<scalar_t>  &nodeFuncs,
        std::vector<char>      &nodeType,
		std::vector<arc_t>     &arcs
        );

/// \brief Given a set of arcs (that survive after cancellation),
///        remove unncessary nodes and reindex
void sqeezeCT(
        std::vector<int64_t>         &nodeIds,
        std::vector<scalar_t>        &nodeFuncs,
        std::vector<char>            &nodeTypes,
		std::vector<arc_t>           &sarcs
        );

/// \brief Simplification variant intended for presimplification using persistence
///        All cancelled nodes will be eliminated
/// \returns a mapping from old arc idxs to new arc idxs
std::vector<int64_t>  preSimplifyPers(
        std::vector<int64_t>         &nodeIds,
        std::vector<scalar_t>        &nodeFuncs,
        std::vector<char>            &nodeType,
		std::vector<arc_t>           &arcs,
        float preSimpThreshNorm = 0.01
        );


/// \brief Simpler implementation of Persistence simplification
void simplifyPers(
		const std::vector<scalar_t>  &nodeFuncs,
		const std::vector<char>      &nodeType,
		const std::vector<arc_t>     &arcs,
		std::vector<arc_t>           &carcs, // Cancellation arcs
		std::vector<float>           &wts,
		std::vector<uint32_t>        &featureHierarchy,
		std::vector<arc_t>           &sarcs, // Surviving arcs
		int reqNumFeatures=-1,
		float reqPersThresh = 1
		);


/// \brief Hypervolume based simplification
void simplifyHyperVolume(
		const std::vector<scalar_t>  &nodeFuncs,
		const std::vector<char>      &nodeType,
		const std::vector<arc_t>   &arcs,
		const std::vector<float>   &arcVols,

		std::vector<arc_t>    &carcs,
		std::vector<float>    &wts,
		std::vector<uint32_t> &featureHierarchy,
		std::vector<arc_t>    &sarcs, // Surviving arcs

		int reqNumFeatures = -1
		);
}
