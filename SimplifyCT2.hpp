/// \brief This file is meant to provide a simpler interface to simplification

#include "ContourTreeData.hpp"
#include "SimFunction.hpp"
#include <queue>
#include <vector>

#if defined (WIN32)
#include <functional>
#endif

namespace contourtree {

/// \brief Simpler implementation of Persistence simplification
void simplifyPers(
        const std::vector<scalar_t>  &nodeFuncs,
        const std::vector<char>      &nodeType,
        const std::vector<int64_t>   &arcs,
        std::vector<uint32_t>        &orderPairs,
        std::vector<float>           &wts,
        std::vector<uint32_t>        &featureHierarchy,
        std::vector<int64_t>         &sarcs, // Surviving arcs
        int reqNumFeatures=-1,
        float reqPersThresh = 1
        );

/// \brief Splitting Monkey saddles and Nazi configurations
/// \returns a mapping from old arc idxs to new arc idxs
std::vector<int64_t> splitMonkeysAndNazis(
        std::vector<int64_t>   &nodeids,
        std::vector<scalar_t>  &nodeFuncs,
        std::vector<char>      &nodeType,
        std::vector<int64_t>   &arcs
        );

/// \brief Given a set of arcs (that survive after cancellation),
///        remove unncessary nodes and reindex
void sqeezeCT(
        std::vector<int64_t>         &nodeIds,
        std::vector<scalar_t>        &nodeFuncs,
        std::vector<char>            &nodeTypes,
        std::vector<int64_t>         &sarcs
        );

/// \brief Simplification variant intended for presimplification using persistence
///        All cancelled nodes will be eliminated
/// \returns a mapping from old arc idxs to new arc idxs
std::vector<int64_t>  preSimplifyPers(
        std::vector<int64_t>         &nodeIds,
        std::vector<scalar_t>        &nodeFuncs,
        std::vector<char>            &nodeType,
        std::vector<int64_t>         &arcs,
        float preSimpThreshNorm = 0.01
        );


}
