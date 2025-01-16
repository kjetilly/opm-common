#pragma once

#include "../LibCore/cvfBase.h"
#include "../LibCore/cvfVector3.h"

#include <array>
#include <cmath>
#include <vector>

namespace external 
{
namespace cvf
{
class Plane;
};

namespace caf
{
//==================================================================================================
//
//
//==================================================================================================
class HexGridIntersectionTools
{
public:
    //--------------------------------------------------------------------------------------------------
    ///
    //--------------------------------------------------------------------------------------------------
    struct ClipVx
    {
        ClipVx()
            : vx( cvf::Vec3d::ZERO )
            , normDistFromEdgeVx1( HUGE_VAL )
            , clippedEdgeVx1Id( -1 )
            , clippedEdgeVx2Id( -1 )
            , isVxIdsNative( true )
            , derivedVxLevel( -1 )
        {
        }

        cvf::Vec3d vx;

        double normDistFromEdgeVx1;
        size_t clippedEdgeVx1Id;
        size_t clippedEdgeVx2Id;

        bool isVxIdsNative; //< Pointing to real vertices, or indices to ClipVx's in the supplied triangle vertices array
        long long derivedVxLevel; //< Helper data to make it possible to track what set of ClipVx's the indices is reffering
                            // to in case of consecutive clips
    };

    static bool planeLineIntersect( const cvf::Plane& plane,
                                    const cvf::Vec3d& a,
                                    const cvf::Vec3d& b,
                                    cvf::Vec3d*       intersection,
                                    double*           normalizedDistFromA,
                                    double            epsilon );

    static bool planeTriangleIntersection( const cvf::Plane& plane,
                                           const cvf::Vec3d& p1,
                                           size_t            p1Id,
                                           const cvf::Vec3d& p2,
                                           size_t            p2Id,
                                           const cvf::Vec3d& p3,
                                           size_t            p3Id,
                                           ClipVx*           newVx1,
                                           ClipVx*           newVx2,
                                           bool*             isMostVxesOnPositiveSide );

    static void clipTrianglesBetweenTwoParallelPlanes( const std::vector<ClipVx>& triangleVxes,
                                                       const std::vector<long long>&    cellFaceForEachTriangleEdge,
                                                       const cvf::Plane&          p1Plane,
                                                       const cvf::Plane&          p2Plane,
                                                       std::vector<ClipVx>*       clippedTriangleVxes,
                                                       std::vector<long long>*          cellFaceForEachClippedTriangleEdge );

    static void clipPlanarTrianglesWithInPlaneTriangle( const std::vector<cvf::Vec3d>& triangleVxes,
                                                        const std::vector<long long>&        cellFaceForEachTriangleEdge,
                                                        const cvf::Vec3d&              tp1,
                                                        const cvf::Vec3d&              tp2,
                                                        const cvf::Vec3d&              tp3,
                                                        std::vector<cvf::Vec3d>*       clippedTriangleVxes,
                                                        std::vector<long long>* cellFaceForEachClippedTriangleEdge );

    static cvf::Vec3d planeLineIntersectionForMC( const cvf::Plane& plane,
                                                  const cvf::Vec3d& p1,
                                                  const cvf::Vec3d& p2,
                                                  double*           normalizedDistFromP1 );

    static long long planeHexIntersectionMC( const cvf::Plane&    plane,
                                       const cvf::Vec3d     cell[8],
                                       const size_t         hexCornersIds[8],
                                       std::vector<ClipVx>* triangleVxes,
                                       std::vector<long long>*    cellFaceForEachTriangleEdge );

    static long long       planeHexIntersectionMCTet( const cvf::Plane&    plane,
                                                const cvf::Vec3d     cell[8],
                                                const size_t         hexCornersIds[8],
                                                std::vector<ClipVx>* triangleVxes,
                                                std::vector<long long>*    cellFaceForEachTriangleEdge );
    static cvf::uint planeMcTetIntersection( const cvf::Plane&         plane,
                                             const cvf::Vec3d          hexCell[8],
                                             const size_t              hexCornersIds[8],
                                             const double              cornerDistToPlane[8],
                                             const std::array<long long, 4>& tetCell,
                                             std::vector<ClipVx>*      triangleVxes,
                                             std::vector<long long>*         cellFaceForEachTriangleEdge );
};

}; // namespace caf
} //namespace external
