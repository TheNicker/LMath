
#ifndef __PlaneBoundedVolume_H_
#define __PlaneBoundedVolume_H_

#include "Plane.h"
#include "Sphere.h"
#include "Ray.h"
#include <vector>

namespace LMath 
{

	template <typename T>
    class PlaneBoundedVolumeBase
    {
    public:
		
		using Plane = PlaneBase<T>;
        using PlaneList = std::vector<Plane> ;
		using AxisAlignedBox = BoundsBase<T, 3>;
		using Vector3 = VectorBase<T, 3>;
		using Ray = RayBase<T>;
		using Sphere = SphereBase<T>;
		using ElementType = T;
        /// Publicly accessible plane list, you can modify this direct
        PlaneList planes;
		
        typename Plane::Side outside;

		PlaneBoundedVolumeBase() :outside(Plane::Side::Negative) {}
        /** Constructor, determines which side is deemed to be 'outside' */
		PlaneBoundedVolumeBase(typename Plane::Side theOutside)
            : outside(theOutside) {}

        /** Intersection test with AABB
        @remarks May return false positives but will never miss an intersection.
        */
        inline bool intersects(const AxisAlignedBox& box) const
        {
            if (box.IsValid()) return false;
            //if (box.isInfinite()) return true;

            // Get centre of the box
            Vector3 centre = box.getCenter();
            // Get the half-size of the box
            Vector3 halfSize = box.getHalfSize();
            
            PlaneList::const_iterator i, iend;
            iend = planes.end();
            for (i = planes.begin(); i != iend; ++i)
            {
                const Plane& plane = *i;

                Plane::Side side = plane.getSide(centre, halfSize);
                if (side == outside)
                {
                    // Found a splitting plane therefore return not intersecting
                    return false;
                }
            }

            // couldn't find a splitting plane, assume intersecting
            return true;

        }
        /** Intersection test with Sphere
        @remarks May return false positives but will never miss an intersection.
        */
        inline bool intersects(const Sphere& sphere) const
        {
            PlaneList::const_iterator i, iend;
            iend = planes.end();
            for (i = planes.begin(); i != iend; ++i)
            {
                const Plane& plane = *i;

                // Test which side of the plane the sphere is
                Real d = plane.getDistance(sphere.getCenter());
                // Negate d if planes point inwards
                if (outside == Plane::NEGATIVE_SIDE) d = -d;

                if ( (d - sphere.getRadius()) > 0)
                    return false;
            }

            return true;

        }

        /** Intersection test with a Ray
        @return std::pair of hit (bool) and distance
        @remarks May return false positives but will never miss an intersection.
        */
        inline std::pair<bool, ElementType> intersects(const Ray& ray)
        {

			std::vector<Plane>::const_iterator planeit, planeitend;
			planeitend = planes.end();
			bool allInside = true;
			std::pair<bool, Real> ret;
			std::pair<bool, Real> end;
			ret.first = false;
			ret.second = 0.0f;
			end.first = false;
			end.second = 0;


			// derive side
			// NB we don't pass directly since that would require Plane::Side in 
			// interface, which results in recursive includes since Math is so fundamental
			Plane::Side outside = normalIsOutside ? Plane::POSITIVE_SIDE : Plane::NEGATIVE_SIDE;

			for (planeit = planes.begin(); planeit != planeitend; ++planeit)
			{
				const Plane& plane = *planeit;
				// is origin outside?
				if (plane.getSide(ray.getOrigin()) == outside)
				{
					allInside = false;
					// Test single plane
					std::pair<bool, Real> planeRes =
						ray.intersects(plane);
					if (planeRes.first)
					{
						// Ok, we intersected
						ret.first = true;
						// Use the most distant result since convex volume
						ret.second = std::max(ret.second, planeRes.second);
					}
					else
					{
						ret.first = false;
						ret.second = 0.0f;
						return ret;
					}
				}
				else
				{
					std::pair<bool, Real> planeRes =
						ray.intersects(plane);
					if (planeRes.first)
					{
						if (!end.first)
						{
							end.first = true;
							end.second = planeRes.second;
						}
						else
						{
							end.second = std::min(planeRes.second, end.second);
						}

					}

				}
			}

			if (allInside)
			{
				// Intersecting at 0 distance since inside the volume!
				ret.first = true;
				ret.second = 0.0f;
				return ret;
			}

			if (end.first)
			{
				if (end.second < ret.second)
				{
					ret.first = false;
					return ret;
				}
			}
			return ret;


        }

    };
}

#endif

