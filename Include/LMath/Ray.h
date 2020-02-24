#ifndef __LMATH_Ray_H__
#define __LMATH_Ray_H__

#include "Vector.h"
#include "Plane.h"
#include "Sphere.h"


//#include "OgrePlaneBoundedVolume.h"

namespace LMath 
{
	template <typename T>
    class RayBase
    {
    
    public:

		using ElementType = T;
		using Vector3 = VectorBase<T, 3>;
		using AxisAlignedBox = BoundsBase<T, 3>;
		using Result = std::pair<bool, T> ;
		using Plane = PlaneBase<T>;
		using Sphere = SphereBase<T>;

		

		RayBase() : RayBase(Vector3::Zero, Vector3::Forward) {}
        
		RayBase(const Vector3& origin, const Vector3& direction)
            :mOrigin(origin), mDirection(direction) {}

        /** Sets the origin of the ray. */
        void setOrigin(const Vector3& origin) {mOrigin = origin;} 
        /** Gets the origin of the ray. */
        const Vector3& getOrigin(void) const {return mOrigin;} 

        /** Sets the direction of the ray. */
        void setDirection(const Vector3& dir) {mDirection = dir;} 
        /** Gets the direction of the ray. */
        const Vector3& getDirection(void) const {return mDirection;} 

        /** Gets the position of a point t units along the ray. */
        Vector3 getPoint(ElementType t) const {
            return Vector3(mOrigin + (mDirection * t));
        }
        
        /** Gets the position of a point t units along the ray. */
        Vector3 operator*(ElementType t) const {
            return getPoint(t);
        }

        /** Tests whether this ray intersects the given plane. */
		Result intersects(const Plane& p) const
        {
			ElementType denom = p.normal.Dot(mDirection);
            if ( std::abs(denom) < std::numeric_limits<ElementType>::epsilon())
            {
                // Parallel
                return Result(false, (ElementType)0);
            }
            else
            {
				ElementType nom = p.normal.Dot(mOrigin) + p.d;
				ElementType t = -(nom / denom);
                return Result(t >= 0, static_cast<ElementType>(t));
            }
        }
  //      /** Tests whether this ray intersects the given plane bounded volume. */
		//Result intersects(const PlaneBoundedVolume& p) const
  //      {
  //          return Math::intersects(*this, p.planes, p.outside == Plane::POSITIVE_SIDE);
  //      }
        /** Tests whether this ray intersects the given sphere. */
        Result intersects(const Sphere& s, bool discardInside = true) const
        {
            // Adjust ray origin relative to sphere center
            Vector3 rayorig = mOrigin - s.getCenter();
            ElementType radius = s.getRadius();

            // Check origin inside first
            if (rayorig.squaredLength() <= radius*radius && discardInside)
            {
                return RayTestResult(true, (ElementType)0);
            }

            // Mmm, quadratics
            // Build coeffs which can be used with std quadratic solver
            // ie t = (-b +/- sqrt(b*b + 4ac)) / 2a
            ElementType a = mDirection.dotProduct(mDirection);
            ElementType b = 2 * rayorig.dotProduct(mDirection);
            ElementType c = rayorig.dotProduct(rayorig) - radius*radius;

            // Calc determinant
            ElementType d = (b*b) - (4 * a * c);
            if (d < 0)
            {
                // No intersection
                return RayTestResult(false, (ElementType)0);
            }
            else
            {
                // BTW, if d=0 there is one intersection, if d > 0 there are 2
                // But we only want the closest one, so that's ok, just use the
                // '-' version of the solver
                ElementType t = ( -b - Math::Sqrt(d) ) / (2 * a);
                if (t < 0)
                    t = ( -b + Math::Sqrt(d) ) / (2 * a);
                return RayTestResult(true, t);
            }
        }
        /** Tests whether this ray intersects the given box. */
        Result intersects(const AxisAlignedBox& box) const
        {
			if (box.isNull()) return std::pair<bool, ElementType>(false, (ElementType)0);
			//if (box.isInfinite()) return std::pair<bool, ElementType>(true, (ElementType)0);

			ElementType lowt = 0.0f;
			ElementType t;
			bool hit = false;
			Vector3 hitpoint;
			const Vector3& min = box.getMinimum();
			const Vector3& max = box.getMaximum();
			const Vector3& rayorig = ray.getOrigin();
			const Vector3& raydir = ray.getDirection();

			// Check origin inside first
			if (rayorig > min&& rayorig < max)
			{
				return std::pair<bool, ElementType>(true, (ElementType)0);
			}

			// Check each face in turn, only check closest 3
			// Min x
			if (rayorig.x <= min.x && raydir.x > 0)
			{
				t = (min.x - rayorig.x) / raydir.x;

				// Substitute t back into ray and check bounds and dist
				hitpoint = rayorig + raydir * t;
				if (hitpoint.y >= min.y && hitpoint.y <= max.y &&
					hitpoint.z >= min.z && hitpoint.z <= max.z &&
					(!hit || t < lowt))
				{
					hit = true;
					lowt = t;
				}
			}
			// Max x
			if (rayorig.x >= max.x && raydir.x < 0)
			{
				t = (max.x - rayorig.x) / raydir.x;

				// Substitute t back into ray and check bounds and dist
				hitpoint = rayorig + raydir * t;
				if (hitpoint.y >= min.y && hitpoint.y <= max.y &&
					hitpoint.z >= min.z && hitpoint.z <= max.z &&
					(!hit || t < lowt))
				{
					hit = true;
					lowt = t;
				}
			}
			// Min y
			if (rayorig.y <= min.y && raydir.y > 0)
			{
				t = (min.y - rayorig.y) / raydir.y;

				// Substitute t back into ray and check bounds and dist
				hitpoint = rayorig + raydir * t;
				if (hitpoint.x >= min.x && hitpoint.x <= max.x &&
					hitpoint.z >= min.z && hitpoint.z <= max.z &&
					(!hit || t < lowt))
				{
					hit = true;
					lowt = t;
				}
			}
			// Max y
			if (rayorig.y >= max.y && raydir.y < 0)
			{
				t = (max.y - rayorig.y) / raydir.y;

				// Substitute t back into ray and check bounds and dist
				hitpoint = rayorig + raydir * t;
				if (hitpoint.x >= min.x && hitpoint.x <= max.x &&
					hitpoint.z >= min.z && hitpoint.z <= max.z &&
					(!hit || t < lowt))
				{
					hit = true;
					lowt = t;
				}
			}
			// Min z
			if (rayorig.z <= min.z && raydir.z > 0)
			{
				t = (min.z - rayorig.z) / raydir.z;

				// Substitute t back into ray and check bounds and dist
				hitpoint = rayorig + raydir * t;
				if (hitpoint.x >= min.x && hitpoint.x <= max.x &&
					hitpoint.y >= min.y && hitpoint.y <= max.y &&
					(!hit || t < lowt))
				{
					hit = true;
					lowt = t;
				}
			}
			// Max z
			if (rayorig.z >= max.z && raydir.z < 0)
			{
				t = (max.z - rayorig.z) / raydir.z;

				// Substitute t back into ray and check bounds and dist
				hitpoint = rayorig + raydir * t;
				if (hitpoint.x >= min.x && hitpoint.x <= max.x &&
					hitpoint.y >= min.y && hitpoint.y <= max.y &&
					(!hit || t < lowt))
				{
					hit = true;
					lowt = t;
				}
			}

			return std::pair<bool, ElementType>(hit, (ElementType)lowt);
        }

    
	private:
		Vector3 mOrigin;
		Vector3 mDirection;
		};

}
#endif
