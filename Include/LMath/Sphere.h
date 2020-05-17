#ifndef __Sphere_H_
#define __Sphere_H_

// Precompiler options


#include "Vector.h"
#include "Plane.h"

namespace LMath
{

	template <typename T>
	class SphereBase
	{
	public:
		using ElementType = T;
		using Vector3 = VectorBase<T, 3>;
		using AxisAlignedBox = BoundsBase<T, 3>;
		using Plane = PlaneBase<T>;

		SphereBase() : SphereBase(static_cast<ElementType>(1.0), Vector3::Zero) {}
		SphereBase(const Vector3& center, ElementType radius)
			: mRadius(radius), mCenter(center) {}

		/** Returns the radius of the sphere. */
		ElementType getRadius(void) const { return mRadius; }

		/** Sets the radius of the sphere. */
		void setRadius(ElementType radius) { mRadius = radius; }

		/** Returns the center point of the sphere. */
		const Vector3& getCenter(void) const { return mCenter; }

		/** Sets the center point of the sphere. */
		void setCenter(const Vector3& center) { mCenter = center; }

		/** Returns whether or not this sphere intersects another sphere. */
		bool intersects(const SphereBase& s) const
		{
			return (s.mCenter - mCenter).squaredLength() <= (s.mRadius + mRadius) * (s.mRadius + mRadius);
		}
		/** Returns whether or not this sphere intersects a box. */
		bool intersects(const AxisAlignedBox& box) const
		{
			if (box.IsValid()) return false;

			// Use splitting planes
			const Vector3& center = getCenter();
			ElementType radius = getRadius();
			const Vector3& min = box.GetMin();
			const Vector3& max = box.GetMax();

			// Arvo's algorithm
			ElementType s, d = 0;
			for (int i = 0; i < 3; ++i)
			{
				if (center.at(i) < min.at(i))
				{
					s = center.at(i) - min.at(i);
					d += s * s;
				}
				else if (center.at(i) > max.at(i))
				{
					s = center.at(i) - max.at(i);
					d += s * s;
				}
			}
			return d <= radius * radius;
		}
		/** Returns whether or not this sphere intersects a plane. */
		bool intersects(const Plane& plane) const
		{
			return std::abs(plane.getDistance(getCenter())) <= getRadius();
		}
		/** Returns whether or not this sphere intersects a point. */
		bool intersects(const Vector3& v) const
		{
			return (v - mCenter).squaredLength() <= mRadius * mRadius;
		}
		/** Merges another Sphere into the current sphere */
		void merge(const SphereBase& oth)
		{
			Vector3 diff = oth.getCenter() - mCenter;
			ElementType lengthSq = diff.squaredLength();
			ElementType radiusDiff = oth.getRadius() - mRadius;

			// Early-out
			if (radiusDiff * radiusDiff >= lengthSq)
			{
				// One fully contains the other
				if (radiusDiff <= 0.0f)
					return; // no change
				else
				{
					mCenter = oth.getCenter();
					mRadius = oth.getRadius();
					return;
				}
			}

			ElementType length = std::sqrt(lengthSq);
			ElementType t = (length + radiusDiff) / (2.0f * length);
			mCenter = mCenter + diff * t;
			mRadius = 0.5f * (length + mRadius + oth.getRadius());
		}




	private:
		ElementType mRadius;
		Vector3 mCenter;


	};
}

#endif

