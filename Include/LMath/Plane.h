
#ifndef __LMATH_PLANE_H__
#define __LMATH_PLANE_H__
#include "Vector.h"
#include "Bounds.h"
#include "Matrix.h"


namespace LMath 
{
	template <typename T>
	class  PlaneBase
	{
	public:
		using ElementType = T;
		using Vector3 = VectorBase<T, 3>;
		using Vector4 = VectorBase<T, 3>;
		using AxisAlignedBox = BoundsBase<T, 3>;

		Vector3 normal;
		ElementType d;

	public:
		/** Default constructor - sets everything to 0.
		*/
		PlaneBase() : normal(Vector3::Zero), d(static_cast<ElementType>(0)) {}
		/** Construct a plane through a normal, and a distance to move the plane along the normal.*/
		PlaneBase(const Vector3& rkNormal, ElementType fConstant)
		{
			normal = rkNormal;
			d = -fConstant;
		}
		/** Construct a plane using the 4 constants directly **/
		PlaneBase(ElementType a, ElementType b, ElementType c, ElementType _d) : normal(a, b, c), d(_d) {}
		/// @overload
		explicit PlaneBase(const Vector4& v) : normal(v.xyz()), d(v.w) {}
		PlaneBase(const Vector3& rkNormal, const Vector3& rkPoint)
		{
			redefine(rkNormal, rkPoint);
		}
		PlaneBase(const Vector3& p0, const Vector3& p1, const Vector3& p2)
		{
			redefine(p0, p1, p2);
		}

		/** The "positive side" of the plane is the half space to which the
			plane normal points. The "negative side" is the other half
			space. The flag "no side" indicates the plane itself.
		*/
		enum class Side
		{
			  None
			, Positive
			, Negative
			, Both
		};

		Side getSide(const Vector3& rkPoint) const
		{
			Real fDistance = getDistance(rkPoint);

			if (fDistance < 0.0)
				return Side::Negative;

			if (fDistance > 0.0)
				return Side::Positive;

			return Side::None;
		}

		/**
		Returns the side where the alignedBox is. The flag BOTH_SIDE indicates an intersecting box.
		One corner ON the plane is sufficient to consider the box and the plane intersecting.
		*/
		Side getSide(const AxisAlignedBox& box) const
		{
			if (box.isNull())
				return NO_SIDE;
			if (box.isInfinite())
				return BOTH_SIDE;

			return getSide(box.getCenter(), box.getHalfSize());
		}

		/** Returns which side of the plane that the given box lies on.
			The box is defined as centre/half-size pairs for effectively.
		@param centre The centre of the box.
		@param halfSize The half-size of the box.
		@return
			POSITIVE_SIDE if the box complete lies on the "positive side" of the plane,
			NEGATIVE_SIDE if the box complete lies on the "negative side" of the plane,
			and BOTH_SIDE if the box intersects the plane.
		*/
		Side getSide(const Vector3& centre, const Vector3& halfSize) const
		{
			// Calculate the distance between box centre and the plane
			ElementType dist = getDistance(centre);

			// Calculate the maximise allows absolute distance for
			// the distance between box centre and plane
			ElementType maxAbsDist = normal.absDotProduct(halfSize);

			if (dist < -maxAbsDist)
				return NEGATIVE_SIDE;

			if (dist > +maxAbsDist)
				return POSITIVE_SIDE;

			return BOTH_SIDE;
		}

		/** This is a pseudodistance. The sign of the return value is
			positive if the point is on the positive side of the plane,
			negative if the point is on the negative side, and zero if the
			point is on the plane.
			@par
			The absolute value of the return value is the true distance only
			when the plane normal is a unit length vector.
		*/
		ElementType getDistance(const Vector3& rkPoint) const
		{
			return normal.dotProduct(rkPoint) + d;
		}

		/** Redefine this plane based on 3 points. */
		void redefine(const Vector3& p0, const Vector3& p1, const Vector3& p2)
		{
			normal = Math::calculateBasicFaceNormal(p0, p1, p2);
			d = -normal.dotProduct(p0);
		}

		/** Redefine this plane based on a normal and a point. */
		void redefine(const Vector3& rkNormal, const Vector3& rkPoint)
		{
			normal = rkNormal;
			d = -rkNormal.Dot(rkPoint);
		}

		/** Project a vector onto the plane.
		@remarks This gives you the element of the input vector that is perpendicular
			to the normal of the plane. You can get the element which is parallel
			to the normal of the plane by subtracting the result of this method
			from the original vector, since parallel + perpendicular = original.
		@param v The input vector
		*/
		Vector3 projectVector(const Vector3& v) const
		{
			// We know plane normal is unit length, so use simple method
			MatrixBase<T,3,3> xform;
			xform[0][0] = 1.0f - normal.x * normal.x;
			xform[0][1] = -normal.x * normal.y;
			xform[0][2] = -normal.x * normal.z;
			xform[1][0] = -normal.y * normal.x;
			xform[1][1] = 1.0f - normal.y * normal.y;
			xform[1][2] = -normal.y * normal.z;
			xform[2][0] = -normal.z * normal.x;
			xform[2][1] = -normal.z * normal.y;
			xform[2][2] = 1.0f - normal.z * normal.z;
			return xform * v;
		}

		/** Normalises the plane.
			@remarks
				This method normalises the plane's normal and the length scale of d
				is as well.
			@note
				This function will not crash for zero-sized vectors, but there
				will be no changes made to their components.
			@return The previous length of the plane's normal.
		*/
		ElementType normalise(void)
		{
			ElementType fLength = normal.length();

			// Will also work for zero-sized vectors, but will change nothing
			// We're not using epsilons because we don't need to.
			// Read http://www.ogre3d.org/forums/viewtopic.php?f=4&t=61259
			if (fLength > ElementType(0.0f))
			{
				ElementType fInvLength = 1.0f / fLength;
				normal *= fInvLength;
				d *= fInvLength;
			}

			return fLength;
		}

		/// Get flipped plane, with same location but reverted orientation
		PlaneBase operator - () const
		{
			return PlaneBase(-(normal.x), -(normal.y), -(normal.z), -d); // not equal to Plane(-normal, -d)
		}

		/// Comparison operator
		bool operator==(const PlaneBase& rhs) const
		{
			return (rhs.d == d && rhs.normal == normal);
		}
		bool operator!=(const PlaneBase& rhs) const
		{
			return (rhs.d != d || rhs.normal != normal);
		}

		


	};

	/*inline Plane operator * (const Matrix4& mat, const Plane& p)
	{
		Plane ret;
		Matrix4 invTrans = mat.inverse().transpose();
		Vector4 v4(p.normal.x, p.normal.y, p.normal.z, p.d);
		v4 = invTrans * v4;
		ret.normal.x = v4.x;
		ret.normal.y = v4.y;
		ret.normal.z = v4.z;
		ret.d = v4.w / ret.normal.normalise();

		return ret;
	}

	inline bool Math::intersects(const Plane& plane, const AxisAlignedBox& box)
	{
		return plane.getSide(box) == Plane::BOTH_SIDE;
	}*/

	//typedef std::vector<Plane> PlaneList;
	/** @} */
	/** @} */

} // namespace Ogre

#endif
