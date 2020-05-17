#ifndef __LMATH_QUATERNOIN_H__
#define __LMATH_QUATERNOIN_H__
//#include "LMathConfig.h"
#include "Vector.h"
#include "Math.h"


namespace LMath
{
	template <typename T>
	class QuaternionBase : public VectorBase<T, 4>
	{
	public:
		using Vector3 = VectorBase<T, 3>;
		using Vector4 = VectorBase<T, 4>;
		using Base = Vector4;
		using Vector4::Vector4;
		using Vector4::at;
		using Vector4::X;
		using Vector4::Y;
		using Vector4::Z;
		using Vector4::W;
		using Vector4::L_Epsilon;
		using Vector4::L_Zero;
		using Vector4::L_One;
		using Vector4::L_Half;
		using Vector4::L_Two;
		using Vector4::L_Four;
		using Vector4::Dot;
		using Vector4::Zero;
		using VectorBase<T, 4>::operator*;


		using ElementType = typename Vector4::ElementType;

		static const QuaternionBase Identity;


		ElementType Angle(const QuaternionBase& rhs) const
		{
			return acos(fmin(fabs(Base::Dot(rhs)), L_One)) * L_Two;
		}

		QuaternionBase Conjugate() const
		{
			return QuaternionBase(-X(), -Y(), -Z(), W());
		}

		static QuaternionBase FromAngleAxis(ElementType angle, const Vector3& axis)
		{
			ElementType m = axis.Length();
			ElementType s = sin(angle / L_Two) / m;
			return { axis.X() * s , axis.Y() * s  , axis.Z() * s , cos(angle / L_Two) };
		}

		static QuaternionBase FromEuler(const Vector3& rotation)
		{
			return FromEuler(rotation.X(), rotation.Y(), rotation.Z());
		}

		static QuaternionBase FromEuler(ElementType x, ElementType y, ElementType z)
		{
			ElementType cx = cos(x * L_Half);
			ElementType cy = cos(y * L_Half);
			ElementType cz = cos(z * L_Half);
			ElementType sx = sin(x * L_Half);
			ElementType sy = sin(y * L_Half);
			ElementType sz = sin(z * L_Half);
			return
			{
				cx * sy * sz + cy * cz * sx,
				cx * cz * sy - cy * sx * sz,
				cx * cy * sz - cz * sx * sy,
				sx * sy * sz + cx * cy * cz
			};

		}

		static QuaternionBase FromToRotation(const Vector3& fromVector, const Vector3& toVector)
		{
			ElementType dot = fromVector.Dot(toVector);
			ElementType k = sqrt(fromVector.SquaredLength() * toVector.SquaredLength());

			if (fabs(dot / k + L_One) < L_Epsilon)
			{
				Vector3 normalizedOrtho = fromVector.Orthogonal().Normalized();
				return  { normalizedOrtho.X(), normalizedOrtho.Y(), normalizedOrtho.Z(), L_Zero };
			}
			Vector3 cross = fromVector.Cross(toVector);

			return QuaternionBase(cross.X(), cross.Y(), cross.Z(), dot + k).Normalized();
		}

		QuaternionBase Inverse() const
		{
			double norm = Base::NormSquared();
			return norm > 0.0 ? Conjugate() / norm : Zero;
		}



		QuaternionBase LerpClamped(const QuaternionBase& b, ElementType t) const
		{
			if (t < L_Zero) return Base::Normalized();
			else if (t > L_One) return b.Normalized();
			return LerpUnclamped(b, t);
		}

		QuaternionBase Lerp(const QuaternionBase& b, ElementType t) const
		{
			QuaternionBase quaternion;
			if (Dot(b) >= L_Zero)
				quaternion = *this * (L_One - t) + b * t;
			else
				quaternion = *this * (L_One - t) - b * t;
			return quaternion.Normalized();
		}


		static QuaternionBase LookRotation(const Vector3& forward)
		{
			return LookRotation(forward, Vector3::Up);
		}

		static QuaternionBase LookRotation(const Vector3& forwardInput, const Vector3& upwardsInput)
		{
			// Normalize inputs
			Vector3 forward = forwardInput.Normalized();
			Vector3 upwards = upwardsInput.Normalized();
			// Don't allow zero vectors
			if (forward.isZeroLength() || upwards.isZeroLength())
				return Identity;

			// Handle alignment with up direction
			if (L_One - fabs(forward.Dot(upwards)) < L_Epsilon)
				return FromToRotation(Vector3::Forward, forward);
			// Get orthogonal vectors
			Vector3 right = upwards.Cross(forward).Normalized();
			upwards = forward.Cross(right);
			// Calculate rotation
			QuaternionBase quaternion;
			ElementType radicand = right.X() + upwards.Y() + forward.Z();
			if (radicand > L_Zero)
			{
				quaternion.W() = sqrt(L_One + radicand) * L_Half;
				ElementType recip = L_One / (L_Four * quaternion.W());
				quaternion.X() = (upwards.Z() - forward.Y()) * recip;
				quaternion.Y() = (forward.X() - right.Z()) * recip;
				quaternion.Z() = (right.Y() - upwards.X()) * recip;
			}
			else if (right.X() >= upwards.Y() && right.X() >= forward.Z())
			{
				quaternion.X() = sqrt(L_One + right.X() - upwards.Y() - forward.Z()) * L_Half;
				ElementType recip = L_One / (L_Four * quaternion.X());
				quaternion.W() = (upwards.Z() - forward.Y()) * recip;
				quaternion.Z() = (forward.X() + right.Z()) * recip;
				quaternion.Y() = (right.Y() + upwards.X()) * recip;
			}
			else if (upwards.Y() > forward.Z())
			{
				quaternion.Y() = sqrt(L_One - right.X() + upwards.Y() - forward.Z()) * L_Half;
				ElementType recip = L_One / (L_Four * quaternion.Y());
				quaternion.Z() = (upwards.Z() + forward.Y()) * recip;
				quaternion.W() = (forward.X() - right.Z()) * recip;
				quaternion.X() = (right.Y() + upwards.X()) * recip;
			}
			else
			{
				quaternion.Z() = sqrt(L_One - right.X() - upwards.Y() + forward.Z()) * L_Half;
				ElementType recip = L_One / (L_Four * quaternion.Z());
				quaternion.Y() = (upwards.Z() + forward.Y()) * recip;
				quaternion.X() = (forward.X() + right.Z()) * recip;
				quaternion.W() = (right.Y() - upwards.X()) * recip;
			}
			return quaternion;
		}

		QuaternionBase RotateTowards(QuaternionBase to,
			ElementType maxRadiansDelta)
		{
			ElementType angle = Angle(to);
			if (angle == L_Zero)
				return to;
			maxRadiansDelta = (std::max)(maxRadiansDelta, angle - LMath::Constans::Pi);
			ElementType t = (std::min)(L_One, maxRadiansDelta / angle);
			return Slerp(to, t);
		}

		QuaternionBase SlerpClamped(QuaternionBase b, ElementType t)
		{
			if (t < L_Zero) return this->Normalized();
			else if (t > L_One) return b.Normalized();
			return Slerp(b, t);
		}

		QuaternionBase Slerp(QuaternionBase b, ElementType t)
		{
			ElementType n1;
			ElementType n2;
			ElementType n3 = Dot(b);
			bool flag = false;
			if (n3 < L_Zero)
			{
				flag = true;
				n3 = -n3;
			}
			if (n3 > 0.999999)
			{
				n2 = L_One - t;
				n1 = flag ? -t : t;
			}
			else
			{
				ElementType n4 = acos(n3);
				ElementType n5 = L_One / sin(n4);
				n2 = sin((L_One - t) * n4) * n5;
				n1 = flag ? -sin(t * n4) * n5 : sin(t * n4) * n5;
			}
			QuaternionBase quaternion;
			quaternion.X() = (n2 * X()) + (n1 * b.X());
			quaternion.Y() = (n2 * Y()) + (n1 * b.Y());
			quaternion.Z() = (n2 * Z()) + (n1 * b.Z());
			quaternion.W() = (n2 * W()) + (n1 * b.W());
			return quaternion.Normalized();
		}

		void ToAngleAxis(ElementType& angle, Vector3& axis)
		{
			QuaternionBase rotation = *this;
			if (rotation.W() > L_One)
				rotation = rotation.Normalized();
			angle = L_Two * acos(rotation.W());
			ElementType s = sqrt(L_One - rotation.W() * rotation.W());
			if (s < L_Epsilon) {
				axis.X() = L_One;
				axis.Y() = L_Zero;
				axis.Z() = L_Zero;
			}
			else {
				axis.X() = rotation.X() / s;
				axis.Y() = rotation.Y() / s;
				axis.Z() = rotation.Z() / s;
			}
		}

		Vector3 ToEuler() const
		{
			// If normalized is one, otherwise is correction factor
			ElementType unit = Base::NormSquared();
			ElementType test = X() * W() - Y() * Z();
			Vector3 v;
			// Singularity at north pole
			if (test > 0.4995f * unit)
			{
				v.Y() = L_Two * atan2(Y(), X());
				v.X() = LMath::Constans::PiOver2;
				v.Z() = L_Zero;
				return v;
			}
			// Singularity at south pole
			if (test < -0.4995f * unit)
			{
				v.Y() = -L_Two * atan2(Y(), X());
				v.X() = -LMath::Constans::PiOver2;
				v.Z() = L_Zero;
				return v;
			}
			// Yaw
			v.Y() = atan2(L_Two * W() * Y() + L_Two * Z() * X(),
				L_One - L_Two * (X() * X() + Y() * Y()));
			// Pitch
			v.X() = asin(L_Two * (W() * X() - Y() * Z()));
			// Roll
			v.Z() = atan2(L_Two * W() * Z() + L_Two * X() * Y(),
				L_One - L_Two * (Z() * Z() + X() * X()));
			return v;
		}


		QuaternionBase& operator*=(const QuaternionBase& rhs)
		{
			W() = W * rhs.W() - X() * rhs.X() - Y() * rhs.Y() - Z() * rhs.Z();
			X() = X * rhs.W() + W() * rhs.X() + Y() * rhs.Z() - Z() * rhs.Y();
			Y() = W * rhs.Y() - X() * rhs.Z() + Y() * rhs.W() + Z() * rhs.X();
			Z() = W * rhs.Z() + X() * rhs.Y() - Y() * rhs.X() + Z() * rhs.W();
			return *this;
		}

		QuaternionBase operator*(const QuaternionBase& rhs) const
		{
			QuaternionBase ret;
			ret.W() = W() * rhs.W() - X() * rhs.X() - Y() * rhs.Y() - Z() * rhs.Z();
			ret.X() = X() * rhs.W() + W() * rhs.X() + Y() * rhs.Z() - Z() * rhs.Y();
			ret.Y() = W() * rhs.Y() - X() * rhs.Z() + Y() * rhs.W() + Z() * rhs.X();
			ret.Z() = W() * rhs.Z() + X() * rhs.Y() - Y() * rhs.X() + Z() * rhs.W();
			return ret;
		}


		Vector4 operator*(const Vector4& rhs) const
		{
			//make sure W component in vector 4 is one.
			
			return Vector4(*this * ( static_cast<Vector3>(rhs) / rhs.at(3)) , L_One);
		}


		Vector3 operator*(const VectorBaseTemplate<T,3>& rhs) const
		{
			return *this * static_cast<const Vector3&>(rhs);
		}


		Vector3 operator*(const Vector3& rhs) const
		{
			Vector3 u = static_cast<Vector3>(*this);
			ElementType s = W();
			return u * (u.Dot(rhs) * L_Two)
				+ rhs * (s * s - u.NormSquared())
				+ u.Cross(rhs) * (L_Two * s);
		}

	};

	template <class T>
	const QuaternionBase<T> QuaternionBase<T>::Identity(L_Zero, L_Zero, L_Zero, L_One);
}

#endif