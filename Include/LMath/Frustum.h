#pragma once
#ifndef __LMATH_FRUSTUM_H__
#define __LMATH_FRUSTUM_H__ 1
#include "Bounds.h"
#include "Vector.h"
#include "Matrix.h"
namespace LMath
{
	template <typename T>
	struct FrustumHalfAngles
	{
		T left;
		T bottom;
		T right;
		T top;
	};

	template <typename T>
	class  FrustumBase
	{
	public:
		using ElementType = T;
		using HalfAnglesType = FrustumHalfAngles<ElementType>;
		using Rect = BoundsBase<ElementType,2>;
		using Vector2 = VectorBase<ElementType, 2>;
		using Vector3 = VectorBase<ElementType, 3>;
		using Quaternion = QuaternionBase<ElementType>;
		using Matrix4 = MatrixBase<ElementType, 4, 4>;

		void SetPosition(const Vector3& position)
		{
			fPosition = position;
		}

		void SetOrientation(const Quaternion& orientation)
		{
			fOrientation = orientation;
		}

		const Vector3& GetPosition() const
		{
			return fPosition;
		}

		const Quaternion& GetOrientation() const
		{
			return fOrientation;
		}

		/*const Vector3* worldSpaceCorners = frustum.getWorldSpaceCorners();
		const Vector3 cameraPosition = frustum.GetPosition();*/

		void SetNearClipPlane(ElementType clipPlane)
		{
			fNearClipPlane = clipPlane;
			fProjectionMatrixDirty = true;
		}

		void SetFarClipPlane(ElementType clipPlane)
		{
			fFarClipPlane = clipPlane;
			fProjectionMatrixDirty = true;
		}

		void SetHalfAngles(const HalfAnglesType& halfAngles)
		{
			fHalfAngles = halfAngles;
			fProjectionMatrixDirty = true;
		}


		void SetHalfFOVX(ElementType halfAngle, Vector2 viewportSize, ElementType referenceWidth)
		{
			//Preserve scale
			if (referenceWidth != 0)
				halfAngle = atan(std::tan(halfAngle) * viewportSize.at(0) / referenceWidth);
				

			fHalfAngles.left = halfAngle;
			fHalfAngles.right = halfAngle;


			//Preserve aspect ratio
			if (viewportSize != Vector2::Zero)
			{
				ElementType fovYHalfAngle = atan(std::tan(halfAngle) * viewportSize.at(1) / viewportSize.at(0));
				fHalfAngles.top = fovYHalfAngle;
				fHalfAngles.bottom = fovYHalfAngle;
			}
			fProjectionMatrixDirty = true;
		}

		void SetHalfFOVY(ElementType halfAngle, Vector2 viewportSize, ElementType referenceHeight)
		{
			//Preserve scale
			if (referenceHeight != 0)
				halfAngle = atan(std::tan(halfAngle) * viewportSize.at(1) / referenceHeight);

			fHalfAngles.bottom = halfAngle;
			fHalfAngles.top = halfAngle;

			//Preserve aspect ratio
			if (viewportSize != Vector2::Zero)
			{
				ElementType fovXHalfAngle = atan(std::tan(halfAngle) * viewportSize.at(0) / viewportSize.at(1));
				fHalfAngles.left = fovXHalfAngle;
				fHalfAngles.right = fovXHalfAngle;

			}
			fProjectionMatrixDirty = true;
		}

		const Matrix4& GetProjectionMatrix()
		{
			if (fProjectionMatrixDirty == true)
			{
				UpdateProjectionMatrix();
				fProjectionMatrixDirty = false;
			}
			return fProjectionMatrix;
		}

		Matrix4 GetViewMatrix()
		{
			return Matrix4::CreateViewMatrix(fPosition, fOrientation);
		}

	private:
		void UpdateProjectionMatrix()
		{
			Rect rect = CalcProjectionParameters();
			ElementType left = rect.GetMin().at(0);
			ElementType bottom = rect.GetMin().at(1);
			ElementType right = rect.GetMax().at(0);
			ElementType top = rect.GetMax().at(1);

			fProjectionMatrix =  Matrix4::CreateProjectionMatrix(left, bottom,top,right, fNearClipPlane, fFarClipPlane, -1, true);
		}

		Rect CalcProjectionParameters()
		{
			ElementType thetaLeft = std::tan(fHalfAngles.left);
			ElementType thetaRight = std::tan(fHalfAngles.right);
			ElementType omegaBottom = std::tan(fHalfAngles.bottom);
			ElementType omegaTop = std::tan(fHalfAngles.top);


			ElementType nearFocal = fNearClipPlane / fFarClipPlane;
			ElementType nearOffsetX = fFrustrumOffset.X()  * nearFocal;
			ElementType nearOffsetY = fFrustrumOffset.X() * nearFocal;
			ElementType w_left  = thetaLeft * fNearClipPlane;
			ElementType w_right = thetaRight * fNearClipPlane;
			ElementType h_bottom = omegaBottom * fNearClipPlane;
			ElementType h_top = omegaTop * fNearClipPlane;
			           
			            // left bottom             //right top
			return { { -w_left + nearOffsetX, -h_bottom + nearOffsetY } , { w_right + nearOffsetX, h_top + nearOffsetY } };
		}


	private:
		Vector3 fPosition = Vector3::Zero;
		Quaternion fOrientation = Quaternion::Identity;
		HalfAnglesType fHalfAngles;
		Vector2 fFrustrumOffset = Vector2::Zero;
		Matrix4 fProjectionMatrix = Matrix4::Zero;
		bool fKeepAspectRatio = true;
		bool fProjectionMatrixDirty = true;
		ElementType fFarClipPlane = 1000.0;
		ElementType fNearClipPlane = 1.0;
	};
}
#endif