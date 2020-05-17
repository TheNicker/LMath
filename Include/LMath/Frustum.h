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

	enum ProjectionType
	{
		Orthographic,
		Perspective
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
		using Matrix3 = MatrixBase<ElementType, 3, 3>;
		using Matrix4 = MatrixBase<ElementType, 4, 4>;
		
		void SetProjectionType()
		{
			throw std::logic_error("Not implemented");
		}


		ProjectionType GetProjectionType() const
		{
			return fProjectionType;
		}

		void SetPosition(const Vector3& position)
		{
			fPosition = position;
			fViewMatrixDirty = true;
			fWorldSpaceCornersDirty = true;
		}

		void SetOrientation(const Quaternion& orientation)
		{
			fOrientation = orientation;
			fViewMatrixDirty = true;
			fWorldSpaceCornersDirty = true;
		}


		const Vector3 GetDirection() const
		{
			return fOrientation * Vector3::Forward;
		}

		const std::array<Vector3,8>& GetWorldSpaceCorners() const
		{

			if (fWorldSpaceCornersDirty == true)
			{
				using NonConstType = std::remove_const_t<std::remove_reference_t<decltype(*this)>>;
				const_cast<NonConstType*>(this)->UpdateWorldSpaceCorners();
				fWorldSpaceCornersDirty = false;
			}
			
			return fWorldSpaceCorners;
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
			fWorldSpaceCornersDirty = true;
		}

		void SetFarClipPlane(ElementType clipPlane)
		{
			fFarClipPlane = clipPlane;
			fProjectionMatrixDirty = true;
			fWorldSpaceCornersDirty = true;
		}

		void SetHalfAngles(const HalfAnglesType& halfAngles)
		{
			fHalfAngles = halfAngles;
			fProjectionMatrixDirty = true;
			fWorldSpaceCornersDirty = true;
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
			fWorldSpaceCornersDirty = true;
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
			fWorldSpaceCornersDirty = true;
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
			if (fViewMatrixDirty == true)
			{
				fViewMatrix = Matrix4::CreateViewMatrix(fPosition, fOrientation);
				fViewMatrixDirty = false;
			}
			return fViewMatrix;
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

		void UpdateWorldSpaceCorners()
		{
			Matrix4 eyeToWorld = fViewMatrix.Inverse();

			// Note: Even though we can dealing with general projection matrix here,
			//       but because it's incompatibly with infinite far plane, thus, we
			//       still need to working with projection parameters.

			// Calc near plane corners
			
			Rect vp = CalcProjectionParameters();


			
			ElementType nearLeft = vp.GetMin().at(0);
			ElementType nearRight = vp.GetMax().at(0);
			ElementType nearBottom = vp.GetMin().at(1);
			ElementType nearTop = vp.GetMax().at(1);


			// Treat infinite fardist as some arbitrary far value
			ElementType farDist = (fFarClipPlane == 0) ? 100000 : fFarClipPlane;

			// Calc far palne corners
			ElementType radio = fProjectionType == ProjectionType::Perspective ? farDist / fNearClipPlane : 1;
			ElementType farLeft = nearLeft * radio;
			ElementType farRight = nearRight * radio;
			ElementType farBottom = nearBottom * radio;
			ElementType farTop = nearTop * radio;


			fWorldSpaceCorners[0] = eyeToWorld * Vector3(nearRight, nearTop, -fNearClipPlane);
			fWorldSpaceCorners[1] = eyeToWorld * Vector3(nearLeft, nearTop, -fNearClipPlane);
			fWorldSpaceCorners[2] = eyeToWorld * Vector3(nearLeft, nearBottom, -fNearClipPlane);
			fWorldSpaceCorners[3] = eyeToWorld * Vector3(nearRight, nearBottom, -fNearClipPlane);
			// far
			fWorldSpaceCorners[4] = eyeToWorld * Vector3(farRight, farTop, -farDist);
			fWorldSpaceCorners[5] = eyeToWorld * Vector3(farLeft, farTop, -farDist);
			fWorldSpaceCorners[6] = eyeToWorld * Vector3(farLeft, farBottom, -farDist);
			fWorldSpaceCorners[7] = eyeToWorld * Vector3(farRight, farBottom, -farDist);


		}


	private:
		Vector3 fPosition = Vector3::Zero;
		Quaternion fOrientation = Quaternion::Identity;
		HalfAnglesType fHalfAngles;
		Vector2 fFrustrumOffset = Vector2::Zero;
		Matrix4 fProjectionMatrix = Matrix4::Zero;
		Matrix4 fViewMatrix = Matrix4::Identity;
		bool fKeepAspectRatio = true;
		bool fProjectionMatrixDirty = true;
		mutable bool fWorldSpaceCornersDirty = true;
		bool fViewMatrixDirty = true;

		ElementType fFarClipPlane = 1000.0;
		ElementType fNearClipPlane = 1.0;
		std::array<Vector3, 8> fWorldSpaceCorners;
		ProjectionType fProjectionType = ProjectionType::Perspective;
	};
}
#endif