
#ifndef  _POSE6D_H_
#define _POSE6D_H_

#include "Vector3.h"
#include "Quaternion.h"

namespace alg {

	/*!
	* \brief This class represents a tree-dimensional pose of an object
	*
	* The tree-dimensional pose is represented by a three-dimensional
	* translation vector representing the position of the object and
	* a Quaternion representing the attitude of the object
	*/
	class Pose6D {
	public:

		Pose6D():
			translation(0.0,0.0,0.0),
			rotation(0.0,0.0,0.0)
		{}
		~Pose6D(){}

		/*
		* Constructs a pose from given translation and rotation.
		*/
		Pose6D(const Vector3& trans, const Quaternion& rot):
			translation(trans),
			rotation(rot)
		{}

		/*
		* Constructs a pose from a translation represented by
		* its x, y, z-values and a rotation represented by its
		* Tait-Bryan angles roll, pitch, and yaw
		*/
		Pose6D(float x, float y, float z, double roll, double pitch, double yaw):
			translation(x, y, z),
			rotation(roll, pitch, yaw)
		{ }
		// Pose6D(float x, float y, float z, double q00, double q01, double q02,
		// 		double q10, double q11, double q12,
		// 		double q20, double q21, double q22):
		// 	translation(x, y, z),
		// 	rotation(q00, q01, q02, q10, q11, q12, q20, q21, q22)
		// { }

		Pose6D& operator= (const Pose6D& other) {
			translation = other.trans();
			rotation = other.rot();
			return *this;
		}
		bool operator==(const Pose6D& other) const {
			return translation == other.translation
				&& rotation == other.rotation;
		}

		bool operator!=(const Pose6D& other) const {
			return !(*this == other);
		}


		/*
		* @return the translational component of this pose
		*/
		inline Vector3& trans() { return translation; }
		/*
		* @return the rotational component of this pose
		*/
		inline Quaternion& rot() { return rotation; }
		/*
		* @return the translational component of this pose
		*/
		const Vector3& trans() const { return translation; }
		/*
		* @return the rotational component of this pose
		*/
		const Quaternion& rot() const { return rotation; }


		inline float& x() { return translation(0); }
		inline float& y() { return translation(1); }
		inline float& z() { return translation(2); }
		inline const float& x() const { return translation(0); }
		inline const float& y() const { return translation(1); }
		inline const float& z() const { return translation(2); }

		inline double roll()  const {return (rotation.toEuler())(0); }
		inline double pitch() const {return (rotation.toEuler())(1); }
		inline double yaw()   const {return (rotation.toEuler())(2); }

		/*
		* Transforms the vector v by the transformation which is
		* specified by this.
		* @return the vector which is translated by the translation of
		* this and afterwards rotated by the rotation of this.
		*/
		Vector3 transform (const Vector3 &v) const {
			Vector3 res = this->rot().rotate(v);
			res = res + this->trans();
			return res;
		}
		/*
		* Inverts the coordinate transformation represented by this pose
		* @return a copy of this pose inverted
		*/
		Pose6D inv() const {
			Pose6D result(*this);
			result.rot() = rot().inv().normalized();
			result.trans() = result.rot().rotate(-trans());
			return result;
		}
		/*
		* Inverts the coordinate transformation represented by this pose
		* @return a reference to this pose
		*/
		
		Pose6D& inv_IP() {
			rot() = rot().inv().normalized();
			trans() = rot().rotate(-trans());
			return *this;
		}

		/*!
		* \brief Concatenation
		*
		* Concatenates the coordinate transformations represented
		* by this and p.
		* @return this * p (applies first this, then p)
		*/
		Pose6D operator* (const Pose6D& other) const {
			Quaternion rot_new   = rotation * other.rot();
			Vector3    trans_new = rotation.rotate(other.trans()) + trans();
			return Pose6D(trans_new, rot_new.normalized());
		}
		/*
		* Concatenates p to this Pose6D.
		* @return this which results from first moving by this and
		* afterwards by p
		*/
		const Pose6D& operator*= (const Pose6D& other) {
			trans() += rotation.rotate(other.trans());
			rot() = rot() * other.rot();
			return *this;
		}

		/*!
		* \brief Translational distance
		*
		* @return the translational (euclidian) distance to p
		*/
		double distance (const Pose6D &other) const {
			double dist_x = x() - other.x();
			double dist_y = y() - other.y();
			double dist_z = z() - other.z();
			return sqrt(dist_x*dist_x + dist_y*dist_y + dist_z*dist_z);
		}

		/*!
		* \brief Translational length
		*
		* @return the translational (euclidian) length of the translation
		* vector of this Pose6D
		*/
		double transLength() const {
			return sqrt(x()*x() + y()*y() + z()*z());
		}

		/*!
		* \brief Output operator
		*
		* Output to stream in a format which can be parsed using read().
		*/
		std::ostream& write(std::ostream &s) const {
			translation.write(s);
			s << " ";
			rotation.write(s);
			return s;
		}
		/*!
		* \brief Input operator
		*
		* Parsing from stream which was written by write().
		*/
		std::istream& read(std::istream &s) {
			translation.read(s);
			rotation.read(s);
			return s;
		}
		/*!
		* \brief Binary output operator
		*
		* Output to stream in a binary format which can be parsed using readBinary().
		*/
		std::ostream& writeBinary(std::ostream &s) const {
			translation.writeBinary(s);
			rotation.writeBinary(s);
			return s;
		}
		/*!
		* \brief Binary input operator
		*
		* Parsing from binary stream which was written by writeBinary().
		*/
		std::istream& readBinary(std::istream &s) {
			translation.readBinary(s);
			rotation.readBinary(s);
			return s;
		}

	protected:
		Vector3 translation;
		Quaternion rotation;
	};

	//! user friendly output in format (x y z, u x y z) which is (translation, rotation)
	std::ostream& operator<<(std::ostream& s, const Pose6D& p) {
		s << "(" << p.x() << " " << p.y() << " " << p.z()
			<< ", " << p.rot().u() << " " << p.rot().x() << " " << p.rot().y() << " " << p.rot().z() << ")";
		return s;
	}
}

#endif
