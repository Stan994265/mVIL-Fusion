
#ifndef  _QUATERNION_H_
#define _QUATERNION_H_

#include "Vector3.h"

#include <iostream>
#include <vector>


namespace alg {

	/*
	* The Unit Quaternion is one possible representation of the
	* attitude of an object in tree-dimensional space.
	*
	* This Quaternion class is implemented according to Diebel,
	* James. Representing Attitude: Euler Angle, Unit Quaternions, and
	* Rotation Vectors. Stanford University. 2006. - Technical Report.
	*/

	class Quaternion {

	public:

		/*!
		* \brief Default constructor
		*
		* Constructs the (1,0,0,0) Unit Quaternion
		* representing the identity rotation.
		*/
		inline Quaternion() { u() = 1;  x() = 0; y() = 0; z() = 0;  }

		/*!
		* \brief Copy constructor
		*/
		Quaternion(const Quaternion& other){
			data[0] = other(0);
			data[1] = other(1);
			data[2] = other(2);
			data[3] = other(3);
		}

		/*!
		* \brief Constructor
		*
		* Constructs a Quaternion from four single
		* values
		*/
		Quaternion(float uu, float xx, float yy, float zz) {
			u() = uu;
			x() = xx;
			y() = yy;
			z() = zz;
		}

		/*!
		* \brief Constructor
		*
		* @param other a vector containing euler angles
		*/
		Quaternion(const Vector3 &other) {
			operator=(Quaternion(other.roll(), other.pitch(), other.yaw()));
		}

		/*!
		* \brief Constructor from Euler angles
		*
		* Constructs a Unit Quaternion from Euler angles / Tait Bryan
		* angles (in radians) according to the 1-2-3 convention.
		* @param roll phi/roll angle (rotation about x-axis)
		* @param pitch theta/pitch angle (rotation about y-axis)
		* @param yaw psi/yaw angle (rotation about z-axis)
		*/
		Quaternion(double roll, double pitch, double yaw) {
			double sroll   = sin(roll);
			double spitch = sin(pitch);
			double syaw   = sin(yaw);
			double croll   = cos(roll);
			double cpitch = cos(pitch);
			double cyaw   = cos(yaw);

			double m[3][3] = { //create rotational Matrix
				{cyaw*cpitch, cyaw*spitch*sroll - syaw*croll, cyaw*spitch*croll + syaw*sroll},
				{syaw*cpitch, syaw*spitch*sroll + cyaw*croll, syaw*spitch*croll - cyaw*sroll},
				{    -spitch,                  cpitch*sroll,                  cpitch*croll}
			};

			float _u = (float) (sqrt(std::max(0., 1 + m[0][0] + m[1][1] + m[2][2]))/2.0);
			float _x = (float) (sqrt(std::max(0., 1 + m[0][0] - m[1][1] - m[2][2]))/2.0);
			float _y = (float) (sqrt(std::max(0., 1 - m[0][0] + m[1][1] - m[2][2]))/2.0);
			float _z = (float) (sqrt(std::max(0., 1 - m[0][0] - m[1][1] + m[2][2]))/2.0);
			u() = _u;
			x() = (m[2][1] - m[1][2])>=0?fabs(_x):-fabs(_x);
			y() = (m[0][2] - m[2][0])>=0?fabs(_y):-fabs(_y);
			z() = (m[1][0] - m[0][1])>=0?fabs(_z):-fabs(_z);
		}
		/*
		* construct from 3*3 rotation maxtrix
		*/
		Quaternion(double q00, double q01, double q02,
				double q10, double q11, double q12,
				double q20, double q21, double q22) 
		{
			float _u = (float) (sqrt(std::max(0., 1 + q00 + q11+ q22))/2.0);
			float _x = (float) (sqrt(std::max(0., 1 + q00 - q11 - q22))/2.0);
			float _y = (float) (sqrt(std::max(0., 1 - q00 + q11 - q22))/2.0);
			float _z = (float) (sqrt(std::max(0., 1 - q00 - q11 + q22))/2.0);
			u() = _u;
			x() = (q21 - q12)>=0?fabs(_x):-fabs(_x);
			y() = (q02- q20)>=0?fabs(_y):-fabs(_y);
			z() = (q10 - q01)>=0?fabs(_z):-fabs(_z);
		}

		//! Constructs a Unit Quaternion from a rotation angle and axis.  
		Quaternion(const Vector3& axis, double angle) {
			double sa = sin(angle/2);
			double ca = cos(angle/2);
			x() = (float) (axis.x()*sa);
			y() = (float) (axis.y()*sa);
			z() = (float) (axis.z()*sa);
			u() = (float) ca;
		}

		/*!
		* \brief Conversion to Euler angles
		*
		* Converts the attitude represented by this to
		* Euler angles (roll, pitch, yaw).
		*/
		Vector3 toEuler() const {
			// create rotational matrix
			double n = norm ();
			double s = n > 0?2./(n*n):0.;

			double xs = x()*s;
			double ys = y()*s;
			double zs = z()*s;

			double ux = u()*xs;
			double uy = u()*ys;
			double uz = u()*zs;

			double xx = x()*xs;
			double xy = x()*ys;
			double xz = x()*zs;

			double yy = y()*ys;
			double yz = y()*zs;
			double zz = z()*zs;

			double m[3][3];

			m[0][0] = 1.0 - (yy + zz);
			m[1][1] = 1.0 - (xx + zz);
			m[2][2] = 1.0 - (xx + yy);

			m[1][0] = xy + uz;
			m[0][1] = xy - uz;

			m[2][0] = xz - uy;
			m[0][2] = xz + uy;
			m[2][1] = yz + ux;
			m[1][2] = yz - ux;

			float roll  = (float) atan2(m[2][1], m[2][2]);
			float pitch = (float) atan2(-m[2][0], sqrt(m[2][1]*m[2][1] + m[2][2]*m[2][2]));
			float yaw   = (float) atan2(m[1][0], m[0][0]);

			return Vector3(roll, pitch, yaw);
		}

		void toRotMatrix(std::vector <double>& rot_matrix_3_3) const {

			// create rotational matrix
			double n = norm ();
			double s = n > 0?2./(n*n):0.;

			double xs = x()*s;
			double ys = y()*s;
			double zs = z()*s;

			double ux = u()*xs;
			double uy = u()*ys;
			double uz = u()*zs;

			double xx = x()*xs;
			double xy = x()*ys;
			double xz = x()*zs;

			double yy = y()*ys;
			double yz = y()*zs;
			double zz = z()*zs;

			double m[3][3];
			m[0][0] = 1.0 - (yy + zz);
			m[1][1] = 1.0 - (xx + zz);
			m[2][2] = 1.0 - (xx + yy);

			m[1][0] = xy + uz;
			m[0][1] = xy - uz;

			m[2][0] = xz - uy;
			m[0][2] = xz + uy;
			m[2][1] = yz + ux;
			m[1][2] = yz - ux;

			rot_matrix_3_3.clear();
			rot_matrix_3_3.resize(9,0.);
			for (unsigned int i=0; i<3; i++) {
				rot_matrix_3_3[i*3] = m[i][0];
				rot_matrix_3_3[i*3+1] = m[i][1];
				rot_matrix_3_3[i*3+2] = m[i][2];
			}
		}

		inline const float& operator() (unsigned int i) const { return data[i]; }
		inline float& operator() (unsigned int i) { return data[i]; }

		float norm () const {
			double n = 0;
			for (unsigned int i=0; i<4; i++) {
				n += operator()(i) * operator()(i);
			}
			return (float) sqrt(n);
		}
		Quaternion& normalize (){
			double len = norm ();
			if (len > 0)
				*this /= (float) len;
			return *this;
		}

		Quaternion normalized () const {
			Quaternion result(*this);
			result.normalize ();
			return result;
		}

		void operator/= (float x) 
		{ 
			for (unsigned int i=0; i<4; ++i)
				operator()(i) /= x;
		}

		Quaternion& operator= (const Quaternion& other) {
			u() = other.u();
			x() = other.x();
			y() = other.y();
			z() = other.z();
			return *this;
		}
		bool operator== (const Quaternion& other) const 
		{
			for (unsigned int i=0; i<4; i++) 
			{
				if (operator()(i) != other(i)) 
					return false;
			}
			return true;
		}

		/*!
		* \brief Quaternion multiplication
		*
		* Standard Quaternion multiplication which is not
		* commutative.
		* @return this * other
		*/
		Quaternion operator* (const Quaternion& other) const {
			return Quaternion(u()*other.u() - x()*other.x() - y()*other.y() - z()*other.z(),
				y()*other.z() - other.y()*z() + u()*other.x() + other.u()*x(),
				z()*other.x() - other.z()*x() + u()*other.y() + other.u()*y(),
				x()*other.y() - other.x()*y() + u()*other.z() + other.u()*z());
		}
		/*!
		* \brief Quaternion multiplication with extended vector
		*
		* @return q * (0, v)
		*/
		Quaternion operator* (const Vector3& v) const {
			return *this * Quaternion(0, v(0), v(1), v(2));
		}

		/*!
		* \brief Inversion
		*
		* @return A copy of this Quaterion inverted
		*/
		inline Quaternion inv() const {  return Quaternion(u(), -x(), -y(), -z()); }

		/*!
		* \brief Inversion
		*
		* Inverts this Quaternion
		* @return a reference to this Quaternion
		*/
		Quaternion& inv_IP() {
			x() = -x();
			y() = -y();
			z() = -z();
			return *this;
		}
		/*!
		* \brief Rotate a vector
		*
		* Rotates a vector to the body fixed coordinate
		* system according to the attitude represented by
		* this Quaternion.
		* @param v a vector represented in world coordinates
		* @return v represented in body-fixed coordinates
		*/
		Vector3 rotate(const Vector3& v) const {
			Quaternion q = *this * v * this->inv();
			return Vector3(q.x(), q.y(), q.z());
		}

		inline float& u() { return data[0]; }
		inline float& x() { return data[1]; }
		inline float& y() { return data[2]; }
		inline float& z() { return data[3]; }

		inline const float& u() const { return data[0]; }
		inline const float& x() const { return data[1]; }
		inline const float& y() const { return data[2]; }
		inline const float& z() const { return data[3]; }

		std::istream& read(std::istream &s) {
			int temp;
			s >> temp; // should be 4
			for (unsigned int i=0; i<4; i++)
				s >> operator()(i);
			return s;
		}


		std::ostream& write(std::ostream &s) const {
			s << 4;
			for (unsigned int i=0; i<4; i++)
				s << " " << operator()(i);
			return s;
		}



		std::istream& readBinary(std::istream &s) {
			int temp;
			s.read((char*)&temp, sizeof(temp));
			double val = 0;
			for (unsigned int i=0; i<4; i++) {
				s.read((char*)&val, sizeof(val));
				operator()(i) = (float) val;
			}
			return s;
		}


		std::ostream& writeBinary(std::ostream &s) const {
			int temp = 4;
			s.write((char*)&temp, sizeof(temp));
			double val = 0;
			for (unsigned int i=0; i<4; i++) {
				val = operator()(i);
				s.write((char*)&val, sizeof(val));
			}
			return s;
		}
	protected:
		float data[4];

	};

	//! user friendly output in format (u x y z)
	std::ostream& operator<<(std::ostream& s, const Quaternion& q) {
		s << "(" << q.u() << " " << q.x() << " " << q.y() << " " << q.z() << ")";
		return s;
	}
}

#endif
