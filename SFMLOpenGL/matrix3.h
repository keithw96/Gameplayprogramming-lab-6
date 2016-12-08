#pragma once
#include <iostream>
#include <vector>
#include <iostream>


struct vector3
{
	float X, Y, Z;

	inline vector3(void) {}
	inline vector3(const float x, const float y, const float z)
	{
		X = x; Y = y; Z = z;
	}

	inline vector3 operator + (const vector3& A) const
	{
		return vector3(X + A.X, Y + A.Y, Z + A.Z);
	}

	inline vector3 operator + (const float A) const
	{
		return vector3(X + A, Y + A, Z + A);
	}

	inline float Dot(const vector3& A) const
	{
		return A.X*X + A.Y*Y + A.Z*Z;
	}
};

class Matrix3
{
	
public:
	double A11;
	double A12;
	double A13;
	double A21;
	double A22;
	double A23;
	double A31;
	double A32;
	double A33;

	vector3 m_row1;
	vector3 m_row2;
	vector3 m_row3;

	Matrix3(vector3 row1, vector3 row2, vector3 row3);

	Matrix3();

	Matrix3(double _A11, double _A12, double _A13, double _A21, double _A22, double _A23, double _A31, double _A32 , double _A33);

	~Matrix3();
private:

};

Matrix3::~Matrix3()
{
}
