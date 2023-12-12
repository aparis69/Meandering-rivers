#pragma once

#include <cmath>
#include <iostream>

/* Forward Declarations */
struct Vector2i;
struct Vector2;
struct Vector3;
struct Vector4;

// Maths utility
namespace Math
{
	const double Pi = 3.14159265358979323846;
	const double TwoPiOverThree = 2.0943951023931954923084;
	const double FourPiOverThree = 4.1887902047863909846168;
	const double HalfPi = Pi / 2.0f;

	inline double Angle(int k, int n)
	{
		return (2.0f * Math::Pi * k) / n;
	}

	inline double Clamp(double x, double a = 0.0, double b = 1.0)
	{
		return x < a ? a : x > b ? b : x;
	}

	inline int Clamp(int x, int a, int b)
	{
		return x < a ? a : x > b ? b : x;
	}

	template<typename T>
	inline T Min(T a, T b)
	{
		return a < b ? a : b;
	}

	template<typename T>
	inline T Max(T a, T b)
	{
		return a > b ? a : b;
	}

	inline double Abs(double a)
	{
		return a < 0 ? -a : a;
	}

	inline double CubicSmoothCompact(double x, double r)
	{
		return (x > r) ? 0.0f : (1.0f - x / r) * (1.0f - x / r) * (1.0f - x / r);
	}

	inline double CubicSmooth(double x, double r)
	{
		return (1.0f - x / r) * (1.0f - x / r) * (1.0f - x / r);
	}

	inline double Pow(double x, double e)
	{
		if (x == 0.0)
		{
			return 0.0;
		}
		else if (x > 0.0)
		{
			return pow(x, e);
		}
		else
		{
			return -pow(-x, e);
		}
	}

	inline double Sqr(double x)
	{
		return x * x;
	}

	inline bool IsNumber(double x)
	{
		return (x == x);
	}

	inline double Bilinear(double a00, double a10, double a11, double a01, double u, double v)
	{
		return (1 - u) * (1 - v) * a00 + (1 - u) * (v) * a01 + (u) * (1 - v) * a10 + (u) * (v) * a11;
	}

	inline double Fract(double x)
	{
		return x - floor(x);
	}

	inline double DegreeToRadian(double a)
	{
		return a * Math::Pi / 180.0;
	}
}

/* Vector3 */
struct Vector3
{
public:
	double x, y, z;

	explicit Vector3() : x(0.0f), y(0.0f), z(0.0f) { }
	explicit Vector3(double n) : x(n), y(n), z(n) {}
	explicit Vector3(double x, double y, double z) : x(x), y(y), z(z) {}

	friend bool operator> (const Vector3&, const Vector3&);
	friend bool operator< (const Vector3&, const Vector3&);
	friend bool operator>= (const Vector3&, const Vector3&);
	friend bool operator<= (const Vector3&, const Vector3&);

	Vector3& operator+= (const Vector3& v)
	{
		x += v.x;
		y += v.y;
		z += v.z;
		return *this;
	}
	Vector3 operator-= (const Vector3& v)
	{
		x -= v.x;
		y -= v.y;
		z -= v.z;
		return *this;
	}
	Vector3 operator*= (double f)
	{
		x *= f;
		y *= f;
		z *= f;
		return *this;
	}
	Vector3 operator/= (double f)
	{
		x /= f;
		y /= f;
		z /= f;
		return *this;
	}
	Vector3 operator*(const Vector3& u) const
	{
		return Vector3(x * u.x, y * u.y, z * u.z);
	}
	Vector3 operator*(double k) const
	{
		return Vector3(x * k, y * k, z * k);
	}
	Vector3 operator/(double k) const
	{
		return Vector3(x / k, y / k, z / k);
	}
	bool operator==(const Vector3& u) const
	{
		return (x == u.x && y == u.y && z == u.z);
	}
	bool operator!=(const Vector3& u) const
	{
		return (x != u.x || y != u.y || z != u.z);
	}
	Vector3 operator-(const Vector3& u) const
	{
		return Vector3(x - u.x, y - u.y, z - u.z);
	}
	Vector3 operator+(const Vector3& u) const
	{
		return Vector3(x + u.x, y + u.y, z + u.z);
	}
	Vector3 operator+(double k) const
	{
		return Vector3(x + k, y + k, z + k);
	}
	double operator[](int i) const
	{
		if (i == 0)
			return x;
		else if (i == 1)
			return y;
		return z;
	}
	double& operator[](int i)
	{
		if (i == 0)
			return x;
		else if (i == 1)
			return y;
		return z;
	}
	friend std::ostream& operator<<(std::ostream& stream, const Vector3& u);
	inline double Max() const
	{
		return Math::Max(Math::Max(x, y), z);
	}
	inline double Min() const
	{
		return Math::Min(Math::Min(x, y), z);
	}
	inline int MaxIndex() const
	{
		if (x >= y)
		{
			if (x >= z)
				return 0;
			else
				return 2;
		}
		else
		{
			if (y >= z)
				return 1;
			else
				return 2;
		}
	}
	Vector3 Orthogonal() const;
	void Orthonormal(Vector3& x, Vector3& y) const;
	static inline Vector3 Min(const Vector3& a, const Vector3& b)
	{
		return Vector3(Math::Min(a.x, b.x), Math::Min(a.y, b.y), Math::Min(a.z, b.z));
	}
	static inline Vector3 Max(const Vector3& a, const Vector3& b)
	{
		return Vector3(Math::Max(a.x, b.x), Math::Max(a.y, b.y), Math::Max(a.z, b.z));
	}
};
inline std::ostream& operator<<(std::ostream& stream, const Vector3& u)
{
	stream << "(" << u.x << ", " << u.y << ", " << u.z << ");";
	return stream;
}
inline Vector3 Cross(const Vector3& u, const Vector3& v)
{
	return Vector3((u.y * v.z) - (u.z * v.y), (u.z * v.x) - (u.x * v.z), (u.x * v.y) - (u.y * v.x));
}
inline double Dot(const Vector3& u, const Vector3& v)
{
	return u.x * v.x + u.y * v.y + u.z * v.z;
}
inline double Magnitude(const Vector3& u)
{
	return sqrt(u.x * u.x + u.y * u.y + u.z * u.z);
}
inline double SquaredMagnitude(const Vector3& u)
{
	return u.x * u.x + u.y * u.y + u.z * u.z;
}
inline Vector3 Normalize(const Vector3& v)
{
	double kk = 1.0f / Magnitude(v);
	return v * kk;
}
inline Vector3 operator-(const Vector3& v)
{
	return Vector3(-v.x, -v.y, -v.z);
}
inline bool operator>(const Vector3& u, const Vector3& v)
{
	return (u.x > v.x) && (u.y > v.y) && (u.z > v.z);
}
inline bool operator<(const Vector3& u, const Vector3& v)
{
	return (u.x < v.x) && (u.y < v.y) && (u.z < v.z);
}
inline bool operator>=(const Vector3& u, const Vector3& v)
{
	return (u.x >= v.x) && (u.y >= v.y) && (u.z >= v.z);
}
inline bool operator<=(const Vector3& u, const Vector3& v)
{
	return (u.x <= v.x) && (u.y <= v.y) && (u.z <= v.z);
}
inline Vector3 operator*(double a, const Vector3& v)
{
	return v * a;
}
inline Vector3 Abs(const Vector3& u)
{
	return Vector3(u[0] > 0.0 ? u[0] : -u[0], u[1] > 0.0 ? u[1] : -u[1], u[2] > 0.0 ? u[2] : -u[2]);
}
inline Vector3 Vector3::Orthogonal() const
{
	Vector3 a = Abs(*this);
	int i = 0;
	int j = 1;
	if (a[0] > a[1])
	{
		if (a[2] > a[1])
		{
			j = 2;
		}
	}
	else
	{
		i = 1;
		j = 2;
		if (a[0] > a[2])
		{
			j = 0;
		}
	}
	a = Vector3(0);
	a[i] = operator[](j);
	a[j] = -operator[](i);
	return a;
}
inline void Vector3::Orthonormal(Vector3& x, Vector3& y) const
{
	x = Normalize(Orthogonal());
	y = Normalize(Cross(*this, x));
}
inline Vector3 Bilinear(const Vector3& a00, const Vector3& a10, const Vector3& a11, const Vector3& a01, double u, double v)
{
	return (1 - u) * (1 - v) * a00 + (1 - u) * (v)*a01 + (u) * (1 - v) * a10 + (u) * (v)*a11;
}

/* Vector2 */
struct Vector2
{
public:
	double x, y;

	explicit Vector2() : x(0.0), y(0.0) { }
	explicit Vector2(double n) : x(n), y(n) { }
	explicit Vector2(double x, double y) : x(x), y(y) { }
	explicit Vector2(const Vector3& v) : x(v.x), y(v.z) { }

	friend bool operator> (const Vector2&, const Vector2&);
	friend bool operator< (const Vector2&, const Vector2&);
	friend bool operator>= (const Vector2&, const Vector2&);
	friend bool operator<= (const Vector2&, const Vector2&);

	Vector2 operator+= (const Vector2& v)
	{
		x += v.x;
		y += v.y;
		return *this;
	}
	Vector2 operator-= (const Vector2& v)
	{
		x -= v.x;
		y -= v.y;
		return *this;
	}
	Vector2 operator*= (double f)
	{
		x *= f;
		y *= f;
		return *this;
	}
	Vector2 operator/= (double f)
	{
		x /= f;
		y /= f;
		return *this;
	}
	Vector2 operator*(const Vector2& v) const
	{
		return Vector2(x * v.x, y * v.y);
	}
	Vector2 operator*(double k) const
	{
		return Vector2(x * k, y * k);
	}
	Vector2 operator/(double k) const
	{
		return Vector2(x / k, y / k);
	}
	Vector2 operator-(const Vector2& u) const
	{
		return Vector2(x - u.x, y - u.y);
	}
	Vector2 operator+(const Vector2& u) const
	{
		return Vector2(x + u.x, y + u.y);
	}
	Vector2 operator+(double k) const
	{
		return Vector2(x + k, y + k);
	}
	Vector2 operator-(double k) const
	{
		return Vector2(x - k, y - k);
	}
	double operator[](int i) const
	{
		if (i == 0)
			return x;
		return y;
	}
	double& operator[](int i)
	{
		if (i == 0)
			return x;
		return y;
	}
	friend std::ostream& operator<< (std::ostream& stream, const Vector2& u);
	inline Vector3 ToVector3(double yy) const
	{
		return Vector3(x, yy, y);
	}
	inline double Max() const
	{
		return Math::Max(x, y);
	}
	inline double Min() const
	{
		return Math::Min(x, y);
	}
	inline Vector2 Orthogonal() const
	{
		return Vector2(-y, x);
	}
};
inline double operator/ (const Vector2& u, const Vector2& v)
{
	return u.x * v.y - u.y * v.x;
}
inline std::ostream& operator<<(std::ostream& stream, const Vector2& u)
{
	stream << "(" << u.x << ", " << u.y << ");";
	return stream;
}
inline double Dot(const Vector2& u, const Vector2& v)
{
	return u.x * v.x + u.y * v.y;
}
inline double Magnitude(const Vector2& u)
{
	return sqrt(u.x * u.x + u.y * u.y);
}
inline double SquaredMagnitude(const Vector2& u)
{
	return u.x * u.x + u.y * u.y;
}
inline Vector2 Normalize(const Vector2& v)
{
	double kk = 1.0f / Magnitude(v);
	return v * kk;
}
inline Vector2 operator-(const Vector2& v)
{
	return Vector2(-v.x, -v.y);
}
inline bool operator>(const Vector2& u, const Vector2& v)
{
	return (u.x > v.x) && (u.y > v.y);
}
inline bool operator<(const Vector2& u, const Vector2& v)
{
	return (u.x < v.x) && (u.y < v.y);
}
inline bool operator>=(const Vector2& u, const Vector2& v)
{
	return (u.x >= v.x) && (u.y >= v.y);
}
inline bool operator<=(const Vector2& u, const Vector2& v)
{
	return (u.x <= v.x) && (u.y <= v.y);
}
inline Vector2 Abs(const Vector2& u)
{
	return Vector2(u[0] > 0.0 ? u[0] : -u[0], u[1] > 0.0 ? u[1] : -u[1]);
}
inline Vector2 operator*(double k, const Vector2& p)
{
	return Vector2(p.x * k, p.y * k);
}
inline bool operator== (const Vector2& u, const Vector2& v)
{
	return ((u.x == v.x) && (u.y == v.y));
}
inline bool operator!= (const Vector2& u, const Vector2& v)
{
	return (!(u == v));
}
inline Vector2 Lerp(const Vector2& a, const Vector2& b, double t)
{
	return a + (b - a) * t;
}
inline Vector2 Bilinear(const Vector2& a00, const Vector2& a10, const Vector2& a11, const Vector2& a01, double u, double v)
{
	return (1 - u) * (1 - v) * a00 + (1 - u) * (v)*a01 + (u) * (1 - v) * a10 + (u) * (v)*a11;
}


/* Matrix2 */
struct Matrix2
{
public:
	double c[4];

	inline Matrix2() { }
	inline Matrix2(double x0, double x1, double x2, double x3)
	{
		c[0] = x0;
		c[1] = x1;
		c[2] = x2;
		c[3] = x3;
	}
	inline Vector2 operator*(const Vector2& v) const
	{
		return Vector2(v[0] * c[0] + v[1] * c[2], v[0] * c[1] + v[1] * c[3]);
	}
};
inline Matrix2 RotationMatrix(double a)
{
	double c = cos(a);
	double s = sin(a);
	return Matrix2(c, s, -s, c);
}
