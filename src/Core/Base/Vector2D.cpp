#include "Core/Base/Vector2D.h"

using namespace std;

Vector2D::Vector2D() {}
Vector2D::Vector2D(double _x, double _y) :x(_x), y(_y) {}
Vector2D::~Vector2D() {}

void Vector2D::Set(const double _x, const double _y)
{
	this->x = _x;
	this->y = _y;
}

Vector2D Vector2D::operator= (const Vector2D& v)
{
	this->x = v.x;
	this->y = v.y;
	return *this;
}

Vector2D Vector2D::operator+ (const Vector2D& v) const
{
	Vector2D sum(this->x + v.x, this->y + v.y);
	return sum;
}

Vector2D Vector2D::operator- (const Vector2D& v) const
{
	Vector2D sum(this->x - v.x, this->y - v.y);
	return sum;
}

Vector2D Vector2D::operator*(const double& k) const
{
	Vector2D result(this->x * k, this->y * k);
	return result;
}

 Vector2D& Vector2D::operator+= (const Vector2D& v)
{
	this->x += v.x;
	this->y += v.y;
	return *this;
}

Vector2D& Vector2D::operator-= (const Vector2D& v)
{
	this->x -= v.x;
	this->y -= v.y;
	return *this;
}

Vector2D& Vector2D::operator*= (const double& k)
{
	this->x *= k;
	this->y *= k;
	return *this;
}

double Vector2D::dotProduct(const Vector2D& v)
{
	return this->x * v.x + this->y * v.y;;
}

double Vector2D::crossProduct(const Vector2D& v)
{
	return this->x * v.y - this->y * v.x;
}

void Vector2D::normalize()
{
	double mod_inv = 1.0 / sqrt(this->x * this->x + this->y * this->y);
	this->x = this->x * mod_inv;
	this->y = this->y * mod_inv;
}

double Vector2D::mod() const
{
	return sqrt(this->x * this->x + this->y * this->y);
}

double Vector2D::mod2() const
{
	return this->x * this->x + this->y * this->y;
}

std::ostream& operator<< (std::ostream& os, const Vector2D& v)
{
	os << "(" << v.x << ", " << v.y << ")";
	return os;
}