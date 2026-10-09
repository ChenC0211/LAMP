#pragma once

#include <iostream>
#include <cmath>

class Vector2D
{
public:
	double x = 0.0, y = 0.0;

	Vector2D();
	Vector2D(double _x, double _y);
	~Vector2D();

	void Set(const double _x, const double _y);

	Vector2D operator= (const Vector2D& v);

	Vector2D operator+ (const Vector2D& v) const;
	Vector2D operator- (const Vector2D& v) const;
	Vector2D operator* (const double& k) const;

	Vector2D& operator+= (const Vector2D& v);
	Vector2D& operator-= (const Vector2D& v);
	Vector2D& operator*= (const double& k);

	// 叉乘
	double crossProduct(const Vector2D& v);
	// 点乘
	double dotProduct(const Vector2D& v);

	// 归一化
	void normalize();
	// 求模
	double mod() const;
	// 求模的平方
	double mod2() const;
};

std::ostream& operator<< (std::ostream& os, const Vector2D& v);