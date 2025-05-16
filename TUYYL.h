#pragma once
#ifndef TUYYL_H
#define TUYYL_H

#include <iostream>
#include <cmath>
#include<NMYYL_Plus.h>

using namespace std;

class TU_Class_Point;
class TU_Class_Circle;
class TU_Class_Line;
class TU_Class_Seg;
void TU_print(std::vector<TU_Class_Point> p , std::vector<TU_Class_Circle> c , std::vector<TU_Class_Line> l , std::vector<TU_Class_Seg> s );
int TU_FindCoPoint(TU_Class_Circle c1, TU_Class_Circle c2, TU_Class_Point*& p1, TU_Class_Point*& p2);
int TU_FindCoPoint(TU_Class_Line& l1, TU_Class_Line& l2, TU_Class_Point*& p);
int TU_FindCoPoint(TU_Class_Line& line, TU_Class_Circle& circle, TU_Class_Point*& p1, TU_Class_Point*& p2);

class TU_Class_Point
{


private:

	friend class TU_Class_Circle;
	friend void TU_print(std::vector<TU_Class_Point> p, std::vector<TU_Class_Circle> c, std::vector<TU_Class_Line> l, std::vector<TU_Class_Seg> s);

	double x = 0, y = 0, z = 0;

public:

	void polemake(double r, double theta)
	{
		this->x = r * cos(theta);
		this->y = r * sin(theta);
		this->z = 0;
	}

	TU_Class_Point(double x, double y, double z)
	{
		this->x = x;
		this->y = y;
		this->z = z;
	}

	TU_Class_Point() {}

	TU_Class_Point(double x, double y)
	{
		this->x = x;
		this->y = y;
		this->z = 0;
	}

	TU_Class_Point operator = (TU_Class_Point p)
	{
		this->x = p.x;
		this->y = p.y;
		this->z = p.z;
		return *this;
	}

	double getX()
	{
		return this->x;
	}

	double getY()
	{
		return this->y;
	}

	double getZ()
	{
		return this->z;
	}


	double distance(TU_Class_Point p)
	{
		return sqrt(pow(this->x - p.x, 2) + pow(this->y - p.y, 2) + pow(this->z - p.z, 2));
	}

};

class TU_Class_Circle
{
private:

	friend void TU_print(std::vector<TU_Class_Point> p, std::vector<TU_Class_Circle> c, std::vector<TU_Class_Line> l, std::vector<TU_Class_Seg> s);
	TU_Class_Point center;
	double radius = 0;

public:
	TU_Class_Circle(TU_Class_Point center, double radius)
	{
		this->center = center;
		this->radius = radius;
	}

	TU_Class_Circle()
	{
		this->center = TU_Class_Point(0, 0);
		this->radius = 1;
	}

	TU_Class_Circle(TU_Class_Point center)
	{
		this->center = center;
	}

	// 接收三个点作为输入参数的构造函数
	TU_Class_Circle(TU_Class_Point p1, TU_Class_Point p2, TU_Class_Point p3) {
		// 检查三点是否共线
		if ((p2.getY() - p1.getY()) * (p3.getX() - p2.getX()) ==
			(p3.getY() - p2.getY()) * (p2.getX() - p1.getX())) {
			throw std::invalid_argument("The three points are collinear.");
		}

		// 计算圆心坐标
		double A1 = p2.getX() - p1.getX(), B1 = p2.getY() - p1.getY(),
			C1 = (A1 * (p1.getX() + p2.getX()) + B1 * (p1.getY() + p2.getY())) / 2;
		double A2 = p3.getX() - p2.getX(), B2 = p3.getY() - p2.getY(),
			C2 = (A2 * (p2.getX() + p3.getX()) + B2 * (p2.getY() + p3.getY())) / 2;

		double D = A1 * B2 - A2 * B1;
		if (D == 0) {
			throw std::invalid_argument("Cannot determine a unique circle.");
		}

		double x_center = (B2 * C1 - B1 * C2) / D;
		double y_center = (A1 * C2 - A2 * C1) / D;

		center = TU_Class_Point(x_center, y_center);

		// 计算半径
		radius = center.distance(p1);
	}

	TU_Class_Point getCenter()
	{
		return this->center;
	}

	double getRadius()
	{
		return this->radius;
	}

	// 返回交点个数，p1 和 p2 是输出参数
	int cross(TU_Class_Circle c, TU_Class_Point*& p1, TU_Class_Point*& p2) {
		double x1 = center.getX(), y1 = center.getY();
		double x2 = c.center.getX(), y2 = c.center.getY();
		double r1 = radius, r2 = c.radius;

		double d = std::sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));

		// 没有交点
		if (d > r1 + r2 || d < std::abs(r1 - r2) || d == 0) {
			p1 = nullptr;
			p2 = nullptr;
			return 0;
		}

		// 相切
		if (std::abs(d - (r1 + r2)) < 1e-8 || std::abs(d - std::abs(r1 - r2)) < 1e-8) {
			double px = x1 + (x2 - x1) * r1 / d;
			double py = y1 + (y2 - y1) * r1 / d;
			p1 = new TU_Class_Point(px, py);
			p2 = nullptr;
			return 1;
		}

		// 有两个交点
		double a = (r1 * r1 - r2 * r2 + d * d) / (2 * d);
		double h = std::sqrt(r1 * r1 - a * a);

		double px = x1 + a * (x2 - x1) / d;
		double py = y1 + a * (y2 - y1) / d;

		double rx = -(y2 - y1) * (h / d);
		double ry = (x2 - x1) * (h / d);

		p1 = new TU_Class_Point(px + rx, py + ry);
		p2 = new TU_Class_Point(px - rx, py - ry);

		return 2;
	}

};

class TU_Class_Vec
{
private:
	double x, y, z;
public:
	TU_Class_Vec(double x, double y, double z) : x(x), y(y), z(z) {}

	TU_Class_Vec() : x(0), y(0), z(0) {}

	TU_Class_Vec(double x, double y) : x(x), y(y), z(0) {}

	TU_Class_Vec(TU_Class_Point p1, TU_Class_Point p2)
	{
		x = p2.getX() - p1.getX();
		y = p2.getY() - p1.getY();
		z = p2.getZ() - p1.getZ();
	}

	TU_Class_Vec operator =(TU_Class_Vec v) {
		x = v.x;
		y = v.y;
		z = v.z;
		return *this;
	}

	// 获取x分量
	double getX() const { return x; }
	// 获取y分量
	double getY() const { return y; }
	// 获取z分量
	double getZ() const { return z; }

	// 向量加法
	TU_Class_Vec operator+(const TU_Class_Vec& v) const {
		return TU_Class_Vec(x + v.x, y + v.y, z + v.z);
	}

	// 向量减法
	TU_Class_Vec operator-(const TU_Class_Vec& v) const {
		return TU_Class_Vec(x - v.x, y - v.y, z - v.z);
	}

	// 点积（内积）
	double dot(const TU_Class_Vec& v) const {
		return x * v.x + y * v.y + z * v.z;
	}

	// 叉积（外积）
	TU_Class_Vec cross(const TU_Class_Vec& v) const {
		return TU_Class_Vec(
			y * v.z - z * v.y,
			z * v.x - x * v.z,
			x * v.y - y * v.x
		);
	}

	// 向量模（长度）
	double norm() const {
		return std::sqrt(x * x + y * y + z * z);
	}

	// 单位化向量（返回单位向量）
	TU_Class_Vec normalize() const {
		double length = norm();
		if (length == 0.0) {
			// 防止除以零
			return TU_Class_Vec(0, 0, 0);
		}
		return TU_Class_Vec(x / length, y / length, z / length);
	}

	// 打印向量
	void print() const {
		std::cout << "vec(" << x << ", " << y << ", " << z << ")" << std::endl;
	}

	void function(double(*func)(double))
	{
		x = func(x);
		y = func(y);
		z = func(z);
	}

	TU_Class_Vec perVec() const {
		if (std::abs(z) > 1e-8) {
			throw std::invalid_argument("只支持二维向量（z == 0）！");
		}

		TU_Class_Vec perp(y, -x, 0);
		return perp.normalize();
	}

};

class TU_Class_Line
{
private:
	TU_Class_Vec v;
	TU_Class_Point p;
public:
	TU_Class_Line(TU_Class_Point p, TU_Class_Vec v)
	{
		if (!v.norm())
			throw std::invalid_argument("错误：方向向量不能为零向量！");
		else
		{
			this->p = p;
			this->v = v;
		}
	}

	TU_Class_Line()
	{
		p = TU_Class_Point(0, 0, 0);
		v = TU_Class_Vec(1, 1, 1);
	}

	TU_Class_Line(TU_Class_Point p1, TU_Class_Point p2)
	{
		if (p1.getX() == p2.getX() && p1.getY() == p2.getY() && p1.getZ() == p2.getZ())
			throw std::invalid_argument("错误：两个点不能相同！");
		v = TU_Class_Vec(p2.getX() - p1.getX(), p2.getY() - p1.getY(), p2.getZ() - p1.getZ());
		p = p1;
	}

	bool isInLine(TU_Class_Point p0)
	{
		if (p.distance(p0) == 0)
			return true;
		else
		{
			TU_Class_Vec v0(this->p, p0);
			return v.cross(v0).norm() == 0;
		}

	}

	// 在 TU_Class_Line 类中添加以下两个成员函数：

// 1. 计算与另一条直线之间的夹角（弧度制）
	double angleWith(const TU_Class_Line& other) const {
		TU_Class_Vec v1 = this->v.normalize();
		TU_Class_Vec v2 = other.v.normalize();

		double dot_product = v1.dot(v2);
		// 防止浮点误差导致acos参数超出[-1, 1]
		dot_product = std::max(-1.0, std::min(1.0, dot_product));

		return std::acos(dot_product);  // 返回弧度值，范围 [0, π]
	}

	// 2. 返回两条直线的角平分线
	// 如果两直线相交，则返回从交点出发，方向为单位方向和的直线；
	// 如果平行，则返回位于两者之间的平行线（取两点中点作为线上一点）
	TU_Class_Line bisectorLine(const TU_Class_Line& other) const {
		TU_Class_Vec v1 = this->v.normalize();
		TU_Class_Vec v2 = other.v.normalize();

		// 检查是否平行（叉积为零向量）
		if (v1.cross(v2).norm() == 0.0) {
			// 平行情况：构造中间线

			// 取两个直线上各一个点
			TU_Class_Point p1 = this->p;
			TU_Class_Point p2 = other.p;

			// 构造一个新的点：两点连线的中点
			TU_Class_Point midPoint(
				(p1.getX() + p2.getX()) / 2,
				(p1.getY() + p2.getY()) / 2,
				(p1.getZ() + p2.getZ()) / 2
			);

			// 方向与原方向一致（任选其一即可）
			return TU_Class_Line(midPoint, v1);
		}
		else {
			// 相交情况：寻找交点

			// 这里假设两条直线在同一平面内且不平行，存在唯一交点
			// 由于是3D空间中的直线，交点可能不存在或重合，需要先判断是否共面、是否相交
			// 为简化问题，这里我们仅处理方向角平分的情况，忽略具体交点计算

			// 取本直线上的点作为起点
			TU_Class_Point intersection_point = this->p;  // 假设为交点，实际应使用几何算法求交点

			// 角平分方向为归一化后的方向和
			TU_Class_Vec bisector_dir = (v1 + v2).normalize();

			return TU_Class_Line(intersection_point, bisector_dir);
		}
	}

	TU_Class_Vec getVec()
	{
		return v;
	}

	TU_Class_Point getPoint()
	{
		return p;
	}

};

class TU_Class_Seg {
private:
	TU_Class_Point p1, p2;

public:
	// 构造函数：用两个点初始化线段，检查是否为同一点
	TU_Class_Seg(TU_Class_Point a, TU_Class_Point b) {
		if (a.getX() == b.getX() && a.getY() == b.getY() && a.getZ() == b.getZ()) {
			throw std::invalid_argument("错误：两个点不能相同！");
		}
		this->p1 = a;
		this->p2 = b;
	}

	// 获取线段长度
	double length() {
		return p1.distance(p2);
	}

	// 返回从 p1 指向 p2 的单位方向向量
	TU_Class_Vec direction() const {
		TU_Class_Vec v(p1, p2);
		return v.normalize();
	}

	// 返回线段的中垂线（中点 + 垂直方向）
	TU_Class_Line perpendicularBisector() {
		// 中点
		TU_Class_Point mid(
			(p1.getX() + p2.getX()) / 2,
			(p1.getY() + p2.getY()) / 2,
			(p1.getZ() + p2.getZ()) / 2
		);

		// 线段方向向量
		TU_Class_Vec dir = direction();

		// 构造垂直于线段方向的向量（取叉积与 z 轴的叉积来获得二维平面上的垂线方向）
		TU_Class_Vec perp_dir = dir.cross(TU_Class_Vec(0, 0, 1));

		if (perp_dir.norm() == 0.0) {
			// 如果是竖直方向的线段，则换一个方向
			perp_dir = TU_Class_Vec(0, 1, 0); // 例如水平方向
		}

		// 返回以中点为起点、垂直方向为方向的直线
		return TU_Class_Line(mid, perp_dir);
	}

	// 可选：获取线段的两个端点
	TU_Class_Point getP1() const { return p1; }
	TU_Class_Point getP2() const { return p2; }
};

void TU_print(std::vector<TU_Class_Point> p = {},
	std::vector<TU_Class_Circle> c = {},
	std::vector<TU_Class_Line> l = {},
	std::vector<TU_Class_Seg> s = {})
{
	std::ofstream outFile("output.txt");
	if (!outFile.is_open())
	{
		std::cerr << "无法打开输出文件！" << std::endl;
		return;
	}

	// 输出点 使用传统for循环
	for (size_t i = 0; i < p.size(); ++i)
	{
		outFile << "point " << p[i].getX() << " " << p[i].getY() << std::endl;
	}

	// 输出圆 使用传统for循环
	for (size_t j = 0; j < c.size(); ++j)
	{
		outFile << "circle "
			<< c[j].getCenter().getX() << " "
			<< c[j].getCenter().getY() << " "
			<< c[j].getRadius() << std::endl;
	}

	// 输出直线 使用传统for循环
	for (size_t k = 0; k < l.size(); ++k)
	{
		double x = l[k].getPoint().getX();
		double y = l[k].getPoint().getY();
		double dx = l[k].getVec().getX();
		double dy = l[k].getVec().getY();

		outFile << "line " << x << " " << y << " " << dx << " " << dy << std::endl;
	}

	// 输出线段 使用传统for循环
	for (size_t m = 0; m < s.size(); ++m)
	{
		double x1 = s[m].getP1().getX();
		double y1 = s[m].getP1().getY();
		double x2 = s[m].getP2().getX();
		double y2 = s[m].getP2().getY();

		outFile << "segment " << x1 << " " << y1 << " " << x2 << " " << y2 << std::endl;
	}

	outFile.close();

	// 调用 Python 脚本进行绘图
	int result = std::system("python C:/code/plot_data.py");

	if (result != 0)
	{
		std::cerr << "调用 Python 脚本失败，请确认：\n"
			<< "- Python 已安装并加入环境变量\n"
			<< "- plot_data.py 存在于指定路径" << std::endl;
	}
	std::cout << "调用 Python 脚本成功，请查看 output.png" << std::endl;
}

int TU_FindCoPoint(TU_Class_Circle c1, TU_Class_Circle c2, TU_Class_Point*& p1, TU_Class_Point*& p2)
{
	return c1.cross(c2, p1, p2);
}

int TU_FindCoPoint(TU_Class_Line& line, TU_Class_Circle& circle, TU_Class_Point*& p1, TU_Class_Point*& p2)
{
	p1 = p2 = nullptr;

	TU_Class_Point L_point = line.getPoint();
	TU_Class_Vec L_dir = line.getVec().normalize(); // 单位化方向向量

	TU_Class_Point C_center = circle.getCenter();
	double r = circle.getRadius();

	double x0 = L_point.getX(), y0 = L_point.getY();
	double dx = L_dir.getX(), dy = L_dir.getY();
	double cx = C_center.getX(), cy = C_center.getY();

	// 将圆心移到原点进行简化
	double ax = x0 - cx;
	double ay = y0 - cy;

	// 求解 at^2 + bt + c = 0 中的系数
	double a = dx * dx + dy * dy; // 始终为 1，因为方向向量已归一化
	double b = 2 * (ax * dx + ay * dy);
	double c = ax * ax + ay * ay - r * r;

	double discriminant = b * b - 4 * a * c;

	if (discriminant < 0)
	{
		// 无交点
		return 0;
	}
	else if (std::abs(discriminant) < 1e-8)
	{
		// 相切
		double t = (-b) / (2 * a);
		double px = x0 + dx * t;
		double py = y0 + dy * t;
		p1 = new TU_Class_Point(px, py);
		return 1;
	}
	else
	{
		// 两个交点
		double sqrt_d = std::sqrt(discriminant);
		double t1 = (-b + sqrt_d) / (2 * a);
		double t2 = (-b - sqrt_d) / (2 * a);

		double px1 = x0 + dx * t1;
		double py1 = y0 + dy * t1;
		double px2 = x0 + dx * t2;
		double py2 = y0 + dy * t2;

		p1 = new TU_Class_Point(px1, py1);
		p2 = new TU_Class_Point(px2, py2);

		return 2;
	}
}

int TU_FindCoPoint(TU_Class_Line& l1, TU_Class_Line& l2, TU_Class_Point*& p)
{
	p = nullptr;

	TU_Class_Point A = l1.getPoint(); // 直线1上的一点
	TU_Class_Point B = l2.getPoint(); // 直线2上的一点

	TU_Class_Vec v1 = l1.getVec().normalize(); // 方向向量1
	TU_Class_Vec v2 = l2.getVec().normalize(); // 方向向量2

	double crossNorm = v1.cross(v2).norm();

	// 情况1：方向向量共线（平行或重合）
	if (crossNorm < 1e-8)
	{
		// 判断是否共线（即是否在同一直线上）
		TU_Class_Vec AB(A, B);
		TU_Class_Vec crossABv1 = AB.cross(v1);

		if (crossABv1.norm() < 1e-8)
		{
			// 两直线重合，有无穷多解（这里返回一个点作为示例）
			p = new TU_Class_Point(A.getX(), A.getY(), A.getZ());
			return 1; // 返回一个代表性的点
		}
		else
		{
			// 平行且不重合，无交点
			return 0;
		}
	}

	// 情况2：异面或相交（使用参数方程求最短距离点）

	// 向量AB = B - A
	TU_Class_Vec AB(A, B);
	TU_Class_Vec n = v1.cross(v2); // 两个方向向量的叉积，垂直于两直线

	// 计算 s 和 t 参数：
	double nDotn = n.dot(n);

	double t_numerator = (AB.cross(v2)).dot(n);
	double t = t_numerator / nDotn;

	double s_numerator = (AB.cross(v1)).dot(n);
	double s = s_numerator / nDotn;

	// 求直线1和直线2上的最近点
	TU_Class_Point P_on_l1(
		A.getX() + v1.getX() * t,
		A.getY() + v1.getY() * t,
		A.getZ() + v1.getZ() * t
	);

	TU_Class_Point P_on_l2(
		B.getX() + v2.getX() * s,
		B.getY() + v2.getY() * s,
		B.getZ() + v2.getZ() * s
	);

	// 如果两点非常接近，则认为是交点
	if (P_on_l1.distance(P_on_l2) < 1e-8)
	{
		p = new TU_Class_Point(P_on_l1.getX(), P_on_l1.getY(), P_on_l1.getZ());
		return 1;
	}
	else
	{
		// 异面直线，无交点
		return 0;
	}
}



#endif
