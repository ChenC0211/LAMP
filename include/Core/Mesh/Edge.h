#pragma once

#include <vector>
#include "Core/Base/Base.h"
#include "Core/Base/Vector2D.h"
#include "Core/Base/Enumeration.h"

class Edge
{
public:
	std::vector<unsigned char> isCut;

	int edgeToVert[2] = {-1, -1};// 0——起点，1——终点
	int masterCell = -1;// 主单元
    int slaveCell = -1;// 从单元

	int parent = -1;// 父边
	int son[2] = {-1, -1};// 子边

	std::vector<int> relatedPreIndex;// 相关的前一时刻的单元编号
	std::vector<Vector2D> relatedPreCoord;// 相关的前一时刻的单元坐标

	Edge();
	Edge(int label, Vector2D midpointCoord, double length, bool isHorizontal);
	~Edge();

	void ensureLevelSetSize(std::size_t n)
	{
		if (isCut.size() != n)
			isCut.assign(n, 0u);
	}

	inline int label() const { return label_; }
	inline int index() const { return index_; }
	inline Vector2D midpointCoord() const { return midpointCoord_; }
	inline double length() const { return length_; }
	inline bool isHorizontal() const { return isHorizontal_; }
	inline int preIndex() const { return preIndex_; }
	inline Boundary boundary() const { return boundary_; }
	inline bool needUpdate() const { return needUpdate_; }
	inline bool cut() const { return cut_; }
	inline double cutLength() const { return cutLength_; }

	inline void setLabel(int label) { label_ = label; }
	inline void setIndex(int index) { index_ = index; }
	inline void setMidpointCoord(Vector2D midpointCoord) { midpointCoord_ = midpointCoord; }
	inline void setLength(double length) { length_ = length; }
	inline void setIsHorizontal(bool isHorizontal) { isHorizontal_ = isHorizontal; }
	inline void setPreIndex(int preIndex) { preIndex_ = preIndex; }
	inline void setBoundary(Boundary boundary) { boundary_ = boundary; }
	inline void setNeedUpdate(bool needUpdate) { needUpdate_ = needUpdate; }
	inline void setCut(bool cut) { cut_ = cut; }
	inline void setCutLength(double cutLength) { cutLength_ = cutLength; }

protected:
	int label_ = -1;// 边标签
	int index_ = -1;// 边编号
	Vector2D midpointCoord_;// 中点坐标
	double length_ = 0.0;// 边长
	bool isHorizontal_ = true;// 是否水平边
	int preIndex_ = -1;

	Boundary boundary_ = Boundary::non;// 边界标签

	bool needUpdate_ = false;// 是否需要更新

	bool cut_ = false;// 是否被界面截断
	double cutLength_ = 0.0;// 截断长度
};
