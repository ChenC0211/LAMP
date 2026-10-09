#pragma once

#include <vector>
#include "Core/Base/Base.h"
#include "Core/Base/Vector2D.h"
#include "Core/Base/Enumeration.h"

class Cell
{
public:
	std::vector<unsigned char> isCut;
	
	int cellToVert[4];// 单元包含的顶点

	int cellToHEdge[2][2];// 单元包含的边
	int cellToVEdge[2][2];// 单元包含的边
	int neiCell[4][2];// 相邻单元

	std::vector<int> relatedPreIndex;// 相关的前一时刻的单元编号
	std::vector<Vector2D> relatedPreCoord;// 相关的前一时刻的单元坐标

	int son[4];// 子网格 0——x-y+, 1——x+y+, 2——x+y-, 3——x-y-
	int parent = -1;// 父网格

	int supCellIndex = -1;

	int coarsenFlag = 0;
	int refineFlag = 0;
	bool needCoarsen = false;// 粗化标志
	bool needRefine = false;// 细化标志

	Cell();
	Cell(int label, Vector2D centerCoord, int level, Vector2D size);
	~Cell();

	void ensureLevelSetSize(std::size_t n)
	{
		if (isCut.size() != n)
			isCut.assign(n, 0u);
	}

	inline int label() const { return label_; }
	inline int index() const { return index_; }
	inline int preIndex() const { return preIndex_; }
	inline Vector2D centerCoord() const { return centerCoord_; }
	inline int level() const { return level_; }
	inline Vector2D size() const { return size_; }
	inline double area() const { return size_.x * size_.y; }
	inline bool needUpdate() const { return needUpdate_; }
	inline Boundary boundary() const { return boundary_; }

	inline void setLabel(int label) { label_ = label; }
	inline void setIndex(int index) { index_ = index; }
	inline void setPreIndex(int preIndex) { preIndex_ = preIndex; }
	inline void setCenterCoord(Vector2D centerCoord) { centerCoord_ = centerCoord; }
	inline void setLevel(int level) { level_ = level; }
	inline void setSize(Vector2D size) { size_ = size; }
	inline void setNeedUpdate(bool needUpdate) { needUpdate_ = needUpdate; }
	inline void setBoundary(Boundary boundary) { boundary_ = boundary; }

private:
	int label_ = -1;// 单元标签
	int index_ = -1;// 单元编号
	int preIndex_ = -1;
	Vector2D centerCoord_;// 中心坐标
	int level_ = 0;// 网格等级
	Vector2D size_;// 网格尺寸

	bool needUpdate_ = false;// 是否需要更新

	Boundary boundary_ = Boundary::non;// 边界标签
};
