#include "Core/Mesh/Cell.h"

using namespace std;

Cell::Cell() 
{
	for (int i = 0; i < 4; i++)
	{
		cellToVert[i] = -1;
		son[i] = -1;

		for (int j = 0; j < 2; j++)
		{
			neiCell[i][j] = -1;
		}
	}

	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			cellToHEdge[i][j] = -1;
			cellToVEdge[i][j] = -1;
		}
	}
}

Cell::Cell
(
	int label,
	Vector2D centerCoord, 
	int level, 
	Vector2D size
) :
	label_(label),
	centerCoord_(centerCoord), 
	level_(level), 
	size_(size)
{
	for (int i = 0; i < 4; i++)
	{
		cellToVert[i] = -1;
		son[i] = -1;

		for (int j = 0; j < 2; j++)
		{
			neiCell[i][j] = -1;
		}
	}

	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			cellToHEdge[i][j] = -1;
			cellToVEdge[i][j] = -1;
		}
	}
}

Cell::~Cell() {}