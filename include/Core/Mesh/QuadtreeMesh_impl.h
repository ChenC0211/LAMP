#pragma once

#include <cstdlib>
#include <limits>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <vector>

#include "Core/IO/ConsoleControl.h"

template<typename Eligible>
bool QuadtreeMesh::coarseningFamilyEligible(const Cell& parent, const Eligible& eligible) const
{
    bool result = true;
    // Do not short circuit: the original submerge callback also consumes all
    // four flags, including siblings of a rejected coarsening request.
    for (int i=0; i<4; ++i)
        if (!eligible(parent.son[i])) result = false;
    return result;
}

template<typename LevelCellFieldT>
bool QuadtreeMesh::mayAdaptToTargets(const LevelCellFieldT& targetLevel) const
{
    if (targetLevel.size() != cellIndexNum_)
        throw std::invalid_argument("adaptation preflight requires current target-cell storage");
    for (int i=0; i<cellIndexNum_; ++i)
    {
        const int label = cellIndexToLabel[i];
        if (label < 0 || label >= maxCellNum_ || cells_[label].label() != label
            || targetLevel.label(i) != label || !std::isfinite(targetLevel.value[i])
            || targetLevel.value[i] < std::numeric_limits<int>::min()
            || targetLevel.value[i] > std::numeric_limits<int>::max())
            throw std::invalid_argument("invalid or stale adaptation target");
    }

    const int minimum = minimumAdaptationLevel();
    std::vector<unsigned char> coarsen(maxCellNum_, 0), leaf(maxCellNum_, 0);
    for (int i=0; i<cellIndexNum_; ++i)
    {
        const auto& cell = cells_[cellIndexToLabel[i]];
        // Pending legacy votes are outside this target-only shortcut. Let the
        // unchanged adapt() consume them rather than incorrectly skipping work.
        if (cell.refineFlag != 0 || cell.coarsenFlag != 0) return true;
        const int aim = boundedTargetLevel(targetLevel.value[i], minimum);
        if (cell.level() < aim && cell.level()+1 < static_cast<int>(nCellsPerLevel.size())
            && targetRefinementAllowed(cell)) return true;
        leaf[cell.label()] = 1;
        coarsen[cell.label()] = cell.level() > 0 && cell.level() > aim && targetCoarseningAllowed(cell);
    }
    for (int i=0; i<cellIndexNum_; ++i)
    {
        const auto& cell = cells_[cellIndexToLabel[i]];
        if (!coarsen[cell.label()]) continue;
        if (cell.parent < 0 || cell.parent >= maxCellNum_) return true;
        const auto& parent = cells_[cell.parent];
        for (int child:parent.son)
            if (child < 0 || child >= maxCellNum_) return true;
        if (coarseningFamilyEligible(parent, [&](int child) {
                // Non-leaf flags are not cleared by the original leaf marker.
                // Preserve that legacy behavior conservatively in this query.
                return leaf[child] ? coarsen[child] && cells_[child].level() == cell.level()
                    : cells_[child].needCoarsen;
            })) return true;
    }
    // No first change at any level: subsequent sweeps cannot acquire a new
    // eligible neighbor/family without a preceding topology change.
    return false;
}

template<typename Field, typename VertField, typename HEdgeField, typename VEdgeField,
         typename CellFieldT, typename LevelCellFieldT>
void QuadtreeMesh::adapt(std::vector<Field>& fields, CellFieldT& magGradU,
                         std::vector<VertField*>& vertFields, std::vector<HEdgeField*>& hEdgeFields,
                         std::vector<VEdgeField*>& vEdgeFields, std::vector<CellFieldT*>& cellFields,
                         const LevelCellFieldT* targetLevel)
{
    if (ConsoleControl::enabled(ConsoleControl::Level::kDebug))
    {
        std::cout << "--------------------------------\n";
        std::cout << "[Mesh] adapt begin\n";
    }

    for (int currentLevel = 0; currentLevel < nCellsPerLevel.size() - 1; currentLevel++)
    {
        prepareTransferIndices();
        refineOrCoarsen0(fields, magGradU, currentLevel, targetLevel);

        subdivideCell(currentLevel);

        updateLabelSpace(verts_, vertLabelSpace_, vertPtr_, maxVertNum_);
        updateLabelSpace(hEdges_, hEdgeLabelSpace_, hEdgePtr_, maxHEdgeNum_);
        updateLabelSpace(vEdges_, vEdgeLabelSpace_, vEdgePtr_, maxVEdgeNum_);
        updateLabelSpace(cells_, cellLabelSpace_, cellPtr_, maxCellNum_);

        for (Field& f : fields)
        {
            f.updateWithMesh(*this);
        }
        magGradU.updateWithMesh(*this);
        for (VertField* f : vertFields)
        {
            f->updateWithMesh(*this);
        }
        for (HEdgeField* f : hEdgeFields)
        {
            f->updateWithMesh(*this);
        }
        for (VEdgeField* f : vEdgeFields)
        {
            f->updateWithMesh(*this);
        }
        for (CellFieldT* f : cellFields)
        {
            f->updateWithMesh(*this);
        }
    }
    // std::cout << "Mesh refinement done." << std::endl;

    for (int currentLevel = static_cast<int>(nCellsPerLevel.size()) - 1; currentLevel > 0; currentLevel--)
    {
        prepareTransferIndices();
        refineOrCoarsen0(fields, magGradU, currentLevel, targetLevel);

        submergeCell(currentLevel);

        updateLabelSpace(verts_, vertLabelSpace_, vertPtr_, maxVertNum_);
        updateLabelSpace(hEdges_, hEdgeLabelSpace_, hEdgePtr_, maxHEdgeNum_);
        updateLabelSpace(vEdges_, vEdgeLabelSpace_, vEdgePtr_, maxVEdgeNum_);
        updateLabelSpace(cells_, cellLabelSpace_, cellPtr_, maxCellNum_);

        for (Field& f : fields)
        {
            f.updateWithMesh(*this);
        }
        magGradU.updateWithMesh(*this);
        for (VertField* f : vertFields)
        {
            f->updateWithMesh(*this);
        }
        for (HEdgeField* f : hEdgeFields)
        {
            f->updateWithMesh(*this);
        }
        for (VEdgeField* f : vEdgeFields)
        {
            f->updateWithMesh(*this);
        }
        for (CellFieldT* f : cellFields)
        {
            f->updateWithMesh(*this);
        }
    }
    // std::cout << "Mesh coarsening done." << std::endl;

    if (ConsoleControl::enabled(ConsoleControl::Level::kDebug))
    {
        std::cout << "--------------------------------\n";
        std::cout << "level  cellNum\n";
        for (int currentLevel = 0; currentLevel < nCellsPerLevel.size(); currentLevel++)
        {
            std::cout << currentLevel << "      "
                << nCellsPerLevel[currentLevel] << "\n";
        }
        std::cout << "--------------------------------\n";
    }
    if (ConsoleControl::enabled(ConsoleControl::Level::kDebug))
    {
        std::cout << "[Mesh] adapt done: cells=" << cellIndexNum_
                  << " edges(H/V)=" << hEdgeIndexNum_ << "/" << vEdgeIndexNum_
                  << " verts=" << vertIndexNum_ << "\n";
        std::cout << "--------------------------------\n";
    }

    defineTypes();
}

template<typename VertField, typename CellFieldT, typename LevelCellFieldT>
void QuadtreeMesh::refineOrCoarsen0(std::vector<VertField>& LSs, CellFieldT& magGradU, const int currentLevel,
                                    const LevelCellFieldT* targetLevel)
{
    [[maybe_unused]] const bool useOmp = (cellIndexNum_ > 4096);
    int fieldsInCurrentLevel = 0;
    for (int levelSetI = 0; levelSetI < LSs.size(); levelSetI++)
    {
        VertField& field = LSs[levelSetI];
        if (currentLevel <= maxLevel_local[levelSetI])
        {
            fieldsInCurrentLevel++;
            // int nRefine = 0;
            // int nCoarsen = 0;
            #ifdef _OPENMP
            #pragma omp parallel for if(useOmp)
            #endif
            for (int c = 0; c < cellIndexNum_; c++)
            {
                Cell& cell = cells_[cellIndexToLabel[c]];
                // cell.needRefine = false;
                // cell.needCoarsen = false;
                if (cell.level() == currentLevel)
                {
                    double value[4];

                    for (int i = 0; i < 4; i++)
                    {
                        const Vert& vert = verts_[cell.cellToVert[i]];
                        value[i] = field.value[vert.index()];
                    }

                    double minValue 
                        = min
                        (
                            min(fabs(value[0]), fabs(value[1])), 
                            min(fabs(value[2]), fabs(value[3]))
                        );
                
                    int sgn = 0;
                    if (value[0] > 0.0 && value[1] > 0.0 && value[2] > 0.0 && value[3] > 0.0)
                    {
                        sgn = 1;
                    }
                    else if (value[0] < 0.0 && value[1] < 0.0 && value[2] < 0.0 && value[3] < 0.0)
                    {
                        sgn = -1;
                    }

                    if 
                    (
                        (sgn <= 0 && cell.level() < minLevel_neg[levelSetI])
                        || (sgn > 0 && cell.level() < minLevel_pos[levelSetI])
                        || (minValue < C_Lipschitz * cell.size().mod())
                    )
                    {
                        cell.refineFlag++;
                        // cell.needRefine = true;
                        // nRefine++;
                        for (int i = 0; i < 4; i++)// 避免相邻单元超过两个
                        {
                            if (cell.neiCell[i][0] != -1)
                            {
                                if (cells_[cell.neiCell[i][0]].level() < currentLevel)
                                {
                                    cell.refineFlag--;
                                    // cell.needRefine = false;
                                    // nRefine--;
                                    break;
                                }
                            }
                        }
                    }
                    else if 
                    (
                        (sgn <= 0 && cell.level() > minLevel_neg[levelSetI])
                        || (sgn > 0 && cell.level() > minLevel_pos[levelSetI])
                    )
                    {
                        cell.coarsenFlag++;
                        // cell.needCoarsen = true;
                        // nCoarsen++;
                        for (int i = 0; i < 4; i++)// 避免相邻单元超过两个
                        {
                            if (cell.neiCell[i][0] != -1)
                            {
                                if (cells_[cell.neiCell[i][0]].level() > currentLevel)
                                {
                                    cell.coarsenFlag--;
                                    // cell.needCoarsen = false;
                                    // nCoarsen--;
                                    break;
                                }
                            }
                        } 
                    }
                }
            }
            // cout << "level: " << currentLevel << ", refine: " << nRefine << ", coarsen: " << nCoarsen << endl;
        }
    }

    const int _minLevel = minimumAdaptationLevel();
    double _dLevel_1 = static_cast<double>(maxLevel_ - _minLevel) + 1.0;

    double _magGradUMax = 0.0;
    double _magGradUMin = MAXVAL;
    double _magGradU_point1 = 0.0;
    double _magGradU_point2 = 0.0;
    if (targetLevel == nullptr)
    {
        for (int cellI = 0; cellI < cellIndexNum_; cellI++)
        {
            Cell& cell = cells_[cellIndexToLabel[cellI]];
            _magGradUMax = max(_magGradUMax, magGradU.value[cellI]);
            _magGradUMin = min(_magGradUMin, magGradU.value[cellI]);
        }
        _magGradU_point1 = 0.75 * (_magGradUMax - _magGradUMin) + _magGradUMin;
        _magGradU_point2 = 0.25 * (_magGradUMax - _magGradUMin) + _magGradUMin;
    }
    
    #ifdef _OPENMP
    #pragma omp parallel for if(useOmp)
    #endif
    for (int cellI = 0; cellI < cellIndexNum_; cellI++)
    {
        Cell& cell = cells_[cellIndexToLabel[cellI]];

        if (cell.level() != currentLevel)
        {
            continue;
        }

        int _aimLevel = _minLevel;
        if (targetLevel != nullptr)
        {
            double v = 0.0;
            if (cellI >= 0 && cellI < targetLevel->size())
                v = targetLevel->value[cellI];
            _aimLevel = boundedTargetLevel(v, _minLevel);
        }
        else
        {
            if (magGradU.value[cellI] > _magGradU_point1)
            {
                _aimLevel = maxLevel_;
            }
            else if (magGradU.value[cellI] < _magGradU_point2)
            {
                _aimLevel = _minLevel;
            }
            else
            {
                _aimLevel = static_cast<int>(std::floor((magGradU.value[cellI] - _magGradU_point2) / (_magGradU_point1 - _magGradU_point2) * _dLevel_1)) + _minLevel;
            }
        }

        if (cell.level() < _aimLevel)
        {
            cell.refineFlag++;

            if (!targetRefinementAllowed(cell)) cell.refineFlag--;
        }
        if (cell.level() > _aimLevel)
        {
            cell.coarsenFlag++;
            if (!targetCoarseningAllowed(cell)) cell.coarsenFlag--;
        }
    }

    int nRefine = 0;
    int nCoarsen = 0;
    #ifdef _OPENMP
    #pragma omp parallel for if(useOmp) reduction(+:nRefine, nCoarsen)
    #endif
    for (int c = 0; c < cellIndexNum_; c++)
    {
        Cell& cell = cells_[cellIndexToLabel[c]];
        cell.needRefine = false;
        cell.needCoarsen = false;
        if (cell.level() == currentLevel)
        {
            if (cell.refineFlag > 0)
            {
                cell.needRefine = true;
                nRefine++;
            }

            if (cell.coarsenFlag == fieldsInCurrentLevel + 1)
            {
                cell.needCoarsen = true;
                nCoarsen++;
            }

            cell.refineFlag = 0;
            cell.coarsenFlag = 0;
        }
    }
    // cout << "level: " << currentLevel << ", refine: " << nRefine << ", coarsen: " << nCoarsen << endl;
}

template<typename T>
void QuadtreeMesh::updateLabelSpace(T*& list, int*& space, int& ptr, const int maxNum)
{
    int _count = 0;
    for (int i = 0; i < maxNum; i++)
    {
        T& t = list[i];
        if (t.label() != -1)
        {
            _count++;
            ptr = i + 1;
        }
    }

    int _remainNum = maxNum - _count;

    int _remainCount = 0;
    for (int i = 0; i < maxNum; i++)
    {
        T& t = list[i];
        if (t.label() == -1)
        {
            space[_remainCount] = i;
            _remainCount++;
        }
    }
}

template<typename T>
void QuadtreeMesh::count(T*& list, const int ptr, int& indexNum, int*& indexToLabel)
{
    indexNum = 0;
    for (int i = 0; i < ptr; i++)
    {
        T& t = list[i];
        if (t.label() != -1 && t.son[0] == -1)
        {
            indexNum++;
        }
    }

    int* tempIndexToLabel = new int[indexNum];
    int count = 0;
    for (int i = 0; i < ptr; i++)
    {
        T& t = list[i];
        if (t.label() != -1 && t.son[0] == -1)
        {
            tempIndexToLabel[count] = i;
            t.setIndex(count);
            count++;
        }
    }

    delete[] indexToLabel;
    indexToLabel = tempIndexToLabel;
}

template<typename VertField>
void QuadtreeMesh::outputFields(const std::string& name, const std::vector<std::string>& valNames, const std::vector<VertField>& fields) const
{
    std::cout << std::endl << "Writing fields..." << std::endl;

    std::ostringstream fileName;
    fileName << name;
    std::string filePath = "./output/" + fileName.str() + ".dat";

    std::ofstream fp(filePath, std::ofstream::out);
	if (!fp.is_open())
	{
        std::cerr << filePath << " write mesh file failed!" << std::endl;
        std::exit(-1);
	}

    fp << "Title= \"Quadtree Mesh" << "\"" << std::endl;
    fp << "Variables= \"X\" \"Y\" ";
    for (int i = 0; i < valNames.size(); i++)
    {
        fp << "\"" << valNames[i] << "\" ";
    }
    fp << std::endl;

	fp << "ZONE N = " << vertIndexNum_ << ", E = " << cellIndexNum_
        << ", F = FEPOINT, ET = QUADRILATERAL" << std::endl;

	for (int i = 0; i < vertIndexNum_; i++)
    {
        const int vertLabel = vertIndexToLabel[i];
        const Vert& vert = verts_[vertLabel];
        fp << verts_[vertLabel].coord().x << " " << verts_[vertLabel].coord().y << " ";
        for (int j = 0; j < fields.size(); j++)
        {
            fp << fields[j].value[i] << " ";
        }
        fp << std::endl;
    }

	for (int i = 0; i < cellIndexNum_; i++)
    {
        const int cellLabel = cellIndexToLabel[i];
        const Cell& cell = cells_[cellLabel];
        for (int j = 0; j < 4; j++)
        {
            const Vert& vert = verts_[cell.cellToVert[j]];
            fp << vert.index() + 1 << " ";
        }
        fp << std::endl;
    }

    fp.close();
    std::cout << "Field output done." << std::endl << std::endl;
}
