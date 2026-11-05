#pragma once

#include <cstdlib>
#include <limits>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <vector>
#include <algorithm>
#include <array>
#include <new>
#include <type_traits>

namespace LAMP::Mesh { class NativeTargetLevelField; }

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
std::size_t QuadtreeMesh::completeRefinementPrerequisites(LevelCellFieldT& targetLevel) const
{
    if (targetLevel.size() != cellIndexNum_)
        throw std::invalid_argument("adaptation prerequisites require current target-cell storage");
    const int minimum = minimumAdaptationLevel();
    std::vector<int> aim(cellIndexNum_), originalAim(cellIndexNum_), queue;
    std::vector<unsigned char> queued(cellIndexNum_, 0);
    queue.reserve(cellIndexNum_);
    // Validate the entire field first. No partially raised result is exposed
    // when a later stale label or nonfinite target rejects.
    for (int i=0; i<cellIndexNum_; ++i)
    {
        const int label = cellIndexToLabel[i];
        if (label < 0 || label >= maxCellNum_ || cells_[label].label() != label
            || cells_[label].index() != i || targetLevel.label(i) != label
            || !std::isfinite(targetLevel.value[i])
            || targetLevel.value[i] < std::numeric_limits<int>::min()
            || targetLevel.value[i] > std::numeric_limits<int>::max())
            throw std::invalid_argument("invalid or stale adaptation prerequisite target");
        const auto& cell = cells_[label];
        if (cell.level() < 0 || cell.level() > maxLevel_)
            throw std::invalid_argument("invalid current level in adaptation prerequisites");
        originalAim[i] = aim[i] = boundedTargetLevel(targetLevel.value[i], minimum);
        if (cell.level() < aim[i]) { queue.push_back(i); queued[i] = 1; }
    }
    for (std::size_t head=0; head<queue.size(); ++head)
    {
        const auto& cell = cells_[cellIndexToLabel[queue[head]]];
        // Each actual leaf needs this examination at most once: the condition
        // depends on its CURRENT level, not on how far above it the desire is.
        for (int side=0; side<4; ++side) for (int half=0; half<2; ++half)
        {
            const int label = cell.neiCell[side][half];
            if (label == -1) continue;
            if (label < 0 || label >= maxCellNum_ || cells_[label].label() != label)
                throw std::invalid_argument("invalid actual face neighbor in adaptation prerequisites");
            const auto& neighbor = cells_[label];
            const int i = neighbor.index();
            if (i < 0 || i >= cellIndexNum_ || cellIndexToLabel[i] != label
                || targetLevel.label(i) != label)
                throw std::invalid_argument("nonleaf or stale face neighbor in adaptation prerequisites");
            if (neighbor.level() >= cell.level()) continue;
            const int required = neighbor.level()+1;
            if (aim[i] < required) aim[i] = required;
            if (!queued[i]) { queue.push_back(i); queued[i] = 1; }
        }
    }
    std::size_t raised = 0;
    for (int i=0; i<cellIndexNum_; ++i) if (aim[i] > originalAim[i])
    {
        targetLevel.value[i] = static_cast<double>(aim[i]);
        ++raised;
    }
    return raised;
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

    // A local receipt is safe only for the built-in nonvirtual target-field
    // updates. Opaque/generic callbacks retain every original capacity scan.
    // It never survives this adapt(), and its first use scans ALL pool slots.
    struct LabelSpaceReceipt {std::vector<int> free;int ptr=0;bool valid=false;};
    std::unique_ptr<std::array<LabelSpaceReceipt,4>> labelSpaces;
    bool useLabelSpaces=std::is_same_v<CellFieldT,LAMP::Mesh::NativeTargetLevelField>
        &&std::is_same_v<LevelCellFieldT,LAMP::Mesh::NativeTargetLevelField>
        &&targetLevel!=nullptr&&fields.empty()&&vertFields.empty()
        &&hEdgeFields.empty()&&vEdgeFields.empty();
    auto nativeWorkspace=beginNativeAdaptationMaintenance(useLabelSpaces&&levelSetCount_==0);
    const auto originalLabelSpaces=[&] {
        updateLabelSpace(verts_, vertLabelSpace_, vertPtr_, maxVertNum_);
        updateLabelSpace(hEdges_, hEdgeLabelSpace_, hEdgePtr_, maxHEdgeNum_);
        updateLabelSpace(vEdges_, vEdgeLabelSpace_, vEdgePtr_, maxVEdgeNum_);
        updateLabelSpace(cells_, cellLabelSpace_, cellPtr_, maxCellNum_);
    };
    const auto captureLabelSpace=[](const auto& list,int*& space,int& ptr,int maximum,LabelSpaceReceipt& receipt) {
        receipt.valid=false;receipt.free.clear();receipt.free.reserve(static_cast<std::size_t>(maximum));
        int nextPtr=ptr;
        for(int i=0;i<maximum;++i) {
            if(list.peek(i).label()!=-1)nextPtr=i+1;
            else receipt.free.push_back(i);
        }
        // The entire original ascending free prefix, never its unused tail.
        std::copy(receipt.free.begin(),receipt.free.end(),space);
        ptr=nextPtr;receipt.ptr=nextPtr;receipt.valid=true;
    };
    const auto restoreLabelSpace=[](int*& space,int& ptr,const LabelSpaceReceipt& receipt) {
        std::copy(receipt.free.begin(),receipt.free.end(),space);ptr=receipt.ptr;
    };
    const auto updateLabelSpaces=[&](bool noPoolWrites) {
        if(updateNativeLabelSpaces(nativeWorkspace.get()))return;
        if(!useLabelSpaces) {originalLabelSpaces();return;}
        if(noPoolWrites&&labelSpaces&&std::all_of(labelSpaces->begin(),labelSpaces->end(),
            [](const auto& receipt) {return receipt.valid;})) {
            // subdivideCell's no-work path still overwrites part of the free
            // prefix. Restore it; simply skipping maintenance is incorrect.
            restoreLabelSpace(vertLabelSpace_,vertPtr_,(*labelSpaces)[0]);
            restoreLabelSpace(hEdgeLabelSpace_,hEdgePtr_,(*labelSpaces)[1]);
            restoreLabelSpace(vEdgeLabelSpace_,vEdgePtr_,(*labelSpaces)[2]);
            restoreLabelSpace(cellLabelSpace_,cellPtr_,(*labelSpaces)[3]);return;
        }
        bool captured=false;
        try {
            if(!labelSpaces)labelSpaces=std::make_unique<std::array<LabelSpaceReceipt,4>>();
            captureLabelSpace(verts_,vertLabelSpace_,vertPtr_,maxVertNum_,(*labelSpaces)[0]);
            captureLabelSpace(hEdges_,hEdgeLabelSpace_,hEdgePtr_,maxHEdgeNum_,(*labelSpaces)[1]);
            captureLabelSpace(vEdges_,vEdgeLabelSpace_,vEdgePtr_,maxVEdgeNum_,(*labelSpaces)[2]);
            captureLabelSpace(cells_,cellLabelSpace_,cellPtr_,maxCellNum_,(*labelSpaces)[3]);
            captured=true;
        }catch(const std::bad_alloc&) {
        }catch(const std::length_error&) {
        }
        if(captured)return;
        // Extra receipt storage is optional execution state. On resource
        // failure, release it and complete the original allocation-free scan.
        useLabelSpaces=false;
        labelSpaces.reset();
        originalLabelSpaces();
    };
    for (int currentLevel = 0; currentLevel < nCellsPerLevel.size() - 1; currentLevel++)
    {
        prepareTransferIndices(nativeWorkspace.get());
        refineOrCoarsen0(fields, magGradU, currentLevel, targetLevel);

        // This is an exact branch-body proof, not a leaf-count/generation
        // guess: with no needRefine, no subdivision mutation body is entered.
        bool noPoolWrites=useLabelSpaces;
        if(noPoolWrites)for(int i=0;i<cellIndexNum_;++i)
            if(cells_[cellIndexToLabel[i]].needRefine) {noPoolWrites=false;break;}
        subdivideCellWithWorkspace(currentLevel,nativeWorkspace.get());

        updateLabelSpaces(noPoolWrites);

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
        prepareTransferIndices(nativeWorkspace.get());
        refineOrCoarsen0(fields, magGradU, currentLevel, targetLevel);

        // Keep even rejected coarsening-family work conservative: any
        // needCoarsen triggers the full original pool-occupancy inspection.
        bool noPoolWrites=useLabelSpaces;
        if(noPoolWrites)for(int i=0;i<cellIndexNum_;++i)
            if(cells_[cellIndexToLabel[i]].needCoarsen) {noPoolWrites=false;break;}
        submergeCellWithWorkspace(currentLevel,nativeWorkspace.get());

        updateLabelSpaces(noPoolWrites);

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
void QuadtreeMesh::updateLabelSpace(LAMP::Mesh::Detail::EntityPool<T>& list, int*& space, int& ptr, const int maxNum)
{
    int _count = 0;
    int _remainCount = 0;
    // Preserve the complete ascending free-label prefix and last live ptr
    // in one capacity scan. Leave the unused suffix of space untouched.
    for (int i = 0; i < maxNum; i++)
    {
        const T& t = list.peek(i);
        if (t.label() != -1)
        {
            _count++;
            ptr = i + 1;
        }
        else
        {
            space[_remainCount] = i;
            _remainCount++;
        }
    }

    int _remainNum = maxNum - _count;
}

template<typename T>
void QuadtreeMesh::count(LAMP::Mesh::Detail::EntityPool<T>& list, const int ptr, int& indexNum, int*& indexToLabel)
{
    indexNum = 0;
    for (int i = 0; i < ptr; i++)
    {
        const T& t = list.peek(i);
        if (t.label() != -1 && t.son[0] == -1)
        {
            indexNum++;
        }
    }

    int* tempIndexToLabel = new int[indexNum];
    int count = 0;
    for (int i = 0; i < ptr; i++)
    {
        const T& t = list.peek(i);
        if (t.label() != -1 && t.son[0] == -1)
        {
            tempIndexToLabel[count] = i;
            list[i].setIndex(count);
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
