#pragma once

#include "Core/Field/CellField.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <iterator>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace LAMP::Mesh
{
// Discrete native adaptation targets, NOT a physical/interpolated CellField.
// Core adapt invokes updateWithMesh after each ONE-level topology operation.
// Refine: actual old parent -> all four actual children. Coarsen: max of all
// four actual old leaf siblings. LS-neighbor transfer stencils are never used.
class NativeTargetLevelField : public CellField<double>
{
    struct Leaf {int label=-1,parent=-1,level=0,slot=-1;};
    std::vector<Leaf> previous_;
    int allocatedLevel_=0;
    static void require(bool ok,const char* why)
    {if(!ok)throw std::invalid_argument(std::string("[NativeTargetLevelField] ")+why);}
    static std::vector<Leaf> leaves(const QuadtreeMesh& mesh)
    {
        std::vector<Leaf> result;result.reserve(mesh.cellIndexNum());
        std::unordered_map<int,int> unique;
        for(int i=0;i<mesh.cellIndexNum();++i){
            const int label=mesh.cellIndexToLabel[i];
            require(label>=0&&label<mesh.cellLabelCapacity(),"invalid current native label");
            const auto& cell=mesh.cells(label);
            require(cell.label()==label&&cell.index()==i&&cell.level()>=0&&cell.level()<=mesh.maxLevel(),"invalid current leaf identity/level");
            require(std::all_of(std::begin(cell.son),std::end(cell.son),[](int q){return q==-1;}),"native target storage contains a nonleaf");
            require(unique.emplace(label,i).second,"duplicate current leaf label");
            int slot=-1;
            if(cell.parent>=0){
                require(cell.parent<mesh.cellLabelCapacity(),"invalid actual parent label");
                const auto& parent=mesh.cells(cell.parent);
                require(parent.label()==cell.parent&&parent.level()+1==cell.level(),"actual parent is not one level above child");
                for(int k=0;k<4;++k)if(parent.son[k]==label){require(slot<0,"duplicate actual parent child slot");slot=k;}
                require(slot>=0,"actual parent does not own current leaf");
            }else require(cell.level()==0,"non-root leaf has no actual parent");
            result.push_back({label,cell.parent,cell.level(),slot});
        }
        return result;
    }
    void validTarget(double target) const
    {require(std::isfinite(target)&&target>=0&&target<=allocatedLevel_&&std::floor(target)==target,"requires a finite allocated integer target");}
public:
    explicit NativeTargetLevelField(const QuadtreeMesh& mesh,double initial=0.)
        :CellField<double>(mesh,initial),previous_(leaves(mesh)),allocatedLevel_(mesh.maxLevel())
    {validTarget(initial);}
    void updateWithMesh(const QuadtreeMesh& mesh)
    {
        require(mesh.maxLevel()==allocatedLevel_&&this->size()==static_cast<int>(previous_.size()),"different allocation or stale prior storage");
        const auto current=leaves(mesh);
        std::unordered_map<int,int> oldByLabel;
        std::unordered_map<int,std::array<int,4>> oldByParent;
        for(int i=0;i<this->size();++i){const auto& old=previous_[i];
            require(this->label(i)==old.label,"stale previous target labels");validTarget(this->value[i]);
            require(oldByLabel.emplace(old.label,i).second,"duplicate previous leaf label");
            if(old.parent>=0){
                auto inserted=oldByParent.emplace(old.parent,std::array<int,4>{{-1,-1,-1,-1}});
                auto& group=inserted.first->second;
                require(old.slot>=0&&old.slot<4&&group[old.slot]<0,"duplicate/missing previous child slot");group[old.slot]=i;
            }
        }
        // usage kinds: surviving leaf=1; refined old parent=2 (four uses);
        // coarsened old child=3. No old leaf may disappear or serve two roles.
        std::vector<int> usage(previous_.size()),count(previous_.size()),refinedSlots(previous_.size());
        std::vector<double> next(current.size());
        for(std::size_t i=0;i<current.size();++i){const auto& now=current[i];
            const auto survivor=oldByLabel.find(now.label);
            if(survivor!=oldByLabel.end()){
                const int old=survivor->second;
                require(previous_[old].parent==now.parent&&previous_[old].level==now.level&&previous_[old].slot==now.slot,"survivor changed actual ancestry");
                require(usage[old]==0,"duplicate/mixed survivor use");usage[old]=1;count[old]=1;next[i]=this->value[old];continue;
            }
            const auto parent=oldByLabel.find(now.parent);
            if(parent!=oldByLabel.end()){
                const int old=parent->second;
                require(previous_[old].level+1==now.level&&now.slot>=0,"refinement is not one actual level");
                require(usage[old]==0||usage[old]==2,"mixed old-parent use");
                const int bit=1<<now.slot;require((refinedSlots[old]&bit)==0,"duplicate refined child slot");
                usage[old]=2;++count[old];refinedSlots[old]|=bit;next[i]=this->value[old];continue;
            }
            const auto family=oldByParent.find(now.label);
            require(family!=oldByParent.end(),"new leaf has neither actual old parent nor old sibling family");
            double target=0;
            for(int k=0;k<4;++k){const int old=family->second[k];
                require(old>=0&&previous_[old].parent==now.label&&previous_[old].slot==k
                    &&previous_[old].level==now.level+1,"coarsening requires four complete one-level old leaf siblings");
                require(usage[old]==0,"duplicate/mixed old sibling use");usage[old]=3;count[old]=1;
                target=std::max(target,this->value[old]);
            }
            next[i]=target;
        }
        for(std::size_t i=0;i<previous_.size();++i){
            require(usage[i]!=0,"old target leaf has no actual successor");
            if(usage[i]==2)require(count[i]==4&&refinedSlots[i]==15,"refinement requires all four actual new leaf children");
            else require(count[i]==1,"old target leaf used more than once");
        }
        // All provenance checks precede replacement. This deliberately avoids
        // CellField::updateWithMesh and all relatedPreIndex LS neighbors.
        CellField<double>::initialize(mesh,0.);
        for(std::size_t i=0;i<next.size();++i)this->value[i]=next[i];
        previous_=current;
    }
};
}
