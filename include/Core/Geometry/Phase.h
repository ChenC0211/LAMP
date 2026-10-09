#pragma once

#include <iostream>
#include <string>
#include "Core/Base/Enumeration.h"
#include <vector>
#include "Core/Base/Structures.h"

class Phase
{
public:
    Phase();
    ~Phase();

    inline std::string name() const { return name_; }
    inline std::vector<LevelSetOperation> operation() const { return operation_; }
    inline std::vector<SignInt> levelSetIndexA() const { return levelSetIndexA_; }
    inline std::vector<SignInt> levelSetIndexB() const { return levelSetIndexB_; }
    inline std::vector<SignInt> boundaryLines() const { return boundaryLines_; }

    inline void setName(const std::string& name) { name_ = name; }
    inline void setOperation(std::vector<LevelSetOperation> operation) { operation_ = operation; }
    inline void setLevelSetIndexA(std::vector<SignInt> levelSetIndexA) { levelSetIndexA_ = levelSetIndexA; }
    inline void setLevelSetIndexB(std::vector<SignInt> levelSetIndexB) { levelSetIndexB_ = levelSetIndexB; }
    inline void setBoundaryLines(std::vector<SignInt> boundaryLines) { boundaryLines_ = boundaryLines; }

private:
    std::string name_;

    std::vector<LevelSetOperation> operation_;
    std::vector<SignInt> levelSetIndexA_;
    std::vector<SignInt> levelSetIndexB_;

    std::vector<SignInt> boundaryLines_;

};