#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <utility>
#include <stack>
#include <regex>
#include <vector>
#include "Core/Base/Vector2D.h"
#include "Core/Base/Structures.h"

// Forward declare to keep Dictionary lightweight.
class UnitsContext;

class Dictionary
{
public:
    Dictionary();
    virtual ~Dictionary();

protected: 
    std::string trim(const std::string& line) const;
    std::string removeInlineComments(const std::string &line) const;
    std::vector<std::string> split(const std::string &str) const;
    std::vector<std::string> split(const std::string& str, const std::string& delimiter) const;
    bool parseCoord(const std::string& line, Vector2D& coord) const;
    // Parse a coordinate where each component may include a unit, e.g.
    //   (3 cm, 0 cm)
    //   (3, 0) cm
    // Values with units are converted to code coordinates using units->LengthSIToCode.
    // If units==nullptr, unit-bearing inputs are interpreted as SI and returned directly.
    bool parseCoordWithUnits(const std::string& line, Vector2D& coord, const UnitsContext* units) const;
    bool parseIndex(const std::string& line, std::vector<SignInt>& index) const;
    bool parseIndex_brace(const std::string& line, std::vector<SignInt>& index) const;
    bool parseCoords(const std::string& line, std::vector<Vector2D>& coords) const;
    // Parse multiple coordinate pairs in one line, supporting optional units.
    bool parseCoordsWithUnits(const std::string& line, std::vector<Vector2D>& coords, const UnitsContext* units) const;

private:

};