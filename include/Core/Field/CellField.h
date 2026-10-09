#pragma once

#include <iostream>
#include "Core/Field/Field.h"

template<typename T>
class VertField;

template<typename T>
class CellField : public Field<T>
{
public:
    CellField();
    CellField(const CellField& field);
    template<typename Mesh>
    CellField(Mesh& mesh, const T& value);
    ~CellField();

    template<typename Mesh>
    void initialize(const Mesh& mesh, const T& value);

    template<typename Mesh>
    void updateWithMesh(const Mesh& mesh);

    template<typename Mesh>
    std::unique_ptr<VertField<T>> toVertField(const Mesh& mesh) const;

    template<typename Mesh>
    void output(const Mesh& mesh, const std::string& _name) const;

private:
    // template<typename Mesh>
    // void subInterpolation(const Mesh& mesh, std::unique_ptr<T[]>& values, unique_ptr<int[]>& labels, int preI, Cell& cell);

};

#include "Core/Field/CellField_impl.h"
