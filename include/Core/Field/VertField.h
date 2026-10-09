#pragma once

#include <iostream>
#include "Core/Field/Field.h"

template<typename T>
class CellField;

template<typename T>
class HEdgeField;

template<typename T>
class VEdgeField;

template<typename T>
class VertField : public Field<T>
{
public:
    VertField();
    VertField(const VertField& field);
    template<typename Mesh>
    VertField(Mesh& mesh, const T& value);
    ~VertField();

    template<typename Mesh>
    void initialize(const Mesh& mesh, const T& value);

    template<typename Mesh>
    void updateWithMesh(const Mesh& mesh);

    template<typename Mesh>
    std::unique_ptr<CellField<T>> toCellField(const Mesh& mesh) const;

    template<typename Mesh>
    std::unique_ptr<HEdgeField<T>> toHEdgeField(const Mesh& mesh) const;

    template<typename Mesh>
    std::unique_ptr<VEdgeField<T>> toVEdgeField(const Mesh& mesh) const;

private:
    template<typename Mesh>
    void subInterpolation(const Mesh& mesh, const Cell& cell);

};

#include "Core/Field/VertField_impl.h"
