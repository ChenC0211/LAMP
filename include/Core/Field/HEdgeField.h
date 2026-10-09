#pragma once

#include <iostream>
#include "Core/Field/Field.h"

template<typename T>
class VertField;

template<typename T>
class HEdgeField : public Field<T>
{
public:
    HEdgeField();
    HEdgeField(const HEdgeField& field);
    template<typename Mesh>
    HEdgeField(Mesh& mesh, const T& value);
    ~HEdgeField();

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
    // void subInterpolation(const Mesh& mesh, std::unique_ptr<T[]>& values, unique_ptr<int[]>& labels, int preI, Edge& hEdge);
    template<typename Mesh>
    void findRelated(const Mesh& mesh, Cell** cells, const Vert& vert, std::vector<double>& values, std::vector<Vector2D>& coords) const;
};

#include "Core/Field/HEdgeField_impl.h"
