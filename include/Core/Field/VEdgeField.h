#pragma once

#include <iostream>
#include "Core/Field/Field.h"

template<typename T>
class VertField;

template<typename T>
class VEdgeField : public Field<T>
{
public:
    VEdgeField();
    VEdgeField(const VEdgeField& field);
    template<typename Mesh>
    VEdgeField(Mesh& mesh, const T& value);
    ~VEdgeField();

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
    // void subInterpolation(const Mesh& mesh, std::unique_ptr<T[]>& values, unique_ptr<int[]>& labels, int preI, Edge& vEdge);
    template<typename Mesh>
    void findRelated(const Mesh& mesh, Cell** cells, const Vert& vert, std::vector<double>& values, std::vector<Vector2D>& coords) const;
};

#include "Core/Field/VEdgeField_impl.h"
