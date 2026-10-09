#pragma once

#include <iostream>
#include <cstddef>
#include <memory>
#include <utility>
#include <string>
#include <fstream>
#include <sstream>
#include "Core/Mesh/Vert.h"
#include "Core/Mesh/Edge.h"
#include "Core/Mesh/Cell.h"
#include "Core/Base/Interpolation.h"
#include "Core/Base/Base.h"

// Preserve unique array ownership while propagating const to the elements.
// A const Field may be borrowed across module boundaries without granting
// write access through unique_ptr::operator[] or unique_ptr::get(). This is
// still a borrow, not a snapshot; the owning Field may change its storage.
template<typename T>
class FieldValues
{
public:
    FieldValues() = default;
    FieldValues(FieldValues&&) noexcept = default;
    FieldValues& operator=(FieldValues&&) noexcept = default;
    FieldValues(const FieldValues&) = delete;
    FieldValues& operator=(const FieldValues&) = delete;
    FieldValues& operator=(std::unique_ptr<T[]> values) noexcept
    { values_ = std::move(values); return *this; }

    T& operator[](std::size_t index) noexcept { return values_[index]; }
    const T& operator[](std::size_t index) const noexcept { return values_[index]; }
    T* get() noexcept { return values_.get(); }
    const T* get() const noexcept { return values_.get(); }
    explicit operator bool() const noexcept { return static_cast<bool>(values_); }
    friend bool operator==(const FieldValues& a, const FieldValues& b) noexcept { return a.get() == b.get(); }
    friend bool operator!=(const FieldValues& a, const FieldValues& b) noexcept { return a.get() != b.get(); }
    friend bool operator==(const FieldValues& values, std::nullptr_t) noexcept { return !values; }
    friend bool operator!=(const FieldValues& values, std::nullptr_t) noexcept { return static_cast<bool>(values); }
    friend bool operator==(std::nullptr_t, const FieldValues& values) noexcept { return !values; }
    friend bool operator!=(std::nullptr_t, const FieldValues& values) noexcept { return static_cast<bool>(values); }
    void reset(T* values = nullptr) noexcept { values_.reset(values); }
private:
    std::unique_ptr<T[]> values_;
};

template<typename T>
class Field
{
public:
    FieldValues<T> value;

    Field();
    Field(const Field<T>& field);
    Field(int size);
    Field(int size, const T& value);
    ~Field();

    inline int size() const { return size_; }
    inline int label(int i) const { return label_[i]; }

    void inline setLabel(int i, int label) { label_[i] = label; }

    Field& operator=(const Field& field);
    // Exchange complete native storage without allocation, including labels.
    // Used when restoring an already prepared field snapshot.
    void swap(Field& other) noexcept {
        std::swap(size_,other.size_);
        label_.swap(other.label_);
        std::swap(value,other.value);
    }
    Field operator+(const Field& field) const;
    Field operator-(const Field& field) const;
    Field operator*(const double& scale) const;

    void initialize(int size, const T& value);

    template<typename Mesh>
    void output(const Mesh& mesh, const std::string& _name) const;

protected:
    int size_ = 0;
    std::unique_ptr<int[]> label_;

};

#include "Core/Field/Field_impl.h"
