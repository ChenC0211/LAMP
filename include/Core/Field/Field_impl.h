#pragma once

using namespace std;

template<typename T>
Field<T>::Field() : size_(0) {}

template<typename T>
Field<T>::Field(const Field<T>& field) : size_(field.size())
{
    this->value = std::make_unique<T[]>(size_);
    this->label_ = std::make_unique<int[]>(size_);
#pragma omp parallel for
    for (int i = 0; i < size_; i++)
    {
        this->value[i] = field.value[i];
        this->label_[i] = field.label(i);
    }
}

template<typename T>
Field<T>::Field(int size)
{
    this->size_ = size;
    this->value = std::make_unique<T[]>(size_);
    this->label_ = std::make_unique<int[]>(size_);
}

template<typename T>
Field<T>::Field(int size, const T& value)
{
    this->size_ = size;
    this->value = std::make_unique<T[]>(size_);
    this->label_ = std::make_unique<int[]>(size_);
#pragma omp parallel for
    for (int i = 0; i < size_; i++)
    {
        this->value[i] = value;
    }
}

template<typename T>
Field<T>::~Field() {}

template<typename T>
Field<T>& Field<T>::operator=(const Field& field)
{
    if (this == &field)
    {
        return *this;
    }

    if (size_ != field.size())
    {
        this->value = std::make_unique<T[]>(field.size());
        this->label_ = std::make_unique<int[]>(field.size());
        this->size_ = field.size();
    }

    std::copy(field.value.get(), field.value.get() + size_, this->value.get());
    std::copy(field.label_.get(), field.label_.get() + size_, this->label_.get());

    return *this;
}

template<typename T>
Field<T> Field<T>::operator+(const Field& field) const
{
    if (size_ != field.size())
    {
        cerr << "[Field] ERROR: '+' operand sizes do not match." << endl;
        exit(1);
    }

    Field<T> result(size_);
#pragma omp parallel for
    for (int i = 0; i < size_; i++)
    {
        if (label_[i] != field.label(i))
        {
            cerr << "[Field] ERROR: field labels do not match." << endl;
            exit(1);
        }
        result.value[i] = value[i] + field.value[i];
    }

    return result;
}

template<typename T>
Field<T> Field<T>::operator-(const Field& field) const
{
    if (size_ != field.size())
    {
        cerr << "[Field] ERROR: '-' operand sizes do not match." << endl;
        exit(1);
    }
    
    Field<T> result(this->size_);
#pragma omp parallel for
    for (int i = 0; i < this->size_; i++)
    {
        if (label_[i] != field.label(i))
        {
            cerr << "[Field] ERROR: field labels do not match." << endl;
            exit(1);
        }
        result.value[i] = this->value[i] - field.value[i];
    }

    return result;
}

template<typename T>
Field<T> Field<T>::operator*(const double& scale) const
{
    Field<T> result(this->size_);
#pragma omp parallel for
    for (int i = 0; i < this->size_; i++)
    {
        result.value[i] = this->value[i] * scale;
    }

    return result;
}

template<typename T>
void Field<T>::initialize(int size, const T& value)
{
    this->size_ = size;
    this->value = std::make_unique<T[]>(size_);
    this->label_ = std::make_unique<int[]>(size_);
#pragma omp parallel for
    for (int i = 0; i < size_; i++)
    {
        this->value[i] = value;
        this->label_[i] = i;
    }
}

template<typename T>
template<typename Mesh>
void Field<T>::output(const Mesh& mesh, const string& _name) const
{
    cout << endl << "[Field] Writing field data: " << _name << "..." << endl;

    if (size_ != mesh.vertIndexNum())
    {
        cerr << "[Field] ERROR: field size does not match mesh vertex count." << endl;
        cerr << "  fieldSize=" << size_ << " meshVertCount=" << mesh.vertIndexNum() << endl;
        exit(-1);
    }

    ostringstream fileName;
	fileName << _name;
	string filePath = "./output/" + fileName.str() + ".dat";

	ofstream fp(filePath, ofstream::out);
	if (!fp.is_open())
	{
        cerr << "[Field] ERROR: failed to write field file: " << filePath << endl;
		exit(-1);
	}

	fp << "Title= \"" << _name << "\"" << endl;
	fp << "Variables= \"X\" \"Y\" \"value\"" << endl;
	fp << "ZONE N = " << mesh.vertIndexNum() << ", E = " << mesh.cellIndexNum()
		<< ", F = FEPOINT, ET = QUADRILATERAL" << endl;

	for (int i = 0; i < mesh.vertIndexNum(); i++)
    {
        const int vertI = mesh.vertIndexToLabel[i];
        Vert& vert = mesh.verts(vertI);
        fp << vert.coord().x << " " << vert.coord().y << " "
            << value[i] << endl;
    }

	for (int i = 0; i < mesh.cellIndexNum(); i++)
    {
        const int cellI = mesh.cellIndexToLabel[i];
        const Cell& cell = mesh.cells(cellI);
        for (int j = 0; j < 4; j++)
        {
            const Vert& vert = mesh.verts(cell.cellToVert[j]);
            fp << vert.index() + 1 << " ";
        }
        fp << endl;
    }

	fp.close();
    cout << "[Field] Field written: " << _name << "." << endl;
}