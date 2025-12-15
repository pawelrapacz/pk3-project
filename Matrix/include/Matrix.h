#include <cstddef>
#include <cstring>
#include <iostream>
#include <stdexcept>

template<typename Tp, std::size_t Rows, std::size_t Cols>
using BasicMatrix = Tp[Rows][Cols];

template<typename Tp>
class CRSMatrix {
public:
    using size_type = std::size_t;

    template<size_type Rows, size_type Cols>
    using basic_matrix = BasicMatrix<Tp, Rows, Cols>;

    struct Dimensions {
        size_type rows;
        size_type cols;

        inline constexpr bool operator==(const Dimensions& o) {
            return rows == o.rows && cols == o.cols;
        }

        inline constexpr bool operator!=(const Dimensions& o) {
            return rows != o.rows || cols != o.cols;
        }
    };


public:
    CRSMatrix() = default;

    template<size_type Rows, size_type Cols>
    CRSMatrix(const basic_matrix<Rows, Cols>& m)
        : _dim(Rows, Cols),
          _nnz(number_of_non_zeros(m)),
          _v(new Tp[_nnz]),
          _col_index(new size_type[_nnz]),
          _row_index(new size_type[Rows + 1]) {
        _row_index[0] = 0;
        size_type nnzTmp {};

        for (size_type i = 0; i < rows(); i++) {
            for (size_type j = 0; j < cols(); j++) {
                if (m[i][j]) {
                    _v[nnzTmp]         = m[i][j];
                    _col_index[nnzTmp] = j;
                    nnzTmp++;
                }
            }
            _row_index[i + 1] = nnzTmp;
        }
    }

    CRSMatrix(const CRSMatrix& other)
        : _dim(other._dim), _nnz(other._nnz) {
        if (not other.empty()) {
            _v         = new Tp[_nnz];
            _col_index = new size_type[_nnz];
            _row_index = new size_type[ridx_size()];
            std::memcpy(_v, other._v, sizeof(Tp) * _nnz);
            std::memcpy(_col_index, other._col_index, sizeof(size_type) * _nnz);
            std::memcpy(_row_index, other._row_index, sizeof(size_type) * (ridx_size()));
        }
    }

    CRSMatrix(CRSMatrix&& other) noexcept
        : _dim(other._dim),
          _nnz(other._nnz),
          _v(other._v),
          _col_index(other._col_index),
          _row_index(other._row_index) {
        other._dim       = Dimensions();
        other._nnz       = size_type();
        other._v         = nullptr;
        other._col_index = nullptr;
        other._row_index = nullptr;
    }

    ~CRSMatrix() {
        delete[] _v;
        delete[] _col_index;
        delete[] _row_index;
    }

    template<size_type Rows, size_type Cols>
    CRSMatrix& operator=(const basic_matrix<Rows, Cols>& m) {
        clear();

        _dim.rows  = Rows;
        _dim.cols  = Cols;
        _nnz       = number_of_non_zeros(m);
        _v         = new Tp[_nnz];
        _col_index = new size_type[_nnz];
        _row_index = new size_type[ridx_size()];

        _row_index[0] = 0;
        size_type nnzTmp {};

        for (size_type i = 0; i < rows(); i++) {
            for (size_type j = 0; j < cols(); j++) {
                if (m[i][j]) {
                    _v[nnzTmp]         = m[i][j];
                    _col_index[nnzTmp] = j;
                    nnzTmp++;
                }
            }
            _row_index[i + 1] = nnzTmp;
        }
        return *this;
    }

    CRSMatrix& operator=(const CRSMatrix& other) {
        clear();
        if (not other.empty()) {
            _dim       = other._dim;
            _nnz       = other._nnz;
            _v         = new Tp[_nnz];
            _col_index = new size_type[_nnz];
            _row_index = new size_type[ridx_size()];
            std::memcpy(_v, other._v, sizeof(Tp) * _nnz);
            std::memcpy(_col_index, other._col_index, sizeof(size_type) * _nnz);
            std::memcpy(_row_index, other._row_index, sizeof(size_type) * (ridx_size()));
        }
        return *this;
    }

    CRSMatrix& operator=(CRSMatrix&& other) noexcept {
        if (this != &other) {
            _dim       = other._dim;
            _nnz       = other._nnz;
            _v         = other._v;
            _col_index = other._col_index;
            _row_index = other._row_index;

            other._dim       = Dimensions();
            other._nnz       = size_type();
            other._v         = nullptr;
            other._col_index = nullptr;
            other._row_index = nullptr;
        }
        return *this;
    }

    void clear() noexcept {
        delete[] _v;
        delete[] _col_index;
        delete[] _row_index;
        _dim       = Dimensions();
        _nnz       = size_type();
        _v         = nullptr;
        _col_index = nullptr;
        _row_index = nullptr;
    }

    inline bool is_zero_matrix() const noexcept {
        return _nnz == size_type();
    }

    inline bool empty() const noexcept {
        return _row_index == nullptr;
    }

    inline Dimensions dim() const noexcept {
        return _dim;
    }

    inline size_type rows() const noexcept {
        return _dim.rows;
    }

    inline size_type cols() const noexcept {
        return _dim.cols;
    }

    // scalar multiplication
    inline CRSMatrix operator*(Tp val) const {
        CRSMatrix out = *this;
        out *= val;
        return out;
    }

    inline CRSMatrix& operator*=(Tp val) noexcept {
        for (size_type i = 0; i < _nnz; i++) _v[i] *= val;
        return *this;
    }

    // matrix multiplication
    CRSMatrix operator*(const CRSMatrix& other) const {
        if (cols() != other.rows())
            throw std::invalid_argument("The number of columns in the first matrix must be equal "
                                        "to the number of rows in the second matrix.");

        CRSMatrix out;
        out._dim          = { rows(), other.cols() };
        out._row_index    = new size_type[ridx_size()];
        out._row_index[0] = 0;

        CRSMatrix m = other;
        m.transpose();

        // wyznaczanie nnz i row_index
        for (size_type i = 0; i < rows(); i++) {
            if (nnz_row(i)) {
                for (size_type j = 0; j < m.rows(); j++) {
                    if (!m.nnz_row(j))
                        continue;
                    Tp val {};
                    size_type posA = _row_index[i], posB = m._row_index[j];
                    while (posA < _row_index[i + 1] && posB < m._row_index[j + 1]) {
                        if (_col_index[posA] == m._col_index[posB]) {
                            val += _v[posA] * m._v[posB];
                            posA++;
                            posB++;
                        }
                        else if (_col_index[posA] < m._col_index[posB]) {
                            posA++;
                        }
                        else {
                            posB++;
                        }
                    }

                    if (val) {
                        out._nnz++;
                    }
                }
            }
            out._row_index[i + 1] = out._nnz;
        }

        out._v         = new Tp[out._nnz];
        out._col_index = new size_type[out._nnz];
        out._nnz       = size_type();


        for (size_type i = 0; i < rows(); i++) {
            if (nnz_row(i)) {
                for (size_type j = 0; j < other.cols(); j++) {
                    if (!m.nnz_row(j))
                        continue;
                    Tp val {};
                    size_type posA = _row_index[i], posB = m._row_index[j];
                    while (posA < _row_index[i + 1] && posB < m._row_index[j + 1]) {
                        if (_col_index[posA] == m._col_index[posB]) {
                            val += _v[posA] * m._v[posB];
                            posA++;
                            posB++;
                        }
                        else if (_col_index[posA] < m._col_index[posB]) {
                            posA++;
                        }
                        else {
                            posB++;
                        }
                    }

                    if (val) {
                        out._v[out._nnz]         = val;
                        out._col_index[out._nnz] = j;
                        out._nnz++;
                    }
                }
            }
        }

        return out;
    }

    inline CRSMatrix& operator*=(const CRSMatrix& other) {
        *this = *this * other;
        return *this;
    }

    // addition
    CRSMatrix operator+(const CRSMatrix& other) const {
        if (dim() != other.dim())
            throw std::invalid_argument("The dimensions of both matricies must be equal.");

        CRSMatrix out;
        out._dim          = dim();
        out._row_index    = new size_type[out.ridx_size()];
        out._row_index[0] = 0;

        // wyznaczanie nnz i _row_index
        for (size_type i = 0; i < rows(); i++) {
            size_type posA = _row_index[i], posB = other._row_index[i];

            while (posA < _row_index[i + 1] && posB < other._row_index[i + 1]) {
                if (_col_index[posA] == other._col_index[posB]) {
                    if (_v[posA] + other._v[posB] != Tp())
                        out._nnz++;

                    posA++;
                    posB++;
                }
                else if (_col_index[posA] < other._col_index[posB]) {
                    posA++;
                    out._nnz++;
                }
                else {
                    posB++;
                    out._nnz++;
                }
            }
            // obliczamy pozostałe nnz w wierszu (jeden z nich zawsze 0)
            out._nnz += (_row_index[i + 1] - posA) + (other._row_index[i + 1] - posB);
            out._row_index[i + 1] = out._nnz;
        }

        out._v         = new Tp[out._nnz];
        out._col_index = new size_type[out._nnz];
        out._nnz       = size_type();

        for (size_type i = 0; i < rows(); i++) {
            size_type posA = _row_index[i], posB = other._row_index[i];

            // zakresy pokrywają się
            while (posA < _row_index[i + 1] && posB < other._row_index[i + 1]) {
                if (_col_index[posA] == other._col_index[posB]) {
                    if (_v[posA] + other._v[posB] != Tp()) {
                        out._v[out._nnz] = _v[posA] + other._v[posB];

                        out._col_index[out._nnz] = _col_index[posA];
                        out._nnz++;
                    }
                    posA++;
                    posB++;
                }
                else if (_col_index[posA] < other._col_index[posB]) {
                    out._v[out._nnz]           = _v[posA];
                    out._col_index[out._nnz++] = _col_index[posA++];
                }
                else {
                    out._v[out._nnz]           = other._v[posB];
                    out._col_index[out._nnz++] = other._col_index[posB++];
                }
            }

            while (posA < _row_index[i + 1]) {
                out._v[out._nnz]           = _v[posA];
                out._col_index[out._nnz++] = _col_index[posA++];
            }

            while (posB < other._row_index[i + 1]) {
                out._v[out._nnz]           = other._v[posB];
                out._col_index[out._nnz++] = other._col_index[posB++];
            }
        }

        return out;
    }

    inline CRSMatrix& operator+=(const CRSMatrix& other) {
        *this = *this + other;
        return *this;
    }

    inline CRSMatrix operator-() const {
        CRSMatrix out(*this);
        for (size_type i = 0; i < out._nnz; i++) out._v[i] = -out._v[i];

        return out;
    }

    inline CRSMatrix operator-(const CRSMatrix& other) const {
        return *this + (-other);
    }

    inline CRSMatrix& operator-=(const CRSMatrix& other) {
        *this = *this - other;
        return *this;
    }

    void transpose() {
        CRSMatrix n;
        n._dim       = { cols(), rows() };
        n._nnz       = _nnz;
        n._v         = new Tp[_nnz];
        n._col_index = new size_type[_nnz];
        n._row_index = new size_type[cols() + 1]();

        // kopia _row_index, pozwala określić odpowiednią pozycję wartości
        auto pos_ptr = new size_type[cols() + 1];
        pos_ptr[0]   = 0;

        for (size_type i = 0; i < _nnz; i++) n._row_index[_col_index[i] + 1]++;

        for (size_type i = 0; i < n.rows(); i++) {
            n._row_index[i + 1] += n._row_index[i];
            pos_ptr[i + 1] = n._row_index[i + 1];  // kopia
        }

        // idziemy po wierszach, czyli nowch kolumnach
        for (size_type i = 0; i < rows(); i++) {
            // teraz wyłuskujemy każdy wiersz i go zmieniamy na kolumnę
            for (size_type j = _row_index[i]; j < _row_index[i + 1]; j++) {
                auto new_row = _col_index[j];  // poprzednia kolumna to nowy wiersz
                auto pos     = pos_ptr
                    [new_row]++;  // zebranie pozycji i przesunięcie na kolejną pozycję w nowym wierszu

                n._v[pos]         = _v[j];
                n._col_index[pos] = i;  // nowa kolumna to stary wiersz
            }
        }

        delete[] pos_ptr;
        *this = std::move(n);
    }

    void print() const {
        if (!empty()) {
            std::cout << "\nV = [ ";
            for (size_type i = 0; i < _nnz; i++) std::cout << _v[i] << ", ";

            std::cout << "]\nCOL_INDEX = [ ";
            for (size_type i = 0; i < _nnz; i++) std::cout << _col_index[i] << ", ";

            std::cout << "]\nROW_INDEX = [ ";
            for (size_type i = 0; i < ridx_size(); i++) std::cout << _row_index[i] << ", ";
            std::cout << "]\n";
        }
        else {
            std::cout << "\nEMPTY!";
        }
    }

    void printm() const {
        if (!empty()) {
            std::cout << '\n';
            for (size_type i = 0; i < rows(); i++) {
                auto pos = _row_index[i];
                for (size_type j = 0; j < cols(); j++) {
                    if (pos < _row_index[i + 1] && j == _col_index[pos])
                        std::cout << _v[pos++] << "\t";
                    else
                        std::cout << "0\t";
                }
                std::cout << std::endl;
            }
        }
        else {
            std::cout << "\nEMPTY!";
        }
    }


protected:
    inline size_type ridx_size() const noexcept {
        return _dim.rows + 1;
    }

    inline size_type nnz_row(size_type idx) const noexcept {
        return _row_index[idx + 1] - _row_index[idx];
    }

    inline size_type nnz_col(size_type idx) const noexcept {
        size_type count {};
        for (size_type i = 0; i < _nnz; i++)
            if (_col_index[i] == idx)
                count++;

        return count;
    }

    template<size_type Rows, size_type Cols>
    static constexpr size_type inline number_of_non_zeros(
        const basic_matrix<Rows, Cols>& m) noexcept {
        size_type nnz {};
        for (size_type i = 0; i < Rows; i++)
            for (size_type j = 0; j < Cols; j++)
                if (m[i][j])
                    nnz++;
        return nnz;
    }


private:
    Dimensions _dim       = Dimensions();
    size_type _nnz        = 0;
    Tp* _v                = nullptr;
    size_type* _col_index = nullptr;
    size_type* _row_index = nullptr;
};

template<typename Tp>
using CSRMatrix = CRSMatrix<Tp>;
