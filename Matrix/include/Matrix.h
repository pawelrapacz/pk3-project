#include <cstddef>
#include <cstring>
#include <iostream>

template<typename Tp, std::size_t Rows, std::size_t Cols>
using BasicMatrix = Tp[ Rows ][ Cols ];

template<typename Tp>
class CRSMatrix {
public:
    using size_type    = std::size_t;

    template<std::size_t Rows, std::size_t Cols>
    using basic_matrix = BasicMatrix<Tp, Rows, Cols>;

    struct Dimensions {
        size_type rows;
        size_type cols;
    };


public:
    CRSMatrix() = default;

    template<std::size_t Rows, std::size_t Cols>
    CRSMatrix(const basic_matrix<Rows, Cols>& m)
        : _dim(Rows, Cols),
          _nnz(number_of_non_zeros(m)),
          _v(new Tp[ _nnz ]),
          _col_index(new size_type[ _nnz ]),
          _row_index(new size_type[ Rows + 1 ]) {
        _row_index[ 0 ]  = 0;
        size_type nnzTmp = 0;

        for (size_type i = 0; i < rows(); i++) {
            for (size_type j = 0; j < cols(); j++) {
                if (m[ i ][ j ]) {
                    _v[ nnzTmp ]         = m[ i ][ j ];
                    _col_index[ nnzTmp ] = j;
                    nnzTmp++;
                }
            }
            _row_index[ i + 1 ] = nnzTmp;
        }
    }

    CRSMatrix(const CRSMatrix& other)
        : _dim(other._dim), _nnz(other._nnz) {
        if (not other.empty()) {
            _v         = new Tp[ _nnz ];
            _col_index = new size_type[ _nnz ];
            _row_index = new size_type[ ridx_count() ];
            std::memcpy(_v, other._v, sizeof(Tp) * _nnz);
            std::memcpy(_col_index, other._col_index, sizeof(size_type) * _nnz);
            std::memcpy(_row_index, other._row_index, sizeof(size_type) * (ridx_count()));
        }
    }

    CRSMatrix(CRSMatrix&& other) noexcept
        : _dim(other._dim),
          _nnz(other._nnz),
          _v(other._v),
          _col_index(other._col_index),
          _row_index(other._row_index) {
        other._dim = Dimensions();
        other._nnz       = 0;
        other._v         = nullptr;
        other._col_index = nullptr;
        other._row_index = nullptr;
    }

    ~CRSMatrix() {
        delete[] _v;
        delete[] _col_index;
        delete[] _row_index;
    }

    template<std::size_t Rows, std::size_t Cols>
    CRSMatrix& operator=(const basic_matrix<Rows, Cols>& m) {
        clear();

        _dim.rows      = Rows;
        _dim.cols      = Cols;
        _nnz       = number_of_non_zeros(m);
        _v         = new Tp[ _nnz ];
        _col_index = new size_type[ _nnz ];
        _row_index = new size_type[ ridx_count() ];

        _row_index[ 0 ]  = 0;
        size_type nnzTmp = 0;

        for (size_type i = 0; i < rows(); i++) {
            for (size_type j = 0; j < cols(); j++) {
                if (m[ i ][ j ]) {
                    _v[ nnzTmp ]         = m[ i ][ j ];
                    _col_index[ nnzTmp ] = j;
                    nnzTmp++;
                }
            }
            _row_index[ i + 1 ] = nnzTmp;
        }
        return *this;
    }

    CRSMatrix& operator=(const CRSMatrix& other) {
        clear();
        if (not other.empty()) {
            _dim = other._dim;
            _nnz       = other._nnz;
            _v         = new Tp[ _nnz ];
            _col_index = new size_type[ _nnz ];
            _row_index = new size_type[ ridx_count() ];
            std::memcpy(_v, other._v, sizeof(Tp) * _nnz);
            std::memcpy(_col_index, other._col_index, sizeof(size_type) * _nnz);
            std::memcpy(_row_index, other._row_index, sizeof(size_type) * (ridx_count()));
        }
        return *this;
    }

    CRSMatrix& operator=(CRSMatrix&& other) noexcept {
        if (this != &other) {
            _dim = other._dim;
            _nnz       = other._nnz;
            _v         = other._v;
            _col_index = other._col_index;
            _row_index = other._row_index;

            other._dim = Dimensions();
            other._nnz       = 0;
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
        _dim = Dimensions();
        _nnz       = 0;
        _v         = nullptr;
        _col_index = nullptr;
        _row_index = nullptr;
    }

    inline bool empty() const noexcept
    { return _row_index == nullptr; }

    inline Dimensions dim() const noexcept
    { return _dim; }

    inline size_type rows() const noexcept
    { return _dim.rows; }

    inline size_type cols() const noexcept
    { return _dim.cols; }

    

    // scalar multiplication
    CRSMatrix operator*(Tp val) const {
        CRSMatrix out = *this;
        out *= val;
        return out;
    }
    
    CRSMatrix& operator*=(Tp val) noexcept {
        for (size_type i = 0; i < _nnz; i++)
            _v[i] *= val;
    }
    
    void transpose() {
        CRSMatrix n;
        n._dim = {cols(), rows()};
        n._nnz = _nnz;
        n._v = new Tp[_nnz];
        n._col_index = new size_type[_nnz];
        n._row_index = new size_type[cols() + 1]();

        auto pos_ptr = new size_type[cols() + 1]; // kopia _row_index, pozwala określić odpowiednią pozycję wartości
        pos_ptr[0] = 0;

        for (size_type i = 0; i < _nnz; i++)
            n._row_index[_col_index[i] + 1]++;

        for (size_type i = 0; i < n.ridx_count(); i++) {
            n._row_index[i + 1] += n._row_index[i];
            pos_ptr[i + 1] = n._row_index[i + 1]; // kopia
        }

        // idziemy po wierszach, czyli nowch kolumnach
        for (size_type i = 0; i < rows(); i++) {
            // teraz wyłuskujemy każdy wiersz i go zmieniamy na kolumnę
            for (size_type j = _row_index[i]; j < _row_index[i + 1]; j++) {
                auto new_row = _col_index[j]; // poprzednia kolumna to nowy wiersz
                auto pos = pos_ptr[new_row]++; // zebranie pozycji i przesunięcie na kolejną pozycję w nowym wierszu

                n._v[pos] = _v[j];
                n._col_index[pos] = i; // nowa kolumna to stary wiersz
            }
        }

        delete[] pos_ptr;
        *this = std::move(n);
    }
    

    void print() const {
        if (_v) {
            std::cout << "\nV = [ ";
            for (size_type i = 0; i < _nnz; i++) std::cout << _v[ i ] << ", ";

            std::cout << "]\nCOL_INDEX = [ ";
            for (size_type i = 0; i < _nnz; i++) std::cout << _col_index[ i ] << ", ";

            std::cout << "]\nROW_INDEX = [ ";
            for (size_type i = 0; i < ridx_count(); i++) std::cout << _row_index[ i ] << ", ";
            std::cout << "]\n";
        }
        else { std::cout << "\nV = [ ]\nCOL_INDEX = [ ]\nROW_INDEX = [ ]\n"; }
    }


protected:
    inline size_type ridx_count() const noexcept
    { return _dim.rows + 1; }

    inline size_type nnz_row(size_type idx) const noexcept {
        return _row_index[idx + 1] - _row_index[idx];
    }

    inline size_type nnz_col(size_type idx) const noexcept {
        size_type count { };
        for (size_type i = 0; i < _nnz; i++)
            if (_col_index[i] == idx) count++;
        
        return count;
    }

    template<std::size_t Rows, std::size_t Cols>
    static constexpr size_type inline number_of_non_zeros(const basic_matrix<Rows, Cols>& m) noexcept {
        size_type nnz {};
        for (size_type i = 0; i < Rows; i++)
            for (size_type j = 0; j < Cols; j++)
                if (m[ i ][ j ]) nnz++;
        return nnz;
    }



private:
    Dimensions _dim = Dimensions();
    std::size_t _nnz                   = 0;
    Tp* _v                             = nullptr;
    std::size_t* _col_index            = nullptr;
    std::size_t* _row_index            = nullptr;
};


template<typename Tp>
using CSRMatrix = CRSMatrix<Tp>;
