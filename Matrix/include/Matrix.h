#include <cstddef>
#include <cstring>
#include <iostream>

template<typename Tp, std::size_t Rows, std::size_t Cols>
using BasicMatrix = Tp[ Rows ][ Cols ];

template<typename Tp>
class CRSMatrix {
    using size_type    = std::size_t;

    template<std::size_t Rows, std::size_t Cols>
    using basic_matrix = BasicMatrix<Tp, Rows, Cols>;

    std::size_t _rows                  = 0;
    std::size_t _cols                  = 0;
    std::size_t _nnz                   = 0;
    Tp* _v                             = nullptr;
    std::size_t* _col_index            = nullptr;
    std::size_t* _row_index            = nullptr;

public:
    CRSMatrix() = default;

    template<std::size_t Rows, std::size_t Cols>
    CRSMatrix(const basic_matrix<Rows, Cols>& m)
        : _rows(Rows),
          _cols(Cols),
          _nnz(number_of_non_zeros(m)),
          _v(new Tp[ _nnz ]),
          _col_index(new size_type[ _nnz ]),
          _row_index(new size_type[ Rows + 1 ]) {
        _row_index[ 0 ]  = 0;
        size_type nnzTmp = 0;

        for (size_type i = 0; i < _rows; i++) {
            for (size_type j = 0; j < _cols; j++) {
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
        : _rows(other._rows), _cols(other._cols), _nnz(other._nnz) {
        if (not other.empty()) {
            _v         = new Tp[ _nnz ];
            _col_index = new size_type[ _nnz ];
            _row_index = new size_type[ ridx_count() ];
            std::memcpy(_v, other._v, sizeof(Tp) * _nnz);
            std::memcpy(_col_index, other._col_index, sizeof(size_type) * _nnz);
            std::memcpy(_row_index, other._row_index, sizeof(size_type) * (ridx_count()));
        }
    }

    CRSMatrix(CRSMatrix&& other)
        : _rows(other._rows),
          _cols(other._cols),
          _nnz(other._nnz),
          _v(other._v),
          _col_index(other._col_index),
          _row_index(other._row_index) {
        other._rows      = 0;
        other._cols      = 0;
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

        _rows      = Rows;
        _cols      = Cols;
        _nnz       = number_of_non_zeros(m);
        _v         = new Tp[ _nnz ];
        _col_index = new size_type[ _nnz ];
        _row_index = new size_type[ ridx_count() ];

        _row_index[ 0 ]  = 0;
        size_type nnzTmp = 0;

        for (size_type i = 0; i < _rows; i++) {
            for (size_type j = 0; j < _cols; j++) {
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

    CRSMatrix& operator=(const CRSMatrix& other) noexcept {
        clear();
        if (not other.empty()) {
            _rows      = other._rows;
            _cols      = other._cols;
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
        _rows      = other._rows;
        _cols      = other._cols;
        _nnz       = other._nnz;
        _v         = other._v;
        _col_index = other._col_index;
        _row_index = other._row_index;

        other._rows      = 0;
        other._cols      = 0;
        other._nnz       = 0;
        other._v         = nullptr;
        other._col_index = nullptr;
        other._row_index = nullptr;
        return *this;
    }

    void clear() noexcept {
        delete[] _v;
        delete[] _col_index;
        delete[] _row_index;
        _rows      = 0;
        _cols      = 0;
        _nnz       = 0;
        _v         = nullptr;
        _col_index = nullptr;
        _row_index = nullptr;
    }

    bool empty() const noexcept { return _row_index == nullptr; }

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
    template<std::size_t Rows, std::size_t Cols>
    static constexpr size_type inline number_of_non_zeros(const basic_matrix<Rows, Cols>& m) {
        size_type nnz {};
        for (size_type i = 0; i < Rows; i++)
            for (size_type j = 0; j < Cols; j++)
                if (m[ i ][ j ]) nnz++;
        return nnz;
    }

    inline size_type ridx_count() const noexcept { return _rows + 1; }
};

template<typename Tp>
using CSRMatrix = CRSMatrix<Tp>;
