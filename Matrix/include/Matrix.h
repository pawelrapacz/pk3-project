#include <cstddef>
#include <cstring>
#include <iostream>

template<typename Tp, std::size_t Rows, std::size_t Colms = Rows>
class CRSMatrix {
    using size_type    = std::size_t;
    using basic_matrix = Tp[ Rows ][ Colms ];

    static inline constexpr size_type RowsI = Rows + 1;

    std::size_t _nnz                   = 0;
    Tp* _v                             = nullptr;
    std::size_t* _col_index            = nullptr;
    std::size_t _row_index[ Rows + 1 ] = { 0 };

public:
    CRSMatrix() = default;

    CRSMatrix(const basic_matrix& m)
        : _nnz(number_of_non_zeros(m)), _v(new Tp[ _nnz ]), _col_index(new size_type[ _nnz ]) {
        _row_index[ 0 ]  = 0;
        size_type nnzTmp = 0;

        for (size_type i = 0; i < Rows; i++) {
            for (size_type j = 0; j < Colms; j++) {
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
        : _nnz(other._nnz), _v(new Tp[ _nnz ]), _col_index(new size_type[ _nnz ]) {
        std::memcpy(_v, other._v, sizeof(Tp) * _nnz);
        std::memcpy(_col_index, other._col_index, sizeof(size_type) * _nnz);
        std::memcpy(_row_index, other._row_index, sizeof(_row_index));
    }

    CRSMatrix(CRSMatrix&& other)
        : _nnz(other._nnz), _v(other._v), _col_index(other._col_index) {
        std::memcpy(_row_index, other._row_index, sizeof(_row_index));

        other._v         = nullptr;
        other._col_index = nullptr;
        std::memset(other._row_index, 0, sizeof(other._row_index));
    }

    ~CRSMatrix() {
        delete[] _v;
        delete[] _col_index;
    }

    CRSMatrix& operator=(const CRSMatrix& other) noexcept {
        _nnz = other._nnz;
        std::memcpy(_v, other._v, sizeof(Tp) * _nnz);
        std::memcpy(_col_index, other._col_index, sizeof(size_type) * _nnz);
        std::memcpy(_row_index, other._row_index, sizeof(_row_index));
    }

    CRSMatrix& operator=(CRSMatrix&& other) noexcept {
        _nnz       = other._nnz;
        _v         = other._v;
        _col_index = other._col_index;
        std::memcpy(_row_index, other._row_index, sizeof(_row_index));

        other._v         = nullptr;
        other._col_index = nullptr;
        std::memset(other._row_index, 0, sizeof(other._row_index));
    }

    void print() const {
        std::cout << "\nV = [ ";
        for (size_type i = 0; i < _nnz; i++) std::cout << _v[ i ] << ", ";

        std::cout << "]\nCOL_INDEX = [ ";
        for (size_type i = 0; i < _nnz; i++) std::cout << _col_index[ i ] << ", ";

        std::cout << "]\nROW_INDEX = [ ";
        for (size_type i = 0; i < Rows + 1; i++) std::cout << _row_index[ i ] << ", ";
        std::cout << "]\n";
    }

protected:
    static constexpr size_type inline number_of_non_zeros(const basic_matrix& m) {
        size_type nnz {};
        for (size_type i = 0; i < Rows; i++)
            for (size_type j = 0; j < Colms; j++)
                if (m[ i ][ j ]) nnz++;
        return nnz;
    }
};

template<typename Tp, std::size_t M, std::size_t N>
using CSRMatrix = CRSMatrix<Tp, M, N>;
