#include <iostream>
#include <vector>
#import <cmath>

using namespace std;

template<typename T>
struct Matrix {
    int column, row;
    vector<vector<T> > model;

    Matrix() {}

    Matrix(vector<vector<T> > coped) {
        model = coped;
    }

    Matrix(int rows, int cols, T value = 0) {
        model.resize(cols);
        for (int i = 0; i < cols; ++i)
            model[i].resize(rows, value);
        column = cols;
        row = rows;
    }

    Matrix operator+(Matrix<T> other) {
        if (column != other.column || row != other.row)
            return Matrix<T>();
        Matrix<T> result(row, column);
        for (int i = 0; i < row; ++i) {
            for (int j = 0; j < column; ++j) {
                result.set(i, j, model[i][j] + other[i][j]);
            }
        }
        return result;
    }

    Matrix operator-(Matrix<T> other) {
        if (column != other.column || row != other.row)
            return Matrix<T>();
        Matrix<T> result(row, column);
        for (int i = 0; i < row; ++i) {
            for (int j = 0; j < column; ++j) {
                result.set(i, j, model[i][j] - other[i][j]);
            }
        }
        return result;
    }

    void fill(T value) {
        for (int i = 0; i < row; ++i) {
            for (int j = 0; j < column; ++j) {
                model[i][j] = value;
            }
        }
    }

    void print() {
        for (int i = 0; i < row; ++i) {
            for (int j = 0; j < column; ++j) {
                cout << model[i][j] << ' ';
            }
            cout << endl;
        }
    }

    const vector<vector<T>> clone() {
        vector<vector<T>> copy;
        for (int i = 0; i < row; ++i) {
            for (int j = 0; j < column; ++j) {
                copy[i].push_back(model[i][j]);
            }
        }
        return copy;
    }

    void set(int i, int j, T value) {
        model[i][j] = value;
    }

    vector<T> operator[](int index) {
        if (model.size() > index)
            return model[index];
        return vector<T>();
    }

    Matrix<T> operator*(Matrix<T> other) {
        if (column != other.row)
            return Matrix();
        Matrix result = Matrix(row, other.column);
        for (int i = 0; i < row; ++i)
            for (int j = 0; j < other.column; ++j)
                for (int k = 0; k < column; ++k)
                    result.model[i][j] += model[i][k] * other[k][j];
        return result;
    }

    Matrix<T> strassen(Matrix<T> a, Matrix<T> b) {
        if (a.row <= 1 || a.column <= 1 || b.row <= 1 || b.column <= 1)
            return a * b;

        Matrix<T> result = Matrix(a.row, b.column);
        Matrix<T> a11(a.row / 2, a.column / 2), a12(a.row / 2, a.column / 2),
                a21(a.row / 2, a.column / 2), a22(a.row / 2, a.column / 2);
        Matrix<T> b11(b.row / 2, b.column / 2), b12(b.row / 2, b.column / 2),
                b21(b.row / 2, b.column / 2), b22(b.row / 2, b.column / 2);

        int la = a.row / 2, ca = a.column / 2;
        int lb = b.row / 2, cb = b.column / 2;
        for (int i = 0; i < la; ++i)
            for (int j = 0; j < ca; ++j) {
                a11.set(i, j, a[i][j]);
                a12.set(i, j, a[i][j + ca]);
                a21.set(i, j, a[i + la][j]);
                a22.set(i, j, a[i + la][j + ca]);
            }
        for (int i = 0; i < lb; ++i)
            for (int j = 0; j < cb; ++j) {
                b11.set(i, j, b[i][j]);
                b12.set(i, j, b[i][j + cb]);
                b21.set(i, j, b[i + lb][j]);
                b22.set(i, j, b[i + lb][j + cb]);
            }
        Matrix<T> m1 = strassen(a11 + a22, b11 + b22);
        Matrix<T> m2 = strassen(a21 + a22, b11);
        Matrix<T> m3 = strassen(a11, b12 - b22);
        Matrix<T> m4 = strassen(a22, b21 - b11);
        Matrix<T> m5 = strassen(a11 + a12, b22);
        Matrix<T> m6 = strassen(a21 - a11, b11 + b12);
        Matrix<T> m7 = strassen(a12 - a22, b21 + b22);

        Matrix<T> c11 = m1 + m4 - m5 + m7;
        Matrix<T> c12 = m3 + m5;
        Matrix<T> c21 = m2 + m4;
        Matrix<T> c22 = m1 + m2 - m3 - m6;

        for (int i = 0; i < la; ++i)
            for (int j = 0; j < cb; ++j) {
                result.set(i, j, c11[i][j]);
                result.set(i + la, j, c21[i][j]);
                result.set(i, j + cb, c12[i][j]);
                result.set(i + la, j + cb, c22[i][j]);

            }
        return result;
    }
};

template<typename N>
Matrix<N> randomMT(int row, int col, int max = 1000) {
    srand(time(NULL));
    Matrix<N> result(row, col);
    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j) {
            result.set(i, j, sqrt(rand() % max));
        }
    }
    return result;
}


int main() {

    Matrix<double> a, b;

    a = randomMT<double>(2, 2);
    cout << "mtz base\n";
    a.print();
    cout << "mtz muilt normal\n";
    (a * a).print();
    cout << "mtz multi stra\n";
    b.strassen(a, a).print();
    return 0;
}