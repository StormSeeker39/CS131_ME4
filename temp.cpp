#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <vector>
#include <cmath>
#define MAX_SIZE = 5;

using namespace std;

class Matrix {
public:
	int size;
	int col;
	int row;
	double val[5][5];

	Matrix() {
		row = col = 0;
	}

	Matrix(int r, int c) {
		row = r;
		col = c;
		init();
	}

	void print() {
		for (int i = 0; i < row; i++) {
			for (int j = 0; j < col; j++) {
				printf("%.7lf\t", val[i][j]);
			}
			cout << endl;
		}
		cout << endl;
	}

private:

	void init() {
		for (int i = 0; i < row; i++) {
			for (int j = 0; j < col; j++) {
				val[i][j] = 0.0;
			}
		}
	}
};

void householder(Matrix A);
int sgn(double a);
Matrix operator*(Matrix A, Matrix B);
Matrix operator*(int A, Matrix B);
Matrix operator-(Matrix A, Matrix B);

int main() {
	Matrix A(2, 2);
	A.val[0][0] = 1;
	A.val[0][1] = 2;
	A.val[1][0] = 2;
	A.val[1][1] = 1;
	A.print();

	Matrix B(2, 2);
	B.val[0][0] = 1;
	B.val[0][1] = 2;
	B.val[1][0] = 2;
	B.val[1][1] = 1;
	B.print();

	Matrix C = A*B;
	C.print();
	
	Matrix D = 2*A;
	D.print();
	
	Matrix E = A - B;
	E.print();
	return 0;
}

void householder(Matrix A) {
	vector < Matrix > H(A.col);
	for (int j = 0; j < A.col; j++) {
		Matrix h(A.row, 1);

		vector<double> v(A.row, 0.0);
		double alpha = 0.0;
		for (int k = j; j < A.row; k++) {
			v[k] = A.val[k][j];
			alpha += pow(v[k], 2);
		}
		alpha = (-1)*sgn(v[j])*sqrt(alpha);
		v[j] -= alpha;

	}
}

Matrix operator*(Matrix A, Matrix B) {
	Matrix m(A.row, B.col);
	if (A.col != B.row) {
		return m;
	}
	for (int i = 0; i < m.row; i++) {
		for (int j = 0; j < m.col; j++) {
			for (int k = 0; k < A.col; k++) {
				m.val[i][j] += A.val[i][k] * B.val[k][j];
			}
		}
	}
	return m;
}

Matrix operator*(int c, Matrix A) {
	Matrix m(A.row, A.col);
	for (int i = 0; i < m.row; i++) {
		for (int j = 0; j < m.col; j++) {
			m.val[i][j] = c * A.val[i][j];
		}
	}
	return m;
}

Matrix operator-(Matrix A, Matrix B) {
	Matrix m(A.row, A.col);
	if (A.row != B.row && A.col != B.col) {
		cout << "Subtracting mismatched matrices" << endl;
		return m;
	}
	for (int i = 0; i < m.row; i++) {
		for (int j = 0; j < m.col; j++) {
			m.val[i][j] = A.val[i][j] - B.val[i][j];
		}
	}
	return m;
}

int sgn(double a) {
	return (a > 0.0) - (a < 0.0);
}