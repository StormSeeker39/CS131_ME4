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

	/*Default Constructor*/
	Matrix() {
		row = col = 0;
	}

	/*Creates identity matrix of size n.*/
	Matrix(int n) {
		row = col = n;
		init();
		for (int i = 0; i < n; i++) {
			val[i][i] = 1;
		}
	}

	/*Creates a matrix with r rows and c columns.*/
	Matrix(int r, int c) {
		row = r;
		col = c;
		init();
	}

	/*Creates a matrix that resembles a column vector*/
	Matrix(vector<double> v) {
		row = v.size();
		col = 1;
		init();
		for (int i = 0; i < row; i++) {
			val[i][0] = v[i];
		}
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

	Matrix transpose() {
		Matrix t(col, row);
		for (int i = 0; i < row; i++) {
			for (int j = 0; j < col; j++) {
				t.val[j][i] = val[i][j];
			}
		}
		return t;
	}

	vector<double> slice(int n) {
		vector<double> a;
		for (int i = 0; i < row; i++) {
			a.push_back(val[i][n]);
		}
		return a;
	}

	double* operator[](int i) {
		return val[i];
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

void householder(Matrix & A);
void givens(Matrix & A);
void gram_schmidt1(Matrix & A);
void gram_schmidt2(Matrix & A);
int sgn(double a);
Matrix operator*(Matrix A, Matrix B);
Matrix operator*(double A, Matrix B);
Matrix operator-(Matrix A, Matrix B);
double norm2(vector<double> v);
double dot(vector<double> v);
double dot(vector<double> v1, vector<double> v2);

int main() {
	Matrix A(2, 2);
	A[0][0] = 1;
	A[0][1] = 2;
	A[1][0] = 2;
	A[1][1] = 1;
	A.print();
	
	Matrix B(2, 2);
	B[0][0] = 1;
	B[0][1] = 2;
	B[1][0] = 2;
	B[1][1] = 1;
	B.print();

	Matrix C = A*B;
	C.print();

	Matrix D = 2 * A;
	D.print();

	Matrix E = A - B;
	E.print();

	Matrix sam(5, 3);
	double c = -1.0;
	for (int i = 0; i < sam.row; i++) {
		sam[i][0] = 1;
		sam[i][1] = c;
		sam[i][2] = pow(c, 2);
		c += 0.5;
	}
	sam.print();

	//householder(sam);
	//givens(sam);
	gram_schmidt1(sam);
	return 0;
}

void householder(Matrix & A) {
	vector < Matrix > H;
	for (int j = 0; j < A.col; j++) {
		Matrix h(A.row);
		vector<double> v(A.row, 0.0);
		double alpha = 0.0;
		for (int k = j; k < A.row; k++) {
			v[k] = A[k][j];
		}
		alpha = (-1) * sgn(v[j]) * norm2(v);
		v[j] -= alpha;

		h = h - (2 / dot(v))*(Matrix(v) * Matrix(v).transpose());
		//h.print();

		A = h*A;
		//A.print();
		H.push_back(h);
	}
	Matrix Q(A.row);
	for (int i = 0; i < H.size(); i++) {
		Q = Q * H[i];
	}
	cout << "Q:" << endl;
	Q.print();
}

void givens(Matrix & A) {
	vector<Matrix> G;
	for (int j = 0; j < A.col; j++) {
		for (int i = A.row - 1; i > j; i--) {
			vector<double> a(2, 0.0);
			a[0] = A[j][j];
			a[1] = A[i][j];
			double hyp = norm2(a);
			double c = a[0] / hyp;
			double s = a[1] / hyp;

			Matrix g(A.row);
			g[j][j] = c;
			g[j][i] = s;
			g[i][j] = -s;
			g[i][i] = c;

			g.print();
			A = g * A;
			A.print();
			G.push_back(g);
		}
	}
}

void gram_schmidt1(Matrix & A) {
	Matrix R(A.col, A.col);
	Matrix Q(A.row, A.col);
	for (int k = 0; k < A.col; k++) {
		R[k][k] = norm2(A.slice(k));
		vector<double> q(A.row);
		for (int i = 0; i < A.row; i++) {
			Q[i][k] = A[i][k] / R[k][k];
			q[i] = Q[i][k];
		}
		for (int j = k + 1; j < A.col; j++) {
			R[k][j] = dot(q, A.slice(j));
			for (int i2 = 0; i2 < A.row; i2++) {
				A[i2][j] -= R[k][j] * q[i2];
			}
		}
	}
	Q.print();
	R.print();
}

Matrix operator*(Matrix A, Matrix B) {
	Matrix m(A.row, B.col);
	if (A.col != B.row) {
		cout << "ERROR: Multiplying incompatible matrices" << endl;
		return m;
	}
	for (int i = 0; i < m.row; i++) {
		for (int j = 0; j < m.col; j++) {
			for (int k = 0; k < A.col; k++) {
				m[i][j] += A[i][k] * B[k][j];
			}
		}
	}
	return m;
}

Matrix operator*(double c, Matrix A) {
	Matrix m(A.row, A.col);
	for (int i = 0; i < m.row; i++) {
		for (int j = 0; j < m.col; j++) {
			m[i][j] = c * A[i][j];
		}
	}
	return m;
}

Matrix operator-(Matrix A, Matrix B) {
	Matrix m(A.row, A.col);
	if (A.row != B.row && A.col != B.col) {
		cout << "Subtracting incompatible matrices" << endl;
		return m;
	}
	for (int i = 0; i < m.row; i++) {
		for (int j = 0; j < m.col; j++) {
			m[i][j] = A[i][j] - B[i][j];
		}
	}
	return m;
}

int sgn(double a) {
	return (a > 0.0) - (a < 0.0);
}

double norm2(vector<double> v) {
	return sqrt(dot(v));
}

double dot(vector<double> v) {
	double dot = 0.0;
	for (int i = 0; i < v.size(); i++) {
		dot += pow(v[i], 2);
	}
	return dot;
}

double dot(vector<double> v1, vector<double> v2) {
	if (v1.size() != v2.size()) {
		cout << "ERROR: Computing dot product of incompatible vectors" << endl;
		return 0;
	}
	double dot = 0.0;
	for (int i = 0; i < v1.size(); i++) {
		dot += v1[i] * v2[i];
	}
	return dot;
}