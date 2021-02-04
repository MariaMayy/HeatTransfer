#include <iostream>
#include <stdio.h>
#include <vector>
#include <math.h>
using namespace std;

void Set(double* f, double* x, double* k, double* q, int N, int i);
vector <double> RealDecision(int N, double x0 = 0.5);
vector <double> MethodRun(int N, double* k, double* q, double* f, double* x, double B1, double M1, double B2, double M2);
double Max(vector <double> u1, vector <double> u2, int N);

void Set(double* f, double* x, double* k, double* q, int N, int i) {
	double h = 1.0 / N;

	for (int j = 0; j < N + 1; j++) {
		x[j] = j*h;

		// если модельная задача ( i = 0 ) 
		if (i == 0) {
			f[j] = sin(0.5);
			k[j] = exp(0.5);
			q[j] = exp(0.5);
		}
		else {
			f[j] = sin(x[j]);
			k[j] = exp(x[j]);
			q[j] = exp(x[j]);
		}
	}
};

vector <double> RealDecision(int N, double Beta1, double M1, double Beta2, double M2, double x0) {
	double h = 1.0 / N;
	double* x;
	x = new double[N + 1];

	int i = 0;

	vector <double> VectorY;

	for (i = 0; i < N + 1; i++) x[i] = i * h;
	double C1 = sin(0.5) / (exp(0.5) * (-exp(1.5) + exp(-0.5) - exp(1.0) - exp(-1.0)) );

	for (i = 0; i < N + 1; i++) VectorY.push_back(C1 * (exp(x[i]) + exp(-x[i])) + sin(0.5) / exp(0.5));
	
	return VectorY;
}

vector <double> RMethod(double* A, double* B, double* C, double* f, int N) {
	vector <double> VectorY;
	double* delta;
	double* alfa;
	double* beta;
	double* y;
	int i;

	y = new double[N + 2];
	delta = new double[N + 2];
	alfa = new double[N + 2];
	beta = new double[N + 2];

	delta[0] = C[0];
	beta[1] = f[0] / C[0];

	alfa[1] = B[0] / C[0];

	for (i = 1; i < N; i++) {
		delta[i] = C[i] - A[i] * alfa[i];
		alfa[i + 1] = B[i] / delta[i];

		beta[i + 1] = (f[i] + A[i] * beta[i]) / delta[i];
	}
	beta[N + 1] = (f[N] + A[N] * beta[N]) / (C[N] - A[N] * alfa[N]);

	y[N] = beta[N + 1];
	for (i = N - 1; i >= 0; i--) {
		y[i] = alfa[i + 1] * y[i + 1] + beta[i + 1];
	}

	for (i = 0; i <= N; i++) {
		VectorY.push_back(y[i]);
	}

	return VectorY;

	delete[]alfa;
	delete[]beta;
	delete[]delta;
}

void SetMatrix(int N, double* k, double* q, double* f, double* x, double B1, double M1, double B2, double M2, double* DiagElemC, double* UpDiagElemB, double* UnderDiagElemA) {
	double h = 1.0 / N;

	UpDiagElemB[0] = k[0] / h;
	DiagElemC[0] = UpDiagElemB[0] + B1;

	f[0] = M1;

	for (int i = 1; i < N; i++) {
		UnderDiagElemA[i] = k[i] / (h * h);
		UpDiagElemB[i] = k[i + 1] / (h * h);
		DiagElemC[i] = UnderDiagElemA[i] + UpDiagElemB[i] + q[i];
	}

	UnderDiagElemA[N] = -k[N] / h;
	DiagElemC[N] = UnderDiagElemA[N] - B2;
	f[N] = -M2;

}


double Max(vector <double> u1, vector <double> u2, int N) {
	double MaxCur = 0;

	for (int i = 0; i < N; i++) {
		if (MaxCur < abs(u1[i] - u2[i])) {
			MaxCur = abs(u1[i] - u2[i]);
		}
	}
	return MaxCur;
}


int main() {
	int N, NPrev;
	int iMode; // 0 - if model task, 1 - if real task
	int iEnd; // end of task

	double* q, * k, * x, * f;
	double x0 = 0.5;
	vector <double> u1; // analytical solution
	vector <double> u2; // numerical solution

	double* A;
	double* B;
	double* C;



	cout << "Enter task mode ( 0 - model task, 1 - real task): ";
	cin >> iMode;

	cout << "Enter Start N: ";
	cin >> N;

	double B1 = 0;
	double B2 = 1;
	double M1 = 0;
	double M2 = 0;
	double Eps = 0.01;
	double dMax = 0; // current error


	if (iMode) {
		// real task
		do {
			dMax = 0;
			// allocate memory for array
			x = new double[N + 1];
			q = new double[N + 1];
			k = new double[N + 1];
			f = new double[N + 1];

			Set(f, x, k, q, N, 1);

			A = new double[N + 1];
			B = new double[N];
			C = new double[N + 1];

			SetMatrix(N, k, q, f, x, B1, M1, B2, M2, C, B, A);
			u1 = RMethod(A, B, C, f, N); // numerical solution for N

			cout << "N = " << N << '\n';

			delete[] x;
			delete[] k;
			delete[] q;
			delete[] f;
			delete[] A;
			delete[] B;
			delete[] C;

			NPrev = N;
			N = 2 * N;
			// allocate memory for array
			x = new double[N + 1];
			q = new double[N + 1];
			k = new double[N + 1];
			f = new double[N + 1];

			Set(f, x, k, q, N, 1);

			A = new double[N + 1];
			B = new double[N];
			C = new double[N + 1];

			SetMatrix(N, k, q, f, x, B1, M1, B2, M2, C, B, A);
			u2 = RMethod(A, B, C, f, N); // numerical solution for 2N

			for (int i = 0; i < NPrev + 1; i++) {
				cout << " y_N[" << i << "] = " << u1[i] << " y_2N[" << 2*i << "] = " << u2[2*i] << " delta = " << abs(u2[2*i] - u1[i]) << '\n';
				if (dMax < abs(u2[2 * i] - u1[i])) dMax = abs(u2[2 * i] - u1[i]);
			};
			iEnd = dMax < Eps;

			if (iEnd) {
				cout << "Max delta = " << dMax << " < " << Eps << '\n';
			}
			else
			{
				cout << "Max delta = " << dMax << " > " << Eps << '\n';
			}
				delete[] x;
				delete[] k;
				delete[] q;
				delete[] f;
				delete[] A;
				delete[] B;
				delete[] C;
				u1.clear();
				u2.clear();
			


		} while ( !iEnd );


	}
	else {
		// model task
		do {
			u1 = RealDecision(N, 0, 1, 0, 0, x0); // analytical solution

// allocate memory for array
			x = new double[N + 1];
			q = new double[N + 1];
			k = new double[N + 1];
			f = new double[N + 1];

			Set(f, x, k, q, N, 0);

			A = new double[N + 1];
			B = new double[N];
			C = new double[N + 1];

			SetMatrix(N, k, q, f, x, B1, M1, B2, M2, C, B, A);
			u2 = RMethod(A, B, C, f, N); // numerical solution

			cout << "N = " << N << '\n';
			for (int i = 0; i < N + 1; i++) {
				cout << " u[" << i << "] = " << u1[i] << " y[" << i << "] = " << u2[i] << " delta = " << abs(u2[i] - u1[i]) << '\n';
			};

			iEnd = Max(u1, u2, N) < Eps;

			if (iEnd) {
				cout << "Max delta = " << Max(u1, u2, N) << " < " << Eps << '\n';
			}
			else
			{
				cout << "Max delta = " << Max(u1, u2, N) << " > " << Eps << '\n';
				N = 2 * N;
			}
				delete[] x; 
				delete[] k;
				delete[] q;
				delete[] f;
				delete[] A;
				delete[] B;
				delete[] C;
				u1.clear();
				u2.clear();
				
			

		} while ( !iEnd );

	}


	return 0;
}

