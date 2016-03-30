#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#define D 4.4
#define V -1.2
#define H -1
#define eps 0.001
using namespace std;
double f(double x, double y) {
	return (0.2*exp(x)*cos(y));
}
double phi(double x, double y) {
	return (exp(x)*cos(y));
}
template <typename T>
void writeMatr(T ** arr, int n) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++)
			cout<<arr[i][j]<<" ";
		printf("\n");
	}
	printf("\n");
}
void writeVec(double * Vec, int n) {
	for (int i = 0; i < n; i++) {
		printf("%lf ", Vec[i]);
	}
	printf("\n");
}
void createMatr(double **arr, int n,double h) {
	int sred = (n - 1) / 2;
	for (int i = 0; i < n; i++) {
		arr[n-1][i] = phi(i*h,0);
		arr[n-1-i][n-1] = phi(1, i*h);
		if (i <= sred) {
			arr[sred][i] = phi(i*h, 0.5);
			arr[n - i - 1][0] = phi(0, i*h);
		}
		else
		{
			arr[n-i-1][i] = phi(i*h, i*h);
		}
	}
}
int QueryMatr(int **Ind,double ** arr,int n) {
	int t = 1;
	for (int i = n - 2; i > 0; i--) {
		for (int j = 1; j < n - 1; j++) {
			if ((i <= (n - 1) / 2) && (j < n - i))
				continue;
			else
			{
				Ind[i][j] = t++;
			}
		}
	}
	return (t-1);
}
void returnToMatr(double **arr, double *x, int n) {
	int t = 0;
	for (int i = n - 2; i > 0; i--) {
		for (int j = 1; j < n - 1; j++) {
			if ((i <= (n - 1) / 2) && (j < n - i))
				continue;
			else
			{
				arr[i][j] = x[t++];
			}
		}
	}
}
void copy(double *d, double *s, int n)
{
	for (int i = 0; i < n; i++)
		d[i] = s[i];
}
double norm(double *r, int n)
{
	double rezult = 0;
	for (int i = 0; i < n; i++)
		rezult = rezult + r[i] * r[i];
	rezult = sqrt(rezult);
	return (rezult);
}
double matric(double **arr, double *z, int n, int i)
{
	double rezult = 0;
	for (int j = 0; j < n; j++)
		rezult = rezult + arr[i][j] * z[j];
	return rezult;
}
double alfa(double *r, double **arr, double *z, int n)
{
	double rezult = 0, rezult2 = 0;
	for (int i = 0; i < n; i++)
		rezult = rezult + r[i] * r[i];
	for (int i = 0; i < n; i++)
	{
		rezult2 = rezult2 + matric(arr, z, n, i)*z[i];
	}
	return (rezult / rezult2);
}
double beta(double *r, double *tempr, int n)
{
	double rezult = 0, rezult2 = 0;
	for (int i = 0; i < n; i++) {
		rezult = rezult + r[i] * r[i];
		rezult2 = rezult2 + tempr[i] * tempr[i];
	}
	return (rezult / rezult2);
}
void fillMatr(int ** Ind, double **A, double **arr, double h, int n, int m,double *F) {
	for (int i = n - 2; i > 0; i--) {
		for (int j = 1; j < n - 1; j++) {
			if ((i <= (n - 1) / 2) && (j < n - i))
				continue;
			else
			{
				A[Ind[i][j]-1][Ind[i][j]-1] = D / (h*h);
				F[Ind[i][j]-1] = f(i*h, j*h);
				if (Ind[i + 1][j] > 0)
					A[Ind[i][j] - 1][Ind[i + 1][j] - 1] = V/(h*h);
				else 
					F[Ind[i][j]-1] += (- V / (h*h))*arr[i+1][j];
				if (Ind[i - 1][j]>0)
					A[Ind[i][j] - 1][Ind[i - 1][j] - 1] = V / (h*h);
				else
					F[Ind[i][j]-1] += (-V / (h*h))*arr[i - 1][j];
				if (Ind[i][j+1]>0)
					A[Ind[i][j] - 1][Ind[i][j+1] - 1] = H / (h*h);
				else
					F[Ind[i][j] - 1] += (-H / (h*h))*arr[i][j+1];
				if (Ind[i][j - 1]>0)
					A[Ind[i][j] - 1][Ind[i][j - 1] - 1] = H / (h*h);
				else
					F[Ind[i][j] - 1] += (-H / (h*h))*arr[i][j - 1];
			}
		}
	}
}
double * sopryazh(double *F, double **A, int n) {
	double * x = (double*)calloc(sizeof(double), n);
	double * r = (double*)calloc(sizeof(double), n);
	double * z = (double*)calloc(sizeof(double), n);
	double * tempr = (double*)calloc(sizeof(double), n);
	double * tempz = (double*)calloc(sizeof(double), n);
	double a, b;
	copy(r, F, n);
	copy(z, r, n);
	while (norm(r, n) > eps)
	{
		copy(tempr, r, n);
		copy(tempz, z, n);
		a = alfa(r, A, z, n);
		for (int i = 0; i < n; i++)
		{
			x[i] = x[i] + a*z[i];
			r[i] = r[i] - a*matric(A, z, n, i);
		}
		b = beta(r, tempr, n);
		for (int i = 0; i < n; i++)
			z[i] = r[i] + b*tempz[i];
	}
	free(z); free(r); free(tempr); free(tempz);
	return x;
}
void writeToFile(double ** arr, int n,double h) {
	FILE * f = fopen("output.txt", "w");
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++)
			fprintf(f, "%lf %lf %lf\n", i*h, j*h, arr[i][j]);
		fprintf(f, "\n");
	}
}
void main() {
	double h, **arr, **A,*F,*x;
	int m,n,**Ind;
	printf("ENTER ODD NUMBER: ");
	scanf("%d", &n);
	h = 1.0 / (n-1);
	F = (double*)calloc(sizeof(double), n);
	arr = (double**)calloc(sizeof(double*), n);
	for (int i = 0; i < n; i++) {
		arr[i] = (double*)calloc(sizeof(double), n);
	}
	Ind = (int**)calloc(sizeof(int*), n);
	for (int i = 0; i < n; i++) {
		Ind[i] = (int*)calloc(sizeof(int), n);
	}
	createMatr(arr, n, h);
	m = QueryMatr(Ind,arr, n);
	A = (double**)calloc(sizeof(double*),m);
	for (int i = 0; i < m; i++) {
		A[i] = (double*)calloc(sizeof(double), m);
	}
	fillMatr(Ind,A, arr, h,n, m,F);
	x = sopryazh(F,A, m);
	writeMatr<double>(arr, n);
	returnToMatr(arr, x, n);
	writeMatr<double>(arr, n);
	writeToFile(arr, n,h);
}