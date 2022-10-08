#include <iostream>
#include <vector>
#include <stdlib.h>
#include <cmath>
#include <stdbool.h>

#define ABS(x) ((x)>0?(x):-(x))
#define MAX(x,y) ((x)>(y)?(x):(y)
#define PI 3.1415926

//кусочно-линейные базисные функции
class Equation {
public:

  double* xk; //сетка
  int n; //размер сетки
  double a, b;
  double* ak; //ответ
  //phi(x) = sum from i=0 to n ak*phi_k(x)

  Equation(double a, double b, int size) {
    this->a = a;
    this->b = b;
    n = size;
    xk = new double[n + 1];
    double h = (b - a) / n;
    for (int i = 0; i <= n; ++i)
      xk[i] = a + h * i;

    ak = nullptr;
  }

  ~Equation() {
    delete[] xk;
    delete[] ak;
  }

  //коэффициенты СЛАУ метода Ритца

   //a_ki = integral from 0 to 1 phi_k(x)' * phi_i(x)' dx
  double A(int k, int i) {

    if (k == 0 || i == 0 || k == n || i == n)
      return 0.0;

    if (i == k - 1) {
      return -1.0 / (xk[k] - xk[k - 1]);
    }
    else if (i == k) {
      return 1.0 / (xk[k] - xk[k - 1]) + 1.0 / (xk[k + 1] - xk[k]);
    }
    else if (i == k + 1) {
      return -1.0 / (xk[k + 1] - xk[k]);
    }
    else
      return 0.0;

  }

  //b_k = integral from 0 to 1 (f * phi_k(x)) dx
  double B(int k) {

    if (k == 0 || k == n) {
      return 0.0;
    }

    //f(x) = 2
    /*
    else
      return 2.0 * (1.0 / (xk[k] - xk[k - 1]) * (0.5 * xk[k] * xk[k] + 0.5 * xk[k - 1] * xk[k - 1] - xk[k - 1] * xk[k])\
        + 1.0 / (xk[k + 1] - xk[k]) * (0.5 * xk[k + 1] * xk[k + 1] + 0.5 * xk[k] * xk[k] - xk[k] * xk[k + 1]));
    */

    //f(x) = x*sin(PI*x)
    else {
      double a = xk[k - 1], b = xk[k], c = xk[k + 1];
      return (PI * (a - 2 * b) * sin(PI * b) + (PI * PI * b * (b - a) - 2) * cos(PI * b) + PI * a * sin(PI * a) + 2 * cos(PI * a)) / (PI * PI * PI * (a - b)) + \
        (PI * (c - 2 * b) * sin(PI * b) + (PI * PI * b * (b - c) - 2) * cos(PI * b) + PI * c * sin(PI * c) + 2 * cos(PI * c)) / (PI * PI * PI * (a - b));
    }
  }

  double p(double x) {
    return -1.0;
  }

  double q(double x) {
    return 0.0;
  }
  double r(double x) {
    return 0.0;
  }

  double f(double x) {
    //-u'' = 2
    //return 2.0;

    //-u'' = x*sin(PI*x)
    return x * sin(PI * x);
  }

  double phi(int i, double x) {
    if (i == 0 || i == n)
      return 0.0;

    if (xk[i - 1] <= x && x >= xk[i])
      return (x - xk[i - 1]) / (xk[i] - xk[i - 1]);
    else if (xk[i] < x && x >= xk[i + 1])
      return (xk[i + 1] - x) / (xk[i + 1] - xk[i]);
    else
      return 0.0;
  }

  double dphidx(int i, double x) {
    if (i == 0 || i == n)
      return 0.0;

    if (xk[i - 1] <= x && x >= xk[i])
      return 1.0 / (xk[i] - xk[i - 1]);
    else if (xk[i] < x && x >= xk[i + 1])
      return -1.0 / (xk[i + 1] - xk[i]);
    else
      return 0.0;
  }

  //solution
  double u(double x) {
    double ans = 0.0;
    for (int i = 1; i <= n - 1; ++i)
      ans += ak[i] * phi(i, x);
    return ans;
  }

  //first derivative of solution
  double dudx(double x) {
    double ans = 0.0;
    for (int i = 1; i <= n - 1; ++i)
      ans += ak[i] * dphidx(i, x);
    return ans;
  }

  double exact(double x) {
    //-u'' = 2
    //return x * (1.0 - x);

    //-u'' = x*sin(PI*x)
    return (4.0 * x + PI * x * sin(PI * x) + 2.0 * cos(PI * x) - 2.0) / (PI * PI * PI);
  }

  double dexactdx(double x) {
    //-u'' = 2
    //return -2.0 * x + 1.0;

    //-u'' = x*sin(PI*x)
    return  (-PI * sin(PI * x) + PI * PI * x * cos(PI * x) + 4.0) / (PI * PI * PI);
  }

  //значение функционала в точке 
  double J(double* a) {
    //phi(x) = sum from i=0 to n a_i * phi_i(x)
    //J(phi) = integral from 0 to 1 (|phi|')^2dx - integral from 0 to 1 (f * phi(x))dx
    double ans = 0.0;

    for (int i = 1; i <= n - 1; ++i) {
      for (int k = 1; k <= n - 1; ++k)
        ans += (a[i] * a[k] * A(i, k));
    }

    for (int i = 1; i <= n - 1; ++i)
      ans -= (a[i] * B(i));

    return ans;
  }

  double L2norm() {
    //right traingles

    int i;
    double n = 1.0, h;
    double sum = 10000000.0, prevsum = 10000000.0;
    double x;
    double f;

    do {
      prevsum = sum;
      sum = 0.0;
      n *= 2.0;
      h = (b - a) / n;
      x = a + h;

      for (i = 0; i < n; i++) {
        //|| u - u*||
        //f = ABS(u(x) - exact(x))*ABS(u(x) - exact(x));

        //||u' - u*'||
        f = ABS(dudx(x) - dexactdx(x)) * ABS(dudx(x) - dexactdx(x));

        sum += f;
        x += h;
      }
      sum *= h;

    } while (ABS((sum - prevsum))/ABS(sum) > 0.05);

    return sqrt(sum);
  }
};


void FiniteDifferencesMethod(Equation* eq, double a, double b, double y0, double yn) {
  int n = eq->n;
  double h, x;
  int k;
  double* F, * C, * D, * B, * alpha, * beta, * y;
  double pk, qk;

  y = new double[n + 1];
  F = new double[n + 1];
  C = new double[n + 1];
  D = new double[n + 1];
  B = new double[n + 1];
  alpha = new double[n + 1];
  beta = new double[n + 1];

  //заполняем трехдиагональную матрицу коэффициентов системы
  B[0] = 0;
  C[0] = 1;
  D[0] = 0;
  F[0] = y0;
  x = a;
  h = (b - a) / n;
  for (k = 1; k < n; k++) {
    x += h;
    pk = eq->p(x);
    qk = eq->q(x);
    B[k] = pk - h * qk / 2.0;
    C[k] = eq->r(x) * h * h - 2.0 * pk;
    D[k] = pk + h * qk / 2.0;
    F[k] = eq->f(x) * h * h;
  }
  B[n] = 0;
  C[n] = 1;
  D[n] = 0;
  F[n] = yn;

  //Метод прогонки

  //прямой ход
  beta[0] = -D[0] / C[0];
  alpha[0] = F[0] / C[0];
  for (k = 1; k < n; k++) {
    beta[k] = -D[k] / (B[k] * beta[k - 1] + C[k]);
    alpha[k] = (F[k] - B[k] * alpha[k - 1]) / (B[k] * beta[k - 1] + C[k]);
  }
  beta[n] = 0;
  alpha[n] = (F[n] - B[n] * alpha[n - 1]) / (B[n] * beta[n - 1] + C[n]);

  //обратный ход
  y[n] = yn;
  for (k = n - 1; k > 0; k--)
    y[k] = beta[k] * y[k + 1] + alpha[k];
  y[0] = y0;

  delete[] F;
  delete[] C;
  delete[] D;
  delete[] B;
  delete[] alpha;
  delete[] beta;

  eq->ak = y;
}

void GalerkinRitzMethod(Equation* eq) {
  double h = (eq->b - eq->a) / eq->n;
  int n = eq->n;
  eq->ak = new double[n + 1];
  // -u'' = f, u(a) = u(b) = 0
  //u_ = sum y_i*phi_i(x), i=1,n
  //phi_i(x) - базисные функции
  //СЛАУ вида sum a_ki * y_i = b_k
  //a_ki = integral from 0 to 1 |phi_k(x)' * phi_i(x)'|dx
  //b_k = integral from 0 to 1 f*phi_k(x) dx
  // 
  //Метод прогонки
  double* R = new double[n + 1];

  double* B = new double[n + 1];
  double* C = new double[n + 1];
  double* D = new double[n + 1];

  for (int j = 0; j <= n; ++j) {
    D[j] = eq->A(j + 1, j);
    C[j] = eq->A(j, j);
    B[j] = eq->A(j - 1, j);
    R[j] = eq->B(j);
  }
  B[0] = 0.0;
  D[n] = 0.0;

  double* d = new double[n + 1], * l = new double[n + 1];

  //Метод прогонки
  d[1] = -D[1] / C[1];
  l[1] = R[1] / C[1];
  for (int j = 2; j < n; ++j) {
    d[j] = -D[j] / (C[j] + B[j] * d[j - 1]);
    l[j] = (R[j] - B[j] * l[j - 1]) / (C[j] + B[j] * d[j - 1]);
  }
  d[n] = 0.0;

  eq->ak[0] = 0.0;
  eq->ak[n] = 0.0;

  for (int j = n - 1; j > 0; --j) {
    eq->ak[j] = d[j] * eq->ak[j + 1] + l[j];
  }

  delete[] B;
  delete[] C;
  delete[] D;
  delete[] R;
  delete[] d;
  delete[] l;
}

void RitzRelaxationMethod(Equation* eq, double* arr) {
  //J(v) = (integral from 0 to 1 (|u'|^2)dx - integral from 0 to 1 (f*u)dx) -> min
  // 
  //u(x) = sum from i=0 to n a_i*phi_i(x)

  if (arr == nullptr) {
    arr = new double[eq->n + 1];
    for (int i = 0; i <= eq->n; ++i)
      arr[i] = 2.0;

    arr[0] = 0.0;
    arr[eq->n] = 0.0;
  }

  double eps = 1.0;
  double* arr1 = new double[eq->n + 1], * arr2 = new double[eq->n + 1], * arr3 = new double[eq->n + 1];

  int N = 3;
  double y1, y2, y3;
  double** matrix = new double* [N];
  for (int j = 0; j < N; ++j)
    matrix[j] = new double[N];

  double* b = new double[N];

  for (int i = 1; i <= eq->n - 1; ++i) {

    for (int j = 0; j <= eq->n; ++j) {
      arr1[j] = arr[j];
      arr2[j] = arr[j];
      arr3[j] = arr[j];
    }

    arr1[i] -= 1;
    arr3[i] += 1;

    y1 = eq->J(arr1);
    y2 = eq->J(arr2);
    y3 = eq->J(arr3);

    //TO DO: find parabola equation using 3 dots and find its point of minimum
  }

  if (eq->ak != nullptr)
    delete[] eq->ak;

  eq->ak = arr;
  delete[] arr1;
  delete[] arr2;
  delete[] arr3;
  delete[] b;

  for (int j = 0; j < N; ++j)
    delete[] matrix[j];
  delete[] matrix;

}

int main(void) {
  double* arr = nullptr;
  int n = 2;
  printf("h\n");

  // ||u-u*|| = O(h^2)
  n = 2;
  for (int i = 0; i < 15; ++i, n *= 2) {
    printf("%lf ", ABS(log(1.0 / n)));
  }
  printf("\n\nnorm\n");

  n = 2;
  for (int i = 0; i < 15; ++i, n *= 2) {
    Equation* eq = new Equation(0.0, 1.0, n);
    GalerkinRitzMethod(eq);
    printf("%lf ", ABS(log(eq->L2norm())));
    delete eq;
  }
  

  //FiniteDifferencesMethod(&eq, a, b, A, B);
  //GalerkinRitzMethod(&eq);

  return 0;
}