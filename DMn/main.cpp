#include "main.h"

real* g;
real* dg;
real** A;
real* b;
int N; // число приемников
real gamma;
int nx, ny, nz; // число ячеек по коордам

typedef real** Matrix;
typedef real* Vector;

struct origin
{
   real x[2], y[2], z[2];
} Orig;

int main()
{
   //reciever Recs;
   //cout << Integrate(-2, 3, 1, 2, 0, 2, [](real x, real y, real z) {return pow(x*y, 5.) + x*y*z*z*z; });
   return 0;
}

// Чтение координат приемников
// Возвращает координаты по Х и по Y, а также количество координат по каждой из этих осей.
void ReadRec(real* RecX, real* RecY, int& nRecX, int& nRecY)
{
   ifstream in("Params.txt");

   real xL, xR, yL, yR, xH, yH;

   in >> xL >> xR >> xH; // Вввод параметров (Левая граница оси, Правая граница оси, Шаг по оси)
   in >> yL >> yR >> yH;

   nRecX = (int)((xR - xL) / xH);
   nRecY = (int)((yR - yL) / yH);

   RecX = new real[nRecX]{};
   RecY = new real[nRecY]{};

   for (int yI = 0; yI < nRecY; yI++)
      RecY[yI] = yL + yI * yH;

   for (int xI = 0; xI < nRecX; xI++)
      RecX[xI] = xL + xI * xH;
}

void SLAEgen(int nX, int nY, int nZ, int nRecX, int nRecY, real* X, real* Y, real* Z, real* RecX, real* RecY)
{
   int K = nX * nY * nZ;   // Количество ячеек
   int nRec = nRecX * nRecY;   // Количество приемников
   Matrix A = new real * [K] {};   // Глобальная матрица СЛАУ
   for (int i = 0; i < K; i++)
      A[i] = new real[K]{};
   real* b = new real[K]{}; // Глобальная правая часть СЛАУ

   function<real(real, real, real, real*, int)> f = [](real xB, real yB, real zB, real* args, int argNum) // подынтегральная функция
   {
      real r = sqrt(pow(args[0] - xB, 2) + pow(args[1] - yB, 2) + pow(args[2] - zB, 2));   // Подсчет нормы вектора
      return zB / pow(r, 3);
   };

   real RecCoords[3]{};
   real gq, gs;
   int iX, iY, iZ, jX, jY, jZ;
   int nXY = nX * nY;

   for (int iRecX = 0; iRecX < nRecX; iRecX++)
   {
      RecCoords[0] = RecX[iRecX];
      for (int iRecY = 0; iRecY < nRecY; iRecY++)
      {
         RecCoords[1] = RecY[iRecY];
         for (int q = 0; q < K; q++)
         {
            iZ = q / nXY;
            iY = (q % nXY) / nX;
            iX = (q % nXY) % nX;
            gq = Integrate(X[iX], X[iX + 1], Y[iY], Y[iY + 1], Z[iZ], Z[iZ + 1], f, RecCoords, 3);
            for (int s = 0; s < K; s++)
            {
               jZ = s / nXY;
               jY = (s % nXY) / nX;
               jX = (s % nXY) % nX;
               gs = Integrate(X[jX], X[jX + 1], Y[jY], Y[jY + 1], Z[jZ], Z[jZ + 1], f, RecCoords, 3);
               A[q][s] += gq * gs;
            }
            b[q] += gq * g[q];

         }
      }
   }

   for (int q = 0; q < K; q++)
   {
      b[q] /= 4 * M_PI;
      for (int s = 0; s < K; s++)
      {
         A[q][s] /= 4 * M_PI * 4 * M_PI;
      }
   }

}

void FindEdgesAndCenter()
{
   real xL, xR;
   int nX;
   ifstream in("Mesh.txt");
   in >> xL >> xR >> nX;

   real* xMesh = new real[nX]{};

   real hX = (xR - xL) / nX;
   int i = 0;
   for (real x = xL ; x < xR; x += hX, i++)
   {
      xMesh[i] = x;
   }

   real yL, yR;
   int nY;
   in >> yL >> yR >> nY;
   real* yMesh = new real[nY]{};

   real hY = (yR - yL) / nY;

   i = 0;
   for (real y = yL; y < yR; y += hY, i++)
   {
      yMesh[i] = y;
   }
   real zL, zR;
   int nZ;
   in >> zL >> zR >> nZ;
   real* zMesh = new real[nZ]{};

   real hZ = (zR - zL) / nZ;
   i = 0;
   for (real z = zL; z < zR; z += hZ, i++)
   {
      zMesh[i] = z;
   }



}

void SLAEsolve(double** M, int N, double* b, double* q)
{
   int i, j;

   /*GAUSSE---------------------------------------*/
   double a;
   int iK, jK;
   for (i = 0; i < N; i++)
   {
      a = M[i][i];
      if (a)
      {
         for (j = i; j < N; j++) // Делим текущую строку на коэффициент
         {
            M[i][j] /= a;
         }
         b[i] /= a;

         for (iK = i + 1; iK < N; iK++) // Вычитаем текущую строку из последующих
         {
            a = M[iK][i];
            for (jK = i + 1; jK < N; jK++)
            {
               M[iK][jK] -= a * M[iK][jK];
            }
            b[iK] -= a * b[i];
         }
      }
      else
         exit(-1);
   }

   for (i = N - 1; i > -1; i--)
   {
      for (j = N - 1; j > i; j--)
         b[i] -= M[i][j] * q[j];

      q[i] = b[i] / M[i][i];
   }
}

real Integrate(real xL, real xR, real yL, real yR, real zL, real zR, function<real(real, real, real, real* args, int argNum)> f, real* args, int argNum)
{
   real result = 0.;
   //const int segs = 5;

   real xj[5] = { -sqrt(5. + 2. * (sqrt(10. / 7.))) / 3., -sqrt(5. - 2. * (sqrt(10. / 7.))) / 3. ,
              0. , sqrt(5. - 2. * (sqrt(10. / 7.))) / 3. , sqrt(5. + 2. * (sqrt(10. / 7.))) / 3. };
   real qj[5] = { (322. - 13. * sqrt(70.)) / 900., (322. + 13. * sqrt(70.)) / 900., 128. / 225.,
                  (322. + 13. * sqrt(70.)) / 900., (322. - 13. * sqrt(70.)) / 900. };
   real hX = (xR - xL) / 2.,
      hY = (yR - yL) / 2.,
      hZ = (zR - zL) / 2.,
      cX = (xR + xL) / 2.,
      cY = (yR + yL) / 2.,
      cZ = (zR + zL) / 2.;


   for (int i = 0; i < 5; i++)               // x
      for (int j = 0; j < 5; j++)            // y 
         for (int k = 0; k < 5; k++)         // z
            result += qj[i] * qj[j] * qj[k] * (f(cX + xj[i] * hX, cY + xj[j] * hY, cZ + xj[k] * hZ, args, argNum));

   return (xR - xL) * (yR - yL) * (zR - zL) * result / 8.;
}