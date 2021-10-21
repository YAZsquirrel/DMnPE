#include "main.h"

/*
real* g;
real* dg;
real** A;
real* b;
int N; // число приемников
real gamma;
int nx, ny, nz; // число ячеек по коордам
*/
typedef real** Matrix;
typedef real* Vector;

struct origin
{
   real x[2], y[2], z[2];
} Orig;

void Syntetics(real* RecX, real* RecY, int nX, int nY, int nZ, int& nRecX, int& nRecY, real* X, real* Y, real* Z)
{
    ofstream in("RecData.txt");
    
    real mes, r[3], sumG = 0.0;
    
    in << nRecX * nRecY << endl;

    auto G = [](real ro, real r, real coord, real mes) {return mes * ro / (4.0 * M_PI * r * r * r) * coord; };
    
    int i;
    int i;
    for (i = 0; i < nX && X[i] < figureCoordL; i++);

    sumG = 0;
    for (int ix = xL; ix < xR; ix++)
    {
        r[0] = X[ix] - RecX[recI];
        for (int iy = yL; iy < yR; iy++)
        {
            r[1] = Y[iy] - RecY[recJ];
            for (int iz = zL; iz < zR; iz++)
            {
                r[2] = Z[iz];
                mes = (X[ix + 1] - X[ix]) * (Y[iy + 1] - Y[iy]) * (Z[iz + 1] - Z[iz]);
                sumG += G(ro[iCell], sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]), Z[iz], mes);
            }
        }
    }
    in << sumG << " ";


    for (int recI = 0; recI < nRecX; recI++)
    {
        for (int recJ = 0; recJ < nRecY; recJ++)
        {
        }
    }
}

int readSyntetics(real** G, int nRecX, int nRecY)
{
    int nG;
    ifstream in("RecData.txt");
    in >> nG;
    if (nG != nRecX * nRecY)
        return -1;
    *G = new real[nG]{};

    real* localG = *G;

    for (int i = 0; i < nG; i++) in >> localG[i];

}

int main()
{

   //reciever Recs;
   //cout << Integrate(-2, 3, 1, 2, 0, 2, [](real x, real y, real z) {return pow(x*y, 5.) + x*y*z*z*z; });
   return 0;
}

int WriteRecMesh(real** recCoord, real left, real right, real step)
{
    int n = (int)((right - left) / step);
    real* temp = *recCoord = new real[n]{};

    for (int i = 0; i < n; i++) temp[i] = left + i * step;

    return n;
}

// Чтение координат приемников
// Возвращает координаты по Х и по Y, а также количество координат по каждой из этих осей.
void ReadRec(real** recX, real** recY, int& nRecX, int& nRecY)
{
    ifstream in("Params.txt");

    real xL, xR, yL, yR, xH, yH;

    in >> xL >> xR >> xH; // Вввод параметров (Левая граница оси, Правая граница оси, Шаг по оси)
    in >> yL >> yR >> yH;

    nRecX = WriteRecMesh(recX, xL, xR, xH);
    nRecY = WriteRecMesh(recY, yL, yR, yH);
}

void SLAEgen(int nX, int nY, int nZ, int nRecX, int nRecY, real* X, real* Y, real* Z, real* RecX, real* RecY, real *G)
{
   int K = nX * nY * nZ;   // Количество ячеек
   
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

            b[q] += gq * G[q];
         }
      }
   }

   for (int q = 0; q < K; q++)
   {
      b[q] /= 4 * M_PI;
      for (int s = 0; s < K; s++)  A[q][s] /= 4 * M_PI * 4 * M_PI;
   }

}

void ReadCoordAxis(real** coordAxis, real left, real right, int n)
{
    real* axis = *coordAxis = new real[n]{};

    real h = (left - right) / n;
    int i = 0;
    for (real coord = left; coord < right; coord += h, i++)  axis[i] = coord;
    if (i == 0) cout << "Error 1: Left coord > right coord";
}

void ReadCoords(real**X, real** Y, real**Z, int& nX, int nY, int nZ)
{
   ifstream in("Mesh.txt");

   real xL, xR, yL, yR, zL, zR;
   in >> xL >> xR >> nX;
   in >> yL >> yR >> nY;
   in >> zL >> zR >> nZ;

   ReadCoordAxis(X, xL, xR, nX);
   ReadCoordAxis(Y, yL, yR, nY);
   ReadCoordAxis(Z, zL, zR, nZ);
}

void SLAEsolve(double** M, int N, double* b, double* q)
{
   //GAUSSE
   double a;
   
   for (int i = 0; i < N; i++)
   {
      a = M[i][i];
      if (a)
      {
          // Делим текущую строку на коэффициент
         for (int j = i; j < N; j++) M[i][j] /= a;

         b[i] /= a;

         for (int iK = i + 1; iK < N; iK++) // Вычитаем текущую строку из последующих
         {
            a = M[iK][i];

            for (int jK = i + 1; jK < N; jK++)  M[iK][jK] -= a * M[iK][jK];

            b[iK] -= a * b[i];
         }
      }
      else
         exit(-1);
   }

   for (int i = N - 1; i > -1; i--)
   {
      for (int j = N - 1; j > i; j--)  b[i] -= M[i][j] * q[j];
      q[i] = b[i] / M[i][i];
   }
}

real Integrate(real xL, real xR, real yL, real yR, real zL, real zR, function<real(real, real, real, real* args, int argNum)> f, real* args, int argNum)
{
   
   const int nKnot = 5;

   real xj[nKnot] = { -sqrt(5. + 2. * (sqrt(10. / 7.))) / 3., -sqrt(5. - 2. * (sqrt(10. / 7.))) / 3.,
              0. , sqrt(5. - 2. * (sqrt(10. / 7.))) / 3. , sqrt(5. + 2. * (sqrt(10. / 7.))) / 3. };
   real qj[nKnot] = { (322. - 13. * sqrt(70.)) / 900., (322. + 13. * sqrt(70.)) / 900., 128. / 225.,
                  (322. + 13. * sqrt(70.)) / 900., (322. - 13. * sqrt(70.)) / 900. };
   real hX = (xR - xL) / 2.,
      hY = (yR - yL) / 2.,
      hZ = (zR - zL) / 2.,
      cX = (xR + xL) / 2.,
      cY = (yR + yL) / 2.,
      cZ = (zR + zL) / 2.;

   real result = 0.;
   for (int i = 0; i < nKnot; i++)               // x
      for (int j = 0; j < nKnot; j++)            // y 
         for (int k = 0; k < nKnot; k++)         // z
            result += qj[i] * qj[j] * qj[k] * (f(cX + xj[i] * hX, cY + xj[j] * hY, cZ + xj[k] * hZ, args, argNum));

   return (xR - xL) * (yR - yL) * (zR - zL) * result / 8.;
}