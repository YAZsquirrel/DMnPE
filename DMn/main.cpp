#include "main.h"

real* g;
real* dg;
real** A;
real* b;
int N; // число приемников
real gamma;
int nx, ny, nz; // число ячеек по коордам

typedef real** Matrix;

struct origin
{
   real x[2], y[2], z[2];
} Orig;

int main()
{
   //reciever Recs;
   cout << Integrate(-2, 3, 1, 2, 0, 2, [](real x, real y, real z) {return pow(x*y, 5.) + x*y*z*z*z; });
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
         
   nRecX = (xR - xL) / xH;
   nRecY = (yR - yL) / yH;

   RecX = new real[nRecX]{};
   RecY = new real[nRecY]{};

   for (int yI = 0; yI < nRecY; yI++)
      RecY[yI] = yL + yI*yH;
   
   for (int xI = 0; xI < nRecX; xI++)
      RecX[xI] = xL + xI*xH;
}

void SLAEgen(int nX, int nY, int nZ, int nRecX, int nRecY, real* X, real *Y, real*Z, real*RecX,)
{
    int K = nX*nY*nZ;
    int nRec = nRecX * nRecY;
    Matrix A = new real * [K] {};

    for (int i = 0; i < K; i++)
        A[i] = new real[K]{};
    
    real sum;
    for(int i = 0; i < K; i++)
        for (int j = 0; j < K; j++)
        {
            sum = 0;
            for (int iRec = 0; iRec < nRec; iRec++)
            {
                sum += Integrate(X[i], X[i + 1], real yL, real yR, real zL, real zR, function<real(real, real, real)> f);
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

    for (real x = xL, int i=0; x < xR; x+=hX, i++)
    {
        xMesh[i] = x;
    }

    real yL, yR;
    int nY;
    in >> yL >> yR >> nY;
    real* yMesh = new real[nY] {};

    real hY = (yR - yL) / nY;

    for (real y = yL, int i = 0; y < yR; y += hY, i++)
    {
        yMesh[i] = y;
    }
    real zL, zR;
    int nZ;
    in >> zL >> zR >> nZ;
    real* zMesh = new real[nZ] {};

    real hZ = (zR - zL) / nZ;

    for (real z = zL, int i = 0; z < zR; z += hZ, i++)
    {
        zMesh[i] = z;
    }

   

}

real Integrate(real xL, real xR, real yL, real yR, real zL, real zR, function<real(real, real, real)> f)
{
   real result = 0.;
   //const int segs = 5;

   real xj[5] = { -sqrt(5. + 2. * (sqrt(10. / 7.)))/3., -sqrt(5. - 2. * (sqrt(10. / 7.)))/3. , 
              0. , sqrt(5. - 2. * (sqrt(10. / 7.)))/3. , sqrt(5. + 2. * (sqrt(10. / 7.)))/3. };
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
            result += qj[i] * qj[j] * qj[k] * (f( cX + xj[i] * hX, cY + xj[j] * hY, cZ + xj[k] * hZ));

   return (xR - xL) * (yR - yL) * (zR - zL) * result / 8.;
}