#include "main.h"
#define OPTIM_SYNTETIC
typedef real** Matrix;

struct recInfo	 // Структура данных приемников
{
   real *RecX, *RecY;		// координаты приемников на поверхности (z = 0)
   int nRecX, nRecY;		// колво приемников по осям
};

struct meshInfo // Структура данных сетки
{
   real *X, *Y, *Z;		// координаты узлов секти по осям
   int nX, nY, nZ;	// количество узлов сетки по осям
};

// Ввод координат тела аномалии (для синтетических данных)
void BodyInput(real& xL, real& xR, real& yL, real& yR, real& zL, real& zR)	
{
   ifstream in("Body.txt");

   in >> xL >> xR >> yL >> yR >> zL >> zR;

   in.close();
}

// лол
real ro(real x, real y, real z)
{
   return 1;
}

// каво-чиво
int findIndCoord(real bodyCoord, real* X, int n, bool right)
{
   int i, iPrev = 0;
   int iR = n - 1, iL = 0;

   for (i = iR; (bodyCoord <= X[i - 1] || X[i] <= bodyCoord ) && i != iPrev; i = (iL + iR) / 2)
   {
      iPrev = i;
      (X[i] > bodyCoord ? iR : iL) = i;
   }

   return right ? i: --i;
}

// Ввод синтетических данных для тестирования программы

void Syntetics(meshInfo mesh, recInfo rec)
{
   real* RecX = rec.RecX, *RecY = rec.RecY;
   int nRecX = rec.nRecX, nRecY = rec.nRecY;

   real* X = mesh.X, *Y = mesh.Y, *Z = mesh.Z;
   int nX = mesh.nX, nY = mesh.nY, nZ = mesh.nZ;

   real xL, xR, yL, yR, zL, zR;
   BodyInput(xL, xR, yL, yR, zL, zR);
   
   ofstream out("RecData.txt");
   out << nRecX * nRecY << endl;
   
   
   auto G = [](real ro, real r, real coord, real mes) {return mes * ro / (4.0 * M_PI * pow(r, 3.)) * coord; };

   xL = 3, xR = 3.1, X = new real[3]{ 0, 2, 4 }, nX = 3;
   int ixL = findIndCoord(xL, X, nX,0);
   int ixR = findIndCoord(xR, X, nX,1);
   int iyL = findIndCoord(yL, Y, nX,0);
   int iyR = findIndCoord(yR, Y, nX,1);
   int izL = findIndCoord(zL, Z, nX,0);
   int izR = findIndCoord(zR, Z, nX,1);

   if (ixL == ixR || iyL == iyR || izL == izR) exit(-100);

   real mes, rSqr[3]{}, sumG;

   for (int recI = 0; recI < nRecX; recI++)
   { 
      for (int recJ = 0; recJ < nRecY; recJ++)
      {
         sumG = 0.;
         for (int ix = ixL; ix < ixR; ix++)
         {
             rSqr[0] = pow(X[ix] - RecX[recI], 2.);
            for (int iy = iyL; iy < iyR; iy++)
            {
                rSqr[1] = pow(Y[iy] - RecY[recJ], 2.);
               for (int iz = izL; iz < izR; iz++)
               {
                  rSqr[2] = pow(Z[iz], 2.);
                  mes =  (X[ix + 1] - X[ix])*(Y[iy + 1] - Y[iy])*(Z[iz + 1] - Z[iz]);
                  sumG += G(ro(X[ix], Y[iy], Z[iz]), sqrt(rSqr[0] + rSqr[1] + rSqr[2]), Z[iz], mes);
               }
            }
         }
         out << sumG << "\n";
      }
   }
}

void readSyntetics(string fileName, real** G)
{
   ifstream in(fileName);
   int nG;
   in >> nG;

   real* localG =  *G = new real[nG]{};

   for (int i = 0; i < nG; i++) in >> localG[i];
}

// Заполнение координат по осям в сетку
void FillCoordAxis(real** coordAxis, real left, real right, int n)
{
   real* axis = *coordAxis = new real[n]{};

   real h = (right - left) / n;
   int i = 0;
   for (real coord = left; coord < right; coord += h, i++)  axis[i] = coord;
   if (i == 0) cout << "Error 1: Left coord > right coord";
}

// Создание приемников
void CreateRec(string fileName, recInfo &rec)
{
   ifstream in(fileName);

   real xL, xR, yL, yR, xH, yH;

   in >> xL >> xR >> xH;
   in >> yL >> yR >> yH;

   rec.nRecX = (int)((xR - xL) / xH);
   rec.nRecY = (yR - yL) == 0 ? 1 : (int)((yR - yL) / yH);

   FillCoordAxis(&(rec.RecX), xL, xR, rec.nRecX);
   FillCoordAxis(&(rec.RecY), yL, yR, rec.nRecY);
}

// Генерация матрицы и вектора правой части для решения МНК-функционала 
void SLAEgen(meshInfo mesh, recInfo rec, real* g, Matrix& A, real** b, int& K )
{
   int nX = mesh.nX, nY = mesh.nY, nZ = mesh.nZ;
   real* X = mesh.X, *Y = mesh.Y, *Z = mesh.Z;
   
   int nRecX = rec.nRecX, nRecY = rec.nRecY;
   real* RecX = rec.RecX, *RecY = rec.RecY;

   K = nX * nY * nZ;   // Количество ячеек

   A = new real * [K] {};   // Глобальная матрица СЛАУ
   for (int i = 0; i < K; i++)
     A[i] = new real[K]{};

   real* bb = *b = new real[K]{}; // Глобальная правая часть СЛАУ

   function<real(real, real, real, real*, int)> f = [](real xB, real yB, real zB, real* args, int argNum) // подынтегральная функция
   {
      real r = sqrt(pow(args[0] - xB, 2) + pow(args[1] - yB, 2) + pow(args[2] - zB, 2));   // Подсчет нормы вектора
      return zB / pow(r, 3);
   };

   real RecCoords[3]{};
   real gq, gs;
   int iX, iY, iZ, jX, jY, jZ;
   int nXY = nX * nY, ir = 0;

   // Заполнение вкладов ячеек в значения на приемниках
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
            bb[q] += gq * g[ir];
         }
         ir++;
      }
   }

   // Дозаполнение вкладов ячеек в значения на приемниках
   for (int q = 0; q < K; q++) 
   {
      bb[q] /= 4. * M_PI;
      for (int s = 0; s < K; s++) A[q][s] /= 4. * M_PI * 4. * M_PI;
      //std::cout << A[q][q] << " ";
      //cout << bb[q] << "\n";
   }


}

// Чтение координат сетки
void CreateMesh(string fileName, meshInfo& mesh)
{
   ifstream in(fileName);

   real xL, xR, yL, yR, zL, zR;
   in >> xL >> xR >> mesh.nX;
   in >> yL >> yR >> mesh.nY;
   in >> zL >> zR >> mesh.nZ;

   FillCoordAxis(&mesh.X, xL, xR, mesh.nX + 1);
   FillCoordAxis(&mesh.Y, yL, yR, mesh.nY + 1);
   FillCoordAxis(&mesh.Z, zL, zR, mesh.nZ + 1);
}

// Решение СЛАУ методом Гаусса без выбора ведущего элемента
void SLAEsolve(Matrix M, int N, real* b, real** q)
{
   //GAUSSE
   real* qq = *q = new real[N]{};
   real a;	// Коэффициент Гаусса

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

            for (int jK = i; jK < N; jK++)  M[iK][jK] -= a * M[i][jK];

            b[iK] -= a * b[i]; 
         }
      }
      else
         exit(-1);
   }

   for (int i = N - 1; i > -1; i--)
   {
      for (int j = N - 1; j > i; j--)  b[i] -= M[i][j] * qq[j];
      qq[i] = b[i] / M[i][i];
   }
}

// Численное интегрирование усложненным методом Гаусса
real Integrate(real xL, real xR, real yL, real yR, real zL, real zR, function<real(real, real, real, real* args, int argNum)> f, real* args, int argNum)
{
   const int nKnot = 5; // Количество узлов Гаусса

   real xj[nKnot] = { -sqrt(5. + 2. * (sqrt(10. / 7.))) / 3., -sqrt(5. - 2. * (sqrt(10. / 7.))) / 3.,	// Значения координат на узлах
              0. , sqrt(5. - 2. * (sqrt(10. / 7.))) / 3. , sqrt(5. + 2. * (sqrt(10. / 7.))) / 3. };

   real qj[nKnot] = { (322. - 13. * sqrt(70.)) / 900., (322. + 13. * sqrt(70.)) / 900., 128. / 225.,	// Веса узлов
                  (322. + 13. * sqrt(70.)) / 900., (322. - 13. * sqrt(70.)) / 900. };
       // Шаги
   real hX = (xR - xL) / 2.,
       hY = (yR - yL) / 2.,
       hZ = (zR - zL) / 2., 
       // Центры
       cX = (xR + xL) / 2.,   
       cY = (yR + yL) / 2.,   
       cZ = (zR + zL) / 2.;   

   real result = 0.;
   for (int ix = 0; ix < nKnot; ix++)               
      for (int iy = 0; iy < nKnot; iy++)           
         for (int iz = 0; iz < nKnot; iz++)       
            result += qj[ix] * qj[iy] * qj[iz] * (f(cX + xj[ix] * hX, cY + xj[iy] * hY, cZ + xj[iz] * hZ, args, argNum));
   //cout << (xR - xL) * (yR - yL) * (zR - zL) * result / 8. << endl;
   return (xR - xL) * (yR - yL) * (zR - zL) * result / 8.; // Масштабирование
}  

void WriteResult(string fileName, real* res, int nY, int nZ, int nX)
{
   
   ofstream out(fileName);
   out.precision(3);
   //out.scientific;
   int i, j , k, ind;
   for(j = 0, ind = 0; j < nY; j++)
   {
        for(k = nZ - 1; k > -1; k--)
        {
            for(i = 0; i < nX; i++, ind++)
                out << res[ind] << " ";
            out << endl;
        }
        out << endl;
   }
   /*
   for (int i = 0; i < K; i++)
   {
      out << res[i] << "\n";
   }
   */
}

int main()
{
   meshInfo* mesh = new meshInfo;
   CreateMesh("Mesh.txt", *mesh);
   
   recInfo* rec =  new recInfo;
   CreateRec("Params.txt", *rec);
   
   Syntetics(*mesh, *rec);
   
   real* g; // вектор значений поля с приемников
   readSyntetics("RecData.txt", &g);
   
   Matrix A;
   real *b;
   int K;

   SLAEgen(*mesh, *rec, g, A, &b, K);

   real* res;
   SLAEsolve(A, K, b, &res);

   WriteResult("SolutionF.txt", res, mesh->nX, mesh->nY, mesh->nZ);

   return 0;
}