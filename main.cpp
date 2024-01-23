#include <iostream>
#include <random>
#include <ctime>
#include <fstream>
using namespace std;

void LUD(double **A, double **L, double **U, double *x, double *b, int N);
void ShowMat(double **A, int n, int m);
void ShowVec(double *V, int n);
void MatVec(double **A, double *x, double *r, int N);

int main()
{
    double**A;
    double**L;
    double**U;
    double *b;
    double *x;
    int r;
    int c;


    fstream fin;

    fin.open("LU_Decomposition_Example.txt",ios::in);
    fin >> r;  

    A = new double *[r];

    for(int i = 0; i<r; i++)            ///Read in the textfile
    {
        A[i] = new double[r];
    }

    for(int i=0; i<r; i++)              ///Matrix generator
    {
        for(int j=0; j<r; j++)
        {
            fin>>A[i][j];
        }

    }

    cout << "A Matrix consists of " << r << " rows by " << r << " columns" << endl << endl;     ///Output the number of rows and columns

    b = new double[r];
    x = new double[r];
    for(int i = 0; i < r; i++)           ///Solution
    {
        x[i] = 0.0;
        fin >> b[i];
    }

    fin.close();

    fstream din;

    L = new double *[r];
    U = new double *[r];

    for(int i = 0; i < r; i++)            ///Read in the textfile
    {
        L[i] = new double[r];
        U[i] = new double[r];
    }

    for(int i = 0; i < r; i++)              ///Matrix generator with all 0s
    {
        for(int j = 0; j < r; j++)
        {
            L[i][j] = 0;
        }
    }

    for(int a=0; a<r; a++)              ///Matrix generator with all 0s
    {
        for(int h=0; h<r; h++)
        {
            U[a][h] = 0;
        }

    }

    cout<<"L and U Matrices consist of " << r << " rows by " << r << " columns"<<endl<<endl;     ///Output the number of rows and columns

    LUD(A, L, U, x,b,r);

    ofstream outfile;
    outfile.open("solutionfile.txt");              ///Print out the solution in the square root number of rows

    int rootOfRows = sqrt(r);

    int k;

    for(int m = 0; m < rootOfRows; m++)
    {
        for(int n = 0; n < rootOfRows; n++)
        {

        k = m * rootOfRows + n;
        outfile << x[k] << " ";
        }

        outfile << endl;
    }

    return 0;
}

void MatVec(double **A, double *x, double *r, int n) ///Multiply the vector.
{
    double sum = 0.0;
    for(int row=0; row<n; row++)
    {
        for( int col=0; col<n; col++)
        {
            sum = sum + A[row][col]*x[col];
        }
        r[row]=sum;

        sum = 0.0;
    }
}

void ShowMat(double **A, int n, int m) ///Display the matrix to the screen to see math that is changing the matrix.
{
    for(int i=0; i<n; i++)
    {
        for(int j=0; j<m; j++)
        {
            cout<<A[i][j]<<" ";
        }
        cout<<endl;
    }
    cout<<endl;
}

void ShowVec(double *V, int n)          ///Display the vector to the screen to see math that is changing the single vector.
{
    for(int j=0; j<n; j++)
    {
        cout<<V[j]<<" ";
    }
    cout<<endl<<endl;
}

void LUD(double **A, double **L, double **U, double *x, double *b, int N)
{
    double maxval;
    int maxcol;
    double temp;
    double currentRow;
    double currentCol;
    double sum;

    for(int col=0; col<N; col++)
    {
        for (int i = 0; i < N; i++)         ///Put 1.0 on the diagonal line of the U Matrix
        {
            for (int j = 0; j < N; j++)
            {
                if (col == col)
                {
                    U[col][col] = 1;
                }
            }
        }

    }

    for(int c = 0; c<N; c++)        ///Copy the first column of L = first column of A
    {
        L[c][0] = A[c][0];
    }


    for(int c = 0; c<N; c++)        ///First row of U = first row of A/A[0][0]
    {
        temp = A[0][0];
        U[0][c] = A[0][c]/temp;
    }

    for(int k=1; k<N; k++)
    {

        for(int row=k; row<N; row++)            ///First loop finds a row of U, loop is to solve a col of U (all the rows k and below)
        {
            sum = 0.0;

            for(int c = 0; c<k; c++)
            {
                sum = sum + L[row][c] * U[c][k];
            }

            L[row][k] = A[row][k] - sum;
        }



        for(int row = k; row<N; row++)        ///Second loop finds a column of U, loop to solve a row of L (all different columns > k)
        {
            sum = 0.0;

            for(int c=0; c<k; c++)                  ///everything above the current row k
            {
                sum = sum + L[k][c] * U[c][row];

            }

            U[k][row] = (A[k][row] - sum) / L[k][k];
        }


    }


    double* y = new double[N];                  ///Ax = b, L * U * x = y


    y[0] = b[0]/L[0][0];                        ///y = L/b
    for(int r = 1; r < N; r++)
    {
        sum = 0.0;
        for(int c = 0; c < r; c++)             ///So, all of the values of y below row k are multiplied by 0, however, all the values of u above have been previously computed.
        {
            sum = sum + L[r][c] * y[c];        ///Compute a sum
        }
        y[r] = (b[r] - sum)/L[r][r];
    }


    x[N-1]=y[N-1];
    for(int r = N-2; r>=0; r--)   ///Do the back solve
    {
        sum = 0;
        for(int c = r+1; c < N; c++)
        {
            sum = sum + U[r][c]*x[c];
        }
        x[r]= y[r]-sum;
    }


    double * error = new double[N];
    double * bpred = new double[N];
    double rel_error;

    MatVec(A,x,bpred,N);

    cout << endl << endl;
    sum=0.0;

    double errorSum=0.0;
    double bSum;  
    for(int i = 0; i<N; i++)          ///Compute the error
    {
        error[i] = b[i] - bpred[i];
        sum = sum + pow(error[i],2);
    }
    
    rel_error = sqrt(sum);

}

