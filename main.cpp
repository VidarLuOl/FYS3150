#include <iostream>
#include <stdio.h>
#include <vector>
#include <fstream>
#include <math.h>
#include <cstdlib>      //abs()
#include <chrono>       //time
#include <armadillo>    //Matrix

using namespace std;
using namespace arma;

void first(unsigned int n)
    {
        unsigned int i;

        double * b;
        double * b_;
        double * g_;
        double * g;
        double * u;
        double h = 1.0/(n+1);
        double hh = h*h;
        double * x;

        int * a;
        int * c;

        a = new int [n-1];
        c = new int [n-1];

        b = new double [n];
        b_ = new double [n];
        g = new double [n];
        g_ = new double [n];
        u = new double [n]();
        x = new double [n];

        for(i = 0; i < n; i++){
            a[i] = -1;
            c[i] = -1;
            b[i+1] = 2;
            x[i] = i*h;
            b_[i] = 0;
            g[i] = hh*(100*exp(-10*(i*h)));
        }

        b_[0] = 0;
        b_[1] = 2;
        g_[0] = g[0];
        g_[1] = g[1];


        auto start = std::chrono::high_resolution_clock::now();

        for(i = 2; i < n-1; i++){
            b_[i] = b[i] - (a[i-1]*c[i-1])/b_[i-1];
            g_[i] = g[i] - (g_[i-1]*a[i-1])/b_[i-1];
        }

        for(i = n-2; i > 0; i--){
            u[i] = (g_[i] - a[i]*u[i+1])/b_[i];
        }

        auto finish = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double> elapsed = finish - start;
        cout << elapsed.count() << endl;

        ofstream myfile;
        myfile.open ("verdier.txt");
        myfile << n << " " << elapsed.count() <<  endl;
        for(i = 0; i < n; i++){
             myfile << x[i] << " " << u[i] << endl;
        }
        myfile.close();
}

void second(unsigned int n)
{
    unsigned int i;

    double * b_;
    double * g_;
    double * g;
    double * u;
    double h = 1.0/(n+1);
    double hh = h*h;
    double * x;

    b_ = new double [n];
    g = new double [n];
    g_ = new double [n];
    u = new double [n]();
    x = new double [n];

    for(i = 0; i < n; i++){
        x[i] = i*h;
        g[i] = hh*(100*exp(-10*(i*h)));
    }

    b_[0] = 0;
    b_[1] = 2;
    g_[0] = g[0];
    g_[1] = g[1];

    auto start = std::chrono::high_resolution_clock::now();

    for(i = 2; i < n-1; i++){
        b_[i] = 2 - 1/b_[i-1];
        g_[i] = g[i] - (g_[i-1]*(-1))/b_[i-1];
    }

    for(i = n-2; i > 0; i--){
        u[i] = (g_[i] - (-1)*u[i+1])/b_[i];
    }

    for(i = 0; i < n; i++){
        cout << u[i] << endl;
    }

    auto finish = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed = finish - start;

    ofstream myfile;
    myfile.open ("verdier.txt");
    myfile << n << " " << elapsed.count() << endl;
    for(i = 0; i < n; i++){
         myfile << x[i] << " " << u[i] << endl;
         cout << u[i] << endl;
    }
    myfile.close();
}

void third(unsigned int n)
{
    unsigned int i;

    double * b_;
    double * g_;
    double * g;
    double * u;
    double h = 1.0/(n+1);
    double hh = h*h;
    double * x;
    double * e;
    double * q;
    double m = 1 - exp(-10);

    b_ = new double [n];
    g = new double [n];
    g_ = new double [n];
    u = new double [n]();
    x = new double [n];
    e = new double [n];
    q = new double [n];

    double p = 1/n;
    for(i = 0; i < n; i++){
        x[i] = p*i;
    }

    for(i = 0; i < n; i++){
        x[i] = i*h;
        g[i] = hh*(100*exp(-10*(i*h)));
    }

    auto start = std::chrono::high_resolution_clock::now();


    for(i = 1; i < n; i++){
        b_[i] = 2 - 1/b_[i-1];
        g_[i] = g[i] - (g_[i-1]*(-1))/b_[i-1];
    }

    for(i = n-1; i > 0; i--){
        u[i] = (g_[i] - (-1)*u[i+1])/b_[i];
    }

    for(i = 0; i < n; i++){
        q[i] = 1 - m*i - exp(-10*(int)i);
        e[i] = log10(abs((u[i]-q[i])/q[i]));
    }

    auto finish = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed = finish - start;

    ofstream myfile;
    myfile.open ("verdier.txt");
    myfile << n << " " << elapsed.count() << endl;
    for(i = 1; i < n; i++){
         myfile << e[i] << " " << u[i] << endl;
    }
    myfile.close();
}

void fourth(unsigned int n){
    {
        double h = 1.0/n;
        double hh = h*h;
        unsigned int i;

        vec x;
        x = vec(n);
        vec y;
        y = vec(n);
        vec b;
        b = vec(n);

        auto start = std::chrono::high_resolution_clock::now();
        x(0) = h;
        x(n-1) = x(0) + (n-1)*h;
        b(0) = hh * 100.0*exp(-10.0*h);
        b(n-1) = hh * 100.0*exp(-10.0*(h + (n-1)*h));


        Mat<double> A(n, n, fill::zeros);
        for(unsigned int i=0; i<n; i++){
            A(i,i) = 2;
        }


        for(unsigned int i=1; i<n-1; i++){
            A(i+1,i) = -1;
            A(i,i+1) = -1;
            x(i) = x(i-1) + h;
            b(i) = hh * 100.0*exp(-10.0*x(i));
        }

        Mat<double> L, U;

        lu(L, U, A);

        vec Ans = solve(A, b);

        auto finish = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double> elapsed = finish - start;

        ofstream myfile;
        myfile.open ("verdier2.txt");
        myfile << n << " " << elapsed.count() << endl;
        for(i = 0; i < n; i++){
             myfile << h*i << endl;
        }
        myfile << Ans << endl;
        myfile.close();

    }
}

int main()
{
    unsigned int n = 10000;

    //first(n);     //Opg 1.b
    //second(n);    //Opg 1.c
    //third(n);     //Opg 1.d
    fourth(n);     //Opg 1.e

    return 0;
}
