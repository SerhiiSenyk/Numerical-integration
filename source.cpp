
#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <limits.h>
using namespace std;
void LeftRestanglesMethod (const double,const double,const double); 
void RightRestanglesMethod (const double,const double,const double); 
void MiddleRestanglesMethod (const double,const double,const double); 
void TrapezoidalRule (const double,const double,const double); 
void SimpsonsRule (const double,const double,const double);
double maxsecondDerivative(const double,const double);
inline double function(const double);
void enterData(double*,double*,double*);
void delimitation(const int,const char);
int main()
{
  double A,B,error;
  char ch[50];
  printf("\t\t\tProgram for calculating the integral\n\n");
  printf("Integral function : (sh(x) + 2)/ln^2*(2 + x^2)\n\n");
  printf("Enter 'A' to automatically enter;else enter  another letter : ");
  scanf ("%s",&ch);
	if ( *ch == 'A' || *ch == 'a' )
	{
		A = 0.0;
        B = 4.0;
		error = 0.001;
	}
	else{
       enterData(&A,&B,&error);
	}
  printf("\nA = %lf\n",A);
  printf("B = %lf\n",B);
  printf("Error = %.10lf\n",error);
  LeftRestanglesMethod (A,B,error);
  RightRestanglesMethod (A,B,error);
  MiddleRestanglesMethod (A,B,error);
  TrapezoidalRule (A,B,error);
  SimpsonsRule (A,B,error);
  printf("\n");
  system("pause");
  return 0;
}
void LeftRestanglesMethod (const double A,const double B,const double Eps)
{
  double X,H,Sum,result;
  int N;
  const int count = 3;
  double I[count];
  H = sqrt(Eps);
  for(int i = 0;i < count;++i)
  {
     N = (B - A)/H;
	 H = (B - A)/N;
	 Sum = 0;
	 X = A;
       for(int j = 0;j < N;++j)
       {
	      Sum += function(X);
	      X += H;
       }
     Sum *= H;
	 I[i] = Sum;
	 H /= 2;
  } 
  result = I[0] - ((I[0] - I[1])*(I[0] - I[1]))/(I[0] - 2*I[1] +I[2]);//схема Ейткена
  delimitation(80,'-');
  printf("\t\t\tLeft restangles method\n");
  printf("Choise step by indirect schemes(by empirical estimates).\n");
  printf("Aitken scheme\n");
  printf("The result with the error %.10lf = %.15lf\n",Eps,result);	
  printf("Last N = %d\n",N);
  printf("Last step = %.10lf\n",H);
  return;
}
void RightRestanglesMethod (const double A,const double B, double Eps)
{
  double X,H,Sum,result;//H - Крок
  int N;
  int i = 1;
  double I1 = 0;
  double I2 = 0;
  H = sqrt(Eps);
  while(true)
  {
     N = (B - A)/H;
	 H = (B - A)/N;
	 Sum = 0; 
	 X = B;
     for(int j = 0;j < N;++j)
       {
	      Sum += function(X);
	      X -= H;
       }
     Sum *= H;
	 i%2 == 0?I1 = Sum:I2 = Sum;
	 if(fabs(I1 - I2) < Eps)
	 {
       result = Sum;
	   break;
	 }	
	 H /= 2 ;
	 ++i;
  }
  delimitation(80,'-');
  printf("\n\t\t\tRight restangles method\n");
  printf("Choise step by indirect schemes(by empirical estimates).\n");
  printf("The result with the error %.10lf = %.15lf\n",Eps,result);	
  printf("N = %d\n",N);
  printf("Step = %.10lf\n",H);
  return;
}
void MiddleRestanglesMethod(const double A,const double B,double Eps)
{
  double X,H,Sum,result;//H - Крок
  int N;
  int i = 1;
  double I1 = 0;
  double I2 = 0;
  H = sqrt(Eps);
  while(true)
  {
     N = (B - A)/H;
	 H = (B - A)/N;
	 Sum = 0; 
	 X = A;
     for(int j = 0;j < N;++j)
       {
	      Sum += function(X + H/2);
	      X += H;
       }
     Sum *= H;
	 i%2 == 0?I1 = Sum:I2 = Sum;
	 if(fabs(I1 - I2) < Eps)
	 {
       result = Sum;
	   break;
	 }	
	 H /= 2 ;
	 ++i;
  }
  delimitation(80,'-');
  printf("\n\t\t\tMiddle restangles method\n");
  printf("Choise step by indirect schemes(by empirical estimates).\n");
  printf("The result with the error %.10lf = %.15lf\n",Eps,result);	
  printf("N = %d\n",N);
  printf("Step = %.10lf\n",H);
  return;
} 

void TrapezoidalRule(const double A,const double B,const double Eps)
{
  double X,H,Sum;//H - Крок
  int N;
  H = sqrt(12*Eps/(maxsecondDerivative(A,B)*(B - A)));
  N = (B - A)/H;
  H = (B - A)/N;
  Sum = (function(A) + function(B))/2.0; 
  X = A + H;
     for(int i = 0;i < N - 1;++i)
       {
	      Sum += function(X);
	      X += H;
       }
  Sum *= H;
  delimitation(80,'-');
  printf("\n\t\t\tTrapezoidal rule \n");
  printf("Selection step behind theoretical estimation error.\n");
  printf("The result with the error %.10lf = %.15lf\n",Eps,Sum);	
  printf("N = %d\n",N);
  printf("Step = %.10lf\n",H);
  return;
}
void SimpsonsRule(const double A,const double B,const double Eps)
{
  double X,H,Sum,result;//H - Крок
  int N,E;
  int i = 1;
  double I1 = 0;
  double I2 = 0;
  H = pow(Eps,0.25);
  while(true)
  {
	N =(B-A)/(2*H);
    H =(B-A)/(2*N);                                                        
	Sum = function(A) - function(B) ;                 
    E = -1;
	  for (int j = 2*N;j > 0;--j)
	  {
           Sum += (3+E)*function(A + j*H);
		   E = -E;
	  }  
	 Sum *= H/3;
	i%2 == 0?I1 = Sum:I2 = Sum;
	if(fabs(I1 - I2) < Eps)
	{
       result = Sum;
	   break;
	}	
	H /= 2 ;
    ++i;
  }
  delimitation(80,'-');
  printf("\n\t\t\tSimpson's rule\n");
  printf("Choise step by indirect schemes(by empirical estimates).\n");
  printf("The result with the error %.10lf = %.15lf\n",Eps,result);	
  printf("N = %d\n",N);
  printf("Step = %.10lf\n",H);
  delimitation(80,'-');
  return;
}
inline double function(double x	)
{
  return ((sinh(x) + 2)/(pow(log(2 + pow(x,2)),2)));
}
inline double secondDerivative(double x)
{
   if(x == 0) x = 0.0001;
   return pow(log(x*x + 2),-2)*( sinh(x) + 4*pow((x*x + 2),-2)*pow(log(x*x + 2),-1)* ((x*x + 2)*(-sinh(x) - 2) + 2*x*x* (sinh(x) - cosh(x)*pow(x,-1)*(x*x + 2) + 3*(sinh(x) + 2)*pow(log(x*x + 2),-1) + 2 ) ) );
}
double maxsecondDerivative(const double A,const double B)
{
   double temp = 0.001;
   double max = fabs(secondDerivative(A));
   for(double X = A;X <= B;X += temp)
   {
        if(fabs(secondDerivative(X)) > max)
		{
           max = fabs(secondDerivative(X));
		}
   }
   return max;
}
void delimitation(const int w,const char ch)
{
	for(int i = 0;i < w;++i)
	{	
		printf("%c",ch);
	}
   return;
}
void enterData(double *A,double *B,double *error)
{ 
   delimitation(80,'*');
   printf("\n\t\t\tEntering data\n");
   while(true)
   {
     printf("Enter left border : ");
     while(scanf("%lf",A) != 1)
     {
       printf("Vvedeno ne chyslo.Vvedit' zanovo : ");
	   cin.clear();
	   cin.ignore(cin.rdbuf()->in_avail());
     }
     if(*A < -100 || *A > 100)
     printf("-100 <= left border <= 100 !!!\n");
	 else break;
   }
   cin.clear();
   cin.ignore(cin.rdbuf()->in_avail());
   while(true)
   {
     printf("Enter right border : ");
     while(scanf("%lf",B) != 1)
     { 
       printf("Vvedeno ne chyslo.Vvedit' zanovo : ");
	   cin.clear();
	   cin.ignore(cin.rdbuf()->in_avail());
     }
     if(*B < -100 || *B > 100)
     printf("-100 <= right border <= 100 !!!\n");
	 else break;
   }
   cin.clear();
   cin.ignore(cin.rdbuf()->in_avail());
   while(true)
   {
     printf("Enter error : ");
     while(scanf("%lf",error) != 1)
     {
       printf("Vvedeno ne chyslo.Vvedit' zanovo : ");
       cin.clear();
	   cin.ignore(cin.rdbuf()->in_avail());
     }
     if(*error > 1 || *error <= 0)
     printf("Error < 1 && Error > 0!!!\n");
	 else break;
   }
   putchar('\n');
   delimitation(80,'*');
   return;
}
