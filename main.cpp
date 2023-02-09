//
// Created by patri on 27/01/2023.
//

#include <cstring>
#include "dependencies.h"
#include "rk4.h"


using namespace std;

void twoBody(double t, double y[], double max[], double tau[], double n[], double k[], double yp[]) {

    yp[0] = ((1 - (pow((y[4] / max[4]), n[0])) / (pow((y[4] / max[0]),  n[0]) + pow(k[0],  n[0]))) - (
            y[0] / max[0])) / tau[0];

    yp[1] = (((pow((y[0] / max[0]), n[1])) / (pow((y[0] / max[0]),  n[1]) + pow(k[1],  n[1]))) - (y[1] / max[1])) / tau[1];

    yp[2] = (((pow((y[1] / max[1]),  n[2])) / (pow((y[1] / max[1]),  n[2]) + pow(k[2],  n[2]))) - (y[2] / max[2])) / tau[2];

    yp[3] = (((pow((y[2] / max[2]),  n[3])) / (pow((y[2] / max[2]),  n[3]) + pow(k[3],  n[3]))) - (y[3] / max[3])) / tau[3];

    yp[4] = (((pow((y[3] / max[3]),  n[4])) / (pow((y[3] / max[3]),  n[4]) + pow(k[4],  n[4]))) - (y[4] / max[4])) / tau[4];

}


void twoBodyFixed(double t, double y[], double yp[]) {
    double max[] = {2.5619, 0.9499, 0.9299, 0.9475, 2.5619};
    double tau[] = {1.25, 4, 1.02, 1.57, 3.43};
    double n[] = {13, 4, 3, 4, 16};
    double k[] = {0.72, 0.50, 0.45, 0.51, 0.52};

    yp[0] = ((1 - (pow((y[4] / max[4]), n[0])) / (pow((y[4] / max[0]),  n[0]) + pow(k[0],  n[0]))) - (
            y[0] / max[0])) / tau[0];

    yp[1] = (((pow((y[0] / max[0]), n[1])) / (pow((y[0] / max[0]),  n[1]) + pow(k[1],  n[1]))) - (y[1] / max[1])) / tau[1];

    yp[2] = (((pow((y[1] / max[1]),  n[2])) / (pow((y[1] / max[1]),  n[2]) + pow(k[2],  n[2]))) - (y[2] / max[2])) / tau[2];

    yp[3] = (((pow((y[2] / max[2]),  n[3])) / (pow((y[2] / max[2]),  n[3]) + pow(k[3],  n[3]))) - (y[3] / max[3])) / tau[3];

    yp[4] = (((pow((y[3] / max[3]),  n[4])) / (pow((y[3] / max[3]),  n[4]) + pow(k[4],  n[4]))) - (y[4] / max[4])) / tau[4];

}

void f01 ( double t, double y[], double yp[] )

/******************************************************************************/
/*
  Purpose:

    F01 supplies the right hand side of the ODE for problem 1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 February 2012

  Author:

    John Burkardt

  Parameters:

    Input, double T, the time.

    Input, double Y[], the dependent variable.

    Output, double YP[], the value of the derivative.
*/
{
    yp[0] =   y[1];
    yp[1] = - y[0];

    return;
}


/******************************************************************************/

void test01 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST01 tests ODE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 February 2012

  Author:

    John Burkardt
*/
{
    double abserr;
    int i;
    int iflag;
    int iwork[5];
    int neqn = 2;
    double pi = 3.141592653589793;
    double relerr;
    int step_num = 12;
    double t;
    double tout;
    double *work;
    double *y;

    printf ( "\n" );
    printf ( "TEST01\n" );
    printf ( "  ODE solves a system of ordinary differential\n" );
    printf ( "  equations.\n" );
    printf ( "\n" );
    printf ( "      T           Y(1)         Y(2)\n" );
    printf ( "\n" );

    abserr = 0.00001;
    relerr = 0.00001;

    iflag = 1;

    t = 0.0;
    y = ( double * ) malloc ( neqn * sizeof ( double ) );
    y[0] = 1.0;
    y[1] = 0.0;

    printf ( "  %8g  %14g  %14g\n", t, y[0], y[1] );

    work = ( double * ) malloc ( ( 100 + 21 * neqn ) * sizeof ( double ) );

    for ( i = 1; i <= step_num; i++ )
    {
        tout = ( double ) ( i ) * 2.0 * pi / ( double ) ( step_num );

        ode ( f01, neqn, y, &t, tout, relerr, abserr, &iflag, work, iwork );

        if ( iflag != 2 )
        {
            printf ( "\n" );
            printf ( "TEST01 - Fatal error!\n" );
            printf ( "  ODE returned IFLAG = %d\n", iflag );
            break;
        }
        printf ( "  %8g  %14g  %14g\n", t, y[0], y[1] );
    }

    free ( work );
    free ( y );

    return;
}
/******************************************************************************/

void test02 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests ODE by integrating in the NEGATIVE time direction.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 February 2012

  Author:

    John Burkardt
*/
{
    double abserr;
    int i;
    int iflag;
    int iwork[5];
    int neqn = 2;
    double pi = 3.141592653589793;
    double relerr;
    int step_num = 12;
    double t;
    double tout;
    double *work;
    double *y;

    printf ( "\n" );
    printf ( "TEST02\n" );
    printf ( "  ODE solves a system of ordinary differential\n" );
    printf ( "  equations.\n" );
    printf ( "\n" );
    printf ( "  In this example, we integrate in the negative\n" );
    printf ( "  time direction.\n" );
    printf ( "\n" );
    printf ( "      T           Y(1)         Y(2)\n" );
    printf ( "\n" );

    abserr = 0.00001;
    relerr = 0.00001;

    iflag = 1;

    t = 0.0;
    y = ( double * ) malloc ( neqn * sizeof ( double ) );
    y[0] = 1.0;
    y[1] = 0.0;

    printf ( "  %8g  %14g  %14g\n", t, y[0], y[1] );

    work = ( double * ) malloc ( ( 100 + 21 * neqn ) * sizeof ( double ) );

    for ( i = 1; i <= step_num; i++ )
    {
        tout = - ( double ) ( i ) * 2.0 * pi / ( double ) ( step_num );

        ode ( f01, neqn, y, &t, tout, relerr, abserr, &iflag, work, iwork );

        if ( iflag != 2 )
        {
            printf ( "\n" );
            printf ( "TEST02 - Fatal error!\n" );
            printf ( "  ODE returned IFLAG = %d\n", iflag );
            break;
        }
        printf ( "  %8g  %14g  %14g\n", t, y[0], y[1] );
    }

    free ( work );
    free ( y );

    return;
}
/******************************************************************************/
int main_test ( )

/******************************************************************************/
/*
  Purpose:

    ode_test() tests ode().

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 February 2012

  Author:

    John Burkardt
*/
{
    timestamp ( );
    printf ( "\n" );
    printf ( "ODE_TEST\n" );
    printf ( "  C version\n" );
    printf ( "  Test the ODE library.\n" );

    test01 ( );
    test02 ( );
/*
  Terminate.
*/
    printf ( "\n" );
    printf ( "ODE_TEST\n" );
    printf ( "  Normal end of execution.\n" );
    printf ( "\n" );
    timestamp ( );

    return 0;
}


void outputToFile(string path, string text, bool append){
    ofstream outputf;

    if(append){
        outputf.open(path, std::ios_base::app);
    }else{
        outputf.open(path);
    }

    outputf << text;
    outputf.close();
}


void predator_phase_plot ( int n, int m, double t[], double y[] )

/******************************************************************************/
/*
  Purpose:

    predator_phase_plot makes a phase plot of the results.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 April 2020

  Author:

    John Burkardt

  Input:

    int n: the number of steps to take.

    int m: the number of variables.

    double t[n+1], y[(n+1)*m]: the times and solution values.
*/
{
    char command_filename[80];
    FILE *command;
    char data_filename[80];
    FILE *data;
    char header[] = "predator";
    int j;

    printf ( "\n" );
    printf ( "predator_phase_plot:\n" );
    printf ( "  Write command and data files that can be used\n" );
    printf ( "  by gnuplot for a predator-prey phase plot.\n" );
/*
  Create the data file.
*/
    strcpy ( data_filename, header );
    strcat ( data_filename, "_data.txt" );

    data = fopen ( data_filename, "wt" );

    for ( j = 0; j <= n; j++ )
    {
        fprintf ( data, "%g", t[j]);

        for(int i=0; i<m; i++){
            fprintf ( data, " %g",y[i+j*m]);
        }

        fprintf ( data, "\n");

        //fprintf ( data, "  %g  %g  %g\n", t[j], y[0+j*m], y[1+j*m] );
    }

    fclose ( data );

    printf ( "\n" );
    printf ( "  predator_phase_plot: data stored in \"%s\".\n", data_filename );
/*
  Create the command file.
*/
    strcpy ( command_filename, header );
    strcat ( command_filename, "_commands.txt" );

    command = fopen ( command_filename, "wt" );

    fprintf ( command, "# %s\n", command_filename );
    fprintf ( command, "#\n" );
    fprintf ( command, "# Usage:\n" );
    fprintf ( command, "#  gnuplot < %s\n", command_filename );
    fprintf ( command, "#\n" );
    fprintf ( command, "set term png\n" );
    fprintf ( command, "set output '%s.png'\n", header );
    fprintf ( command, "set xlabel '<-- Prey -->'\n" );
    fprintf ( command, "set ylabel '<-- Predator -->'\n" );
    fprintf ( command, "set title 'Predator-prey solution by rk4'\n" );
    fprintf ( command, "set grid\n" );
    fprintf ( command, "set style data lines\n" );
    fprintf ( command, "plot '%s' using 2:3 with lines\n", data_filename );
    fprintf ( command, "quit\n" );

    fclose ( command );

    printf ( "  predator_phase_plot: plot commands stored in \"%s\".\n",
             command_filename );

    return;
}

void predator_deriv ( double t, double y[], double f[] )

/******************************************************************************/
/*
  Purpose:

    predator_deriv returns the right hand side of the predator ODE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 April 2020

  Author:

    John Burkardt

  Input:

    double T, the current time.

    double Y[M], the current solution value.

  Output:

    double F[M], the value of the derivative, dU/dT.
*/
{
    double dfdt;
    double drdt;
    double fox;
    double rab;

    rab = y[0];
    fox = y[1];

    drdt =   2.0 * rab - 0.001 * rab * fox;
    dfdt = -10.0 * fox + 0.002 * rab * fox;

    f[0] = drdt;
    f[1] = dfdt;

    return;
}

void rk4_predator_test ( )

/*
  Purpose:

    rk4_predator_test tests RK4 on the predator prey ODE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 April 2020

  Author:

    John Burkardt
*/
{
    int m;
    int n = 1000;
    double *t;
    double tspan[2];
    double *y;
    double *y0;

    m = 2;
    n = 1000;

    t = ( double * ) malloc ( ( n + 1 ) * sizeof ( double ) );
    y = ( double * ) malloc ( ( n + 1 ) * m * sizeof ( double ) );
    y0 = ( double * ) malloc ( m * sizeof ( double ) );

    printf ( "\n" );
    printf ( "rk4_predator_test\n" );
    printf ( "  Use rk4() to solve the predator prey ODE.\n" );

    tspan[0] = 0.0;
    tspan[1] = 5.0;
    y0[0] = 5000.0;
    y0[1] = 100.0;

    rk4 ( predator_deriv, tspan, y0, n, m, t, y );

    predator_phase_plot ( n, m, t, y );
/*
  Free memory.
*/
    free ( t );
    free ( y );
    free ( y0 );

    return;
}

void grn_test ( )

{
    int m=5;
    int n = 49;
    double *t;
    double tspan[2];
    double *y;
    double *y0;
    double **vectors = new double*[m+1];

    for(int i=0; i<m+1; i++){
        vectors[i] = new double [50];
    }

    y0 = new double [m];

    for(int i=1; i< m; i++){
        y0[i] = vectors[i][0];
    }

    readFileToVectors("../GRN5.txt", m+1, vectors);


   // t = ( double * ) malloc ( ( n + 1 ) * sizeof ( double ) );
    y = ( double * ) malloc ( ( n + 1 ) * m * sizeof ( double ) );
   //y0 = ( double * ) malloc ( m * sizeof ( double ) );

    printf ( "\n" );
    printf ( "twoBodyFixed\n" );
    printf ( "  Use rk4() to solve the twoBody ODE.\n" );

    tspan[0] = 0.0;
    tspan[1] = 72.0;
    //y0[0] = 5000.0;
    //y0[1] = 100.0;

    rk4 ( twoBodyFixed, tspan, y0, n, m, vectors[0], y );

    predator_phase_plot ( n, m, vectors[0], y );
/*
  Free memory.
*/
    //free ( t );
    free ( y );
    delete y0;
    for(int i=0; i<m; i++){
        delete vectors[i];
    }

    return;
}
int main_test2 ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for rk4_test.

  Discussion:

    rk4_test tests rk4().

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 April 2020

  Author:

    John Burkardt
*/
{
    timestamp ( );
    printf ( "\n" );
    printf ( "rk4_test:\n" );
    printf ( "  C version\n" );
    printf ( "  Test rk4() .\n" );

    rk4_predator_test ( );
/*
  Terminate.
*/
    printf ( "\n" );
    printf ( "rk4_test:\n" );
    printf ( "  Normal end of execution.\n" );
    printf ( "\n" );
    timestamp ( );

    return 0;
}



double testFunc(Individual* ind){
    double eval =0;
    double x = ind->getDimension(0);
    double y = ind->getDimension(1);
    eval = pow(x + 2*y -7, 2) + pow(2*x + y - 5, 2);

    return eval;
}

int main(){

//    main_test2();

    grn_test();

    double *vectors[6];
    for(int i=0; i<6; i++){
        vectors[i] = new double [50];
    }

    readFileToVectors("../GRN5.txt", 6, vectors);
    readFile("../GRN5.txt", 5);
    return 0;

    cout << "\nHello\n";
    int numDim = 2;
    ESAlgorithm esAlgorithm = ESAlgorithm(numDim, 5, 10);

    for(int i=0; i<numDim; i++){
       esAlgorithm.setBounds(i, -10, 10, ESAlgorithm::LOWER_CLOSED, ESAlgorithm::UPPER_CLOSED);
    }
   // esAlgorithm.createPopulation(0, 10);

   esAlgorithm.setEvaluationFunction(testFunc);
   esAlgorithm.run1Plus1ES(1, 1.0, 0.817, 10, 250);
   cout << esAlgorithm.populationToCSVString() + "\n";

    esAlgorithm.setSigmaBounds(0.1, 10);
    esAlgorithm.runPopulationalIsotropicES(1, 0.5, 100, 5, 10);
    cout << esAlgorithm.populationToCSVString() + "\n";

    esAlgorithm.setSigmaBounds(0.1, 10);
    esAlgorithm.runPopulationalNonIsotropicES(1, 0.5, 100, 5, 10);
    cout << esAlgorithm.populationToCSVString() + "\n";

    cout << "\nFinished\n";
}

