//  (C) Copyright Gennadiy Rozental 2005.
//  Distributed under the Boost Software License, Version 1.0.
//  (See accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org/libs/test for the library home page.

// Boost.Test

// each test module could contain no more then one 'main' file with init function defined
// alternatively you could define init function yourself
#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "GRNEDOHelpers.h"
#include "algModes.h"
#include "GRNSeries.h"
using namespace algModes;
using namespace GRNEDOHelpers;

//____________________________________________________________________________//
/*
// most frequently you implement test cases as a free functions with automatic registration
BOOST_AUTO_TEST_CASE( test1 )
{
        // reports 'error in "test1": test 2 == 1 failed'
        BOOST_CHECK( 2 == 1 );
}

//____________________________________________________________________________//

// each test file may contain any number of test cases; each test case has to have unique name
BOOST_AUTO_TEST_CASE( test2 )
{
        int i = 0;

                // reports 'error in "test2": check i == 2 failed [0 != 2]'
        BOOST_CHECK_EQUAL( i, 2 );

        BOOST_CHECK_EQUAL( i, 0 );
}
*/
BOOST_AUTO_TEST_CASE( test_ind_5_var_lsoda )
{
    appContext ctx{};

    //esse indivíduo deveria ter fitness ~26
    double ind0[19] = {1.2163355099083872, 1.1264485098219865, 2.973714367061704,
                       2.952143123315177, 2.998260518457365, 0.5687249950503857,
                       0.4580723119903261, 0.46214892372246563, 0.6182568295500336,
                       0.5213082492659304, 0.7708877748759901, 0.1497642024548283,
                       4.254757908429968, 3.759370669969996, 4.784173526119725,
                       10.935884810737809, 24.595975874929724, 2.8109199678182635,
                       4.922623602327875};

    initializeGRN5Context(&ctx, SINGLE_SET_MODE, 1);
    double eval = grn5EvaluationLSODA(ind0, &ctx);

    clearContext(&ctx);

    BOOST_CHECK_CLOSE_FRACTION( eval, 26.92, 0.001 );

}


BOOST_AUTO_TEST_CASE( test_ind_5_var_rk4 )
{
    appContext ctx{};

    //esse indivíduo deveria ter fitness ~26
    double ind0[19] = {1.2163355099083872, 1.1264485098219865, 2.973714367061704,
                       2.952143123315177, 2.998260518457365, 0.5687249950503857,
                       0.4580723119903261, 0.46214892372246563, 0.6182568295500336,
                       0.5213082492659304, 0.7708877748759901, 0.1497642024548283,
                       4.254757908429968, 3.759370669969996, 4.784173526119725,
                       10.935884810737809, 24.595975874929724, 2.8109199678182635,
                       4.922623602327875};

    initializeGRN5Context(&ctx, SINGLE_SET_MODE, 1);
    double eval = grn5EvaluationRK4(ind0, &ctx);
    clearContext(&ctx);

    BOOST_CHECK_CLOSE_FRACTION( eval, 26.92, 0.001 );

}

BOOST_AUTO_TEST_CASE( test_ind_10_var_lsoda )
{
    appContext ctx{};

    //esse indivíduo deveria ter fitness ~56
    double ind0[40] = {1.73,2,0.81,0.11, 1.23, 1.78,
                       1.14, 1.04, 3.47, 3.21, 0.45,
                       0.56, 0.99, 0.77, 0.71, 0.66,
                       0.46, 0.48, 0.66, 0.99, 0.85,
                       0.61, 0.55, 0.46, 0.17, 20,
                       9, 24, 12, 2, 2, 6, 4, 7,
                       24, 2, 7, 21, 20, 3};

    initializeGRN10Context(&ctx, SINGLE_SET_MODE, 1);
    double eval = grn10EvaluationLSODA(ind0, &ctx);
    clearContext(&ctx);

    BOOST_CHECK_CLOSE_FRACTION( eval, 56.71, 0.001 );

}


BOOST_AUTO_TEST_CASE( test_ind_10_var_rk4 )
{
    appContext ctx{};

    //esse indivíduo deveria ter fitness ~56
    double ind0[40] = {1.73,2,0.81,0.11, 1.23, 1.78,
                       1.14, 1.04, 3.47, 3.21, 0.45,
                       0.56, 0.99, 0.77, 0.71, 0.66,
                       0.46, 0.48, 0.66, 0.99, 0.85,
                       0.61, 0.55, 0.46, 0.17, 20,
                       9, 24, 12, 2, 2, 6, 4, 7,
                       24, 2, 7, 21, 20, 3};

    initializeGRN10Context(&ctx, SINGLE_SET_MODE, 10);
    double eval = grn10EvaluationRK4(ind0, &ctx);
    clearContext(&ctx);

    BOOST_CHECK_CLOSE_FRACTION( eval, 56.71, 0.001 );

}

BOOST_AUTO_TEST_CASE( test_GRNSeries_should_return_correct_number_of_rows_and_columns )
{
    GRNSeries series = GRNSeries();
    series.initializeMatrix("testGRN.txt");
    BOOST_CHECK_EQUAL(series.getNumTimeSteps(), 2 );
    BOOST_CHECK_EQUAL(series.getNumColumns(), 3 );
}

BOOST_AUTO_TEST_CASE(test_GRNSeries_should_return_correct_vector){
    GRNSeries series = GRNSeries();
    series.loadFromFile("testGRN.txt");
    double **vectors = series.getVectors();

    BOOST_CHECK_CLOSE_FRACTION( vectors[0][0], 0.0, 0.001 );
    BOOST_CHECK_CLOSE_FRACTION( vectors[1][0], 0.7095, 0.001 );
    BOOST_CHECK_CLOSE_FRACTION( vectors[2][0], 0.1767, 0.001 );
    BOOST_CHECK_CLOSE_FRACTION( vectors[0][1], 1.4694, 0.001 );
    BOOST_CHECK_CLOSE_FRACTION( vectors[1][1], 1.1517, 0.001 );
    BOOST_CHECK_CLOSE_FRACTION( vectors[2][1], 0.3415, 0.001 );

}

BOOST_AUTO_TEST_CASE(test_GRNSeries_should_split_correctly_start){
    GRNSeries series = GRNSeries();
    series.loadFromFile("testGRN2.txt");
    GRNSeries splitSeries = GRNSeries(series, 0, 1);
    double **vectors = splitSeries.getVectors();

    BOOST_CHECK_CLOSE_FRACTION( vectors[0][0], 0.0, 0.001 );
    BOOST_CHECK_CLOSE_FRACTION( vectors[1][0], 0.7095, 0.001 );
    BOOST_CHECK_CLOSE_FRACTION( vectors[2][0], 0.1767, 0.001 );
    BOOST_CHECK_CLOSE_FRACTION( vectors[0][1], 1.4694, 0.001 );
    BOOST_CHECK_CLOSE_FRACTION( vectors[1][1], 1.1517, 0.001 );
    BOOST_CHECK_CLOSE_FRACTION( vectors[2][1], 0.3415, 0.001 );

}

BOOST_AUTO_TEST_CASE(test_GRNSeries_should_split_correctly_end){
    GRNSeries series = GRNSeries();
    series.loadFromFile("testGRN2.txt");
    GRNSeries splitSeries = GRNSeries(series, 2, 2);
    double **vectors = splitSeries.getVectors();

    BOOST_CHECK_CLOSE_FRACTION( vectors[0][0], 2.0, 0.001 );
    BOOST_CHECK_CLOSE_FRACTION( vectors[1][0], 3.5, 0.001 );
    BOOST_CHECK_CLOSE_FRACTION( vectors[2][0], 4.9, 0.001 );

}

BOOST_AUTO_TEST_CASE(test_GRNSeries_should_return_correct_max_values){
    GRNSeries series = GRNSeries("testGRN.txt");
    double *maxValues = series.getMaxValues();

    BOOST_CHECK_CLOSE_FRACTION( maxValues[0], 1.1517, 0.001 );
    BOOST_CHECK_CLOSE_FRACTION( maxValues[1], 0.3415, 0.001 );

}

//todo: fazer o teste abaixo com diferentes conjuntos
BOOST_AUTO_TEST_CASE(test_initialize_GRN5_context_LSODA){
    //void GRNEDOHelpers::initializeGRNContext(appContext* ctx, int granularity, int numVariables, int numTau, int numN, int numK, int setStart, int setEnd, double** vectors, double * maxValues)
    appContext ctx{};
    GRNSeries series = GRNSeries("GRN5.txt");
    GRNSeries trainingSeries = GRNSeries(series, 0, 49);
    double **vectors = trainingSeries.getVectors();
    double *maxValues = trainingSeries.getMaxValues();
    int numVariables = trainingSeries.getNumColumns() -1;

    double ind0[19] = {1.2163355099083872, 1.1264485098219865, 2.973714367061704,
                       2.952143123315177, 2.998260518457365, 0.5687249950503857,
                       0.4580723119903261, 0.46214892372246563, 0.6182568295500336,
                       0.5213082492659304, 0.7708877748759901, 0.1497642024548283,
                       4.254757908429968, 3.759370669969996, 4.784173526119725,
                       10.935884810737809, 24.595975874929724, 2.8109199678182635,
                       4.922623602327875};

    initializeGRNContext(&ctx, 1, numVariables, 5, 7, 7, 0, 49, vectors, maxValues);

    double eval = grn5EvaluationLSODA(ind0, &ctx);
    clearContext2Test(&ctx);

    BOOST_CHECK_CLOSE_FRACTION( eval, 26.92, 0.001 );

}

//todo: fazer o teste abaixo com diferentes conjuntos
BOOST_AUTO_TEST_CASE(test_evaluate_GRN5_LSODA){
    //void GRNEDOHelpers::initializeGRNContext(appContext* ctx, int granularity, int numVariables, int numTau, int numN, int numK, int setStart, int setEnd, double** vectors, double * maxValues)
   // appContext ctx{};
    GRNSeries series = GRNSeries("GRN5.txt");
    GRNSeries testSeries = GRNSeries(series, 35, 49);
    double **vectors = testSeries.getVectors();
    double *maxValues = testSeries.getMaxValues();
    int numVariables = testSeries.getNumColumns() -1;

    double ind0[19] = {1.2163355099083872, 1.1264485098219865, 2.973714367061704,
                       2.952143123315177, 2.998260518457365, 0.5687249950503857,
                       0.4580723119903261, 0.46214892372246563, 0.6182568295500336,
                       0.5213082492659304, 0.7708877748759901, 0.1497642024548283,
                       4.254757908429968, 3.759370669969996, 4.784173526119725,
                       10.935884810737809, 24.595975874929724, 2.8109199678182635,
                       4.922623602327875};
    double ind1[19] = {1.501389186814979,1.942449012417511,2.274636941429680,
                       1.992580970631320,3.042183132272457,0.653825553189122,
                       0.513122411862982,0.426546191045782,0.489374978181903,
                       0.551634460261957,0.376952220630747,0.640936236050839,
                       13.585063650892280,5.589296776771892,3.715287690756511,
                       3.239909510258489,17.220015477006509,1.317604387880688,
                       4.528470958973156};

    Individual ind = Individual(19);
    ind.setParameters(ind0);
    ind.setMaxValues(maxValues);

    appContext ctx{.series = &series};

    //initializeGRNContext(&ctx, 1, numVariables, 5, 7, 7, 0, 49, vectors, maxValues);

    double eval = grnEvaluationLSODATest(ind.getParameters(), &ctx);
    //clearContext2Test(&ctx);

    BOOST_CHECK_CLOSE_FRACTION( eval, 26.92, 0.001 );

}

BOOST_AUTO_TEST_CASE(test_evaluate_GRN5_RK4){
    //void GRNEDOHelpers::initializeGRNContext(appContext* ctx, int granularity, int numVariables, int numTau, int numN, int numK, int setStart, int setEnd, double** vectors, double * maxValues)
    // appContext ctx{};
    GRNSeries series = GRNSeries("GRN5.txt");
    GRNSeries testSeries = GRNSeries(series, 35, 49);
    double **vectors = testSeries.getVectors();
    double *maxValues = testSeries.getMaxValues();
    int numVariables = testSeries.getNumColumns() -1;

    double ind0[19] = {1.2163355099083872, 1.1264485098219865, 2.973714367061704,
                       2.952143123315177, 2.998260518457365, 0.5687249950503857,
                       0.4580723119903261, 0.46214892372246563, 0.6182568295500336,
                       0.5213082492659304, 0.7708877748759901, 0.1497642024548283,
                       4.254757908429968, 3.759370669969996, 4.784173526119725,
                       10.935884810737809, 24.595975874929724, 2.8109199678182635,
                       4.922623602327875};

    Individual ind = Individual(19);
    ind.setParameters(ind0);
    ind.setMaxValues(maxValues);

    appContext ctx{.series = &series};

    //initializeGRNContext(&ctx, 1, numVariables, 5, 7, 7, 0, 49, vectors, maxValues);

    double eval = grnEvaluationRK4Test(ind.getParameters(), &ctx);
    //clearContext2Test(&ctx);

    BOOST_CHECK_CLOSE_FRACTION( eval, 26.92, 0.001 );

}

//todo: fazer o teste abaixo com diferentes conjuntos
BOOST_AUTO_TEST_CASE(test_initialize_GRN10_context_LSODA){
    //void GRNEDOHelpers::initializeGRNContext(appContext* ctx, int granularity, int numVariables, int numTau, int numN, int numK, int setStart, int setEnd, double** vectors, double * maxValues)
    appContext ctx{};
    GRNSeries series = GRNSeries("GRN10.txt");
    double **vectors = series.getVectors();
    double *maxValues = series.getMaxValues();
    int numVariables = series.getNumColumns() -1;


    initializeGRNContext(&ctx, 1, numVariables, 10, 15, 15, 0, 49, vectors, maxValues);

    //esse indivíduo deveria ter fitness ~56
    double ind0[40] = {1.73,2,0.81,0.11, 1.23, 1.78,
                       1.14, 1.04, 3.47, 3.21, 0.45,
                       0.56, 0.99, 0.77, 0.71, 0.66,
                       0.46, 0.48, 0.66, 0.99, 0.85,
                       0.61, 0.55, 0.46, 0.17, 20,
                       9, 24, 12, 2, 2, 6, 4, 7,
                       24, 2, 7, 21, 20, 3};

    double eval = grn10EvaluationLSODA(ind0, &ctx);
    clearContext(&ctx);

    BOOST_CHECK_CLOSE_FRACTION( eval, 56.71, 0.001 );

}



