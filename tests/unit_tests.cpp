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
using namespace GRNEDOHelpers;

//____________________________________________________________________________//

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

    initializeGRN5Context(&ctx, ctx.TRAINING_MODE, 1);
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

    initializeGRN5Context(&ctx, ctx.TRAINING_MODE, 1);
    double eval = grn5EvaluationRK4(ind0, &ctx);
    clearContext(&ctx);

    BOOST_CHECK_CLOSE_FRACTION( eval, 26.92, 0.001 );

}