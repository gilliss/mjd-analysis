/***************************************************
* A script to generate batch scripts from a run list
****************************************************/

#include <iostream>
#include <fstream>
#include <stdlib.h>
using namespace std;

int main () {

    ////////////////////////////////////
    // Open data
    ////////////////////////////////////
    // system("./2nuBB_Systematics 0 1 0");
    // system("./2nuBB_Systematics 0 1 1");
    // system("./2nuBB_Systematics 0 1 2");
    // system("./2nuBB_Systematics 0 1 3");
    //
    // system("./2nuBB_Systematics 1 1 0");
    // system("./2nuBB_Systematics 1 1 1");
    // system("./2nuBB_Systematics 1 1 2");
    // system("./2nuBB_Systematics 1 1 3");
    //
    // system("./2nuBB_Systematics 2 1 0");
    // system("./2nuBB_Systematics 2 1 1");
    // system("./2nuBB_Systematics 2 1 2");
    // system("./2nuBB_Systematics 2 1 3");
    //
    // system("./2nuBB_Systematics 3 1 0");
    // system("./2nuBB_Systematics 3 1 1");
    // system("./2nuBB_Systematics 3 1 2");
    // system("./2nuBB_Systematics 3 1 3");
    //
    // system("./2nuBB_Systematics 4 2 0");
    // system("./2nuBB_Systematics 4 2 1");
    // system("./2nuBB_Systematics 4 2 2");
    // system("./2nuBB_Systematics 4 2 3");
    //
    // system("./AutoGenBatch 51 1 0");
    // system("./AutoGenBatch 51 1 1");
    // system("./AutoGenBatch 51 1 2");
    // system("./AutoGenBatch 51 1 3");
    //
    // system("./AutoGenBatch 51 2 0");
    // system("./AutoGenBatch 51 2 1");
    // system("./AutoGenBatch 51 2 2");
    // system("./AutoGenBatch 51 2 3");
    //
    // system("./AutoGenBatch 52 1 0");
    // system("./AutoGenBatch 52 1 1");
    // system("./AutoGenBatch 52 1 2");
    // system("./AutoGenBatch 52 1 3");
    //
    // system("./AutoGenBatch 52 2 0");
    // system("./AutoGenBatch 52 2 1");
    // system("./AutoGenBatch 52 2 2");
    // system("./AutoGenBatch 52 2 3");
    //
    // system("./2nuBB_Systematics 53 1 0");
    // system("./2nuBB_Systematics 53 1 1");
    // system("./2nuBB_Systematics 53 1 2");
    // system("./2nuBB_Systematics 53 1 3");
    //
    // system("./2nuBB_Systematics 53 2 0");
    // system("./2nuBB_Systematics 53 2 1");
    // system("./2nuBB_Systematics 53 2 2");
    // system("./2nuBB_Systematics 53 2 3");
    //
    // system("./2nuBB_Systematics 6 1 0");
    // system("./2nuBB_Systematics 6 1 1");
    // system("./2nuBB_Systematics 6 1 2");
    // system("./2nuBB_Systematics 6 1 3");
    //
    // system("./2nuBB_Systematics 6 2 0");
    // system("./2nuBB_Systematics 6 2 1");
    // system("./2nuBB_Systematics 6 2 2");
    // system("./2nuBB_Systematics 6 2 3");

    ////////////////////////////////////
    // Open data: det-by-det
    ////////////////////////////////////
    // system("./2nuBB_RateByDet_BothMs 0 2 1");
    // system("./2nuBB_RateByDet_BothMs 1 2 1");
    // system("./2nuBB_RateByDet_BothMs 2 2 1");
    // system("./2nuBB_RateByDet_BothMs 3 2 1");
    // system("./2nuBB_RateByDet_BothMs 4 2 1");
    // // system("./2nuBB_RateByDet_BothMs 51 2 1");
    // // system("./2nuBB_RateByDet_BothMs 52 2 1");
    // system("./2nuBB_RateByDet_BothMs 53 2 1");
    // system("./2nuBB_RateByDet_BothMs 6 2 1");

    ////////////////////////////////////
    // Blind data
    ////////////////////////////////////
    system("./2nuBB_Systematics 1 1 0");
    system("./2nuBB_Systematics 1 1 1");
    system("./2nuBB_Systematics 1 1 2");
    system("./2nuBB_Systematics 1 1 3");

    system("./2nuBB_Systematics 2 1 0");
    system("./2nuBB_Systematics 2 1 1");
    system("./2nuBB_Systematics 2 1 2");
    system("./2nuBB_Systematics 2 1 3");

    system("./2nuBB_Systematics 53 1 0");
    system("./2nuBB_Systematics 53 1 1");
    system("./2nuBB_Systematics 53 1 2");
    system("./2nuBB_Systematics 53 1 3");

    system("./2nuBB_Systematics 53 2 0");
    system("./2nuBB_Systematics 53 2 1");
    system("./2nuBB_Systematics 53 2 2");
    system("./2nuBB_Systematics 53 2 3");

    system("./2nuBB_Systematics 6 1 0");
    system("./2nuBB_Systematics 6 1 1");
    system("./2nuBB_Systematics 6 1 2");
    system("./2nuBB_Systematics 6 1 3");

    system("./2nuBB_Systematics 6 2 0");
    system("./2nuBB_Systematics 6 2 1");
    system("./2nuBB_Systematics 6 2 2");
    system("./2nuBB_Systematics 6 2 3");

    //////
    // move files
    //////
    system("mv *.txt output_data/blind/");
    system("mv *.root output_data/blind/");

    ////////////////////////////////////
    // Blind data: det-by-det
    ////////////////////////////////////
    system("./2nuBB_RateByDet_BothMs 1 2 1");
    system("./2nuBB_RateByDet_BothMs 2 2 1");
    system("./2nuBB_RateByDet_BothMs 53 2 1");
    system("./2nuBB_RateByDet_BothMs 6 2 1");

    return 0;
}
