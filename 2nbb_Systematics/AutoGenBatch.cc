/***************************************************
* A script to generate batch scripts from a run list
****************************************************/

#include <iostream>
#include <fstream>
#include <cstdlib>
using namespace std;

int main (int argc, char* argv[]) {
    // READ IN ARGUMENTS/SETTINGS FROM COMMAND (the actual arguments start at index 1)
    if (argc < 2)
    {
        cout << "Error: Too few arguments" << argv[0] << endl;
        cout << "Usage is: ./2nuBB_Systematics DS C CutScheme" << endl;
        return 1;
    }
    int DS = atoi(argv[1]); // DS#
    int c = atoi(argv[2]); // M#
    int cutScheme = atoi(argv[3]); // [raw,raw+DC,raw+DC+AE,raw+DC+AE+DCR]
    int gain = 0; //[0=HG,1=LG]

    char batchFileName[200];
    ofstream myfile;
    char command[250];

    sprintf(batchFileName,"batch2nbbSys_DS%d_M%d_CutScheme%d_Gain%d.sh",DS,c,cutScheme,gain);
    myfile.open(batchFileName);

        myfile << "#!/bin/bash -l\n";

        myfile << "#SBATCH --partition=long\n"; // shared-chos
        myfile << "#SBATCH --time=24:00:00\n";
        myfile << "#SBATCH --ntasks=1\n";
        myfile << "#SBATCH --account=majorana\n";
        myfile << "#SBATCH --mail-type=FAIL\n";
        myfile << "#SBATCH --mail-user=gilliss@unc.edu\n";
        myfile << "#SBATCH --job-name="<<DS<<"_"<<c<<"_"<<cutScheme<<"\n";

        myfile << "./2nuBB_Systematics "<<DS<<" "<<c<<" "<<cutScheme<<endl;

    myfile.close();

    sprintf(command,"sbatch %s",batchFileName);
    system(command);

    return 0;
}
