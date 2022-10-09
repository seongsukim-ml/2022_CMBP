#include "AA_Metropolis_program_header2.hpp"

// arguments list that helps to pass the args to model
vector<string> result_to_file = vector<string>();

int main(int argn, char *argv[]){ // Input argument: argv[0]--> file name / argv[1]--> Input parameter
    signal(SIGSEGV, &handler);
    signal(SIGINT, &handler);
    if(argn >= 2){
        string Input_file = argv[1];
        vector<FLOAT1> input = Writer::Argument_reader(Input_file,13);
        kLx = (INT1)input[0]; kLy = (INT1)input[1]; kN = kLx*kLy;
        kBin = (INT1)input[2]; kB = input[3]; kJx = input[4]; kJy = input[5];
        alpha = input[6]; Tsrt = input[7]; Tfin = input[8];
        isTinf = input[9]; Random = input[10];
        equil_time_base = input[11]; mcs = input[12];
    }
    // arguments list that helps to pass the args to model
    vector<FLOAT1> args = {kLx,kLy,kBin,kB,kJx,kJy,alpha,Tsrt,Tfin,isTinf,Random,equil_time_base,mcs};

    // Filename Base: '\Result\(Model Name)_c_(kL)_int[erval]_(kBin) + (blahblah)
    #ifdef _WIN32
    static string kFilename = ".\\Result\\"+Model::Name()+"_c_"+to_string(kLx)+"_"+to_string(kLy)+"_int"+to_string(kBin)+"_mcs"+to_string(mcs)+"_a"+to_string(alpha);
    #endif
    #ifdef linux
    static string kFilename = "./Result/"+Model::Name()+"_c_"+to_string(kLx)+"_"+to_string(kLy)+"_int"+to_string(kBin)+"_mcs"+to_string(mcs)+"_a"+to_string(alpha);
    #endif

    Greetings();

    Model model = Model(args);
    // CalcProbLookUp();

    FLOAT1 MM, HH;
    FLOAT1 mcs_i = 1/FLOAT1(mcs);
    FLOAT1 kNi   = 1/FLOAT1(kN);
    FLOAT1 kNi2  = kNi*kNi;
    FLOAT1 kNir  = 1/pow(kN,0.5);
    // cout << kNir << '\n';

    for(int cBin = 0; cBin < kBin; cBin++){
        model.Initialize(model.BetaV[cBin]);

        // Output data
        // 0: <m>, 1: <m^2>, 2: <m^4>, 3: <E>/sqrt(N), 4: <E^2>/N
        model.res = vector<FLOAT1>(5,0);

        equil_time = equil_time_base;
        // if(model.TV[i]<=2.4 || model.TV[i]>=2.0) equil_time =1000;

        model.IterateUntilEquilibrium(equil_time);

        model.Measure();
        HH = model.HH;
        MM = model.staggered;

        cout <<"idx: " << left << setw(4) << cBin << "|| " << left << setw(10) << model.TV[cBin];
        cout << "|| "  << left << setw(9) << MM/(FLOAT1)kN << "  " << left << setw(12) << HH;
        cout << "|| "  << left << setw(9) << model.sigma/(FLOAT1)kN << "|| ";
        // for(int i = 0; i < kN; i++){
        //     cout << model.sc[i];
        //     if(i%kLx == kLx-1) cout << '/';
        // }
        cout << '\n';

        /***********Monte Carlo Step and Caculate the data***********/
        int block_size = 10000;
        int blocks = mcs/block_size;
        FLOAT1 bsi = 1/(long double) block_size;

        for(int j = 0; j  < blocks; j++){
            INT2 blocksum_MM    = 0;
            INT2 blocksum_MM2   = 0;
            INT2 blocksum_MM4   = 0;
            FLOAT2 blocksum_HH  = 0;
            FLOAT2 blocksum_HH2 = 0;

            for(int k = 0; k < block_size; k++){
                model.Calculate();
                model.Measure_fast();

                HH = model.HH;
                MM = abs(model.staggered);

                blocksum_MM   += MM;
                blocksum_MM2  += MM*MM;
                blocksum_MM4  += MM*MM*MM*MM;
                blocksum_HH   += HH;
                blocksum_HH2  += HH*HH;
            }
            model.res[0] += blocksum_MM *bsi*kNi;          // = <m>
            model.res[1] += blocksum_MM2*bsi*kNi2;         // = <m^2>
            model.res[2] += blocksum_MM4*bsi*kNi2*kNi2;    // = <m^4>
            model.res[3] += blocksum_HH *bsi*kNir;         // = <E>/sqrt(N)
            model.res[4] += blocksum_HH2*bsi*kNi;          // = <E^2>/N
        }
        for(int j = 0; j < 5; j++)
            model.res[j] /= (FLOAT1) blocks;  // average of blocks
        /***********************************************************/

        /*******Calculate Magnetizaition and Specific Heat.*********/
        model.MV[cBin] = model.res[0];
        model.CV[cBin] = (model.BetaV[cBin]*model.BetaV[cBin])*(model.res[4]-model.res[3]*model.res[3]);
        model.BV[cBin] = 0.5*(3-model.res[2]/(model.res[1]*model.res[1]));
        /***********************************************************/

        cout <<"                         ";
        cout << left << setw(8) << model.MV[cBin] << "  " << left << setw(12) << model.CV[cBin] << "|| ";
        cout << left << setw(14) << model.Fliped_Step << "  " << left << setw(10) << model.Total_Step << endl;

        string result = to_string(cBin) + "," + to_string(model.TV[cBin]) + "," + to_string(model.MV[cBin]) + "," + to_string(model.CV[cBin]) + ",";
        result = result + to_string(model.res[0]) + "," + to_string(model.res[1]) + "," + to_string(model.res[2]) + ",";
        result = result + to_string(model.res[3]) + "," + to_string(model.res[4]) + "\n";

        result_to_file.push_back(result);
    }
    /***********Save the result of the Calculation**********/
    Writer modelW = Writer(kFilename+"_Test_");
    modelW.WriteLine("idx,temperture,(staggered)magnetization,specific heat,abs(mm),mm**2,mm**4,HH/L,HH**2/L\n");
    for(int i = 0; i < kBin; i++)
        modelW.WriteLine(result_to_file.at(i));
    modelW.CloseNewFile();
    /******************************************************/
    Farewell();
}