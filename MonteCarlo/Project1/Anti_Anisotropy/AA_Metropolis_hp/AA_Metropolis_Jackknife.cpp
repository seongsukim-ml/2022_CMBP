#include "AA_Metropolis_program_header_hb.hpp"

// arguments list that helps to pass the args to model
vector<vector<FLOAT2>> result_to_file = vector<vector<FLOAT2>>();
// vector<string> result_to_file = vector<string>();
vector<string> result_to_file_cor = vector<string>();
vector<string> result_to_file_cor_err = vector<string>();

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

    FLOAT1 MM, HH, MM_noAbs;
    FLOAT1 mcs_i = 1/FLOAT1(mcs);
    FLOAT1 kNi   = 1/FLOAT1(kN);
    FLOAT1 kNi2  = kNi*kNi;
    FLOAT1 kNir  = 1/pow(kN,0.5);

    for(int cBin = 0; cBin < kBin; cBin++){
        model.Initialize(model.BetaV[cBin]);

        // Output data
        // 0: <m>, 1: <m^2>, 2: <m^4>, 3: <E>/sqrt(N), 4: <E^2>/N
        model.res = vector<FLOAT1>(5,0);
        vector<FLOAT2> cor_short_sum(kN);
        vector<FLOAT2> cor_short_avg(kN);
        vector<FLOAT2> cor_long_sum(kN);
        vector<FLOAT2> cor_long_avg(kN);

        equil_time = equil_time_base;
        // if(model.TV[i]<=2.4 || model.TV[i]>=2.0) equil_time =1000;

        model.IterateUntilEquilibrium(equil_time);

        model.Measure_fast();
        HH = model.HH;
        MM = model.staggered;

        cout <<"idx: " << left << setw(4) << cBin << "|| " << left << setw(10) << model.TV[cBin];
        cout << "|| "  << left << setw(9) << MM/(FLOAT1)kN << "  " << left << setw(12) << HH;
        cout << "|| "  << left << setw(9) << model.sigma/(FLOAT1)kN << "|| ";
        cout << '\n';

        /***********Monte Carlo Step and Caculate the data***********/
        // int block_size = 10000;
        // int blocks = mcs/block_size;
        int blocks = 100;
        int block_size = mcs/blocks;
        FLOAT1 bsi = 1/(long double) block_size;

        MM_noAbs = 0;

        vector<FLOAT2> Block_MM(blocks);
        vector<FLOAT2> Block_MM2(blocks);
        vector<FLOAT2> Block_MM4(blocks);
        vector<FLOAT2> Block_HH(blocks);
        vector<FLOAT2> Block_HH2(blocks);
        vector<FLOAT2> Block_MM_noAbs(blocks);

        for(int bidx = 0; bidx < blocks; bidx++){
            INT2 blocksum_MM    = 0;
            INT2 blocksum_MM2   = 0;
            INT2 blocksum_MM4   = 0;
            FLOAT2 blocksum_HH  = 0;
            FLOAT2 blocksum_HH2 = 0;
            INT2 blocksum_MM_noAbs = 0;

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
                blocksum_MM_noAbs += model.sigma;
            }
            Block_MM[bidx]  = blocksum_MM *bsi*kNi;
            Block_MM2[bidx] = blocksum_MM2*bsi*kNi2;
            Block_MM4[bidx] = blocksum_MM4*bsi*kNi2*kNi2;
            Block_HH[bidx]  = blocksum_HH *bsi*kNir;
            Block_HH2[bidx] = blocksum_HH2*bsi*kNi;

            Block_MM_noAbs[bidx] = blocksum_MM_noAbs*bsi*kNi;

            model.res[0] += Block_MM[bidx];          // = <m>
            model.res[1] += Block_MM2[bidx];         // = <m^2>
            model.res[2] += Block_MM4[bidx];         // = <m^4>
            model.res[3] += Block_HH[bidx];          // = <E>/sqrt(N)
            model.res[4] += Block_HH2[bidx];         // = <E^2>/N
            MM_noAbs += Block_MM_noAbs[bidx];         // = <E^2>/N
        }

        for(int j = 0; j < 5; j++)
            model.res[j] /= (FLOAT1) blocks;  // average of blocks
        FLOAT1 MM_noAbsf =  MM_noAbs/(FLOAT1) blocks;
        /***********************************************************/

        /*******Calculate Magnetizaition and Specific Heat.*********/
        model.MV[cBin] = model.res[0];
        model.CV[cBin] = (model.BetaV[cBin]*model.BetaV[cBin])*(model.res[4]-model.res[3]*model.res[3]);
        model.BV[cBin] = 0.5*(3-model.res[2]/(model.res[1]*model.res[1]));
        /***********************************************************/

        /*******Calculate Jackknife Error.*********/
        FLOAT2 MMerr = 0;
        FLOAT2 CCerr = 0;
        FLOAT2 BBerr = 0;
        FLOAT2 MM2err = 0;
        FLOAT2 MM4err = 0;

        for(int bidx = 0; bidx < blocks; bidx++){
            FLOAT2 MMtemp  = (model.MV[cBin]*blocks-Block_MM[bidx])/(blocks-1);
            MMerr += (model.MV[cBin]-MMtemp)*(model.MV[cBin]-MMtemp);

            FLOAT2 MM2temp = (model.res[1]*blocks-Block_MM2[bidx])/(blocks-1);
            MM2err += (MM2temp-model.res[1])*(MM2temp-model.res[1]);

            FLOAT2 MM4temp = (model.res[2]*blocks-Block_MM4[bidx])/(blocks-1);
            MM4err += (MM4temp-model.res[2])*(MM4temp-model.res[2]);

            FLOAT2 HHtemp  = (model.res[3]*blocks-Block_HH[bidx])/(blocks-1);
            FLOAT2 HH2temp = (model.res[4]*blocks-Block_HH2[bidx])/(blocks-1);

            FLOAT2 CCtemp = model.BetaV[cBin]*model.BetaV[cBin]*(HH2temp-HHtemp*HHtemp);
            CCerr += (CCtemp-model.CV[cBin])*(CCtemp-model.CV[cBin]);

            FLOAT2 BBtemp = 0.5*(3-MM4temp/(MM2temp*MM2temp));
            BBerr += (BBtemp-model.BV[cBin])*(BBtemp-model.BV[cBin]);
        }
        MMerr  = pow(MMerr,0.5);
        MM2err = pow(MM2err,0.5);
        MM4err = pow(MM4err,0.5);
        CCerr  = pow(CCerr,0.5);
        BBerr  = pow(BBerr,0.5);
        /***********************************************************/

        cout <<"                         ";
        cout << left << setw(8) << model.MV[cBin] << "  " << left << setw(12) << model.CV[cBin] << "|| ";
        cout << left << setw(14) << model.Fliped_Step << "  " << left << setw(10) << model.Total_Step << endl;

        // Standard result save
        // string result = to_string(cBin) + "," + to_string(model.TV[cBin]) + "," + to_string(model.MV[cBin]) + "," + to_string(model.CV[cBin]) + ",";
        // result = result + to_string(model.res[0]) + "," + to_string(model.res[1]) + "," + to_string(model.res[2]) + ",";
        // result = result + to_string(model.res[3]) + "," + to_string(model.res[4]) + "," + to_string(model.BV[cBin]) + ",";
        // result = result + to_string(MMerr) + "," + to_string(CCerr) + "," + to_string(BBerr) + ",";
        // result = result + to_string(MM2err) + "," + to_string(MM4err) + "\n";

        vector<FLOAT2> result = {
            FLOAT1(cBin),model.TV[cBin],model.MV[cBin],model.CV[cBin],
            model.res[0],model.res[1],model.res[2],model.res[3],model.res[4],
            model.BV[cBin], MMerr, CCerr, BBerr, MM2err, MM4err
            };

        result_to_file.push_back(result);
    }

    kFilename += "_Jackknife";
    /***********Save the result of the Calculation**********/
    Writer modelW = Writer(kFilename);
    // vector<string> ={"idx","temperature","stag_mag",}
    modelW.WriteLine("idx,temperture,(staggered)magnetization,specific heat,abs(mm),mm**2,mm**4,HH/L,HH**2/L,Binder,MMerr,CCerr,BBerr,MM2err,MM4err\n");
    for(int i = 0; i < kBin; i++)
        modelW.WriteLine(result_to_file.at(i));
    modelW.CloseNewFile();
    /******************************************************/
    Farewell();
}