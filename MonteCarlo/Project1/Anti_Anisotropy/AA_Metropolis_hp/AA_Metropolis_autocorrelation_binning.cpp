#include "AA_Metropolis_program_header_hb.hpp"

// arguments list that helps to pass the args to model
vector<string> result_to_file = vector<string>();
vector<string> result_to_file_auto_MM = vector<string>();
vector<string> result_to_file_auto_CC = vector<string>();

template<typename T1>
FLOAT2 compact_vector(T1 &v, FLOAT2 A, long size, long total_length){
    if(size == 1){
        return 0;
    }
    else{
        FLOAT2 res = 0;
        for(int i = 0; i < size/2; i++){
            res += (A-v[2*i])*(A-v[2*i]);
            res += (A-v[2*i+1])*(A-v[2*i+1]);
            v[i] = (v[2*i] + v[2*i + 1])/2;
        }
        return res/size/(size-1);
    }
}

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
    // equil_time = 0;
    // kBin = 1;
    // Tsrt = 2.269;
    // mcs = 16777216;

    // arguments list that helps to pass the args to model
    vector<FLOAT1> args = {kLx,kLy,kBin,kB,kJx,kJy,alpha,Tsrt,Tfin,isTinf,Random,equil_time_base,mcs};

    // Filename Base: '\Result\(Model Name)_c_(kL)_int[erval]_(kBin) + (blahblah)
    #ifdef _WIN32
    static string kFilename = ".\\Result\\auto\\"+Model::Name()+"_c_"+to_string(kLx)+"_"+to_string(kLy)+"_int"+to_string(kBin)+"_mcs"+to_string(mcs)+"_a"+to_string(alpha);
    #endif
    #ifdef linux
    static string kFilename = "./Result/auto/"+Model::Name()+"_c_"+to_string(kLx)+"_"+to_string(kLy)+"_int"+to_string(kBin)+"_mcs"+to_string(mcs)+"_a"+to_string(alpha);
    #endif

    Greetings();

    Model model = Model(args);

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

        model.Measure_fast();
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
        int block_size = mcs;
        int blocks = mcs/block_size;
        FLOAT1 bsi = 1/(long double) block_size;

        // Autocorrelation data
        vector<FLOAT2> MMauto(mcs);
        vector<FLOAT2> CCauto(mcs);
        vector<FLOAT2> HHauto(mcs);
        vector<FLOAT2> HH2auto(mcs);
        int bLayer = log2(mcs);
        vector<FLOAT2> MMauto_res(bLayer+1);
        vector<FLOAT2> CCauto_res(bLayer+1);

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
                MMauto[(j*block_size + k)]  = MM;
                HHauto[(j*block_size + k)]  = HH;
                HH2auto[(j*block_size + k)] = HH*HH;
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
        int binning_size = mcs;
        for(int j = 0; j <= bLayer; j++){
            FLOAT2 MMtemp = compact_vector<vector<FLOAT2>>(MMauto,model.res[0]*kN,binning_size,mcs);
            MMtemp *= kNi2;
            MMauto_res[j] = pow(MMtemp,0.5);

            for(int i = 0; i < binning_size; i++){
                CCauto[i] = model.BetaV[cBin]*model.BetaV[cBin]*(HH2auto[i]-HHauto[i]*HHauto[i])*kNi;
                CCauto_res[j] += (model.CV[cBin]-CCauto[i])*(model.CV[cBin]-CCauto[i]);
            }
            FLOAT2 HHtemp = compact_vector<vector<FLOAT2>>(HHauto,model.res[3],binning_size,mcs);
            FLOAT2 HH2temp = compact_vector<vector<FLOAT2>>(HH2auto,model.res[4],binning_size,mcs);
            CCauto_res[j] /= binning_size*(binning_size-1);
            CCauto_res[j] = pow(CCauto_res[j],0.5);

            binning_size /= 2;
        }

        cout <<"                         ";
        cout << left << setw(8) << model.MV[cBin] << "  " << left << setw(12) << model.CV[cBin] << "|| ";
        cout << left << setw(14) << model.Fliped_Step << "  " << left << setw(10) << model.Total_Step << endl;

        string result = to_string(cBin) + "," + to_string(model.TV[cBin]) + "," + to_string(model.MV[cBin]) + "," + to_string(model.CV[cBin]) + ",";
        result = result + to_string(model.res[0]) + "," + to_string(model.res[1]) + "," + to_string(model.res[2]) + ",";
        result = result + to_string(model.res[3]) + "," + to_string(model.res[4]) + "\n";

        result_to_file.push_back(result);

        string result2 = to_string(cBin) + "," + to_string(model.TV[cBin]) + "," + to_string(model.MV[cBin]) + "," + to_string(model.CV[cBin]) + ",";
        for(int j = 0; j < bLayer; j++){
            result2 += to_string(MMauto_res[j]) + ",";
        }
        result2.pop_back();
        result_to_file_auto_MM.push_back(result2 + "\n");

        string result3 = to_string(cBin) + "," + to_string(model.TV[cBin]) + "," + to_string(model.MV[cBin]) + "," + to_string(model.CV[cBin]) + ",";
        for(int j = 0; j < bLayer; j++){
            result3 += to_string(CCauto_res[j]) + ",";
        }
        result3.pop_back();
        result_to_file_auto_CC.push_back(result3 + "\n");
    }
    kFilename += "_auto";
    /***********Save the result of the Calculation**********/
    Writer modelW = Writer(kFilename);
    modelW.WriteLine("idx,temperture,(staggered)magnetization,specific heat,abs(mm),mm**2,mm**4,HH/L,HH**2/L\n");
    for(int i = 0; i < kBin; i++)
        modelW.WriteLine(result_to_file.at(i));
    modelW.CloseNewFile();
    /******************************************************/

    /***********Save the result of the Autocorrelation**********/
    Writer modelW2 = Writer(kFilename+"_MM");
    modelW2.WriteLine("idx,temperture,MM,CC,\n");
    for(int i = 0; i < kBin; i++)
        modelW2.WriteLine(result_to_file_auto_MM.at(i));
    modelW2.CloseNewFile();
    /******************************************************/

    /***********Save the result of the Autocorrelation**********/
    Writer modelW3 = Writer(kFilename+"_CC");
    modelW3.WriteLine("idx,temperture,MM,CC,\n");
    for(int i = 0; i < kBin; i++)
        modelW3.WriteLine(result_to_file_auto_CC.at(i));
    modelW3.CloseNewFile();
    /******************************************************/
    Farewell();
}