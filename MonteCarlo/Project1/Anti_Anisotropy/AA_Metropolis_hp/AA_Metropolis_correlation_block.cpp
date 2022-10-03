#include "AA_Metropolis_program_header.hpp"

// arguments list that helps to pass the args to model
vector<string> result_to_file = vector<string>();
vector<string> result_to_file_cor = vector<string>();

int main(int argn, char *argv[]){ // Input argument: argv[0]--> file name / argv[1]--> Input parameter
    signal(SIGSEGV, &handler);
    signal(SIGINT, &handler);
    if(argn >= 2){
        string Input_file = argv[1];
        vector<FLOAT1> input = Writer::Argument_readerL(Input_file,13);
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

    for(int gg = 0; gg < 1; gg++){
        Model model = Model(args);

        FLOAT1 MM, HH;
        FLOAT1 mcs_i = 1/FLOAT1(mcs);
        FLOAT1 kNi  = 1/FLOAT1(kN);
        FLOAT1 kNir = pow(kN,0.5);
        cout << kNir << '\n';

        for(int i = 0; i < kBin; i++){
            model.Initialize(model.BetaV[i]);

            // Output data
            // 0: <m>, 1: <m^2>, 2: <m^4>, 3: <E>/sqrt(N), 4: <E^2>/N
            model.res = vector<FLOAT1>(5,0);
            vector<INT2> cor_sum(kN);
            vector<FLOAT2> cor_avg(kN);

            equil_time = equil_time_base;
            // if(model.TV[i]<=2.4 || model.TV[i]>=2.0) equil_time =1000;

            model.IterateUntilEquilibrium(equil_time);

            model.Measure_fast();
            HH = model.HH;
            MM = model.staggered;

            cout <<"idx: " << left << setw(4) << i << "|| " << left << setw(10) << model.TV[i];
            cout << "|| "  << left << setw(9) << MM/(FLOAT1)kN << "  " << left << setw(12) << HH;
            cout << "|| "  << left << setw(9) << model.sigma/(FLOAT1)kN << "|| ";
            // for(int i = 0; i < kN; i++){
            //     cout << model.sc[i];
            //     if(i%kLx == kLx-1) cout << '/';
            // }
            cout << '\n';

            // Blocking method
            int block_size = 1000;
            int blocks = mcs/block_size;
            vector<FLOAT2> Block_MM(blocks);
            vector<FLOAT2> Block_MM2(blocks);
            vector<FLOAT2> Block_MM4(blocks);
            vector<FLOAT2> Block_HH(blocks);
            vector<FLOAT2> Block_HH2(blocks);
            // vector<FLOAT2> Block_CC(blocks);
            // vector<FLOAT2> Block_BB(blocks);

            /***********Monte Carlo Step and Caculate the data***********/
            for(int j = 0; j < mcs; j++){
                model.Calculate();                     //O(N^2)

                model.Measure_fast();
                HH = model.HH/(FLOAT1) kNir;        // = E
                MM = abs(model.staggered)/(FLOAT1)kN;    // = |M|

                model.res[0] += MM*mcs_i;              // = <m>
                model.res[1] += (MM*mcs_i*MM);         // = <m^2>
                model.res[2] += (MM*mcs_i*MM)*(MM*MM); // = <m^4>
                model.res[3] += HH*mcs_i;              // = <E>/sqrt(N)
                model.res[4] += HH*mcs_i*HH;           // = <E^2>/N

                int bidx = j/block_size;
                long double bs_i = 1/(long double) block_size;
                Block_MM[bidx] += MM*bs_i;
                Block_MM2[bidx] += MM*bs_i*MM;
                Block_MM4[bidx] += MM*bs_i*MM*MM*MM;
                Block_HH[bidx] += HH*bs_i;
                Block_HH2[bidx] += HH*bs_i*HH;

                model.CalculateCorrelation();
                // model.CalculateCorrelationStaggered();
                for(int ii = 0; ii < kN; ii++)
                    cor_sum[ii] += model.cor[ii];
                    // cor_sum[ii] += model.cor_stag[ii];
            }
            /***********************************************************/

            /*******Calculate Magnetizaition and Specific Heat.*********/
            model.MV[i] = model.res[0];
            model.CV[i] = (model.BetaV[i]*model.BetaV[i])*(model.res[4]-model.res[3]*model.res[3]);
            model.BV[i] = 0.5*(3-model.res[4]/(model.res[3]*model.res[3]));
            /***********************************************************/

            /*******Calculate Magnetizaition and Specific Heat.*********/
            for(int j = 0; j < kN; j++)
                cor_avg[j] = (FLOAT2)(2*cor_sum[j]-mcs)/mcs;
            /***********************************************************/

            FLOAT2 MMerr = 0;
            FLOAT2 CCerr = 0;
            FLOAT2 BBerr = 0;
            for(int j = 0; j < blocks; j++){
                MMerr += (Block_MM[j]-model.MV[i])*(Block_MM[j]-model.MV[i]);
                FLOAT2 CCtemp = model.BetaV[i]*model.BetaV[i]*(Block_HH2[j]-Block_HH[j]*Block_HH[j]);
                CCerr += (CCtemp-model.CV[i])*(CCtemp-model.CV[i]);
                FLOAT2 BBtemp = 0.5*(3-Block_MM4[j]/(Block_MM2[j]*Block_MM2[j]));
                BBerr += BBtemp;
            }
            cout <<"                         ";
            cout << left << setw(8) << model.MV[i] << "  " << left << setw(12) << model.CV[i] << "|| ";
            cout << left << setw(14) << model.Fliped_Step << "  " << left << setw(10) << model.Total_Step << endl;

            string result = to_string(i) + "," + to_string(model.TV[i]) + "," + to_string(model.MV[i]) + "," + to_string(model.CV[i]) + ",";
            result = result + to_string(model.res[0]) + "," + to_string(model.res[1]) + "," + to_string(model.res[2]) + ",";
            result = result + to_string(model.res[3]) + "," + to_string(model.res[4]) + ",";
            result = result + to_string(MMerr) + "," + to_string(CCerr) + "," + to_string(BBerr) + "\n";

            result_to_file.push_back(result);

            // correlation result save
            string result2 = to_string(i) + "," + to_string(model.TV[i]) + ",";
            for(int j = 0; j < kN; j++){
                result2 = result2 + to_string(cor_avg[j]) + ",";
            }
            result2.pop_back();
            result_to_file_cor.push_back(result2);
        }
        kFilename += "_Block_";
        /***********Save the result of the Calculation**********/
        Writer modelW = Writer(kFilename);
        modelW.WriteLine("idx,temperture,(staggered)magnetization,specific heat,abs(mm),mm**2,mm**4,HH/L,HH**2/L,MMerr,CCerr,BBerr\n");
        for(int i = 0; i < kBin; i++)
            modelW.WriteLine(result_to_file.at(i));
        modelW.CloseNewFile();
        /******************************************************/

        /***********Save the result of the Correlation**********/
        Writer modelW2 = Writer(kFilename+"_Cor");
        modelW2.WriteLine("idx,temperture,(staggered)magnetization,specific heat,abs(mm),mm**2,mm**4,HH/L,HH**2/L\n");
        for(int i = 0; i < kBin; i++)
            modelW2.WriteLine(result_to_file_cor.at(i));
        modelW2.CloseNewFile();
        /******************************************************/
    }
    Farewell();
}