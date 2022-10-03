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
    vector<FLOAT2> args = {kLx,kLy,kBin,kB,kJx,kJy,alpha,Tsrt,Tfin,isTinf,Random,equil_time_base,mcs};

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

        vector<FLOAT1> Binder = vector<FLOAT1>(kBin,0);

        for(int i = 0; i < kBin; i++){
            model.Initialize(model.BetaV[i]);

            // Output data
            // 0: <m>, 1: <m^2>, 2: <m^4>, 3: <E>/sqrt(N), 4: <E^2>/N
            model.res = vector<FLOAT1>(5,0);

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

            INT1 jB = 1000; // Jackknife blcok
            vector<FLOAT1> Jackknife_1 = vector<FLOAT1>(mcs/jB,0);
            vector<FLOAT1> Jackknife_HH = vector<FLOAT1>(mcs/jB,0);
            vector<FLOAT1> Jackknife_HH2 = vector<FLOAT1>(mcs/jB,0);
            vector<FLOAT1> Jackknife_2 = vector<FLOAT1>(mcs/jB,0);
            vector<FLOAT1> Jackknife_MM2 = vector<FLOAT1>(mcs/jB,0);
            vector<FLOAT1> Jackknife_MM4 = vector<FLOAT1>(mcs/jB,0);
            FLOAT1 c_error = 0;
            FLOAT1 b_error = 0;


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

                Jackknife_HH[j/jB] += HH;
                Jackknife_HH2[j/jB] += HH*HH;
                Jackknife_MM2[j/jB] += MM*MM;
                Jackknife_MM4[j/jB] += MM*MM*MM*MM;
            }
            /***********************************************************/

            /*******Calculate Magnetizaition and Specific Heat.*********/
            model.MV[i] = model.res[0];
            model.CV[i] = (model.BetaV[i]*model.BetaV[i])*(model.res[4]-model.res[3]*model.res[3]);
            Binder[i] = 0.5*(3-model.res[2]/model.res[1]/model.res[1]);
            /***********************************************************/

            /*******Calculate Magnetizaition and Specific Heat.*********/
            for(int j = 0; j < mcs/jB; j++){
                Jackknife_1[j] = (model.BetaV[i]*model.BetaV[i])*((mcs*model.res[4]-Jackknife_HH2[j])/(mcs-jB) \
                            - (mcs*model.res[3]-Jackknife_HH[j])/(mcs-jB)*(mcs*model.res[3]-Jackknife_HH[j])/(mcs-jB));
                c_error += (Jackknife_1[j]-model.CV[i])*(Jackknife_1[j]-model.CV[i]);
                FLOAT1 MM2_dif = (mcs*model.res[1]-Jackknife_MM2[j])/(mcs-jB);
                FLOAT1 MM4_dif = (mcs*model.res[2]-Jackknife_MM4[j])/(mcs-jB);
                Jackknife_2[j] = 0.5*(3-(MM4_dif/MM2_dif/MM2_dif));
                b_error += (Jackknife_2[j]-Binder[i])*(Jackknife_2[j]-Binder[i]);
            }
            c_error = sqrt(c_error/(mcs/jB));
            b_error = sqrt(b_error/(mcs/jB));
            /***********************************************************/


            cout <<"                         ";
            cout << left << setw(8) << model.MV[i] << "  " << left << setw(12) << model.CV[i] << "|| ";
            cout << left << setw(14) << model.Fliped_Step << "  " << left << setw(10) << model.Total_Step << endl;

            string result = to_string(i) + "," + to_string(model.TV[i]) + "," + to_string(model.MV[i]) + "," + to_string(model.CV[i]) + ",";
            result = result + to_string(model.res[0]) + "," + to_string(model.res[1]) + "," + to_string(model.res[2]) + ",";
            result = result + to_string(model.res[3]) + "," + to_string(model.res[4]) + ",";
            result = result + to_string(0) + "," to_string(c_error) +  "," + to_string(b_error) + "\n";

            result_to_file.push_back(result);
        }

        /***********Save the result of the Calculation**********/
        Writer modelW = Writer(kFilename+"_Test_");
        modelW.WriteLine("idx,temperture,(staggered)magnetization,specific heat,abs(mm),mm**2,mm**4,HH/L,HH**2/L,MMerr,CCerr,BBerr\n");
        for(int i = 0; i < kBin; i++)
            modelW.WriteLine(result_to_file.at(i));
        modelW.CloseNewFile();
        /******************************************************/
    }
}