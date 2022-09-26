#pragma once
#include <iostream>
#include <time.h>

namespace model{
    template <typename CONF, typename MODEL, typename RND>
    class Model_Driver{
        public:
            CONF PROFILE;
            MODEL MOD;
            RND rnd;
            clock_t __start__, __finish__;

            // Model_Driver intiailizer and destructor;
            Model_Driver(CONF &PROFILE, MODEL &MOD, RND &rnd)
            :PROFILE(PROFILE), MOD(MOD), rnd(rnd)
            {
                this -> __start__ = clock();
                WELCOME();
            };
            ~Model_Driver(){
                FAREWELL(0);
            };

            // Model_Driver run function
            // virtual void run();

            // function to print the model configuration and calculation time when program starts and ends
            void WELCOME(){
                std::cout << "Radnomness test(seed): " << rnd.seed << '\n';
                PROFILE.print_config();
                // MOD.prinf_config();
                std::cout << "------------------------------------------------------------------------------------------------------------------" << "\n";
                std::cout << "--index--||---Temp----||EQ:sig------HH----------||magnetization---specific heat||Fliped Step------Total Step------" << "\n";
                std::cout << "------------------------------------------------------------------------------------------------------------------" << std::endl;
            };
            void FAREWELL(int ERROR_CODE){
                __finish__ = clock();
                if(ERROR_CODE)
                    std::cout << "\nProgram Abonormally Exit. Spent time: " << (double)(__finish__-__start__)/CLOCKS_PER_SEC << "\n";
                else
                    std::cout << "Program Exit. Spent time: " << (double)(__finish__-__start__)/CLOCKS_PER_SEC << "\n";
                std::cout << "-------------------------------------------------------------------------------------------\n";
            };
    };
}