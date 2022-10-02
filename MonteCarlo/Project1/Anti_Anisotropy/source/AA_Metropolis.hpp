#include "AA_CONf.hpp"
#include "AA_MODEL.hpp"
// #include ""

namespace model::AA{
    class Metropolis : public BIT_MODEL{
        public:
            INT2 total_step;
            INT2 fliped_step;
            void zerostep();
            Metropolis(const model::AA::MODEL_CONF &PROFILE, const ewald_ND &ewd, const myrnd &rand);
            Metropolis(const model::AA::MODEL_CONF &PROFILE);
            void SetTemp();
            void prob();
            void do_step();
            void virtual Measure_fast();
    }
}