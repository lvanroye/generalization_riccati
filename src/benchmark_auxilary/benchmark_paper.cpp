#include "random_problem/random_ocp.hpp"
#include "cart_pendulum_problem/cart_pendulum.hpp"
#include "fatrop_problem_wrap.hpp"
#include <memory>
using namespace genriccati_benchmark;
using namespace std;
int main()
{
    // create a random OCP
    // auto ocp = make_shared<RandomOCP>(100, 10, 4, 10, 0, 10);
    // creat a cart pendulum OCP
    shared_ptr<OCPAbstract> ocp = make_shared<StageOCP>(CartPendulumProblem());
    FatropProblemWrap fatrop_problem_wrap(ocp);
    return 0;
}