
#include "multicut/multicut.h"
#include "visitors/standard_visitor.hxx"
using namespace LP_MP;
int main(int argc, char* argv[])

{
ProblemConstructorRoundingSolver<Solver<FMC_MULTICUT<MessageSendingType::SRMP>,LP,StandardTighteningVisitor>> solver(argc,argv);
solver.ReadProblem(MulticutTextInput::ParseProblem<Solver<FMC_MULTICUT<MessageSendingType::SRMP>,LP,StandardTighteningVisitor>>);
return solver.Solve();

}
