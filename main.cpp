// matrix and vector classes defined inside their own header files

#include "vector.h" 
#include "matrix.h" 
#include <iostream>
#include <vector>
#include <string>
#include <fstream>

/* The N-point mass problem can be put in to the model problem form:
 * "dot" denotes differentation w.r.t time
 * and x is the displacement of a mass w.r.t its equilibrium position (positive to the right)
 *
 * Y' = AY, where Y is the vector [Y1, Y2]^T -> Y1 = [x_1, x_2, x_3,..., x_N]^T and Y2 = [xdot_1, xdot_2,...,xdot_N]^T
 * then Y becomes a 2N-dimensional vector which holds the displacement and velocity information of the i'th mass
 * 1 <= i <= N
 * From forming the equations of motion for the exterior masses (i = 1, N) and interior masses ( 2 <= i <= N)
 * It is clear that the problem to solve is an extension of the tri-diagonal matrix problem, now with damping.
 * Note that Y' = [Ydot_1, Ydot_2]^T such that Ydot1 = [xdot_1, xdot_2,...,xdot_N]^T
 * and Ydot2 = [xdotdot_1, xdotdot_2, ... , xdotdot_N]^T
 * The equation of motions in matrix form reveal that A is a 2N x 2N square matrix.
 * Splitting A into quarters further reveals that:
 * - the top-left quadrant is a zero-matrix
 * - the top-right quadrant is an identity matrix
 * - the bottom-left quadrant is a banded tri-diagonal matrix -> B[1/m_i, -2/m_i, 1/m_i] where m_i is i'th mass' mass.
 * - the bottom-right quadrant is a diagonal matrix with the terms -d/m_i where d is the damping coefficient.
 * Then the problem is solveable with a matrix and vector class.
*/
using namespace std;

// Function declarations not belonging to any class
// Reads parameters from parameters.txt into variables declared in main
void setParameters(int& timeIntegrationScheme, double& upperTimeBound, double& timeStep, double& springConstant, double& dampingCoefficient);

class massSpringSystem {
	// 
private:
	int N{}; // number of masses
	int nS{}; // number of springs
	int k{}; // global spring constant
	int d{}; // damping factor applied individually to each mass

	/* a matrix which holds the initial conditions of each mass in each vector
	* the i'th vector in initialConditions will hold the mass, x(0) and xdot(0) of the i'th mass
	*/

	//Vector<Vector<double>> initialConditions;
public: // will need to adapt the Vector class to suit this 
	void rungeKutta4(const Matrix systemMatrix, const int h, int X0, const int T); // timestep, initial condition vector, upper time bound as parameters
	void explicitEuler(const Matrix system, const int h, const int X0, const int T);

};

void setParameters(int& s, double& tF, double& tStep, double& k, double& d) {
	// Parameters: upper time bound, time step, spring constant, damping ratio
	// This function takes 4 system parameters passed by reference and assigns them the correct values from parameter.txt

	ifstream parameters; parameters.open("parameters.txt", ios_base::in);

	if (parameters.is_open()) { // file validity
		while (parameters.good()) {
			// reading the first line 
			parameters >> s >> tF >> tStep >> k >> d;
			break;
			//if (s != 0 || s != 1 || tF <= 0 || tStep <= 0 || k <= 0 || d < 0){
				//cout << "Invalid parameters in file." << endl;
				//exit(EXIT_FAILURE); // exits if the given parameters are invalid
			//} else { break; }
			}
		
		}
	}
	

int main()
{
	// Mass-spring system parameters
	int scheme{};
	double timeFinal{};
	double timeStep{};
	double springConstant{};
	double damping{};
	
	// Assigning the parameters
	setParameters(scheme, timeFinal, timeStep, springConstant, damping);

	cout << scheme << endl  << timeFinal << endl << timeStep << endl << springConstant << endl << damping;

	Matrix topLeft(6);
	
	
	

	return 0;
}