// matrix and vector classes defined inside their own header files, written by me

#include "vector.h" 
#include "matrix.h" 
#include <iostream>
#include <utility>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>

/* The N-point mass problem can be put in to the model problem form:
 * "dot" denotes differentiation w.r.t time
 * and x is the displacement of a mass w.r.t its equilibrium position (positive to the right)
 *
 * Y' = AY, where Y is the vector [Y1, Y2]^T -> Y1 = [x_1, x_2, x_3,..., x_N]^T and Y2 = [xdot_1, xdot_2,...,xdot_N]^T
 * then Y becomes a 2N-dimensional vector which holds the displacement and velocity information of the i'th mass
 * 1 <= i <= N
 * From forming the equations of motion for the exterior masses (i = 1, N) and interior masses ( 2 <= i <= N)
 * It is clear that the problem to solve is an extension of the tridiagonal matrix problem, now with damping.
 * Note that Y' = [Ydot_1, Ydot_2]^T such that Ydot1 = [xdot_1, xdot_2,...,xdot_N]^T
 * and Ydot2 = [xdotdot_1, xdotdot_2, ... , xdotdot_N]^T
 * The equation of motions in matrix form reveal that A is a 2N x 2N square matrix.
 * Splitting A into quarters further reveals that:
 * - the top-left quadrant is a zero-matrix
 * - the top-right quadrant is an identity matrix
 * - the bottom-left quadrant is a banded tri-diagonal matrix -> B[1/m_i, -2/m_i, 1/m_i] where m_i is i'th mass' mass.
 * - the bottom-right quadrant is a diagonal matrix with the terms -d/m_i where d is the damping coefficient.
 * Then the problem is solvable with a matrix and vector class.
*/
using namespace std;

// Function declarations 
void setParameters(int& s, double& tF, double& tStep, double& k, double& d);
void writeOutput(const Vector solution, double time, bool newFile); 
int getMassCount();

class massSpringSystem {

private:
	int massCount{};
	double springConstant{};
	double dampingCoefficient{};

public:
	massSpringSystem(int N, double k, double d) : massCount(N), springConstant(k), dampingCoefficient(d) {
	}

	void rungeKutta4(const Matrix& systemMatrix, const double h,
	Vector& X0, const double T) const; // timestep, initial condition vector, upper time bound as parameters

	void explicitEuler(const Matrix& system, const double h,
		Vector& X0, const double T) const;

	Vector setInitialState(Vector&, Vector&) const; // takes two empty vectors, returns the initial displacement & velocity vector

	Matrix systemMatrix(Vector&, Vector&) const; // form the associated system matrix

};

Matrix massSpringSystem::systemMatrix(Vector& solutionVector, Vector& massVector) const{ 

	// See lines 12-30 for matrix formation process

	setInitialState(solutionVector, massVector); // initialises the vector of initial conditions

	Matrix topLeft(massCount, massCount); 
	topLeft.zeros(); 

	Matrix topRight(massCount, massCount);
	topRight.Identity();

	Matrix bottomLeft(massCount,  massCount);
	bottomLeft.tridiagonal(springConstant, -2*springConstant, springConstant); // banded matrix coefficients 

	Matrix bottomRight(massCount, massCount);
	
	try {
		bottomRight.diagonal(dampingCoefficient, 0, false); // damping terms
	}
	catch (const logic_error& e) {
		cerr << e.what();
	}

	// dividing each row by the corresponding mass value
	for (int r{ 0 }; r < bottomRight.getSize(0); r++) {
		for (int c{ 0 }; c < bottomRight.getSize(1); c++) {
			bottomLeft.set(r, c, bottomLeft(r, c) / massVector[r]);
			if (!(bottomRight(r, c) == 0.0)) // avoiding negative 0
			{
				bottomRight.set(r, c, bottomRight(r, c) * (-1 / massVector[r])); // damping terms are all negative
			}
		}
	}

	// concatenating the sub-matrices
	try {
		Matrix systemTopHalf(horzCat(topLeft, topRight));
		Matrix systemBottomHalf(horzCat(bottomLeft, bottomRight));
		Matrix system(vertCat(systemTopHalf, systemBottomHalf));

		return system;
	}
	catch (const logic_error& e) {
		cerr << e.what();
	}
	
}

int getMassCount(){

	ifstream parameters; parameters.open("parameters.txt", ios_base::in);

	int massCount{ 0 };
	string lineCounter; 
	if (parameters.is_open()) {
		getline(parameters, lineCounter); // skip the first line
		while (getline(parameters, lineCounter) && !lineCounter.empty()) {
			massCount++;
		}
	}
	parameters.close();
	return massCount;
}

Vector massSpringSystem::setInitialState(Vector& emptyX0, Vector& emptyMass) const {

	ifstream parameters; parameters.open("parameters.txt", ios_base::in);

	string Line;

	if (parameters.is_open()) {
		int j{ 0 }; // initial condition vector tracking
		getline(parameters, Line); // skip the first line
		while (getline(parameters, Line)) { // go through every line
			int i{ 0 }; // tracking the values on each line
			istringstream iss(Line);
			string value;
			
			while (getline(iss, value, ' ')) { // space delimited
				if (i == 0) { 
					if (stod(value) < 0) {
						cout << "Read mass is negative! Check parameters.txt." << endl;
						exit(EXIT_FAILURE);
					} else
					{
						emptyMass.set(j, stod(value));
						i++;
						continue; // skip if value read is mass
					}
					
				}
				else if (i == 1) {
					emptyX0.set(j, stod(value));
					i++;
				}
				else {
					emptyX0.set(j + massCount, stod(value));
					j++;

				}

			}
			
		}
	}
	cout << endl;
	parameters.close();
	return emptyX0;
}

void setParameters(int& s, double& tF, double& tStep, double& k, double& d) {

	// parameters.txt first line: upper time bound, time step, spring constant, damping ratio

	ifstream parameters; parameters.open("parameters.txt", ios_base::in);

	// taking only the first line for parameters
	if (parameters.is_open()) { // file validity
		parameters >> s >> tF >> tStep >> k >> d;
		if (tF <= 0.0) {
			throw invalid_argument("Invalid upper time bound. Check parameters file.");
		}
		else if (tStep <= 0.0) {
			throw invalid_argument("Invalid time step. Check parameters file.");
		}
		else if (k < 0 || d < 0) {
			throw invalid_argument("Negative damping or stiffness not allowed. Check parameters file.");
		}
	}
}

void writeOutput(const Vector solution, double time, bool newFile)
{
	int N{solution.getSize()};

	if (newFile){ // if no file has been made, create one 
		ofstream output; output.open("output.txt"); // new output file
		output.close();
		writeOutput(solution, time, false);
	}
	else {
		ofstream output; output.open("output.txt", ios_base::app | ios_base::out);
		if (output.is_open())
		{
			output << time; // output time every new line
			
			for (int i{0}; i < N/2; i++)
			{
				// first N/2 elements of solution are the displacements of each mass
				// last N/2 elements of the solution are the corresponding velocities 
				
				output << " " << solution[i] << " " << solution[i + N/2]; 
			}
		}
		output << "\n";
		output.close();
	}
}

void massSpringSystem::rungeKutta4(const Matrix& systemMatrix, const double h, Vector& X0, const double T) const {

	double t{ 0 };
	bool newFile {true}; // assuming no output file has been created yet

	writeOutput(X0, t, newFile); // writing the initial conditions first
	newFile = false;

	while (t <= T)
	{
		t += h; // step in time
		Vector k1(systemMatrix * X0 * h);
		Vector k2(systemMatrix * X0 * h + k1/2);
		Vector k3(systemMatrix * X0 * h + k2/2);
		Vector k4(systemMatrix * X0 * h + k3);

		X0 = X0 + k1/6 + (k2+k3)/3 + k4/6;
		
		writeOutput(X0, t, newFile);
		
		
	}
}

void massSpringSystem::explicitEuler(const Matrix& systemMatrix, const double h, Vector& X0, const double T) const {
	double t{ 0 };
	bool newFile{ true };

	writeOutput(X0, t, newFile); // writing the initial conditions first
	newFile = false;

	while (t <= T) {
		t += h;
		X0 = X0 + systemMatrix * X0 * h;
		writeOutput(X0, t, newFile);
		
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
	
	// Assigning the global parameters
	cout << "Assigning parameters to the system..." << endl;
	try {
		setParameters(scheme, timeFinal, timeStep, springConstant, damping);
		cout << "Parameters assigned." << endl;
	}
	catch (const invalid_argument& e) {
		cerr << e.what() << endl;
		exit(EXIT_FAILURE);
	}
	
	int N{ getMassCount() };
	
	
	// creating the system
	cout << "Generating system..." << endl;
	massSpringSystem system(N, springConstant, damping);
	Vector solutionVector(2*N);
	Vector massVector(N); // used in tridiagonal matrix

	// creating the system matrix
	Matrix sysMatrix = system.systemMatrix(solutionVector, massVector);
	cout << solutionVector << endl;


	switch (scheme){
	case 1:
		cout << "Solving with 4th order Runge-Kutta..." << endl;
		system.rungeKutta4(sysMatrix, timeStep, solutionVector, timeFinal);
		cout << "Solution written successfully." << endl;
		exit(EXIT_SUCCESS);
	case 0:
		cout << "Solving with explicit Euler..." << endl;
		system.explicitEuler(sysMatrix, timeStep, solutionVector, timeFinal);
		cout << "Solution written successfully." << endl;
		exit(EXIT_SUCCESS);
	default:
		cout << "Time integration identifier not recognised. Check parameters file." << endl;
		exit(EXIT_FAILURE);
	}
	return 0;
}