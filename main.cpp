#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <math.h>
#include <iostream>


double clamp(double value, double minVal, double maxVal) {
	if (value < minVal) return minVal;
	if (value > maxVal) return maxVal;
	return value;
}

/// DIVERGENCE DEBUGGING
bool showDivergence = false;
GLuint fluidProgram, divergenceProgram;
GLuint fluidTexture, divergenceTexture;


GLuint createShaderProgram(const char* vertSrc, const char* fragSrc) {
	GLuint vertShader = glCreateShader(GL_VERTEX_SHADER);
	glShaderSource(vertShader, 1, &vertSrc, NULL);
	glCompileShader(vertShader);

	GLuint fragShader = glCreateShader(GL_FRAGMENT_SHADER);
	glShaderSource(fragShader, 1, &fragSrc, NULL);
	glCompileShader(fragShader);

	GLuint program = glCreateProgram();
	glAttachShader(program, vertShader);
	glAttachShader(program, fragShader);
	glLinkProgram(program);

	glDeleteShader(vertShader);
	glDeleteShader(fragShader);

	return program;
}



double gravity = -0.9;

double prevTime = 0.0f;
int frames = 0;
int watercells = 0;

const float omega = 1.9f;

const double fldthreld = 0.1;

const double targetFrameTime = 1.0 / 30.0;  // 30 FPS = 0.033 seconds per frame

const int cells = 140;
const double dx = 2.0 / cells;
//float center = cells * stepsize / 2 - 0.75f;
double p[cells][cells] = { 0.0 };
double r[cells][cells] = { 0.0 };
double(*u)[cells] = new double[cells + 1][cells];

double(*v)[cells + 1] = new double[cells][cells + 1];

float(*fluid)[cells] = new float[cells][cells];
double divergence[cells][cells] = { 0.0 };
double vorticity[cells][cells] = { 0.0 };
double(*newV)[cells + 1] = new double[cells][cells + 1];
double(*newU)[cells] = new double[cells + 1][cells];

double precon[cells][cells] = { 0.0 };
double z[cells][cells] = { 0.0 };
double s[cells][cells] = { 0.0 };
double q[cells][cells] = { 0.0 };

double PrevTime = 0.f;

// v^2 = 0^2 + 2(-g)(-S) where S = - 2cells*dx/3
//float UMax = pow(2 * fabs(gravity) * 2.f * (screensize) / 3.f, 0.5f);

/// SOME CALCULATIONS: UMax approx 54.7, so calculating dx at dt = 0.001, apparently dx approx 1!
/// tried g = 1000, UMax 44.7, so with dt = 0.001, dx ~ 0.8
/// SO dt =  is required... 

bool is_fluid(int i, int j) {
	return fluid[i][j] == 1.f; // Assuming 1=fluid, 0=air
}

void mic0();

// Kind of beginplay
void initialize() {

	for (int a = 0; a < cells; a++)
		for (int b = 0; b < cells; b++) {
			if (fabs(a-cells/2)<cells/6  and b < cells / 2 and b != 0 and a != 0 and a != cells - 1 and b != cells - 1) {
				fluid[a][b] = 1.f;
			}
			else
			{
				fluid[a][b] = 0.f;
			}
		}

	for (int a = 0; a <= cells; a++)
		for (int b = 0; b <= cells; b++) {
			if (a != cells) {
				v[a][b] = 0.0;
				newV[a][b] = 0.0;
			}
			if (b != cells) {
				u[a][b] = 0.0;
				newU[a][b] = 0.0;
			}
		}

//	mic0();
}

void CalculateDivergence(int i, int j) {
//	if (fluid[i][j] == 1.0) {
		if (i == 0 or i == cells - 1 or j == 0 or j == cells - 1) {
			if (i != j and i + j != cells - 1) {
				if (i == 0) divergence[i][j] = u[1][j];
				else if (j == 0) divergence[i][j] = v[i][1];
				else if (i == cells - 1) divergence[i][j] = -u[i][j];
				else if (j == cells - 1) divergence[i][j] = -v[i][j];
			}
		}

		else if (i == 1 or i == cells - 2 or j == 1 or j == cells - 2) {
			double du_dx = (u[i + 1][j] * (i + 1 == cells - 1 ? 0.0 : 1.0)
				- u[i][j] * (i == 1 ? 0.0 : 1.0));

			double dv_dy = (v[i][j + 1] * (j + 1 == cells - 1 ? 0.0 : 1.0)
				- v[i][j]  *(j == 1 ? 0.0 : 1.0));

			divergence[i][j] = du_dx + dv_dy;
		}

		else {
			double du_dx = (u[i + 1][j] - u[i][j]);
			double dv_dy = (v[i][j + 1] - v[i][j]);
			divergence[i][j] = du_dx + dv_dy;
		}
//	}
//	else
	//	divergence[i][j] = 0.0;
}



double CalculateNewPressure(int i, int j) {

	int NonSolidCells = 0;
	double new_p = 0.0;

	if (i + 1 != cells - 1) {
	
		NonSolidCells++;
		if(fluid[i + 1][j] == 1.f) new_p += p[i + 1][j];
	}
	if (j + 1 != cells - 1) {
		NonSolidCells++;
		if (fluid[i][j + 1] == 1.f) new_p += p[i][j + 1];
	}
	if (i - 1 != 0) {
		NonSolidCells++;
		if (fluid[i - 1][j] == 1.f) new_p += p[i - 1][j];
	}
	if (j - 1 != 0) {
		NonSolidCells++;
		if (fluid[i][j - 1] == 1.f) new_p += p[i][j - 1];
	}

	new_p -= divergence[i][j];

	return new_p / NonSolidCells;
}

void samplevelocity(double px, double py, double* pu, double* pv) {

	// Clamp and interpolate
	px = clamp(px, 1.0, cells - 1.0); // Avoid sampling boundaries
	py = clamp(py, 1.0, cells - 1.0);

	int i0 = clamp((int)px, 1, cells - 1);
	int j0 = clamp((int)py, 1, cells - 1);
	int i1 = i0 + 1 > cells - 1 ? cells - 1 : i0 + 1;
	int j1 = j0 + 1 > cells - 1 ? cells - 1 : j0 + 1;
	int im1 = (i0 - 1) < 1 ? 1 : i0 - 1;
	int jm1 = (j0 - 1) < 1 ? 1 : j0 - 1;

	double a = px - i0;
	double b = py - j0;

	
		*pu = (1 - a) * (1 - b) * (u[i0][j0] + u[i0][jm1]) / 2 +
			(1 - a) * b * (u[i0][j1] + u[i0][j0]) / 2 +
			a * (1 - b) * (u[i1][j0] + u[i1][jm1]) / 2 +
			a * b * (u[i1][j1] + u[i1][j0]) / 2;
	

	
		*pv = (1 - a) * (1 - b) * (v[i0][j0] + v[im1][j0]) / 2 +
			(1 - a) * b * (v[i0][j1] + v[im1][j1]) / 2 +
			a * (1 - b) * (v[i0][j0] + v[im1][j0]) / 2 +
			a * b * (v[im1][jm1] + v[i0][jm1]) / 2;
}

int Adiag(int i, int j) {
	if (i > 1 and i < cells - 2 and j > 1 and j < cells - 2) {
		return 4;
	}
	else if(i == 0 or j == 0 or i == cells - 1 or j == cells - 1) {
		return 1;
	} 
	else{

		if (i == j or i + j == cells - 1) 
				return 2;
			
			else return 3;
		
	}
}

int Aplusi(int i, int j) {
	if (j == 0 or j == cells - 1)	return 0;

	if (i < cells - 2)	 return -1;
	
	else if (i == cells - 2) return 0;

	else return 0;
}

int Aplusj(int i, int j) {
	if (i == 0 or i == cells - 1)	return 0;

	if (j < cells - 2)	 return -1;

	else if (j == cells - 2) return 0;

	else return 0;
}

void mic0() {
	double tuning = 0.97;
	for(int i = 1; i < cells; i++)
		for (int j = 1; j < cells; j++) {
			double e = Adiag(i, j) 
				- pow(Aplusi(i - 1, j) * precon[i - 1][j], 2.0)
				- pow(Aplusj(i, j - 1) * precon[i][j - 1], 2.0)
				- tuning * (Aplusi(i - 1, j) * Aplusj(i - 1, j) * precon[i - 1][j] * precon[i - 1][j]
					+ Aplusj(i, j - 1) * Aplusi(i, j - 1) * precon[i][j - 1] * precon[i][j - 1]);
			precon[i][j] = 1 / sqrt(e + 1e-30);
		}
	printf("Calculated the precon\n");
}

void applyprecon() {
	// Solving Az = r
	// Holding A = LL^T
	// First solve Lq = r
	for(int i = 1; i < cells; i++)
		for (int j = 1; j < cells; j++) {
			double t = r[i][j] - Aplusi(i - 1, j) * precon[i - 1][j] * q[i - 1][j]
							   - Aplusj(i, j - 1) * precon[i][j - 1] * q[i][j - 1];
			q[i][j] = t * precon[i][j];
		}
	// Next Solve L^T z = q
	for(int i = cells - 2; i >= 0; i--)
		for (int j = cells - 2; j >= 0; j--) {
			double t = q[i][j] - Aplusi(i, j) * precon[i][j] * z[i + 1][j]
							   - Aplusj(i, j) * precon[i][j] * z[i][j + 1];
			z[i][j] = t * precon[i][j];
		}
//	printf("Applied the precon\n");
}

double dotproduct(char rOrs) {
	double dp = 0.0;
	if (rOrs == 'r') {
		for (int i = 0; i < cells; i++)
			for (int j = 0; j < cells; j++) {
				dp += z[i][j] * r[i][j];
			}
	}
	else if (rOrs == 's') {
		for (int i = 0; i < cells; i++)
			for (int j = 0; j < cells; j++) {
				dp += z[i][j] * s[i][j];
			}
	}
	return dp;
}

double A(int i, int j, int k, int l) {
	if (i == k and j == l)	return Adiag(i, j);
	else if (i + 1 == k and j == l) return Aplusi(i, j);
	else if (i - 1 == k and j == l) return Aplusi(i - 1, j);
	else if (i == k and j + 1 == l) return Aplusj(i, j);
	else if (i == k and j - 1 == l) return Aplusj(i, j - 1);
	else return 0.0;
}

void applyA() {
	// Iterate over each cell (i,j) - equation index
	for (int i = 0; i < cells; i++) {
		for (int j = 0; j < cells; j++) {
			// Process center point (i,j)
			z[i][j] += A(i, j, i, j) * s[i][j];

			// Process neighbor points (if within bounds)
			// Top neighbor (i-1, j)
			if (i > 0) z[i - 1][j] += A(i, j, i - 1, j) * s[i - 1][j];
			// Bottom neighbor (i+1, j)
			if (i < cells - 1) z[i + 1][j] += A(i, j, i + 1, j) * s[i + 1][j];
			// Left neighbor (i, j-1)
			if (j > 0) z[i][j - 1] += A(i, j, i, j - 1) * s[i][j - 1];
			// Right neighbor (i, j+1)
			if (j < cells - 1) z[i][j + 1] += A(i, j, i, j + 1) * s[i][j + 1];
		}
	}
}


void Tick(double dt) {

	double dvG = gravity * dt;
	double dtbydx = dt / dx;

//	printf("Time = %lf\n", 1/dt);
//	watercells = 0;
//	frames++;

	// Water loss calculations
/*	for (int i = 0; i < cells; i++)
		for (int j = 0; j < cells; j++){
			if (fluid[i][j] == 1.0) watercells++;
		}
	printf("x=%d, y=%d\n", frames, watercells);*/




	float(*temp_fluid)[cells] = new float[cells][cells];
	for (int a = 0; a < cells; a++)
		for (int b = 0; b < cells; b++)
			temp_fluid[a][b] = 0.f;

	// Advecting the fluid itself
#pragma omp parallel for 
	for (int idx = 0; idx < cells * cells; idx++)
	{
		int i = idx / cells;
		int j = idx % cells;



		if (i > 0 and j > 0 and i < cells - 1 and j < cells - 1) {

			double x = i + 0.50;
			double y = j + 0.50;

			// Forward Euler
		/*	float vx = (u[i][j] + u[i + 1][j]) * 0.5f;
			float vy = (v[i][j] + v[i][j + 1]) * 0.5f;

			float prevX = x - vx * dtbydx;
			float prevY = y - vy * dtbydx;*/


			// RK4
			// Stage 1 : Velocity at starting point (x, y)
			double k1x = (u[i][j] + u[i + 1][j]) * 0.5f;
			double k1y = (v[i][j] + v[i][j + 1]) * 0.5f;

			// Stage 2: Velocity at first midpoint (t - dt/2)
			double x2 = x - 0.5f * dtbydx * k1x;
			double y2 = y - 0.5f * dtbydx * k1y;
			double k2x, k2y;
			samplevelocity(x2, y2, &k2x, &k2y);

			// Stage 3: Velocity at second midpoint (t - dt/2)
			double x3 = x - 0.5f * dtbydx * k2x;
			double y3 = y - 0.5f * dtbydx * k2y;
			double k3x, k3y;
			samplevelocity(x3, y3, &k3x, &k3y);

			// Stage 4: Velocity at endpoint (t - dt)
			double x4 = x - dtbydx * k3x;
			double y4 = y - dtbydx * k3y;
			double k4x, k4y;
			samplevelocity(x4, y4, &k4x, &k4y);

			// Combine velocities with RK4 weights
			double vx = (k1x + 2 * k2x + 2 * k3x + k4x) / 6.0f;
			double vy = (k1y + 2 * k2y + 2 * k3y + k4y) / 6.0f;

			// Trace back full step using combined velocity
			double prevX = x - dtbydx * vx;
			double prevY = y - dtbydx * vy;

			prevX = clamp(prevX, 1.0, cells - 1); // Clamping the prevx between 0 and cells - 1
			prevY = clamp(prevY, 1.0, cells - 1);

			int pi = (int)(prevX);
			int pj = (int)(prevY);

		//	if(fluid[pi][pj] == 1.0 and fluid[i][j] == 0.0)
			temp_fluid[i][j] = fluid[pi][pj];
		}
	}

	auto f = temp_fluid;
	temp_fluid = fluid;
	fluid = f;
	delete[] temp_fluid;

//	printf("Advected the fluid!\n");

			// Advect u and v
#pragma omp parallel for 
	for (int idx = 0; idx <= cells * cells; ++idx) {
		int i = idx / cells;
		int j = idx % cells;

		// --- u‐advection
		if (i > 1 && i < cells - 1 and j > 0 and j < cells - 1) {

			///  why fluid[i-1][j] != 0 ?
		//	if (fluid[i][j] != 0.f or fluid[i - 1][j] != 0.f) {
			double x = i;
			double y = j + 0.5;

			// Stage 1: Velocity at starting point (x, y)
			double k1x = u[i][j];
			double k1y = 0.25 * (v[i][j] + v[i - 1][j] + v[i][j + 1] + v[i - 1][j + 1]); // Average v at u's location

			// Stage 2: Velocity at first midpoint (t - dt/2)
			double x2 = x - 0.5 * k1x * dtbydx;
			double y2 = y - 0.5 * k1y * dtbydx;
			double k2x, k2y;
			samplevelocity(x2, y2, &k2x, &k2y);

			// Stage 3: Velocity at second midpoint (t - dt/2)
			double x3 = x - 0.5 * dtbydx * k2x;
			double y3 = y - 0.5 * dtbydx * k2y;
			double k3x, k3y;
			samplevelocity(x3, y3, &k3x, &k3y);

			// Stage 4: Velocity at endpoint (t - dt)
			double x4 = x - dtbydx * k3x;
			double y4 = y - dtbydx * k3y;
			double k4x, k4y;
			samplevelocity(x4, y4, &k4x, &k4y);

			// Combine velocities with RK4 weights
			double vx = (k1x + 2 * k2x + 2 * k3x + k4x) / 6.0;
			double vy = (k1y + 2 * k2y + 2 * k3y + k4y) / 6.0;

			// Trace back full step using combined velocity
			double prevX = x - dtbydx * vx;
			double prevY = y - dtbydx * vy;

			double NewVel, UselessVel;

			samplevelocity(prevX, prevY, &NewVel, &UselessVel);

			newU[i][j] = NewVel;
			//	}
		}

		// --- v‐advection
		if (j > 1 and j < cells - 1 and i > 0 and i < cells - 1) {
			//	if (fluid[i][j] != 0.f or fluid[i][j - 1] != 0.f) {
			double x = i + 0.5;
			double y = j;

			// Stage 1
			double k1x = 0.25 * (u[i][j] + u[i + 1][j] + u[i][j - 1] + u[i + 1][j - 1]); // Average u at v's location
			double k1y = v[i][j];

			// Stage 2: Velocity at first midpoint (t - dt/2)
			double x2 = x - 0.5 * k1x * dtbydx;
			double y2 = y - 0.5 * k1y * dtbydx;
			double k2x, k2y;
			samplevelocity(x2, y2, &k2x, &k2y);

			// Stage 3: Velocity at second midpoint (t - dt/2)
			double x3 = x - 0.5 * dtbydx * k2x;
			double y3 = y - 0.5 * dtbydx * k2y;
			double k3x, k3y;
			samplevelocity(x3, y3, &k3x, &k3y);

			// Stage 4: Velocity at endpoint (t - dt)
			double x4 = x - dtbydx * k3x;
			double y4 = y - dtbydx * k3y;
			double k4x, k4y;
			samplevelocity(x4, y4, &k4x, &k4y);

			// Combine velocities with RK4 weights
			double vx = (k1x + 2 * k2x + 2 * k3x + k4x) / 6.0;
			double vy = (k1y + 2 * k2y + 2 * k3y + k4y) / 6.0;

			// Trace back full step using combined velocity
			double prevX = x - dtbydx * vx;
			double prevY = y - dtbydx * vy;

			double NewVel, UselessVel;

			samplevelocity(prevX, prevY, &UselessVel, &NewVel);

			newV[i][j] = NewVel;
			//	}
		}
	}

	//	printf("Advected u and v!\n");

	auto tmp = u;
	u = newU;
	newU = tmp;

	auto tmp1 = v;
	v = newV;
	newV = tmp1;

	//	delete[] tmp;
	//	delete[] tmp1;




	// apply gravity
//#pragma omp parallel for 
	for (int i = 1; i <= cells - 2; ++i)
		for (int j = 2; j <= cells - 2; ++j) {
			/// Should I check whether the fluid is there ?
		 //	if (fluid[i][j] == 1.f and fluid[i][j - 1] == 1.f)
			v[i][j] += dvG;


		}

//	printf("Gravity applied!\n");

#if 0
	double(*fx)[cells] = new double[cells][cells];
	double(*fy)[cells] = new double[cells][cells];

	// VORTICITY CONFINEMENT
	for (int i = 0; i < cells; i++)
		for (int j = 0; j < cells; j++) {
			if (i != 0 and j != 0 and i != cells - 1 and j != cells - 1)
				vorticity[i][j] = fabs(
					((v[i + 1][j] + v[i + 1][j + 1]) / 2.f - (v[i - 1][j] + v[i - 1][j + 1]) / 2.f)
					- ((u[i][j + 1] + u[i + 1][j + 1]) / 2.f - (u[i][j - 1] + u[i + 1][j + 1]) / 2.f) / (2 * dx));

			fx[i][j] = 0.0;
			fy[i][j] = 0.0;
		}



	for (int i = 1; i < cells - 1; i++)
		for (int j = 1; j < cells - 1; j++) {
			float grad_x = (vorticity[i + 1][j] - vorticity[i - 1][j]) / (2 * dx);
			float grad_y = (vorticity[i][j + 1] - vorticity[i][j - 1]) / (2 * dx);

			float vecLen = sqrt(grad_x * grad_x + grad_y * grad_y);
			if (vecLen != 0.0) {
				grad_x /= vecLen;
				grad_y /= vecLen;

				float epsilon = 10.0;
				fx[i][j] = epsilon * dx * (grad_y * vorticity[i][j]);
				fy[i][j] = epsilon * dx * (-grad_x * vorticity[i][j]);

				if (i != 1) {
					u[i][j] += (fx[i][j] + fx[i - 1][j]) / 2 * dt;
				}
				if (j != 1) {
					v[i][j] += (fy[i][j] + fy[i][j - 1]) / 2 * dt;
				}
			}
		}

	delete[] fx;
	delete[] fy;
	
#endif 

	/*for (int a = 0; a < cells; a++) {
		// left
		u[1][a] = 0.0;
		v[0][a] = v[1][a]; // free slip wall
		u[cells - 1][a] = 0.0;	// right
		v[cells - 1][a] = v[cells - 2][a];
		v[a][1] = 0.0;	// down
		u[a][0] = u[a][1];
		v[a][cells - 1] = 0.0;	// up
		u[a][cells - 1] = u[a][cells - 2];
	}*/

	// Calculate divergence
	for (int idx = 0; idx < cells * cells; ++idx) {
		int i = idx / cells;
		int j = idx % cells;
		CalculateDivergence(i, j);
	}

/*	for (int i = 1; i < cells - 1; i++)
		for (int j = 1; j < cells - 1; j++) {
			if (i > 1 and i < cells - 2 and j > 1 and j < cells - 2) {
				u[i][j] += divergence[i][j] / 4.0;
				u[i + 1][j] -= divergence[i][j] / 4.0;
				v[i][j] += divergence[i][j] / 4.0;
				v[i][j + 1] -= divergence[i][j] / 4.0;
			}
			else {
				if ((i != j) or (i + j != cells - 1)) {
					u[i][j] += (i == 1 ? 0.0 : 1.0) * divergence[i][j] / 3.0;
					u[i + 1][j] -= (i == cells - 1 ? 0.0 : 1.0) * divergence[i][j] / 3.0;
					v[i][j] += (j == 1 ? 0.0 : 1.0) * divergence[i][j] / 3.0;
					v[i][j + 1] -= (j == cells - 1 ? 0.0 : 1.0) * divergence[i][j] / 3.0;
				}
			}
		}*/

	

//	printf("Divergence calculated!\n");

	// MIC-PCG 
	const double max_error = 1e-6;
		const int max_iter = 100;
	double max_delta = 0.f;
	float old_p = 0.f;
	int iter;
	for (iter = 0; iter < max_iter; ++iter) {

		// Gauss Seidel Method
		for (int i = 0; i < cells; ++i)
			for (int j = 0; j < cells; ++j) {

				// Non Wall cells
				if (i > 0 and j > 0 and i < cells - 1 and j < cells - 1) {
					if (fluid[i][j] == 1.f) {
						old_p = p[i][j];

						float new_p = CalculateNewPressure(i, j);

						p[i][j] = old_p + (1.8) * (new_p - old_p);
						max_delta = fmax(fabs(p[i][j] - old_p), max_delta);
					}
					else if (fluid[i][j] == 0.f) {
						p[i][j] = 0.0; 
					}
				}

				// "Wall Pressures"
			/*	else
				{
					if (i != j and i + j != cells - 1) {
						if (i == 0) {
							p[i][j] = p[i + 1][j] - divergence[i][j] / dtbydx;
						}
						else if (i == cells - 1) {
							p[i][j] = p[i - 1][j] - divergence[i][j] / dtbydx;
						}
						else if (j == 0) {
							p[i][j] = p[i][j + 1] - divergence[i][j] / dtbydx;
						}
						else if (j == cells - 1) {
							p[i][j] = p[i][j - 1] - divergence[i][j] / dtbydx;
						}
					}
				}*/
			}

		if (max_delta < max_error ) break;
	}

	for (int z = 1; z < cells - 1; z++) {
		// left bc 
		p[0][z] = p[1][z] - divergence[0][z];
		// right bc
		p[cells - 1][z] = p[cells - 2][z] - divergence[cells - 1][z];
		// up bc
		p[z][cells - 1] = p[z][cells - 2] - divergence[z][cells - 1];
		// down bc
		p[z][0] = p[z][1] - divergence[z][0];
	}

//		printf("error = %0.30f\n", max_delta);

//	if(iter == max_iter)
	//	printf("report iteration limit exceeded\n");

//	printf("Solved the pressure with max_delta = %lf\n", max_delta);




		// Applying the pressure gradient

		// Horizontal velocity (u)
	//#pragma omp parallel for
	for (int i = 1; i <= cells - 1; ++i) {
		for (int j = 1; j < cells - 1; ++j) {
		//	if (fluid[i][j] or fluid[i - 1][j])
				u[i][j] -= (p[i][j] - p[i - 1][j]);
		}
	}

	// Vertical velocity (v)
//#pragma omp parallel for
	for (int i = 1; i < cells - 1; ++i) {
		for (int j = 1; j <= cells - 1; ++j) {
		//		if (fluid[i][j] or  fluid[i][j-1])
			v[i][j] -= (p[i][j] - p[i][j - 1]);
		}
	}

//	printf("Applied the gradient!\n");


	// Stability Condition work
/*	float Umax = 0;
	for (int a = 0; a < cells; a++)
		for (int b = 0; b < cells; b++) {
			if (fabs(u[a][b]) > Umax)
				Umax = fabs(u[a][b]);
			if (fabs(v[a][b]) > Umax)
				Umax = fabs(v[a][b]);
		}
	Umax += sqrt(5 * dx * -gravity);


	printf("const = %f\n", Umax * dt / dx);*/

//	printf("\n");

}


void processInput(GLFWwindow* window)
{
	if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
	{
		glfwSetWindowShouldClose(window, true);
	}
}

void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
	glViewport(0, 0, width, height);
}


int main()
{
	if (!glfwInit())
		return -1;

	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

	GLFWwindow* window = glfwCreateWindow(700, 700, "Buoyancy Simulation", NULL, NULL);
	if (!window)
	{
		std::cout << "Failed to create the window!" << std::endl;
		glfwTerminate();
		return -1;
	}

	glfwMakeContextCurrent(window);
	glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);

	// Disable VSync (immediately after window creation)
	glfwSwapInterval(0);


	if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
	{
		std::cout << "Failed to load opengl funtion pointers!" << std::endl;
		glfwTerminate();
		return -1;
	}

	// Enable basic OpenGL optimizations
	glEnable(GL_CULL_FACE);
	glCullFace(GL_BACK);

	const char* vertexShaderSrc =
		"#version 330 core\n"
		"layout (location = 0) in vec2 aPos;\n"
		"layout (location = 1) in vec2 aTexCoord;\n"
		"out vec2 texCoord;\n"
		"void main() {\n"
		"    // Rotate 90 degrees clockwise\n"
		"    gl_Position = vec4(-aPos.y, aPos.x, 0.0, 1.0);\n"
		"    texCoord = aTexCoord;\n"
		"}\0";

	const char* fragmentShaderSrc = R"(
#version 330 core
in vec2 texCoord;
out vec4 fragColor;
uniform sampler2D tex;

void main() {
    // Grid parameters
    const float gridSize = 256;
    const float cellSize = 1.0 / gridSize;
    
    // Calculate grid position
    vec2 gridPos = texCoord * gridSize;
    
    // Boundary check (first and last columns/rows)
    bool isBoundary = (gridPos.x < 1.0) || (gridPos.x >= gridSize - 1.0) ||
                     (gridPos.y < 1.0) || (gridPos.y >= gridSize - 1.0);
    
    // DEBUG: Show entire boundary area in red
    if (isBoundary) {
        fragColor = vec4(1.0, 0.0, 0.0, 1.0);
    } else {
        float density = texture(tex, texCoord).r;
        vec3 bgColor = vec3(0.0, 0.0, 0.0);
        vec3 fluidColor = vec3(0.0, 0.0, 1.0);
        fragColor = vec4(mix(bgColor, fluidColor, density), 1.0);
    }
}
)";

	glPixelStorei(GL_UNPACK_ROW_LENGTH, cells);
	glPixelStorei(GL_UNPACK_SWAP_BYTES, GL_FALSE);
	glPixelStorei(GL_UNPACK_SKIP_ROWS, 0);
	glPixelStorei(GL_UNPACK_SKIP_PIXELS, 0);

	unsigned int vertexShader = glCreateShader(GL_VERTEX_SHADER);
	glShaderSource(vertexShader, 1, &vertexShaderSrc, 0);
	glCompileShader(vertexShader);
	int success;
	char infoLog[512];
	glGetShaderiv(vertexShader, GL_COMPILE_STATUS, &success);
	if (!success)
	{
		glGetShaderInfoLog(vertexShader, 512, 0, infoLog);
		std::cout << "Failed to compile the vertex shader! ERR: " << infoLog << std::endl;
	}

	unsigned int fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
	glShaderSource(fragmentShader, 1, &fragmentShaderSrc, 0);
	glCompileShader(fragmentShader);
	glGetShaderiv(fragmentShader, GL_COMPILE_STATUS, &success);
	if (!success)
	{
		glGetShaderInfoLog(fragmentShader, 512, 0, infoLog);
		std::cout << "Failed to compile the fragment shader! ERR: " << infoLog << std::endl;
	}

	unsigned int shaderProgram = glCreateProgram();
	glAttachShader(shaderProgram, vertexShader);
	glAttachShader(shaderProgram, fragmentShader);
	glLinkProgram(shaderProgram);
	glGetProgramiv(shaderProgram, GL_LINK_STATUS, &success);
	if (!success)
	{
		glGetProgramInfoLog(shaderProgram, 512, 0, infoLog);
		std::cout << "Failed to link the shader program! ERR: " << infoLog << std::endl;
	}



	glDeleteShader(vertexShader);
	glDeleteShader(fragmentShader);

	// ===== ADD TEXTURE SETUP CODE HERE =====
	// Create texture
	unsigned int texture;
	glGenTextures(1, &texture);
	// Divergence Debugging
	glGenTextures(1, &divergenceTexture);

	// Initialize both textures
	for (GLuint tex : {fluidTexture, divergenceTexture}) {
		glBindTexture(GL_TEXTURE_2D, tex);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, cells, cells, 0, GL_RED, GL_FLOAT, nullptr);
	}

	// Create fluid shader program (your original) CHECK WHETHER fragmentShaderSrc or what
	fluidProgram = createShaderProgram(vertexShaderSrc, fragmentShaderSrc);

	// Create divergence shader program
	const char* divergenceFragShaderSrc = R"(
	#version 330 core
in vec2 texCoord;
out vec4 fragColor;
uniform sampler2D divergenceTexture;

void main() {
    // Sample raw divergence value
    float div = texture(divergenceTexture, texCoord).r;
    
    // Clamp and visualize
    if (div > 0.0) {
        // Positive divergence (red)
        fragColor = vec4(clamp(div, 0.0, 1.0), 0.0, 0.0, 1.0);
    } else {
        // Negative divergence (blue)
        fragColor = vec4(0.0, 0.0, clamp(-div, 0.0, 1.0), 1.0);
    }
}
	)";

	divergenceProgram = createShaderProgram(vertexShaderSrc, divergenceFragShaderSrc);


	glBindTexture(GL_TEXTURE_2D, texture);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RED, cells, cells, 0, GL_RED, GL_FLOAT, nullptr);

	// Set texture uniform
	glUseProgram(shaderProgram);
	glUniform1i(glGetUniformLocation(shaderProgram, "screenTexture"), 0);

	// ... continue with VAO/VBO setup ...
	// Update vertex data for full-screen quad
	float vertices[] = {
		// Positions   // Texture Coords
		-1.0f, -1.0f,  0.0f, 0.0f,  // Bottom-left
		 1.0f, -1.0f,  1.0f, 0.0f,  // Bottom-right
		-1.0f,  1.0f,  0.0f, 1.0f,  // Top-left
		 1.0f,  1.0f,  1.0f, 1.0f   // Top-right
	};
	unsigned int indices[] = {
		0, 1, 2,  // Triangle 1
		1, 3, 2   // Triangle 2
	};

	unsigned int VAO, VBO, EBO;
	glGenVertexArrays(1, &VAO);
	glGenBuffers(1, &VBO);
	glGenBuffers(1, &EBO);

	glBindVertexArray(VAO);

	// 1. ADD MISSING BUFFER DATA UPLOADS
	glBindBuffer(GL_ARRAY_BUFFER, VBO);
	glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);

	// 2. CORRECTED VERTEX ATTRIBUTE POINTERS
	glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)0);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)(2 * sizeof(float)));
	glEnableVertexAttribArray(1);

	// 3. ADD TEXTURE SWIZZLE TO SHOW GRAYSCALE
	glBindTexture(GL_TEXTURE_2D, texture);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_SWIZZLE_G, GL_RED);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_SWIZZLE_B, GL_RED);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_SWIZZLE_A, GL_ONE);

	glBindVertexArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	initialize();


	while (!glfwWindowShouldClose(window)) {
		processInput(window);

		// Toggle visualization when D is pressed
		static bool dKeyPressed = false;
		if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS) {
			if (!dKeyPressed) {
				showDivergence = !showDivergence;
				dKeyPressed = true;
			}
		}
		else {
			dKeyPressed = false;
		}

		glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT);


		double CurrTime = glfwGetTime();
		double dt = CurrTime - PrevTime;
		Tick(dt);
		PrevTime = CurrTime;

		//	printf("FPS = %f\n", 1/dt);

			// 4. BIND TEXTURE BEFORE UPDATING
			// Update textures
		glBindTexture(GL_TEXTURE_2D, fluidTexture);
		glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, cells, cells, GL_RED, GL_FLOAT, fluid);

		glBindTexture(GL_TEXTURE_2D, divergenceTexture);
		glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, cells, cells, GL_RED, GL_FLOAT, divergence);

		// Choose program and texture based on mode
		GLuint currentProgram;
		GLuint currentTexture;
		float maxValue = 1.0f;  // Default for fluid

		if (showDivergence) {
			currentProgram = divergenceProgram;
			currentTexture = divergenceTexture;

			// Calculate max divergence for normalization
			maxValue = 0.0f;
			for (int i = 0; i < cells; i++)
				for (int j = 0; j < cells; j++) {
					maxValue = fmax(maxValue, fabs(divergence[i][j]));
				}
			if (maxValue < 0.001f) maxValue = 1.0f;  // Avoid division by zero
		}
		else {
			currentProgram = fluidProgram;
			currentTexture = fluidTexture;
		}


		// 5. ACTIVATE TEXTURE UNIT BEFORE RENDERING
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, texture);

		// Render
		glUseProgram(currentProgram);
		glUniform1f(glGetUniformLocation(currentProgram, "maxValue"), maxValue);

		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, currentTexture);
		glBindVertexArray(VAO);
		glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);


		glfwSwapBuffers(window);
		glfwPollEvents();

		// Frame rate limiting to 30 FPS

		double frameEnd = glfwGetTime();
		double elapsedTime = frameEnd - CurrTime;  // Time taken for this frame

		if (elapsedTime < targetFrameTime) {
			double sleepTime = targetFrameTime - elapsedTime;
			// Busy-wait for precision (avoids including extra headers)
			double waitEnd = frameEnd + sleepTime;
			while (glfwGetTime() < waitEnd) {
				// Empty loop - waits until target time is reached
			}
		}
	}

	glDeleteTextures(1, &texture);
	glDeleteProgram(shaderProgram);
	glDeleteBuffers(1, &VBO);
	glDeleteVertexArrays(1, &VAO);

	glfwTerminate();
	return 0;
}