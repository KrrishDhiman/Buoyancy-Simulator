#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <math.h>
#include <iostream>



float clamp(float value, float minVal, float maxVal) {
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



double gravity = -0.8f;

double prevTime = 0.0f;

const float omega = 1.9f;

// Just to center the fluid
float screensize = 1.5f;

const int cells = 256;
const float stepsize = 1.5f / cells;
//float center = cells * stepsize / 2 - 0.75f;
float p[cells][cells] = { 0.f };
float(*u)[cells] = new float[cells + 1][cells];

float(*v)[cells + 1] = new float[cells][cells + 1];

float(*fluid)[cells] = new float[cells][cells];
float divergence[cells][cells] = { 0.f };
float vorticity[cells][cells] = { 0.f };
float cfl = 0.5f;
float(*newV)[cells + 1] = new float[cells][cells + 1];
float(*newU)[cells] = new float[cells + 1][cells];


// grid spacing (same as your stepsize)
const float dx = stepsize;
// time step
//float dt;	//cfl * dx / UMax;
//float dt = 0.001;



float PrevTime = 0.f;

// v^2 = 0^2 + 2(-g)(-S) where S = - 2cells*dx/3
float UMax = pow(2 * fabs(gravity) * 2.f * (screensize) / 3.f, 0.5f);

/// SOME CALCULATIONS: UMax approx 54.7, so calculating dx at dt = 0.001, apparently dx approx 1!
/// tried g = 1000, UMax 44.7, so with dt = 0.001, dx ~ 0.8
/// SO dt =  is required... 

// For faster times, precomputed data
const float inv_dx = 1.0f / dx;
//const float dt_dx = dt * inv_dx;
//const float G_incrementalV = gravity * dt;

bool is_fluid(int i, int j) {
	return fluid[i][j] == 1.f; // Assuming 1=fluid, 0=air
}


// Kind of beginplay
void initialize() {

	for (int a = 0; a < cells; a++)
		for (int b = 0; b < cells; b++) {
			if (a > cells / 6 and a < 5 * cells / 6 and b > cells / 6  and b < 5* cells /6) {
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
				v[a][b] = 0.f;
				newV[a][b] = 0.f;
			}
			if (b != cells) {
				u[a][b] = 0.f;
				newU[a][b] = 0.f;
			}
		}
}

void CalculateDivergence(int i, int j) {
	if (fluid[i][j] != 0.f) {
	if (i == 0 or i == cells - 1 or j == 0 or j == cells - 1) {
	
	}
		
	else if (i == 1 or i == cells - 2 or j == 1 or j == cells - 2) {
		float du_dx = (u[i + 1][j] * (i + 1 == cells - 1 ? 0.f : 1.f) - u[i][j] * (i == 1 ? 0.f : 1.f));

		float dv_dy = (v[i][j + 1] * (j + 1 == cells - 1 ? 0.f : 1.f) - v[i][j] * (j == 1 ? 0.f : 1.f));

		divergence[i][j] = du_dx + dv_dy;
	}

	else {
		float du_dx = (u[i + 1][j]- u[i][j]);
		float dv_dy = (v[i][j + 1]- v[i][j]);
		divergence[i][j] = du_dx + dv_dy;
	}


	}
else
divergence[i][j] = 0.f;
}

void CalculatePressure(int i, int j) {
	// Using Poisson Equation, del^2(P) = del(u*)/dt
	float max_error = 0.001f;
	float LastPressure;
	if (fabs(i - j) + fabs(i + j - cells + 1) < cells - 1) {
		for (int iter = 0; iter < 10; iter++) {
			p[i][j] = (p[i + 1][j]
				+ p[i - 1][j]
				+ p[i][j + 1]
				+ p[i][j - 1]
				- divergence[i][j] * dx * dx /* dt*/) / 4.0;
			LastPressure = p[i][j];
			if (fabs(LastPressure - p[i][j]) <= max_error) break;
		}
	}

	// Apply Neumann boundary conditions (e.g., dp/dn = 0)

}

void ForceIncompressibility(int i, int j) {

	if (fabs(i - j) + fabs(i + j - cells + 1) < cells - 1) // checking for non boundary cells
	{
		float d = divergence[i][j] / 4.f;
		u[i][j] += d;
		u[i + 1][j] -= d;
		v[i][j] += d;
		v[i][j + 1] -= d;

	}

	else {

		if (i == j or i + j == cells - 1)
		{
			if (i == 0 and j == 0) {
				float d = divergence[i][j] / 2;
				//	u[i][j] += d;
				u[i + 1][j] -= d;
				//	v[i][j] += d;
				v[i][j + 1] -= d;
			}
			else if (i == cells - 1 and j == cells - 1) {
				float d = divergence[i][j] / 2;
				u[i][j] += d;
				//	u[i + 1][j] -= d;
				v[i][j] += d;
				//	v[i][j + 1] -= d;
			}
			else if (i == cells - 1 and j == 0) {
				float d = divergence[i][j] / 2;
				u[i][j] += d;
				//	u[i + 1][j] -= d;
				//	v[i][j] += d;
				v[i][j + 1] -= d;
			}
			else if (i == 0 and j == cells - 1) {
				float d = divergence[i][j] / 2;
				//	u[i][j] += d;
				u[i + 1][j] -= d;
				v[i][j] += d;
				//	v[i][j + 1] -= d;
			}
		}

		else {
			if (i == 0)
			{
				float d = divergence[i][j] / 3;
				//	u[i][j] += d;
				u[i + 1][j] -= d;
				v[i][j] += d;
				v[i][j + 1] -= d;
			}
			else if (j == 0)
			{
				float d = divergence[i][j] / 3;
				u[i][j] += d;
				u[i + 1][j] -= d;
				//	v[i][j] += d;
				v[i][j + 1] -= d;
			}
			else if (i == cells - 1)
			{
				float d = divergence[i][j] / 3;
				u[i][j] += d;
				//	u[i + 1][j] -= d;
				v[i][j] += d;
				v[i][j + 1] -= d;
			}
			else if (j == cells - 1)
			{
				float d = divergence[i][j] / 3;
				u[i][j] += d;
				u[i + 1][j] -= d;
				v[i][j] += d;
				//	v[i][j + 1] -= d;
			}
		}

	}
}

void BoundaryConditions(int i) {

	// Apply boundary conditions for all walls
		// Left Wall (i=0)
	u[0][i] = 0.0f;       // Normal velocity: no-penetration
	v[0][i] = v[1][i];   // No-slip for tangential velocity (v)

	// Right Wall (i=cells)
	u[cells - 1][i] = 0.0f;           // Normal velocity
	v[cells - 1][i] = v[cells - 2][i]; // No-slip for v

	// Bottom Wall (j=0)
	v[i][0] = 0.0f;       // Normal velocity
	u[i][0] = -u[i][1];   // No-slip for u

	// Top Wall (j=cells)
	v[i][cells - 1] = 0.0f;           // Normal velocity
	u[i][cells - 1] = -u[i][cells - 2]; // No-slip for u
}

void AdvectUV(int i, int j, double dt) {
	
}

void AdvectFluid(int i, int j) {

}

float CalculateNewPressure(int i, int j, double dxbydt) {

	int NonSolidCells = 0;
	float new_p = 0.f;

	if (i + 1 != cells - 1) {
	//	p[i + 1][j] *= fluid[i + 1][j];
		NonSolidCells++;
		new_p += p[i + 1][j];
	}
	if (j + 1 != cells - 1) {
	//	p[i][j + 1] *= fluid[i][j + 1];
		NonSolidCells++;
		new_p += p[i][j + 1];
	}
	if (i - 1 != 0) {
	//	p[i - 1][j] *= fluid[i - 1][j];
		NonSolidCells++;
		new_p += p[i - 1][j];
	}
	if (j - 1 != 0) {
	//	p[i][j - 1] *= fluid[i][j - 1];
		NonSolidCells++;
		new_p += p[i][j - 1];
	}

	new_p -= divergence[i][j] * dxbydt;
	
	return new_p / NonSolidCells;
}



void Tick(double dt) {

	double dvG = gravity * dt;
	double dtbydx = dt / dx;
//	double dxbydt = dx / dt;

	// velocity at solid boundaries = 0
//#pragma omp parallel for
/*	for (int i = 0; i < cells; i++) {
		v[i][0] = 0.0f;  // bottiom wall
		v[i][cells] = 0.0f; // Top wall
		u[0][i] = 0.0f; // Left wall
		u[cells][i] = 0.0f; // Right wall
	}*/


	// Advect u and v
//#pragma omp parallel for 
	for (int idx = 0; idx <= cells * cells; ++idx) {
		int i = idx / cells;
		int j = idx % cells;

		// --- u‐advection
		if (i > 0 && i < cells and j > 0 and j < cells - 1) {

			///  why fluid[i-1][j] != 0 ?
			if (fluid[i][j] != 0.f and fluid[i - 1][j] != 0.f) {
				float x = i;
				float y = j + 0.5f;

				// Velocity at u's location
				float vx = u[i][j];
				float vy = 0.25f * (v[i][j] + v[i - 1][j] + v[i][j + 1] + v[i - 1][j + 1]); // Average v at u's location

				// Trace back
				float px = x - vx * dtbydx;
				float py = y - vy * dtbydx;

				// Clamp and interpolate
				px = clamp(px, 1.f, cells - 1.f); // Avoid sampling boundaries
				py = clamp(py, 1.f, cells - 1.f);

				// Bi-linear Interpolation Code:
				float a = px - (int)px;
				float b = py - (int)py;
				int x0 = (int)px;
				int y0 = (int)py;

				float val = (1 - a) * (1 - b) * u[x0][y0] +
					(1 - a) * b * u[x0][y0 + 1] +
					a * (1 - b) * u[x0 + 1][y0] +
					a * b * u[x0 + 1][y0 + 1];

				newU[i][j] = val;
			}
		}

		// --- v‐advection
		if (j > 0 and j < cells and i > 0 and i < cells - 1) {
			if (fluid[i][j] != 0.f and fluid[i][j - 1] != 0.f) {
				float x = i + 0.5f;
				float y = j;


				float vx = 0.25f * (u[i][j] + u[i + 1][j] + u[i][j - 1] + u[i + 1][j - 1]); // Average u at v's location
				float vy = v[i][j];

				// Trace back
				float px = x - vx * dtbydx;
				float py = y - vy * dtbydx;

				// Clamp and interpolate
				px = clamp(px, 1.f, cells - 1.f);
				py = clamp(py, 1.f, cells - 1.f); // Avoid sampling boundaries

				// Bi-linear Interpolation Code:
				float a = px - (int)px;
				float b = py - (int)py;
				int x0 = (int)px;
				int y0 = (int)py;


				newV[i][j] = (1 - a) * (1 - b) * v[x0][y0] +
					(1 - a) * b * v[x0][y0 + 1] +
					a * (1 - b) * v[x0 + 1][y0] +
					a * b * v[x0 + 1][y0 + 1];
			}
		}
	}

	auto tmp = u;
	u = newU;
	newU = tmp;

	auto tmp1 = v;
	v = newV;
	newV = tmp1;

	//	delete[] tmp;
	//	delete[] tmp1;





	float(*temp_fluid)[cells] = new float[cells][cells];
	for (int a = 0; a < cells; a++)
		for (int b = 0; b < cells; b++)
			temp_fluid[a][b] = 0.f;

	// Advecting the fluid itself
//#pragma omp parallel for 
	for (int idx = 0; idx < cells * cells; idx++)
	{
		int i = idx / cells;
		int j = idx % cells;



		if (i > 0 and j > 0 and i < cells - 1 and j < cells - 1) {

			float x = i + 0.5f;
			float y = j + 0.5f;

			// Forward Euler
		/*	float vx = (u[i][j] + u[i + 1][j]) * 0.5f;
			float vy = (v[i][j] + v[i][j + 1]) * 0.5f;

			float prevX = x - vx * dt_dx;
			float prevY = y - vy * dt_dx;*/


			// RK2
			// Current velocity at starting point
			float vx0 = (u[i][j] + u[i + 1][j]) * 0.5f;
			float vy0 = (v[i][j] + v[i][j + 1]) * 0.5f;

			// Half-step position (midpoint)
			float midX = x - 0.5f * dtbydx * vx0;
			float midY = y - 0.5f * dtbydx * vy0;

			// Clamp midpoint to grid boundaries
			midX = clamp(midX, 1.f, cells - 1);
			midY = clamp(midY, 1.f, cells - 1);;

			// Calculate velocity at midpoint
			int mi0 = (int)midX;
			int mj0 = (int)midY;
			int mi1 = mi0 + 1;
			int mj1 = mj0 + 1;

			// Safeguard against boundary overflows
			mi1 = (mi1 >= cells) ? cells - 1 : mi1;
			mj1 = (mj1 >= cells) ? cells - 1 : mj1;

			// Bilinear interpolation for velocity at midpoint
			float s = midX - mi0;
			float t = midY - mj0;

			float vxMid = (1 - s) * (1 - t) * u[mi0][mj0]
				+ s * (1 - t) * u[mi1][mj0]
				+ (1 - s) * t * u[mi0][mj1]
				+ s * t * u[mi1][mj1];

			float vyMid = (1 - s) * (1 - t) * v[mi0][mj0]
				+ s * (1 - t) * v[mi1][mj0]
				+ (1 - s) * t * v[mi0][mj1]
				+ s * t * v[mi1][mj1];

			// Full step using midpoint velocity
			float prevX = x - dtbydx * vxMid;
			float prevY = y - dtbydx * vyMid;

			prevX = clamp(prevX, 1.f, cells - 1); // Clamping the prevx between 0 and cells - 1
			prevY = clamp(prevY, 1.f, cells - 1);

			int pi = (int)(prevX);
			int pj = (int)(prevY);

			temp_fluid[i][j] = fluid[pi][pj];
		}
	}

	// Apply Neumann BCs after each full grid sweep
/*	for (int j = 1; j < cells - 1; ++j) {
			p[0][j] = p[1][j] + u[1][j] / dtbydx;
	//	p[0][j] = p[1][j] + 50000.f;

			p[cells - 1][j] = p[cells - 2][j] + u[cells - 1][j] / dtbydx;
	//	p[cells - 1][j] = p[cells - 2][j] + 50000.f;
	}
	for (int i = 1; i < cells - 1; ++i) {
			p[i][0] = p[i][1] + v[i][0] / dtbydx;
	//	p[i][0] = p[i][1] + 50000.f;

			p[i][cells - 1] = p[i][cells - 2] + v[i][cells - 1] / dtbydx;
	//	p[i][cells - 1] = p[i][cells - 2] + 50000.f;
	}*/

	auto f = temp_fluid;
	temp_fluid = fluid;
	fluid = f;
	delete[] temp_fluid;

	// apply gravity
//#pragma omp parallel for 
	for (int i = 0; i <= cells - 1; ++i)
		for (int j = 1; j <= cells - 1; ++j) {
			/// Should I check whether the fluid is there ?
		//	if (fluid[i][j] == 1.f and fluid[i][j - 1] == 1.f)
				v[i][j] += dvG;

			
		}


	float(*fx)[cells] = new float[cells][cells];
	float(*fy)[cells] = new float[cells][cells];

	// VORTICITY CONFINEMENT
	for(int i = 0; i < cells; i++)
		for (int j = 0; j < cells; j++) {
			if (i != 0 and j != 0 and i != cells - 1 and j != cells - 1)
			vorticity[i][j] = fabs(
					((v[i + 1][j] + v[i + 1][j + 1]) / 2.f - (v[i - 1][j] + v[i - 1][j + 1]) / 2.f)
				  - ((u[i][j + 1] + u[i + 1][j + 1]) / 2.f - (u[i][j - 1] + u[i + 1][j + 1]) / 2.f) / (2 * dx));

			fx[i][j] = 0.f;
			fy[i][j] = 0.f;
		}



	for (int i = 1; i < cells - 1; i++)
		for (int j = 1; j < cells - 1; j++) {
			float grad_x = (vorticity[i + 1][j] - vorticity[i - 1][j]) / (2 * dx);
			float grad_y = (vorticity[i][j + 1] - vorticity[i][j - 1]) / (2 * dx);
			
			float vecLen = sqrt(grad_x * grad_x + grad_y * grad_y);
			if (vecLen != 0.f) {
				grad_x /= vecLen;
				grad_y /= vecLen;

				float epsilon = 1.f;
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



	// Calculate divergence
//#pragma omp parallel for 
	for (int idx = 0; idx < cells * cells; ++idx) {
		int i = idx / cells;
		int j = idx % cells;
		CalculateDivergence(i, j);

	//	if (fluid[i][j] != 1 and fluid[i][j] != 0)
		//	printf("fluid = %f\n", fluid[i][j]);
		//	printf("pressure = %f\n", p[i][j]);
	}

	const double max_error = 1e-8;
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

						float new_p = CalculateNewPressure(i, j, 1/dtbydx);

						p[i][j] = old_p + (1.9) * (new_p - old_p);
						max_delta = fmax(fabs(p[i][j] - old_p), max_delta);
					}
					else if (fluid[i][j] == 0.f) {
						p[i][j] = 0.f; continue;
					}
				}

				// "Wall Pressures"
				else
				{
					if (i != j and i + j != cells - 1) {
						if (i == 0) {
							p[i][j] = p[i+1][j] - u[i+1][j] / dtbydx;
						}
						else if (i == cells - 1) {
							p[i][j] = p[i-1][j] + u[i][j] / dtbydx;
						}
						else if (j == 0) {
							p[i][j] = p[i][j + 1] - v[i][1] / dtbydx;
						}
						else if (j == cells - 1) {
							p[i][j] = p[i][j - 1] + v[i][j] / dtbydx;
						}
					}
				}
			}

	
				
			


		if (max_delta < max_error and iter > 8) break;
	}

//	printf("error = %0.15f\n", max_delta);
//	printf("iter = %d\n", iter);




	// Applying the pressure gradient

	// Horizontal velocity (u)
//#pragma omp parallel for
	for (int i = 1; i < cells; ++i) {
		for (int j = 1; j < cells - 1; ++j) {
			if (fluid[i][j] or fluid[i-1][j])
				u[i][j] -= (p[i][j] - p[i - 1][j]) * dtbydx;
		}
	}

	// Vertical velocity (v)
//#pragma omp parallel for
	for (int i = 0; i < cells; ++i) {
		for (int j = 1; j < cells; ++j) {
		//	if (fluid[i][j] or  fluid[i][j-1])
				v[i][j] -= (p[i][j] - p[i][j - 1]) * dtbydx;
		}
	}



	// Stability Condition work
	float Umax = 0;
	for(int a = 0; a < cells; a++)
		for (int b = 0; b < cells; b++) {
			if (fabs(u[a][b]) > Umax)
				Umax = fabs(u[a][b]);
			if (fabs(v[a][b]) > Umax)
				Umax = fabs(v[a][b]);
		}
	Umax += sqrt(5 * dx * -gravity);


	printf("const = %f\n", Umax * dt / dx);

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
    float density = texture(tex, texCoord).r;
    
    // Background color (dark blue/black)
    vec3 bgColor = vec3(0.0, 0.0, 0.0);
    
    // Fluid color (vibrant blue)
    vec3 fluidColor = vec3(0.0, 0.0, 1.0);
    
    // Add glow effect at higher densities
    float glow = smoothstep(0.3, 1.0, density);
    
    // Blend between background and fluid with glow
    vec3 color = mix(bgColor, fluidColor, density + glow*0.5);
    
    fragColor = vec4(color, 1.0);
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




	const double targetFrameTime = 1.0 / 60.f;  // 30 FPS = 0.033 seconds per frame

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
				for(int j = 0; j < cells; j++){
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