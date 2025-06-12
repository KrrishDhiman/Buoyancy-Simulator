#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <math.h>
#include <iostream>



float clamp(float value, float minVal, float maxVal) {
	if (value < minVal) return minVal;
	if (value > maxVal) return maxVal;
	return value;
}



double gravity = -1.f;

double prevTime = 0.0f;

const float omega = 1.8f;

// Just to center the fluid
float screensize = 1.5f;

const int cells = 256;
const float stepsize = 1.5f / cells;
float center = cells * stepsize / 2 - 0.75f;
float p[cells][cells] = { 0.f };
float(*u)[cells] = new float[cells + 1][cells];

float(*v)[cells + 1] = new float[cells][cells + 1];

float(*fluid)[cells] = new float[cells][cells];
float divergence[cells][cells] = { 0.f };
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
	return fluid[i][j] > 0.5f; // Assuming 1=fluid, 0=air
}


// Kind of beginplay
void initializeColours() {

	for (int a = 0; a < cells; a++)
		for (int b = 0; b < cells; b++) {
			float x = -0.75f + a * stepsize;
			float y = -0.75f + b * stepsize;
			if (a > cells / 3 and a < 2 * cells / 3 and b > cells / 3 and b < 2 * cells / 3) {
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
		
		
			float du_dx = (u[i + 1][j] - u[i][j]) / dx;
		
			float dv_dy = (v[i][j + 1] - v[i][j]) / dx;

		divergence[i][j] = du_dx + dv_dy;
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



void Tick(double dt) {

	double dvG = gravity * dt;
	double dtbydx = dt / dx;
	double dx2bydt = dx * dx / dt;

	// apply gravity
//#pragma omp parallel for 
	for (int i = 0; i <= cells - 1; ++i)
		for (int j = 1; j <= cells - 1; ++j) {
			/// Should I check whether the fluid is there ?
			if (fluid[i][j] != 0.f  and fluid[i][j - 1] != 0.f)
				v[i][j] += dvG;
			//	printf("fluid = %f\n", fluid[i][j]);
		}

	// Calculate divergence
//#pragma omp parallel for 
	for (int idx = 0; idx < cells * cells; ++idx) {
		int i = idx / cells;
		int j = idx % cells;
		CalculateDivergence(i, j);
		//	printf("pressure = %f\n", p[i][j]);
	}

	const double max_error = (float)1e-3;
	const int max_iter = 50;
	double max_delta = 0.f;
	float old_p = 0.f;
	int iter;
	for (iter = 0; iter < max_iter; ++iter) {


		// Update interior cells (Gauss-Seidel)
		for (int i = 1; i < cells - 1; ++i) 
			for (int j = 1; j < cells - 1; ++j) {

				if (is_fluid(i, j) == 1.f) {
					old_p = p[i][j];
					float new_p = (p[i + 1][j] + p[i - 1][j] + p[i][j + 1] + p[i][j - 1]
						- divergence[i][j] * dx2bydt) / 4.f;
					p[i][j] = old_p + (1.8) * (new_p - old_p);
					//	p[i][j] *= (float)is_fluid(i, j);
					max_delta = fmax(max_delta, fabs(p[i][j] - old_p));
				}
				else {
					p[i][j] = 0.f; continue;
				}
				
			}


		if (max_delta < max_error) break;
	}

//	printf("error = %f\n", max_delta);

	// Apply Neumann BCs after each full grid sweep
	for (int j = 0; j < cells; ++j) {
		p[0][j] = p[1][j];           // Left
		p[cells - 1][j] = p[cells - 2][j]; // Right
	}
	for (int i = 0; i < cells; ++i) {
		p[i][0] = p[i][1];           // Bottom
		p[i][cells - 1] = p[i][cells - 2]; // Top
	}


	// Applying the pressure gradient

	// Horizontal velocity (u)
//#pragma omp parallel for
	for (int i = 1; i < cells; ++i) {
		for (int j = 0; j < cells; ++j) {
			if (is_fluid(i, j))
				u[i][j] -= (p[i][j] - p[i - 1][j]) * dtbydx;
		}
	}

	// Vertical velocity (v)
//#pragma omp parallel for
	for (int i = 0; i < cells; ++i) {
		for (int j = 0; j <= cells; ++j) {
			if (is_fluid(i, j))
				v[i][j] -= (p[i][j] - p[i][j - 1]) * dtbydx;
		}
	}

	// velocity at solid boundaries = 0
//#pragma omp parallel for
	for (int i = 0; i < cells; i++) {
		v[i][0] = 0.0f;  // bottiom wall
		v[i][cells] = 0.0f; // Top wall
		u[0][i] = 0.0f; // Left wall
		u[cells][i] = 0.0f; // Right wall
	}

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

	auto f = temp_fluid;
	temp_fluid = fluid;
	fluid = f;
	delete[] temp_fluid;

	// Stability Condition work
	float Umax = 0;
	for(int a = 0; a < cells; a++)
		for (int b = 0; b < cells; b++) {
			if (fabs(u[a][b]) > Umax)
				Umax = u[a][b];
			if (fabs(v[a][b]) > Umax)
				Umax = v[a][b];
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

	

	if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
	{
		std::cout << "Failed to load opengl funtion pointers!" << std::endl;
		glfwTerminate();
		return -1;
	}

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

	const char* fragmentShaderSrc =
		"#version 330 core\n"
		"in vec2 texCoord;\n"
		"out vec4 fragColor;\n"
		"uniform sampler2D screenTexture;\n"
		"void main() {\n"
		"    fragColor = texture(screenTexture, texCoord);\n"
		"}\0";

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

	initializeColours();

	// Disable VSync (immediately after window creation)
	glfwSwapInterval(0);

	// Enable basic optimizations
	glEnable(GL_CULL_FACE);
	glCullFace(GL_BACK);

	while (!glfwWindowShouldClose(window)) {
		processInput(window);
		glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT);


		double CurrTime = glfwGetTime();
		double dt = CurrTime - PrevTime;
		Tick(dt);
		PrevTime = CurrTime;

	//	printf("FPS = %f\n", 1/dt);

		// 4. BIND TEXTURE BEFORE UPDATING
		glBindTexture(GL_TEXTURE_2D, texture);
		glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, cells, cells, GL_RED, GL_FLOAT, divergence);

		// 5. ACTIVATE TEXTURE UNIT BEFORE RENDERING
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, texture);

		glUseProgram(shaderProgram);
		glBindVertexArray(VAO);
		glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);

		glfwSwapBuffers(window);
		glfwPollEvents();
	}

	glDeleteTextures(1, &texture);
	glDeleteProgram(shaderProgram);
	glDeleteBuffers(1, &VBO);
	glDeleteVertexArrays(1, &VAO);

	glfwTerminate();
	return 0;
}