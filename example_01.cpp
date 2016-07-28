#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <string>
#include <cstring>

#ifdef _WIN32
#include <windows.h>
#else
#include <sys/time.h>
#endif

#ifdef OSX
#include <GLUT/glut.h>k
#include <OpenGL/glu.h>
#else
#include <GL/glut.h>
#include <GL/glu.h>
#endif

#include <time.h>
#include <math.h>

using namespace std;

#define PI 3.14159265
inline float sqr(float x) { return x*x; }

/*
WARNING: this program creates a lot of objects via the "new" keword, but deallocates very few of them
i.e, it has a very huge memory waste, but modern computers should reclaim the memory once the program is closed
*/

class Viewport;
class Body;
class Camera;
class Point;
class Path;

void draw();
void drawControl();
void drawControlSphere(Point*, float);
void printMatrix(float**, float, float, char*);
void testInvertMethod();

class Viewport {
public:
	int w, h; // width and height
};

class Point {
public:
	float x, y, z;

	Point(float x, float y, float z, bool normalize = false) {
		if (normalize) {
			float mag = sqrt(x*x+y*y+z*z);
			this->x = x / mag;
			this->y = y / mag;
			this->z = z / mag;
		} else {
			this->x = x;
			this->y = y;
			this->z = z;
		}
	};

	Point* cross(Point* p, bool normal = false){
		return new Point(y*p->z - z*p->y, z*p->x - x*p->z, x*p->y - y*p->x, normal);	
	};

	// normal of triangle with vertices: this + p1 + p2
	Point* normal(Point* p1, Point* p2) {
		float dx1 = p1->x - x;
		float dx2 = p2->x - x;
		float dy1 = p1->y - y;
		float dy2 = p2->y - y;
		float dz1 = p1->z - z;
		float dz2 = p2->z - z;
		return new Point(dy1*dz2 - dz1*dy2, dz1*dx2 - dx1*dz2, dx1*dy2 - dy1*dx2);	
	};

	// distance^2
	float distance2(Point* p) {
		float dx = p->x - x;
		float dy = p->y - y;
		float dz = p->z - z;
		return dx*dx + dy*dy + dz*dz;
	};

	// this - p
	Point* subtract(Point* p){
		return new Point(x - p->x, y - p->y, z - p->z);
	}

	// this = this / f
	void divide(float f) {
		x /= f;
		y /= f;
		z /= f;
	}

	Point* merge(Point* p, float u) {
		float v = 1 - u;
		float x = this->x * v + u * p->x;
		float y = this->y * v + u * p->y;
		float z = this->z * v + u * p->z;
		return new Point(x, y, z);
	}
};

class Body {
public:
	float theta, phi;
	float length;

	Body(float theta, float phi, float length) {
		this->theta = theta;
		this->phi = phi;
		this->length = length;
	};

	void move(float dtheta, float dphi) {
		this->theta += dtheta;
		this->phi += dphi;
	};

	void moveTo(float theta, float phi) {
		this->theta = theta;
		this->phi = phi;
	};

	Point getPos(Point base, float thetaB, float phiB) {
		float phCosnLength = cos(phi + phiB) * length;
		base.x += phCosnLength * cos(theta + thetaB);
		base.z += phCosnLength * sin(theta + thetaB);
		base.y += sin(phi + phiB) * length;
		return base;
	}

};

float maxErr2; // maximum error allowed (squared)
bool infinite;
float totalBodyLength;
int subdivide;

class Path {
public:
	Point** p;
	int c, n;
	int adding;

	Path(int n) {
		this->n = n * subdivide - subdivide + 1;
		p = new Point*[this->n];
		c = 0;
		adding = 0;
	};

	void add(Point* p) {
		if (adding == 0) {
			this->p[0] = p;
			adding = 1;
		} else {
			Point *last = this->p[adding - 1];
			for (int i = 1; i <= subdivide; i++) {
				this->p[adding++] = last->merge(p, i * 1.0 / subdivide);
			}
		}
	}

	//bool isGood(Point* p) {
	//	return p->distance2(this->p[c]) <= maxErr2;
	//	//return getErr2(p) <= err2;
	//};

	void next() {
		c++;
		printf("c is now %d\n", c);
	}

	bool isDone() { 
		return (c == n && !infinite);
	};

	Point* get() {
		if (c == n) {
			if (!infinite) 
				return NULL;
			else {
				int k = totalBodyLength * 1.5f;
				p[--c] = new Point(rand() % k - k / 2, rand() % k - k / 2, rand() % k - k / 2);
			}
		}
		return p[c];
	}

	Point* getNext(int i = 1) {
		if (c + i >= n)
			return NULL;
		return p[c + i];
	}

	float getErr2(Point* p) {
		return p->distance2(this->p[c]);
	};

	void move(float x, float y, float z) {
		if (c == n)
			c--;
		p[n - 1] = new Point(p[c]->x + x, p[c]->y + y, p[c]->z + z);
		c = n - 1;
	}
};

class Camera {
public:
	float ox , oy, oz;
	float tx, ty, tz;
	float r, flatr;
	float theta, thetac, thetas;
	float phi, phic, phis;

	bool wireframe, smooth, control, rotate;

	Camera(float tx, float ty, float tz, float ox, float oy, float oz){
		smooth = true;
		control = true;

		this->tx = tx;
		this->ty = ty;
		this->tz = tz;
		this->ox = ox;
		this->oy = oy;
		this->oz = oz;
		float dx = ox - tx;
		float dy = oy - ty;
		float dz = oz - tz;
		r = sqrt(dx*dx + dy*dy + dz*dz);
		flatr = sqrt(dx*dx + dz*dz);
		theta = atan2(dz, dx);
		thetac = cos(theta);
		thetas = sin(theta);
		phi = atan2(dy, flatr);
		phic = cos(phi);
		phis = sin(phi);
	}
};

Viewport viewport;
char *inputFile;
Camera camera(0, 0, 0, 7, 7, 7);
int mouse, mx, my;
bool shift;
int displayListIndex;
Path* path;
int bodyn;
Body** body;
float step;
float speed;
bool updating;

float update(float locaStep = step);

void initScene(){
	//This tells glut to use a double-buffered window with red, green, and blue channels 
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH); // | GLUT_DEPTH

	// Initalize theviewport size
	viewport.w = 800;
	viewport.h = 800;

	// The size and position of the window
	glutInitWindowSize(viewport.w, viewport.h);
	glutInitWindowPosition(0, 0);
	glutCreateWindow("IOnverse KOTLinematics");

	glEnable(GL_DEPTH_TEST);
	glEnable(GL_COLOR_MATERIAL);

	// lighting
	glEnable(GL_LIGHTING);
	glEnable(GL_NORMALIZE);

	glEnable(GL_LIGHT0);
	GLfloat pos0[] = {10000, 5000, 4000, 1.0f};
	GLfloat amb0[] = {.4f, .4f, .4f, 1.0f};
	GLfloat dif0[] = {1, 1, 1, 1.0f};
	GLfloat spec0[] = {.0f, .0f, .0f, 1.0f};
	glLightfv(GL_LIGHT0, GL_POSITION, pos0);
	glLightfv(GL_LIGHT0, GL_AMBIENT, amb0);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, dif0);
	glLightfv(GL_LIGHT0, GL_SPECULAR, spec0);

	/*glEnable(GL_LIGHT1);
	GLfloat pos1[] = {-30000, 100000, 40000, 1.0f};
	GLfloat amb1[] = {.0f, .0f, .0f, 1.0f};
	GLfloat dif1[] = {1, 1, 1, 1.0f};
	GLfloat spec1[] = {.0f, .8f, .0f, 1.0f};
	glLightfv(GL_LIGHT1, GL_POSITION, pos1);
	glLightfv(GL_LIGHT1, GL_AMBIENT, amb1);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, dif1);
	glLightfv(GL_LIGHT1, GL_SPECULAR, spec1);*/

	// default shading settings
	if (camera.wireframe)
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	else
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	if (camera.smooth) {
		glShadeModel(GL_SMOOTH);
		glPolygonMode(GL_FRONT_AND_BACK, GL_SMOOTH);
	} else {
		glShadeModel(GL_FLAT);
		glPolygonMode(GL_FRONT_AND_BACK, GL_FLAT);
	}
}

void myReshape(int w, int h) {
	int s = min(w, h);

	viewport.w = s;
	viewport.h = s;

	glViewport(0, 0, viewport.w, viewport.h);
}

void myDisplay() {
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// inverse kinemaics
	if (updating) {
		float totalDone = 0;
		float c = 0;
		while (totalDone < speed && c >= 0) {
			c = update();
			totalDone += c;
		}
	}
	//for (int i = 0; i < speed; i++) 
	//update();

	// camera
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(60, 1, 0.1, 100);
	gluLookAt(camera.ox, camera.oy, camera.oz, // look from
		camera.tx, camera.ty, camera.tz, // look at
		0, 1 , 0); // up axis

	//printf("%f, %f, %f, %f, %f, %f \n", camera.ox, camera.oy, camera.oz, camera.tx, camera.ty, camera.tz);

	if (camera.rotate) {
		camera.theta -= .004f;
		camera.thetac = cos(camera.theta);
		camera.thetas = sin(camera.theta);
		camera.phic = cos(camera.phi);
		camera.phis = sin(camera.phi);
		camera.flatr = camera.r * camera.phic;
		camera.ox = camera.tx + camera.flatr * camera.thetac;
		camera.oy = camera.ty + camera.r * camera.phis;
		camera.oz = camera.tz + camera.flatr * camera.thetas;
	}

	if (camera.control) 
		drawControl();

	// Start drawing
	draw();

	glFlush();
	glutSwapBuffers();
}

void myFrameMove() {
	//nothing here for now
#ifdef _WIN32
	Sleep(10); //give ~10ms back to OS (so as not to waste the CPU)
#endif
	glutPostRedisplay(); // forces glut to call the display function (myDisplay())
}

Point* computeEndpoint(float *changes) {
	Point base(0, 0, 0);
	float theta = 0, phi = 0;
	for(int i = 0; i < bodyn; i++) {
		base = body[i]->getPos(base, theta + changes[i * 2], phi + changes[i * 2 + 1]);
		theta += body[i]->theta + changes[i*2];
		phi += body[i]->phi + changes[i * 2 + 1];
	}
	return new Point(base.x, base.y, base.z);
}

float** invert3x3(float** m) {
	// code edited from http://stackoverflow.com/questions/983999/simple-3x3-matrix-inverse-code-c
	float det = m[0][0]*(m[1][1]*m[2][2]-m[2][1]*m[1][2])
		-m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0])
		+m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0]);
	float invdet = 1 / det;
	float **result = new float*[3];
	result[0] = new float[3];
	result[1] = new float[3];
	result[2] = new float[3];
	result[0][0] =  (m[1][1] * m[2][2] - m[2][1] * m[1][2]) * invdet;
	result[0][1] = -(m[0][1] * m[2][2] - m[0][2] * m[2][1]) * invdet;
	result[0][2] =  (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * invdet;
	result[1][0] = -(m[1][0] * m[2][2] - m[1][2] * m[2][0]) * invdet;
	result[1][1] =  (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * invdet;
	result[1][2] = -(m[0][0] * m[1][2] - m[1][0] * m[0][2]) * invdet;
	result[2][0] =  (m[1][0] * m[2][1] - m[2][0] * m[1][1]) * invdet;
	result[2][1] = -(m[0][0] * m[2][1] - m[2][0] * m[0][1]) * invdet;
	result[2][2] =  (m[0][0] * m[1][1] - m[1][0] * m[0][1]) * invdet;
	return result;
}

float update(float stepLocal) {
	//if (stepLocal != step)
	//printf("step: %f\n", stepLocal);

	if (path->isDone())
		return -1;

	// initialize matrices and variables

	float delta = 0.001f * PI / 180;

	float **jacob = new float*[3];
	jacob[0] = new float[bodyn * 2];
	jacob[1] = new float[bodyn * 2];
	jacob[2] = new float[bodyn * 2];

	float **jacobProd =  new float*[3];
	jacobProd[0] = new float[3];
	jacobProd[1] = new float[3];
	jacobProd[2] = new float[3];

	float **jacobInv = new float*[3];
	jacobInv[0] = new float[3];
	jacobInv[1] = new float[3];
	jacobInv[2] = new float[3];

	float **jacobProd2 = new float*[bodyn * 2];
	for (int i = 0; i < bodyn * 2; i++) {
		jacobProd2[i] = new float[3];
	}

	float *changes = new float[bodyn * 2];
	for (int i = 0; i < bodyn * 2; i++) {
		changes[i] = 0;
	}
	Point* end = computeEndpoint(changes);

	// goal endpoint - current endpoint
	Point *distance = path->get()->subtract(end);

	// compute jacobian
	Point* t;
	for (int i = 0; i < bodyn * 2; i++) {
		changes[i] = delta;
		t = computeEndpoint(changes)->subtract(end);
		t->divide(delta);
		jacob[0][i] = t->x;
		jacob[1][i] = t->y;
		jacob[2][i] = t->z;
		changes[i] = 0;
	}


	// compute jacobProd = jacob * jacob^T
	for (int r = 0; r < 3; r++) {
		for (int c = 0; c < 3; c++) {
			float prod = 0;
			for (int i = 0 ; i < bodyn * 2; i++) {
				prod +=	jacob[r][i] * jacob[c][i];
			}
			jacobProd[r][c] = prod;
		}
	}


	// compute jacobInv = (jacob * jacob^T)^-1
	jacobInv = invert3x3(jacobProd);


	// compute jacobProd2 = jacob^T * (jacob * jacob^T)^-1
	for (int r = 0; r < bodyn * 2; r++) {
		for (int c = 0; c < 3; c++) {
			float prod = 0;
			for (int i = 0 ; i < 3; i++) {
				prod +=	jacob[i][r] * jacobProd[i][c];
			}
			jacobProd2[r][c] = prod;
		}
	}


	// compute ratio of changes in parameters (and normalize to step size)
	float mag = 0;
	for (int i = 0; i < bodyn * 2; i++) {
		changes[i] += jacobProd2[i][0] * distance->x;
		changes[i] += jacobProd2[i][1] * distance->y;
		changes[i] += jacobProd2[i][2] * distance->z;
		mag += changes[i] * changes[i];
	}
	mag = stepLocal / sqrt(mag);
	for (int i = 0; i < bodyn * 2; i++) {
		changes[i] *= mag;
	}


	// changes	
	float err2prev = path->getErr2(end);
	float err2post = path->getErr2(computeEndpoint(changes));
	if (err2post < err2prev) {
		// apply changes		
		for(int i = 0; i < bodyn; i++) {
			body[i]->move(changes[i * 2], changes[i * 2 + 1]);

			// re-draw
			displayListIndex = 0;
		}

		// next goal if reached current
		if (err2post <= maxErr2) {
			printf("error %f\n", sqrt(err2prev));
			printf("reached goal\n");
			path->next();
			return -1;
		}
		return stepLocal;
	} else if (stepLocal / 2 > 0) { // try again with smaller step
		return update(stepLocal / 2);
	} else { // unreachable
		printf("unreachable goal\n");
		path->next();
		return -1;
	}

	/*printf("\n");
	printMatrix(jacob, 3, bodyn * 2, "jacob");
	printMatrix(jacobProd, 3, 3, "jacobProd");
	printMatrix(jacobInv, 3, 3, "jacobInv");
	printMatrix(jacobProd2, bodyn * 2, 3, "jacobProd2");*/
}

void draw() {
	if (displayListIndex == 0) {
		displayListIndex = glGenLists(1);
		glNewList(displayListIndex, GL_COMPILE_AND_EXECUTE);
		/*glBegin(GL_LINE_STRIP);
		glColor3f(0, 1, 0);

		glVertex3f(0, 0 , 0);
		Point p(0, 0, 0);
		float theta = 0, phi = 0;
		for (int i = 0; i < bodyn; i++) {
		p = body[i]->getPos(p, theta, phi);
		theta += body[i]->theta;
		phi += body[i]->phi;
		glVertex3f(p.x, p.y, p.z);
		glVertex3f(p.x, p.y + .2f, p.z);
		}

		glEnd();*/

		glDisable(GL_LIGHTING);

		glColor3f(1, 1, 0);
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();

		glBegin(GL_LINES);
		glVertex3f(-5, 0, 0);
		glVertex3f(5, 0, 0);
		glVertex3f(0, 5, 0);
		glVertex3f(0, -5, 0);
		glVertex3f(0, 0, 5);
		glVertex3f(0, 0, -5);
		glEnd();

		glEnable(GL_LIGHTING);

		glColor3f(0, .5f, 0);

		Point p(0, 0, 0);
		Point p2(0, 0, 0);
		float theta = 0, phi = 0;
		for (int i = 0; i < bodyn; i++) {
			p2 = body[i]->getPos(p, theta, phi);
			theta += body[i]->theta;
			phi += body[i]->phi;

			glMatrixMode(GL_MODELVIEW);
			glLoadIdentity();
			glTranslatef(p.x, p.y, p.z);
			glRotatef(-theta/PI*180, 0, 1, 0);
			glRotatef(phi/PI*180, 0, 0, 1);
			glRotatef(90, 0, 1, 0);
			/*glBegin(GL_TRIANGLE_FAN);
			glVertex3f(0, 0, body[i]->length);
			Point *tip = new Point(0, 0, body[i]->length);
			glVertex3f(.3 , 0, 0);
			Point *prev = new Point(.3 , 0, 0);
			for (int i = 15; i <= 360; i += 15) {
			Point* cur = new Point(.3 * cos(i * PI / 180), .3 * sin(i * PI / 180), 0);
			Point* normal = tip->normal(prev, cur);
			glNormal3f(normal->x, normal->y, normal->z);
			glVertex3f(.3 * cos(i * PI / 180), .3 * sin(i * PI / 180), 0);
			prev = cur;
			}
			glEnd();*/
			glutSolidCone(0.3, body[i]->length, 20, 20);

			p = p2;
		}

		glEndList();
	} else
		glCallList(displayListIndex);
}

void drawControl() {
	if (path->isDone())
		return;

	glColor3f(.7f, .3f, .3f);

	drawControlSphere(path->get(), .3f);
	for (int i = 1; i < subdivide + 1; i++) {
		drawControlSphere(path->getNext(i), .15f);
	}
}

void drawControlSphere(Point* goal, float r) {
	if (goal == NULL)
		return;

	// line
	glDisable(GL_LIGHTING);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glBegin(GL_LINES);
	glVertex3f(goal->x, goal->y, goal->z);
	glVertex3f(goal->x, 0, goal->z);
	glVertex3f(goal->x, 0, goal->z);
	glVertex3f(0, 0, goal->z);
	glVertex3f(goal->x, 0, goal->z);
	glVertex3f(goal->x, 0, 0);
	glEnd();

	// sphere
	glEnable(GL_LIGHTING);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glTranslatef(goal->x, goal->y, goal->z);
	glutSolidSphere(r, 30, 30);
}

void initCommands(int argc, char *argv[]) {
	if (argc == 5) {		
		inputFile = argv[1];
		maxErr2 = atof(argv[2]);
		step = atof(argv[3]);
		speed = atof(argv[4]);
		printf("file: %s, tolerence: %f, default step: %f, speed-up: %f\n", inputFile, maxErr2, step, speed);
		maxErr2 = maxErr2 * maxErr2;
	} else {
		printf("bad parmameters\n");
		printf("expected: input_file_name tolerence default_step speed_up\n");
		printf("good sample values: input.txt .1 .0005 .015\n");
		exit(0);
	}
}

void readInputFile() {
	//long t1 = GetTickCount();

	string data;
	Point *p;

	ifstream file;
	file.open (inputFile);


	// read body data

	file >> data;
	bodyn = atoi(data.c_str());
	body = new Body*[bodyn];
	totalBodyLength = 0;
	for (int i = 0; i < bodyn; i++) {
		file >> data;
		float len = atof(data.c_str());
		body[i] = new Body(0, 0, len);
		totalBodyLength += len;
		//printf("body %d has length %f\n", i, body[i]->length);
	}


	// read path data

	file >> data;
	int n = atof(data.c_str());
	file >> data;
	subdivide = atof(data.c_str());
	path = new Path(n);

	for (int i  = 0; i < n; i++) {
		file >> data;
		float x = atof(data.c_str());
		file >> data;
		float y = atof(data.c_str());
		file >> data;
		float z = atof(data.c_str());
		p = new Point(x, y, z);
		path->add(p);
	}

	//long t2 = GetTickCount();
	//printf("time: %d\n", t2 - t1);

	file.close();
}

void randomizeBody() {
	for (int i = 0; i < bodyn; i++) {
		body[i]->theta = (rand() % 360) * PI / 180;
		body[i]->phi = (rand() % 360) * PI / 180;
	}
}

void restart() {
	updating = false;
	readInputFile();
	displayListIndex = 0;
}

void myInput(unsigned char key, int x, int y) {
	switch(key) {
	case 'q':
		exit(0);
		break;
	case 'w':
		camera.wireframe = !camera.wireframe;
		if (camera.wireframe)
			glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		else
			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		break;
	case 's':
		camera.smooth = !camera.smooth;
		if (camera.smooth) {
			glShadeModel(GL_SMOOTH);
			glPolygonMode(GL_FRONT_AND_BACK, GL_SMOOTH);
		} else {
			glShadeModel(GL_FLAT);
			glPolygonMode(GL_FRONT_AND_BACK, GL_FLAT);
		}
		break;
	case 'c':
		camera.control = !camera.control;
		break;
	case ' ':
		updating = !updating;
		break;
	case 'r':
		restart();
		break;
	case 't':
		camera.rotate = !camera.rotate;
		break;
	case 'o':
		updating = false;
		randomizeBody();
		displayListIndex = 0;
		break;
	case 'i':
		infinite = !infinite;
		break;
	case '=':
	case '+':
		speed += .015f;
		printf("speed: %f\n", speed);
		break;
	case '-':
		speed -= .015f;
		if (speed < 0)
			speed = 0;
		printf("speed: %f\n", speed);
		break;
	case ',':
		path->move(0, -1, 0);
		break;
	case '.':
		path->move(0, 1, 0);
		break;
	}
}

void mySpecialInput(int key, int x, int y) {
	switch(key)	{
	case GLUT_KEY_UP:
		path->move(0, 0, -1);
		break;	
	case GLUT_KEY_DOWN:
		path->move(0, 0, 1);
		break;
	case GLUT_KEY_LEFT:
		path->move(-1, 0, 0);
		break;
	case GLUT_KEY_RIGHT:
		path->move(1, 0, 0);
		break;
	}
}

void myMouseInput(int button, int state, int x, int y) {
	shift = glutGetModifiers() == GLUT_ACTIVE_SHIFT;

	if (state == 1)
		mouse = -1;
	else {
		mouse = button;
		mx = x;
		my = y;
	}
}

void myMotionInput(int x, int y) {
	float dy = (y - my) / 10.;
	float dx = (x - mx) / 10.;
	my = y;
	mx = x;

	if (shift) {
		switch (mouse) {
		case 0:
			path->move(dx / 5, 0, dy / 5);
			break;
		case 2:
			path->move(0, -dy / 5, 0);
			break;
		}
	}
	else {
		switch (mouse) {
		case 0:
			//printf("theta %f, phi %f, r %f, flatr %f, ox %f, oy %f, oz %f, tx %f, ty %f, tz %f\n", camera.theta, camera.phi, camera.r, camera. flatr, camera.ox, camera.oy, camera.oz, camera.tx, camera.ty, camera.tz);
			camera.theta += dx / 10.;
			camera.phi += dy / 10.;
			if (camera.phi > PI / 2 - .15)
				camera.phi = PI / 2 - .15;
			if (camera.phi < -PI / 2 + .15)
				camera.phi = -PI / 2 + .15;
			camera.thetac = cos(camera.theta);
			camera.thetas = sin(camera.theta);
			camera.phic = cos(camera.phi);
			camera.phis = sin(camera.phi);
			camera.flatr = camera.r * camera.phic;
			camera.ox = camera.tx + camera.flatr * camera.thetac;
			camera.oy = camera.ty + camera.r * camera.phis;
			camera.oz = camera.tz + camera.flatr * camera.thetas;
			//printf("theta %f, phi %f, r %f, flatr %f, ox %f, oy %f, oz %f, tx %f, ty %f, tz %f\n", camera.theta, camera.phi, camera.r, camera. flatr, camera.ox, camera.oy, camera.oz, camera.tx, camera.ty, camera.tz);
			break;
		case 1:
			camera.r += dy;
			camera.flatr = camera.r * camera.phic;
			camera.ox = camera.tx + camera.flatr * camera.thetac;
			camera.oy = camera.ty + camera.r * camera.phis;
			camera.oz = camera.tz + camera.flatr * camera.thetas;
			break;
		case 2:
			float moveh = dy * camera.phis;
			float movev = dy * camera.phic;
			camera.ty += movev;
			camera.oy += movev;
			float movex = dx * camera.thetas + moveh * camera.thetac;
			camera.tx -= movex;
			camera.ox -= movex;
			float movez = dx * camera.thetac - moveh * camera.thetas;
			camera.tz += movez;
			camera.oz += movez;
			break;
		}
	}
}

int main(int argc, char *argv[]) {
	//testInvertMethod();

	initCommands(argc, argv);
	readInputFile();
	//initRandBody();

	glutInit(&argc, argv);

	initScene(); // quick function to set up scene

	glutDisplayFunc(myDisplay); // function to run when its time to draw something
	glutReshapeFunc(myReshape); // function to run when the window gets resized
	glutIdleFunc(myFrameMove); // function to run when not handling any other task
	glutKeyboardFunc(myInput);
	glutSpecialFunc(mySpecialInput);
	glutMouseFunc(myMouseInput);
	glutMotionFunc(myMotionInput);
	glutMainLoop(); // infinite loop that will keep drawing and resizing and whatever else

	return 0;
}

void testInvertMethod() {
	float **result = new float*[3];
	result[0] = new float[3];
	result[1] = new float[3];
	result[2] = new float[3];

	for (int r = 0; r < 3; r++) {
		for (int c = 0; c < 3; c++) {
			result[r][c] = rand() % 10;
		}
	}

	float t0[] = {0.0135529, -0.776603, -0.776595};
	float t1[] = {-0.776603, 89., 0.};
	float t2[] = {-0.776595, 0., 89.};
	result[0] = t0;
	result[1] = t1;
	result[2] = t2;

	printMatrix(result, 3, 3, "matrix");

	result = invert3x3(result);

	printMatrix(result, 3, 3, "inverse");

}

void printMatrix(float** a, float rows, float cols, char* name) {
	printf("\n%s :\n", name);
	for (int r = 0; r < rows; r++) {
		for (int c = 0; c < cols; c++) {
			printf(" %f ", a[r][c]);
		}
		printf("\n");
	}
	printf("\n");
}