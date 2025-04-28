/*

Functions for learnply

Eugene Zhang, 2005
*/



#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#if _WIN32 || __APPLE__
#define NOMINMAX
#include <glut.h>
#else
#include <GL/freeglut.h>
#endif
#include <string.h>
#include <fstream>
#include <algorithm>
#include "icVector.H"
#include "icMatrix.H"
#include "polyhedron.h"
#include "trackball.h"
#include "tmatrix.h"

FILE *this_file;
const int win_width=1024;
const int win_height=1024;
 
int display_mode = 0;
int view_mode = 0;  // 0 = othogonal, 1=perspective
double radius_factor = 0.9; // zoom factor orthogonal view
float s_old, t_old;
float rotmat[4][4];
static Quaternion rvec;
int orientation = 0; // 0=ccw, 1=cw

int mouse_mode = -2;  // -2=no action, -1 = down
int mouse_button = -1; // -1=no button, 0=left, 1=middle, 2=right
int last_x, last_y;

/* Project 1, Problem 1d */
double L = 1.0;

/* Project 2, Problem 1 */
int silhouette_mode = 0; // 0=no silhouette, 1=edge silhouette, 2=face silhouette

/* Project 2, Problem 2 */
int tensor_display_mode = 0; // 0=no tensor, 1=major, 2=minor, 3=cross
int smoothing_scheme = 0; // 0 = uniform weights, 1 = cord weights, 2 = flow weights, 3 = mean weights
double smoothing_step = 0.5; // smoothing step size
int smoothing_iterations = 10; // number of smoothing iterations

/* Project 2, Problem 3 */
int render_mode = 0; // 0=surface, 1=pen ink
int hatch_tracing_mode = 0; // 0=current triangle curvature, 1=project from previous triangle

/* Project 2, Debugging */
int wireframe_mode = 0; // 0=no wireframe, 1=wireframe

int ACSIZE = 1; // Passes for anti-aliasing
struct jitter_struct{
	double x;
	double y;
} jitter_para;
jitter_struct ji1[1] = {{0.0, 0.0}};
jitter_struct ji16[16] = {{0.125, 0.125}, {0.375, 0.125}, {0.625, 0.125}, {0.875, 0.125}, 
						  {0.125, 0.375}, {0.375, 0.375}, {0.625, 0.375}, {0.875, 0.375}, 
						  {0.125, 0.625}, {0.375, 0.625}, {0.625, 0.625}, {0.875, 0.625}, 
						  {0.125, 0.875}, {0.375, 0.875}, {0.625, 0.875}, {0.875, 0.875}, };

Polyhedron *poly;

// callback function forward declaration
void init(void);
void keyboard(unsigned char key, int x, int y);
void motion(int x, int y);
void display(void);
void mouse(int button, int state, int x, int y);
void display_shape(GLenum mode, Polyhedron *poly);

/******************************************************************************
Main program.
******************************************************************************/

int main(int argc, char *argv[])
{
	char *progname;
	int num = 1;
	FILE *this_file;

  	progname = argv[0];

	this_file = fopen("../tempmodels/torus.ply", "r");
	poly = new Polyhedron (this_file);
	fclose(this_file);
	mat_ident( rotmat );	

	poly->initialize(); // initialize everything

	/* Project 1, problem 2 */
	poly->create_corners();

	poly->calc_bounding_sphere();
	poly->calc_face_normals_and_area();
	poly->average_normals();

	/* Project 1, problem 3a-d */
	poly->compute_euler_characteristic();
	poly->compute_gaussian_curvature_angle_deficit();
	poly->compute_gaussian_curvature_valence_deficit();
	poly->compute_handles();

	/* Project 2, Problem */
	poly->compute_vert_voronoi_areas();
	poly->compute_vert_mean_curvature();
	poly->compute_vert_gaussian_curvature();
	poly->compute_vert_curvature_tensor();
	poly->smooth_vert_curvature_tensors(smoothing_scheme, smoothing_step, smoothing_iterations);
	poly->compute_vert_principal_curvatures();

	glutInit(&argc, argv);
	glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowPosition (20, 20);
	glutInitWindowSize (win_width, win_height); 
	glutCreateWindow ("Geometric Modeling");
	init ();
	glutKeyboardFunc (keyboard);
	glutDisplayFunc(display); 
	glutMotionFunc (motion);
	glutMouseFunc (mouse);
	glutMainLoop(); 
	poly->finalize();  // finalize everything

  return 0;    /* ANSI C requires main to return int. */
}

void init(void) 
{
  /* select clearing color */ 
  glClearColor (0.0, 0.0, 0.0, 0.0);  // background
  glShadeModel (GL_FLAT);
  glPolygonMode(GL_FRONT, GL_FILL);

  glDisable(GL_DITHER);
	glEnable(GL_DEPTH_TEST);
  glDepthFunc(GL_LESS);
	// may need it
  glPixelStorei(GL_PACK_ALIGNMENT,1);
	glEnable(GL_NORMALIZE);
	if (orientation == 0) 
		glFrontFace(GL_CW);
	else 
		glFrontFace(GL_CCW);
}

void multmatrix(const Matrix m)
{ 
  int i,j, index = 0;

  GLfloat mat[16];

  for ( i = 0; i < 4; i++)
    for (j = 0; j < 4; j++)
      mat[index++] = m[i][j];

  glMultMatrixf (mat);
}

void set_view(GLenum mode, Polyhedron *poly)
{
	icVector3 up, ray, view;
	GLfloat light_ambient0[] = { 0.3, 0.3, 0.3, 1.0 };
	GLfloat light_diffuse0[] = { 0.7, 0.7, 0.7, 1.0 };
	GLfloat light_specular0[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat light_ambient1[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat light_diffuse1[] = { 0.5, 0.5, 0.5, 1.0 };
	GLfloat light_specular1[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat light_ambient2[] = { 1.0, 1.0, 1.0, 1.0 };
	GLfloat light_diffuse2[] = { 1.0, 1.0, 1.0, 1.0 };
	GLfloat light_specular2[] = { 1.0, 1.0, 1.0, 1.0 };
	GLfloat light_position[] = { 0.0, 0.0, 0.0, 1.0 };

	glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient0);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse0);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular0);
	glLightfv(GL_LIGHT1, GL_AMBIENT, light_ambient1);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, light_diffuse1);
	glLightfv(GL_LIGHT1, GL_SPECULAR, light_specular1);


	glMatrixMode(GL_PROJECTION);
	if (mode == GL_RENDER)
		glLoadIdentity();

	if (view_mode == 0)
		glOrtho(-radius_factor, radius_factor, -radius_factor, radius_factor, 0.0, 40.0);
	else
		gluPerspective(45.0, 1.0, 0.1, 40.0);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	light_position[0] = 5.5;
	light_position[1] = 0.0;
	light_position[2] = 0.0;
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);
	light_position[0] = -0.1;
	light_position[1] = 0.0;
	light_position[2] = 0.0;
	glLightfv(GL_LIGHT2, GL_POSITION, light_position);
}

void set_scene(GLenum mode, Polyhedron *poly)
{
	glTranslatef(0.0, 0.0, -3.0);
	multmatrix( rotmat );

	glScalef(1.0/poly->radius, 1.0/poly->radius, 1.0/poly->radius);
	glTranslatef(-poly->center.entry[0], -poly->center.entry[1], -poly->center.entry[2]);
}

void motion(int x, int y) {
	float r[4];
	float xsize, ysize, s, t;

	switch(mouse_mode){
	case -1:

		xsize = (float) win_width;
		ysize = (float) win_height;
	
		s = (2.0 * x - win_width) / win_width;
		t = (2.0 * (win_height - y) - win_height) / win_height;

		if ((s == s_old) && (t == t_old))
			return;

		mat_to_quat( rotmat, rvec );
		trackball( r, s_old, t_old, s, t );
		add_quats( r, rvec, rvec );
		quat_to_mat( rvec, rotmat );

		s_old = s;
		t_old = t;

		display();
		break;
	}
}

int processHits(GLint hits, GLuint buffer[])
{
	unsigned int i, j;
	GLuint names, *ptr;
	double smallest_depth=1.0e+20, current_depth;
	int seed_id=-1; 
	unsigned char need_to_update;

	printf("hits = %d\n", hits);
	ptr = (GLuint *) buffer;
	for (i = 0; i < hits; i++) {  /* for each hit  */
		need_to_update = 0;
		names = *ptr;
		ptr++;
		
		current_depth = (double) *ptr/0x7fffffff;
		if (current_depth < smallest_depth) {
			smallest_depth = current_depth;
			need_to_update = 1;
		}
		ptr++;
		current_depth = (double) *ptr/0x7fffffff;
		if (current_depth < smallest_depth) {
			smallest_depth = current_depth;
			need_to_update = 1;
		}
		ptr++;
		for (j = 0; j < names; j++) {  /* for each name */
			if (need_to_update == 1)
				seed_id = *ptr - 1;
			ptr++;
		}
	}
	printf("triangle id = %d\n", seed_id);
	return seed_id;
}

void mouse(int button, int state, int x, int y) {
	if (button == GLUT_LEFT_BUTTON || button == GLUT_RIGHT_BUTTON) {
		switch(mouse_mode) {
		case -2:  // no action
			if (state == GLUT_DOWN) {
				float xsize = (float) win_width;
				float ysize = (float) win_height;

				float s = (2.0 * x - win_width) / win_width;
				float t = (2.0 * (win_height - y) - win_height) / win_height;

				s_old = s;
				t_old = t;

				mouse_mode = -1;  // down
				mouse_button = button;
				last_x = x;
				last_y = y;
			}
			break;

		default:
			if (state == GLUT_UP) {
				button = -1;
				mouse_mode = -2;
			}
			break;
		}
	} else if (button == GLUT_MIDDLE_BUTTON) {
		if (state == GLUT_DOWN) {  // build up the selection feedback mode

			GLuint selectBuf[win_width];
		  GLint hits;
		  GLint viewport[4];

		  glGetIntegerv(GL_VIEWPORT, viewport);

			glSelectBuffer(win_width, selectBuf);
		  (void) glRenderMode(GL_SELECT);

		  glInitNames();
		  glPushName(0);

		  glMatrixMode(GL_PROJECTION);
	    glPushMatrix();
			glLoadIdentity();
/*  create 5x5 pixel picking region near cursor location */
	    gluPickMatrix((GLdouble) x, (GLdouble) (viewport[3] - y),
                 1.0, 1.0, viewport);

			set_view(GL_SELECT, poly);
			glPushMatrix ();
			set_scene(GL_SELECT, poly);
			display_shape(GL_SELECT, poly);
	    glPopMatrix();
		  glFlush();

	    hits = glRenderMode(GL_RENDER);
		  poly->seed = processHits(hits, selectBuf);
			display();
		}
	}
}


/******************************************************************************
Process a keyboard action.  In particular, exit the program when an
"escape" is pressed in the window.
******************************************************************************/

void keyboard(unsigned char key, int x, int y) {
	int i;

	/* set escape key to exit */
	switch (key)
	{
	case 27:
		poly->finalize();  // finalize_everything
		exit(0);
		break;

	case '0':
		display_mode = 0;
		display();
		break;

	case '1':
		display_mode = 1;
		display();
		break;

	case '2':
		display_mode = 2;
		display();
		break;

	case '3':
		display_mode = 3;
		display();
		break;

	case '4':
		display_mode = 4;
		display();
		break;

	case '5':
		display_mode = 5;
		display();
		break;

	case '6':
		display_mode = 6;
		display();
		break;

	case '7':
		display_mode = 7;
		display();
		break;

	case '8':
		display_mode = 8;
		display();
		break;

	case '9':
		display_mode = 9;
		display();
		break;

	/* toggle view mode */
	case 'o':
		view_mode = (view_mode + 1) % 2;
		printf("view_mode=%d\n", view_mode);
		display();
		break;

	/* Anti-aliasing */
	case 'x':
		switch (ACSIZE) {
		case 1:
			ACSIZE = 16;
			break;

		case 16:
			ACSIZE = 1;
			break;

		default:
			ACSIZE = 1;
			break;
		}
		printf("ACSIZE=%d\n", ACSIZE);
		display();
		break;

	/* Project 1, Problem 1d */
	case 'p':
		L *= 2.0;
		printf("L=%g\n", L);
		display();
		break;
	case 'l':
		L *= 0.5;
		printf("L=%g\n", L);
		display();
		break;

	/* Project 2, Problem 1 */
	case 's':
		silhouette_mode = (silhouette_mode + 1) % 3;
		printf("silhouette_mode=%d\n", silhouette_mode);
		display();
		break;

	/* Project 2, Problem 2 */
	case 't':
		tensor_display_mode = (tensor_display_mode + 1) % 4;
		printf("tensor_display_mode=%d\n", tensor_display_mode);
		display();
		break;

	case 'k':
		render_mode = (render_mode + 1) % 2;
		display_mode = 0;
		printf("render_mode=%d\n", render_mode);
		display();
		break;

	/* Project2, Debugging */
	case 'w':
		wireframe_mode = (wireframe_mode + 1) % 2;
		printf("wireframe_mode=%d\n", wireframe_mode);
		display();
		break;

	default:
		fprintf(stderr, "Unknown key pressed\n");
		break;
	}
}

static void color_from_scalar(GLdouble *color, GLdouble min, GLdouble max, GLdouble value)
{
	color[0] = color[1] = color[2] = 1.0;
	GLdouble max_mag = std::max(std::abs(min), std::abs(max));
	GLdouble mapped_value = value / max_mag;
	if (mapped_value < 0.0) // make color blue
	{
		color[0] += mapped_value;
		color[1] += mapped_value;
	}
	else // make color red
	{
		color[1] -= mapped_value;
		color[2] -= mapped_value;
	}
}

void display_shape(GLenum mode, Polyhedron *this_poly)
{
	unsigned int i, j;
	GLfloat mat_diffuse[4];

	glEnable (GL_POLYGON_OFFSET_FILL);
	glPolygonOffset (1., 1.);

	glEnable(GL_DEPTH_TEST);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glShadeModel(GL_SMOOTH);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHT1);
	glEnable(GL_LINE_SMOOTH);

	// draw the polyhedron mesh
	for (i = 0; i < this_poly->ntris; i++)
	{
		if (mode == GL_SELECT)
			glLoadName(i + 1);

		Triangle* temp_t = this_poly->tlist[i];

		switch (display_mode)
		{

		/* default vis */
		case 0:
		{
			if (i == this_poly->seed) {
				mat_diffuse[0] = 0.0;
				mat_diffuse[1] = 0.0;
				mat_diffuse[2] = 1.0;
				mat_diffuse[3] = 1.0;
			}
			else {
				mat_diffuse[0] = 1.0;
				mat_diffuse[1] = 1.0;
				mat_diffuse[2] = 0.0;
				mat_diffuse[3] = 1.0;
			}
			glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
			glBegin(GL_TRIANGLES);
			for (j = 0; j < 3; j++) {

				Vertex* v = temp_t->verts[j];
				glNormal3dv(v->normal.entry);
				glVertex3d(v->x, v->y, v->z);
			}
			glEnd();
			break;
		}

		/* Project 1, Problem 1a */
		case 1:
		{
			static const GLfloat color_map[][4] = {
				{ 1.0f, 0.0f, 0.0f, 1.0f },
				{ 0.0f, 1.0f, 0.0f, 1.0f },
				{ 0.0f, 0.0f, 1.0f, 1.0f },
				{ 1.0f, 1.0f, 0.0f, 1.0f },
				{ 0.0f, 1.0f, 1.0f, 1.0f },
				{ 1.0f, 0.0f, 1.0f, 1.0f }
			};
			static const int ncolors = sizeof(color_map) / sizeof(color_map[0]);

			glMaterialfv(GL_FRONT, GL_DIFFUSE, color_map[temp_t->index % ncolors]);
			glBegin(GL_TRIANGLES);
			for (j = 0; j < 3; ++j)
			{
				glNormal3dv(temp_t->verts[j]->normal.entry);
				glVertex3d(temp_t->verts[j]->x, temp_t->verts[j]->y, temp_t->verts[j]->z);
			}
			glEnd();

			break;
		}

		/* Project 1, Problem 1b */
		case 2:
		{
			static const GLfloat color_map[][4] = {
				{ 1.0f, 0.0f, 0.0f, 1.0f },
				{ 0.0f, 1.0f, 0.0f, 1.0f },
				{ 0.0f, 0.0f, 1.0f, 1.0f }
			};

			glBegin(GL_TRIANGLES);
			for (j = 0; j < 3; ++j)
			{
				glMaterialfv(GL_FRONT, GL_DIFFUSE, color_map[j]);
				glNormal3dv(temp_t->verts[j]->normal.entry);
				glVertex3d(temp_t->verts[j]->x, temp_t->verts[j]->y, temp_t->verts[j]->z);
			}
			glEnd();
			break;
		}

		/* Project 1, Problem 1c */
		case 3:
		{
			GLfloat color[4] = {
				(GLfloat)abs(temp_t->normal.x),
				(GLfloat)abs(temp_t->normal.y),
				(GLfloat)abs(temp_t->normal.z),
				1.0f
			};

			glMaterialfv(GL_FRONT, GL_DIFFUSE, color);
			glDisable(GL_LIGHTING);
			glColor3fv(color);

			glBegin(GL_TRIANGLES);
			for (j = 0; j < 3; ++j)
			{
				glNormal3dv(temp_t->verts[j]->normal.entry);
				glVertex3d(temp_t->verts[j]->x, temp_t->verts[j]->y, temp_t->verts[j]->z);
			}
			glEnd();

			break;
		}

		/* Project 1, Problem 1d (554 only) */
		case 4:
		{
			glBegin(GL_TRIANGLES);
			for (j = 0; j < 3; ++j)
			{
				Vertex* v = temp_t->verts[j];
				GLfloat color[4] = {
					1.0f - (GLfloat)abs(fmod(floor(v->x / L), 2)),
					1.0f - (GLfloat)abs(fmod(floor(v->y / L), 2)),
					1.0f - (GLfloat)abs(fmod(floor(v->z / L), 2)),
					1.0f
				};

				glMaterialfv(GL_FRONT, GL_DIFFUSE, color);
				glNormal3dv(v->normal.entry);
				glVertex3d(v->x, v->y, v->z);
			}
			glEnd();
			break;
		}

		/* Project 1, Problem 3 (Gaussian curvature) */
		case 5:
		{
			glDisable(GL_LIGHTING);

			glBegin(GL_TRIANGLES);
			for (j = 0; j < 3; ++j)
			{
				Vertex* v = temp_t->verts[j];
				GLdouble color[3];

				color_from_scalar(color, poly->K_min, poly->K_max, v->K);
				glColor3dv(color);
				glVertex3d(v->x, v->y, v->z);
			}
			glEnd();
			break;
		}

		/* Project 1, Problem 3 (valence deficit) */
		case 6:
		{
			glDisable(GL_LIGHTING);

			glBegin(GL_TRIANGLES);
			for (j = 0; j < 3; ++j)
			{
				Vertex* v = temp_t->verts[j];
				GLdouble color[3];

				color_from_scalar(color, poly->min_valence_deficit, 6.0, v->valence_deficit);
				glColor3dv(color);
				glVertex3d(v->x, v->y, v->z);
			}
			glEnd();
			break;
		}

		/* Project 2, Problem 2 (mean curvature) */
		case 7:
		{


			break;
		}

		/* Project 2, Problem 2 (gaussian curvature) */
		case 8:
		{
			

			break;
		}

		/* Project 2, Problem 2-3 (curvature tensor) */
		case 9:
		{
			

			break;
		}

		}
	}

	// draw the silhouette if enabled
	if (silhouette_mode != 0)
	{
		
	}

	// draw the curvature tensor crosses if enabled
	if (tensor_display_mode > 0)
	{
		
	}

	// draw wireframe if enabled
	if (wireframe_mode)
	{
		// draw the polyhedron mesh edges
		glDisable(GL_LIGHTING);
		glLineWidth(1.0);
		glColor3f(0.0, 0.0, 0.0);
		glBegin(GL_LINES);
		for (int i = 0; i < poly->nedges; i++)
		{
			Edge* e = poly->elist[i];
			glVertex3d(e->verts[0]->x, e->verts[0]->y, e->verts[0]->z);
			glVertex3d(e->verts[1]->x, e->verts[1]->y, e->verts[1]->z);
		}
		glEnd();
	}
}

void display_pen_ink(Polyhedron* this_poly)
{

}

void display(void)
{
	GLint viewport[4];
	int jitter;

	glClearColor (1.0, 1.0, 1.0, 1.0);  // background for rendering color coding and lighting
	glGetIntegerv (GL_VIEWPORT, viewport);
 
	glClear(GL_ACCUM_BUFFER_BIT);
	for (jitter = 0; jitter < ACSIZE; jitter++) 
	{
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		set_view(GL_RENDER, poly);
		glPushMatrix ();
		switch(ACSIZE)
		{
		case 1:
			glTranslatef (ji1[jitter].x*2.0/viewport[2], ji1[jitter].y*2.0/viewport[3], 0.0);
			break;

		case 16:
			glTranslatef (ji16[jitter].x*2.0/viewport[2], ji16[jitter].y*2.0/viewport[3], 0.0);
			break;

		default:
			glTranslatef (ji1[jitter].x*2.0/viewport[2], ji1[jitter].y*2.0/viewport[3], 0.0);
			break;
		}
		set_scene(GL_RENDER, poly);
		if (render_mode == 0)
			display_shape(GL_RENDER, poly);
		else if (render_mode == 1)
			display_pen_ink(poly);
		glPopMatrix ();
		glAccum(GL_ACCUM, 1.0/ACSIZE);
	}
	glAccum (GL_RETURN, 1.0);
	glFlush();
	glutSwapBuffers();
	glFinish();
}