/*

Functions for learnply

Eugene Zhang, 2005
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <vector>
#include <algorithm>
#if _WIN32 || __APPLE__
#include <glut.h>
#else
#include <GL/freeglut.h>
#endif
#include <string.h>
#include <fstream>
#include "ply.h"
#include "icVector.H"
#include "icMatrix.H"
#include "learnply.h"
#include "learnply_io.h"
#include "trackball.h"
#include "tmatrix.h"
#include "noise.h"

static PlyFile *in_ply;

unsigned char orientation;  // 0=ccw, 1=cw

FILE *this_file;
const int win_width=1024;
const int win_height=1024;

double radius_factor = 0.9;

int display_mode = 0; 
double error_threshold = 1.0e-13;
char reg_model_name[128];
FILE *f;	
int ACSIZE = 1; // for antialiasing
int view_mode=0;  // 0 = othogonal, 1=perspective
float s_old, t_old;
float rotmat[4][4];
static Quaternion rvec;

//globals I added: 
float r, g, b;     // colors lol 
float checker_L = 0.2f; // You can adjust this value for different checker sizes

int mouse_mode = -2;  // -2=no action, -1 = down, 0 = zoom, 1 = rotate x, 2 = rotate y, 3 = tranlate x, 4 = translate y, 5 = cull near 6 = cull far
int mouse_button = -1; // -1=no button, 0=left, 1=middle, 2=right
int last_x, last_y;

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

// Watercolor NPR Parameters
int wcolor_mode = 0;    // 0=off, 1=on
float color_intensity = 1.2; // Controls the intensity of the watercolor effect
float edge_darkening = 0.8;  // Controls how much edges are darkened
float noise_scale = 0.1;     // Scale of the noise effect
float granulation_scale = 0.2; // Scale of the granulation effect

// Add at the top with other global variables
bool watercolor_mode = false;  // Toggle for watercolor rendering

// Change const parameters to regular variables
float DT = 0.15f;              // Time step for fluid movement
float VISCOSITY = 0.05f;       // Fluid viscosity
float EVAPORATION_RATE = 0.03f; // Water evaporation rate

//function declarations: 
void initialize_pigments_from_geometry();
float noise(float x, float y, float z);

// Fluid simulation parameters
const int GRID_SIZE = 256;

// Simple vec4 structure for color compositing
struct vec4 {
    float r, g, b, a;
    vec4(float _r, float _g, float _b, float _a) : r(_r), g(_g), b(_b), a(_a) {}
};

struct FluidCell {
    float height;        // Water height
    float pigment[3];    // RGB pigment concentration
    float velocity[2];   // Velocity field (u,v)
    bool wet_mask;       // Wet area mask
};

class WatercolorSimulation {
public:
    FluidCell* grid;  // Made public for direct access
    int size;

    WatercolorSimulation(int grid_size) : size(grid_size) {
        grid = new FluidCell[grid_size * grid_size];
        temp_grid = new FluidCell[grid_size * grid_size];
        reset();
    }

    ~WatercolorSimulation() {
        delete[] grid;
        delete[] temp_grid;
    }

    void simulate() {
        // 1. Solve shallow water (fluid movement)
        solve_shallow_water();
        
        // 2. Advect pigment
        advect_pigment();
        
        // 3. Apply edge darkening based on water height gradient
        apply_edge_darkening();
        
        // 4. Apply evaporation
        apply_evaporation();
    }

    vec4 km_composite(vec4 bg, vec4 fg) {
        float Sa = fg.a * (1.0f - bg.a);
        float Sr = fg.r * Sa + bg.r * bg.a;
        float Sg = fg.g * Sa + bg.g * bg.a;
        float Sb = fg.b * Sa + bg.b * bg.a;
        return vec4(Sr, Sg, Sb, Sa + bg.a);
    }

    // Add pigment at a specific location
    void add_pigment(int x, int y, float r, float g, float b, float amount) {
        if (x >= 0 && x < size && y >= 0 && y < size) {
            int idx = y * size + x;
            grid[idx].pigment[0] += r * amount;
            grid[idx].pigment[1] += g * amount;
            grid[idx].pigment[2] += b * amount;
            grid[idx].height += amount;
            grid[idx].wet_mask = true;  // Mark as wet when adding pigment
        }
    }

    void fill_empty_grid_cells() {
        // Fill empty cells with interpolated values from neighbors
        for (int y = 1; y < size-1; y++) {
            for (int x = 1; x < size-1; x++) {
                int idx = y * size + x;
                
                // If cell is empty, interpolate from neighbors
                if (grid[idx].pigment[0] == 0.0f && grid[idx].pigment[1] == 0.0f && grid[idx].pigment[2] == 0.0f) {
                    int count = 0;
                    float total_r = 0, total_g = 0, total_b = 0;
                    
                    // Check 8-connected neighbors
                    for (int dy = -1; dy <= 1; dy++) {
                        for (int dx = -1; dx <= 1; dx++) {
                            if (dx == 0 && dy == 0) continue;
                            int nidx = (y+dy) * size + (x+dx);
                            if (grid[nidx].pigment[0] > 0 || grid[nidx].pigment[1] > 0 || grid[nidx].pigment[2] > 0) {
                                total_r += grid[nidx].pigment[0];
                                total_g += grid[nidx].pigment[1];
                                total_b += grid[nidx].pigment[2];
                                count++;
                            }
                        }
                    }
                    
                    if (count > 0) {
                        grid[idx].pigment[0] = total_r / count * 0.3f;  // Reduced for smooth transition
                        grid[idx].pigment[1] = total_g / count * 0.3f;
                        grid[idx].pigment[2] = total_b / count * 0.3f;
                    }
                }
            }
        }
    }

    void set_parameters(float viscosity, float evaporation) {
        VISCOSITY = viscosity;
        EVAPORATION_RATE = evaporation;
    }

private:
    FluidCell* temp_grid;

    void reset() {
        for (int i = 0; i < size * size; i++) {
            grid[i].height = 0.0f;
            grid[i].pigment[0] = grid[i].pigment[1] = grid[i].pigment[2] = 0.0f;
            grid[i].velocity[0] = grid[i].velocity[1] = 0.0f;
            grid[i].wet_mask = false;
        }
    }

    // Enhanced edge darkening
    void apply_edge_darkening() {
        // Calculate water height gradients for edge detection
        std::vector<float> height_gradients(size * size, 0.0f);
        
        // Compute gradients
        for (int y = 1; y < size-1; y++) {
            for (int x = 1; x < size-1; x++) {
                int idx = y * size + x;
                if (grid[idx].wet_mask) {
                    // Calculate gradient magnitude
                    float dx = grid[idx+1].height - grid[idx-1].height;
                    float dy = grid[idx+size].height - grid[idx-size].height;
                    height_gradients[idx] = sqrt(dx*dx + dy*dy);
                }
            }
        }
        
        // Apply edge darkening based on gradients
        for (int y = 0; y < size; y++) {
            for (int x = 0; x < size; x++) {
                int idx = y * size + x;
                if (grid[idx].wet_mask) {
                    float gradient = height_gradients[idx];
                    float darken_factor = 1.0f + (edge_darkening * gradient);
                    
                    // Apply darkening to pigment
                    grid[idx].pigment[0] *= darken_factor;
                    grid[idx].pigment[1] *= darken_factor;
                    grid[idx].pigment[2] *= darken_factor;
                }
            }
        }
    }

    // Add evaporation
    void apply_evaporation() {
        for (int i = 0; i < size * size; i++) {
            grid[i].height *= (1.0f - EVAPORATION_RATE * DT);
            grid[i].height = std::max(0.0f, grid[i].height);
            
            // Clear wet mask when fully dry
            if (grid[i].height < 0.001f) {
                grid[i].wet_mask = false;
            }
        }
    }

    // Shallow water equations solver
    void solve_shallow_water() {
        // Copy current state
        memcpy(temp_grid, grid, size * size * sizeof(FluidCell));

        // Solve height field
        for (int y = 1; y < size-1; y++) {
            for (int x = 1; x < size-1; x++) {
                int idx = y * size + x;
                float h = grid[idx].height;
                float u = grid[idx].velocity[0];
                float v = grid[idx].velocity[1];

                // Update height based on velocity divergence
                float div = (grid[idx+1].velocity[0] - grid[idx-1].velocity[0] +
                           grid[idx+size].velocity[1] - grid[idx-size].velocity[1]) * 0.5f;
                temp_grid[idx].height = h - DT * div;

                // Update velocity with viscosity
                temp_grid[idx].velocity[0] = u + VISCOSITY * (
                    grid[idx+1].velocity[0] + grid[idx-1].velocity[0] +
                    grid[idx+size].velocity[0] + grid[idx-size].velocity[0] - 4*u);
                temp_grid[idx].velocity[1] = v + VISCOSITY * (
                    grid[idx+1].velocity[1] + grid[idx-1].velocity[1] +
                    grid[idx+size].velocity[1] + grid[idx-size].velocity[1] - 4*v);
             }
        }

        // Swap grids
        std::swap(grid, temp_grid);
    }

    // Advect pigment based on velocity field
    void advect_pigment() {
        for (int y = 1; y < size-1; y++) {
            for (int x = 1; x < size-1; x++) {
                int idx = y * size + x;
                float u = grid[idx].velocity[0];
                float v = grid[idx].velocity[1];

                // Backtrack position
                float px = x - u * DT;
                float py = y - v * DT;

                // Bilinear interpolation of pigment
                int x0 = (int)px;
                int y0 = (int)py;
                float fx = px - x0;
                float fy = py - y0;

                for (int c = 0; c < 3; c++) {
                    float p00 = grid[y0 * size + x0].pigment[c];
                    float p10 = grid[y0 * size + std::min(x0+1, size-1)].pigment[c];
                    float p01 = grid[std::min(y0+1, size-1) * size + x0].pigment[c];
                    float p11 = grid[std::min(y0+1, size-1) * size + std::min(x0+1, size-1)].pigment[c];

                    temp_grid[idx].pigment[c] = (1-fx)*(1-fy)*p00 + fx*(1-fy)*p10 +
                                              (1-fx)*fy*p01 + fx*fy*p11;
                }
            }
        }
        std::swap(grid, temp_grid);
    }
};

// Global simulation instance
WatercolorSimulation* watercolor_sim = nullptr;

// Function declarations
void display_watercolor();
void print_help();

//for checker board:
int checker_f(int n) { return (n % 2 == 0) ? 1 : 0; }

void init(void);
void keyboard(unsigned char key, int x, int y);
void motion(int x, int y);
void display(void);
void mouse(int button, int state, int x, int y);
void display_shape(GLenum mode, Polyhedron *poly);

// Function to generate a unique color based on polygon ID
// Followed suggestions of using golden ratio from this blog post: 
// https://martin.ankerl.com/2009/12/09/how-to-create-random-colors-programmatically/
void generate_polygon_id_color(int id, float& r, float& g, float& b) {
    // Use golden ratio for good color distribution
    const float golden_ratio = 0.618033988749895f;
    const float h = fmod(id * golden_ratio, 1.0f);  // hue
    const float s = 0.7f;  // saturation
    const float v = 0.95f;  // value/brightness

    // Convert HSV to RGB
    const int hi = static_cast<int>(h * 6);
    const float f = h * 6 - hi;
    const float p = v * (1 - s);
    const float q = v * (1 - f * s);
    const float t = v * (1 - (1 - f) * s);

    switch(hi) {
        case 0: r = v; g = t; b = p; break;
        case 1: r = q; g = v; b = p; break;
        case 2: r = p; g = v; b = t; break;
        case 3: r = p; g = q; b = v; break;
        case 4: r = t; g = p; b = v; break;
        default: r = v; g = p; b = q; break;
    }
}

/******************************************************************************
Main program.
******************************************************************************/

int main(int argc, char *argv[])
{
  char *progname;
  int num = 1;
	FILE *this_file;

  progname = argv[0];

	this_file = fopen("../tempmodels/bunny.ply", "r");
	poly = new Polyhedron (this_file);
	fclose(this_file);
	mat_ident( rotmat );	

	poly->initialize(); // initialize everything

	// Compute curvatures
	poly->compute_vert_mean_curvature();
	poly->compute_vert_gaussian_curvature();

	poly->calc_bounding_sphere();
	poly->calc_face_normals_and_area();
	poly->average_normals();

	// Initialize watercolor simulation
	watercolor_sim = new WatercolorSimulation(GRID_SIZE);

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

	// Cleanup watercolor simulation
	delete watercolor_sim;
	
  return 0;    /* ANSI C requires main to return int. */
}

void color_mapping(double percentage, double col[3])
{
	if (percentage == 0.0){
		col[0] = 1.0;
		col[1] = 1.0;
		col[2] = 1.0;
	}
	else if (percentage <= 1.0/3){
		col[0] = 1.0;
		col[1] = 1.0-percentage*3.0;
		col[2] = 1.0-percentage*3.0;
	}
	else if (percentage <= 2.0/3){
		col[0] = 1.0;
		col[1] = percentage*3.0-1.0;
		col[2] = 0.0;
	}
	else if (percentage <= 3.0/3){
		col[0] = 3.0-percentage*3.0;
		col[1] = 1.0;
		col[2] = 0.0;
	}
	else {
		col[0] = 1.0;
		col[1] = 1.0;
		col[2] = 0.0;
	}
}

/******************************************************************************
Read in a polyhedron from a file.
******************************************************************************/

Polyhedron::Polyhedron(FILE *file)
{
  int i,j;
  int elem_count;
  char *elem_name;

  /*** Read in the original PLY object ***/
  in_ply = read_ply (file);

  for (i = 0; i < in_ply->num_elem_types; i++) {

    /* prepare to read the i'th list of elements */
    elem_name = setup_element_read_ply (in_ply, i, &elem_count);

    if (equal_strings ("vertex", elem_name)) {

      /* create a vertex list to hold all the vertices */
      nverts = max_verts = elem_count;
      vlist = new Vertex *[nverts];

      /* set up for getting vertex elements */

      setup_property_ply (in_ply, &vert_props[0]);
      setup_property_ply (in_ply, &vert_props[1]);
      setup_property_ply (in_ply, &vert_props[2]);
      vert_other = get_other_properties_ply (in_ply, 
					     offsetof(Vertex_io,other_props));

      /* grab all the vertex elements */
      for (j = 0; j < nverts; j++) {
        Vertex_io vert;
        get_element_ply (in_ply, (void *) &vert);

        /* copy info from the "vert" structure */
        vlist[j] = new Vertex (vert.x, vert.y, vert.z);
        vlist[j]->other_props = vert.other_props;
      }
    }
    else if (equal_strings ("face", elem_name)) {

      /* create a list to hold all the face elements */
      ntris = max_tris = elem_count;
      tlist = new Triangle *[ntris];

      /* set up for getting face elements */
      setup_property_ply (in_ply, &face_props[0]);
      face_other = get_other_properties_ply (in_ply, offsetof(Face_io,other_props));

      /* grab all the face elements */
      for (j = 0; j < elem_count; j++) {
        Face_io face;
        get_element_ply (in_ply, (void *) &face);

        if (face.nverts != 3) {
          fprintf (stderr, "Face has %d vertices (should be three).\n",
                   face.nverts);
          exit (-1);
        }

        /* copy info from the "face" structure */
        tlist[j] = new Triangle;
        tlist[j]->nverts = 3;
        tlist[j]->verts[0] = (Vertex *)(intptr_t)face.verts[0];
        tlist[j]->verts[1] = (Vertex *)(intptr_t)face.verts[1];
        tlist[j]->verts[2] = (Vertex *)(intptr_t)face.verts[2];
        tlist[j]->other_props = face.other_props;
      }
    }
    else
      get_other_element_ply (in_ply);
  }

  /* close the file */
  close_ply (in_ply);

  /* fix up vertex pointers in triangles */
  for (i = 0; i < ntris; i++) {
    tlist[i]->verts[0] = vlist[(intptr_t) tlist[i]->verts[0]];
    tlist[i]->verts[1] = vlist[(intptr_t) tlist[i]->verts[1]];
    tlist[i]->verts[2] = vlist[(intptr_t) tlist[i]->verts[2]];
  }

  /* get rid of triangles that use the same vertex more than once */

  for (i = ntris-1; i >= 0; i--) {

    Triangle *tri = tlist[i];
    Vertex *v0 = tri->verts[0];
    Vertex *v1 = tri->verts[1];
    Vertex *v2 = tri->verts[2];

    if (v0 == v1 || v1 == v2 || v2 == v0) {
      free (tlist[i]);
      ntris--;
      tlist[i] = tlist[ntris];
    }
  }
}


/******************************************************************************
Write out a polyhedron to a file.
******************************************************************************/

void Polyhedron::write_file(FILE *file)
{
  int i;
  PlyFile *ply;
  char **elist;
  int num_elem_types;

  /*** Write out the transformed PLY object ***/

  elist = get_element_list_ply (in_ply, &num_elem_types);
  ply = write_ply (file, num_elem_types, elist, in_ply->file_type);

  /* describe what properties go into the vertex elements */
  /* TODO: casting literals to char * to surpress warnings */

  describe_element_ply (ply, (char *)"vertex", nverts);
  describe_property_ply (ply, &vert_props[0]);
  describe_property_ply (ply, &vert_props[1]);
  describe_property_ply (ply, &vert_props[2]);
//  describe_other_properties_ply (ply, vert_other, offsetof(Vertex_io,other_props));

  describe_element_ply (ply, (char *)"face", ntris);
  describe_property_ply (ply, &face_props[0]);

//  describe_other_properties_ply (ply, face_other,
//                                offsetof(Face_io,other_props));

//  describe_other_elements_ply (ply, in_ply->other_elems);

  copy_comments_ply (ply, in_ply);
	char mm[1024];
	snprintf(mm, sizeof(mm), "modified by learnply");
//  append_comment_ply (ply, "modified by simvizply %f");
	  append_comment_ply (ply, mm);
  copy_obj_info_ply (ply, in_ply);

  header_complete_ply (ply);

  /* set up and write the vertex elements */
  put_element_setup_ply (ply, (char *)"vertex");
  for (i = 0; i < nverts; i++) {
    Vertex_io vert;

    /* copy info to the "vert" structure */
    vert.x = vlist[i]->x;
    vert.y = vlist[i]->y;
    vert.z = vlist[i]->z;
    vert.other_props = vlist[i]->other_props;

    put_element_ply (ply, (void *) &vert);
  }

  /* index all the vertices */
  for (i = 0; i < nverts; i++)
    vlist[i]->index = i;

  /* set up and write the face elements */
  put_element_setup_ply (ply, (char *)"face");

  Face_io face;
  face.verts = new int[3];
  
  for (i = 0; i < ntris; i++) {

    /* copy info to the "face" structure */
    face.nverts = 3;
    face.verts[0] = tlist[i]->verts[0]->index;
    face.verts[1] = tlist[i]->verts[1]->index;
    face.verts[2] = tlist[i]->verts[2]->index;
    face.other_props = tlist[i]->other_props;

    put_element_ply (ply, (void *) &face);
  }
  put_other_elements_ply (ply);

  close_ply (ply);
  free_ply (ply);
}

void Polyhedron::initialize(){
	icVector3 v1, v2;

	create_pointers();
	calc_edge_length();
	seed = -1;
}

void Polyhedron::finalize(){
	int i;

	for (i=0; i<ntris; i++){
		free(tlist[i]->other_props);
		free(tlist[i]);
	}
	for (i=0; i<nedges; i++) {
		free(elist[i]->tris);
		free(elist[i]);
	}
	for (i=0; i<nverts; i++) {
		free(vlist[i]->tris);
		free(vlist[i]->other_props);
		free(vlist[i]);
	}

	free(tlist);
	free(elist);
	free(vlist);
	if (!vert_other)
		free(vert_other);
	if (!face_other)
		free(face_other);
}

/******************************************************************************
Find out if there is another face that shares an edge with a given face.

Entry:
  f1    - face that we're looking to share with
  v1,v2 - two vertices of f1 that define edge

Exit:
  return the matching face, or NULL if there is no such face
******************************************************************************/

Triangle *Polyhedron::find_common_edge(Triangle *f1, Vertex *v1, Vertex *v2)
{
  int i,j;
  Triangle *f2;
  Triangle *adjacent = NULL;

  /* look through all faces of the first vertex */

  for (i = 0; i < v1->ntris; i++) {
    f2 = v1->tris[i];
    if (f2 == f1)
      continue;
    /* examine the vertices of the face for a match with the second vertex */
    for (j = 0; j < f2->nverts; j++) {

      /* look for a match */
      if (f2->verts[j] == v2) {

#if 0
	/* watch out for triple edges */

        if (adjacent != NULL) {

	  fprintf (stderr, "model has triple edges\n");

	  fprintf (stderr, "face 1: ");
	  for (k = 0; k < f1->nverts; k++)
	    fprintf (stderr, "%d ", f1->iverts[k]);
	  fprintf (stderr, "\nface 2: ");
	  for (k = 0; k < f2->nverts; k++)
	    fprintf (stderr, "%d ", f2->iverts[k]);
	  fprintf (stderr, "\nface 3: ");
	  for (k = 0; k < adjacent->nverts; k++)
	    fprintf (stderr, "%d ", adjacent->iverts[k]);
	  fprintf (stderr, "\n");

	}

	/* if we've got a match, remember this face */
        adjacent = f2;
#endif

#if 1
	/* if we've got a match, return this face */
        return (f2);
#endif

      }
    }
  }

  return (adjacent);
}


/******************************************************************************
Create an edge.

Entry:
  v1,v2 - two vertices of f1 that define edge
******************************************************************************/

void Polyhedron::create_edge(Vertex *v1, Vertex *v2)
{
  int i,j;
  Triangle *f;

  /* make sure there is enough room for a new edge */

  if (nedges >= max_edges) {

    max_edges += 100;
    Edge **list = new Edge *[max_edges];

    /* copy the old list to the new one */
    for (i = 0; i < nedges; i++)
      list[i] = elist[i];

    /* replace list */
    free (elist);
    elist = list;
  }

  /* create the edge */

  elist[nedges] = new Edge;
  Edge *e = elist[nedges];
  e->index = nedges;
  e->verts[0] = v1;
  e->verts[1] = v2;
  nedges++;

  /* count all triangles that will share the edge, and do this */
  /* by looking through all faces of the first vertex */

  e->ntris = 0;

  for (i = 0; i < v1->ntris; i++) {
    f = v1->tris[i];
    /* examine the vertices of the face for a match with the second vertex */
    for (j = 0; j < 3; j++) {
      /* look for a match */
      if (f->verts[j] == v2) {
        e->ntris++;
        break;
      }
    }
  }

  /* make room for the face pointers (at least two) */
  if (e->ntris < 2)
    e->tris = new Triangle *[2];
  else
    e->tris = new Triangle *[e->ntris];

  /* create pointers from edges to faces and vice-versa */

  e->ntris = 0; /* start this out at zero again for creating ptrs to tris */

  for (i = 0; i < v1->ntris; i++) {

    f = v1->tris[i];

    /* examine the vertices of the face for a match with the second vertex */
    for (j = 0; j < 3; j++)
      if (f->verts[j] == v2) {

        e->tris[e->ntris] = f;
        e->ntris++;

        if (f->verts[(j+1)%3] == v1)
          f->edges[j] = e;
        else if (f->verts[(j+2)%3] == v1)
          f->edges[(j+2)%3] = e;
        else {
          fprintf (stderr, "Non-recoverable inconsistancy in create_edge()\n");
          exit (-1);
        }

        break;  /* we'll only find one instance of v2 */
      }

  }
}


/******************************************************************************
Create edges.
******************************************************************************/

void Polyhedron::create_edges()
{
  int i,j;
  Triangle *f;
  Vertex *v1,*v2;
  double count = 0;

  /* count up how many edges we may require */

  for (i = 0; i < ntris; i++) {
    f = tlist[i];
    for (j = 0; j < f->nverts; j++) {
      v1 = f->verts[j];
      v2 = f->verts[(j+1) % f->nverts];
      Triangle *result = find_common_edge (f, v1, v2);
      if (result)
        count += 0.5;
      else
        count += 1;
    }
  }

  /*
  printf ("counted %f edges\n", count);
  */

  /* create space for edge list */

  max_edges = (int) (count + 10);  /* leave some room for expansion */
  elist = new Edge *[max_edges];
  nedges = 0;

  /* zero out all the pointers from faces to edges */

  for (i = 0; i < ntris; i++)
    for (j = 0; j < 3; j++)
      tlist[i]->edges[j] = NULL;

  /* create all the edges by examining all the triangles */

  for (i = 0; i < ntris; i++) {
    f = tlist[i];
    for (j = 0; j < 3; j++) {
      /* skip over edges that we've already created */
      if (f->edges[j])
        continue;
      v1 = f->verts[j];
      v2 = f->verts[(j+1) % f->nverts];
      create_edge (v1, v2);
    }
  }
}


/******************************************************************************
Create pointers from vertices to faces.
******************************************************************************/

void Polyhedron::vertex_to_tri_ptrs()
{
  int i,j;
  Triangle *f;
  Vertex *v;

  /* zero the count of number of pointers to faces */

  for (i = 0; i < nverts; i++)
    vlist[i]->max_tris = 0;

  /* first just count all the face pointers needed for each vertex */

  for (i = 0; i < ntris; i++) {
    f = tlist[i];
    for (j = 0; j < f->nverts; j++)
      f->verts[j]->max_tris++;
  }

  /* allocate memory for face pointers of vertices */

  for (i = 0; i < nverts; i++) {
    vlist[i]->tris = (Triangle **)
		      malloc (sizeof (Triangle *) * vlist[i]->max_tris);
    vlist[i]->ntris = 0;
  }

  /* now actually create the face pointers */

  for (i = 0; i < ntris; i++) {
    f = tlist[i];
    for (j = 0; j < f->nverts; j++) {
      v = f->verts[j];
      v->tris[v->ntris] = f;
      v->ntris++;
    }
  }
}


/******************************************************************************
Find the other triangle that is incident on an edge, or NULL if there is
no other.
******************************************************************************/

Triangle *Polyhedron::other_triangle(Edge *edge, Triangle *tri)
{
  /* search for any other triangle */

  for (int i = 0; i < edge->ntris; i++)
    if (edge->tris[i] != tri)
      return (edge->tris[i]);

  /* there is no such other triangle if we get here */
  return (NULL);
}


/******************************************************************************
Order the pointers to faces that are around a given vertex.

Entry:
  v - vertex whose face list is to be ordered
******************************************************************************/

void Polyhedron::order_vertex_to_tri_ptrs(Vertex *v)
{
  int i,j;
  Triangle *f;
  Triangle *fnext;
  int nf;
  int vindex;
  int boundary;
  int count;

  nf = v->ntris;
  f = v->tris[0];

  /* go backwards (clockwise) around faces that surround a vertex */
  /* to find out if we reach a boundary */

  boundary = 0;

  for (i = 1; i <= nf; i++) {

    /* find reference to v in f */
    vindex = -1;
    for (j = 0; j < f->nverts; j++)
      if (f->verts[j] == v) {
	vindex = j;
	break;
      }

    /* error check */
    if (vindex == -1) {
      fprintf (stderr, "can't find vertex #1\n");
      exit (-1);
    }

    /* corresponding face is the previous one around v */
    fnext = other_triangle (f->edges[vindex], f);

    /* see if we've reached a boundary, and if so then place the */
    /* current face in the first position of the vertice's face list */

    if (fnext == NULL) {
      /* find reference to f in v */
      for (j = 0; j < v->ntris; j++)
        if (v->tris[j] == f) {
	  v->tris[j] = v->tris[0];
	  v->tris[0] = f;
	  break;
	}
      boundary = 1;
      break;
    }

    f = fnext;
  }

  /* now walk around the faces in the forward direction and place */
  /* them in order */

  f = v->tris[0];
  count = 0;

  for (i = 1; i < nf; i++) {

    /* find reference to vertex in f */
    vindex = -1;
    for (j = 0; j < f->nverts; j++)
      if (f->verts[(j+1) % f->nverts] == v) {
	vindex = j;
	break;
      }

    /* error check */
    if (vindex == -1) {
      fprintf (stderr, "can't find vertex #2\n");
      exit (-1);
    }

    /* corresponding face is next one around v */
    fnext = other_triangle (f->edges[vindex], f);

    /* break out of loop if we've reached a boundary */
    count = i;
    if (fnext == NULL) {
      break;
    }

    /* swap the next face into its proper place in the face list */
    for (j = 0; j < v->ntris; j++)
      if (v->tris[j] == fnext) {
	v->tris[j] = v->tris[i];
	v->tris[i] = fnext;
	break;
      }

    f = fnext;
  }
}


/******************************************************************************
Find the index to a given vertex in the list of vertices of a given face.

Entry:
  f - face whose vertex list is to be searched
  v - vertex to return reference to

Exit:
  returns index in face's list, or -1 if vertex not found
******************************************************************************/

int Polyhedron::face_to_vertex_ref(Triangle *f, Vertex *v)
{
  int j;
  int vindex = -1;

  for (j = 0; j < f->nverts; j++)
    if (f->verts[j] == v) {
      vindex = j;
      break;
    }

  return (vindex);
}

/******************************************************************************
Create various face and vertex pointers.
******************************************************************************/

void Polyhedron::create_pointers()
{
  int i;

  /* index the vertices and triangles */

  for (i = 0; i < nverts; i++)
    vlist[i]->index = i;

  for (i = 0; i < ntris; i++) 
    tlist[i]->index = i;

  /* create pointers from vertices to triangles */
  vertex_to_tri_ptrs();

  /* make edges */
  create_edges();


  /* order the pointers from vertices to faces */
	for (i = 0; i < nverts; i++){
//		if (i %1000 == 0)
//			fprintf(stderr, "ordering %d of %d vertices\n", i, nverts);
    order_vertex_to_tri_ptrs(vlist[i]);
		
	}
  /* index the edges */

  for (i = 0; i < nedges; i++){
//		if (i %1000 == 0)
//			fprintf(stderr, "indexing %d of %d edges\n", i, nedges);
    elist[i]->index = i;
	}

}

void Polyhedron::calc_bounding_sphere()
{
  unsigned int i;
  icVector3 min, max;

  for (i=0; i<nverts; i++) {
    if (i==0)  {
			min.set(vlist[i]->x, vlist[i]->y, vlist[i]->z);
			max.set(vlist[i]->x, vlist[i]->y, vlist[i]->z);
    }
    else {
			if (vlist[i]->x < min.entry[0])
			  min.entry[0] = vlist[i]->x;
			if (vlist[i]->x > max.entry[0])
			  max.entry[0] = vlist[i]->x;
			if (vlist[i]->y < min.entry[1])
			  min.entry[1] = vlist[i]->y;
			if (vlist[i]->y > max.entry[1])
			  max.entry[1] = vlist[i]->y;
			if (vlist[i]->z < min.entry[2])
			  min.entry[2] = vlist[i]->z;
			if (vlist[i]->z > max.entry[2])
			  max.entry[2] = vlist[i]->z;
		}
  }
  center = (min + max) * 0.5;
  radius = length(center - min);
}

void Polyhedron::calc_edge_length()
{
	int i;
	icVector3 v1, v2;

	for (i=0; i<nedges; i++) {
		v1.set(elist[i]->verts[0]->x, elist[i]->verts[0]->y, elist[i]->verts[0]->z);
		v2.set(elist[i]->verts[1]->x, elist[i]->verts[1]->y, elist[i]->verts[1]->z);
		elist[i]->length = length(v1-v2);
	}
}

void Polyhedron::calc_face_normals_and_area()
{
	unsigned int i, j;
	icVector3 v0, v1, v2;
  Triangle *temp_t;
	double length[3];

	area = 0.0;
	for (i=0; i<ntris; i++){
		for (j=0; j<3; j++)
			length[j] = tlist[i]->edges[j]->length;
		double temp_s = (length[0] + length[1] + length[2])/2.0;
		tlist[i]->area = sqrt(temp_s*(temp_s-length[0])*(temp_s-length[1])*(temp_s-length[2]));

		area += tlist[i]->area;
		temp_t = tlist[i];
		v1.set(vlist[tlist[i]->verts[0]->index]->x, vlist[tlist[i]->verts[0]->index]->y, vlist[tlist[i]->verts[0]->index]->z);
		v2.set(vlist[tlist[i]->verts[1]->index]->x, vlist[tlist[i]->verts[1]->index]->y, vlist[tlist[i]->verts[1]->index]->z);
		v0.set(vlist[tlist[i]->verts[2]->index]->x, vlist[tlist[i]->verts[2]->index]->y, vlist[tlist[i]->verts[2]->index]->z);
		tlist[i]->normal = cross(v0-v1, v2-v1);
		normalize(tlist[i]->normal);
	}

	double signedvolume = 0.0;
	icVector3 test = center;
	for (i=0; i<ntris; i++){
		icVector3 cent(vlist[tlist[i]->verts[0]->index]->x, vlist[tlist[i]->verts[0]->index]->y, vlist[tlist[i]->verts[0]->index]->z);
		signedvolume += dot(test-cent, tlist[i]->normal)*tlist[i]->area;
	}
	signedvolume /= area;
	if (signedvolume<0) 
		orientation = 0;
	else {
		orientation = 1;
		for (i=0; i<ntris; i++)
			tlist[i]->normal *= -1.0;
	}
}

void sort(unsigned int *A, unsigned int *B, unsigned int *C, unsigned int sid, unsigned int eid){
  unsigned int i;
	unsigned int *tempA, *tempB, *tempC;
	unsigned int current1, current2, current0;

  if (sid>=eid)
		return;
	sort(A, B, C, sid, (sid+eid)/2);
	sort(A, B, C, (sid+eid)/2+1, eid);
	tempA = (unsigned int *)malloc(sizeof(unsigned int)*(eid-sid+1));
	tempB = (unsigned int *)malloc(sizeof(unsigned int)*(eid-sid+1));
	tempC = (unsigned int *)malloc(sizeof(unsigned int)*(eid-sid+1));
	for (i=0; i<eid-sid+1; i++){
		tempA[i] = A[i+sid];
		tempB[i] = B[i+sid];
		tempC[i] = C[i+sid];
	}
	current1 = sid;
	current2 = (sid+eid)/2+1;
	current0 = sid;
	while ((current1<=(sid+eid)/2) && (current2<=eid)){
		if (tempA[current1-sid] < tempA[current2-sid]) {
			A[current0] = tempA[current1-sid];
			B[current0] = tempB[current1-sid];
			C[current0] = tempC[current1-sid];
			current1++;		
		}
		else if (tempA[current1-sid] > tempA[current2-sid]){
			A[current0] = tempA[current2-sid];
			B[current0] = tempB[current2-sid];
			C[current0] = tempC[current2-sid];
			current2++;		
		}
		else {
			if (tempB[current1-sid] < tempB[current2-sid]) {
				A[current0] = tempA[current1-sid];
				B[current0] = tempB[current1-sid];
				C[current0] = tempC[current1-sid];
				current1++;		
			} else {
				A[current0] = tempA[current2-sid];
				B[current0] = tempB[current2-sid];
				C[current0] = tempC[current2-sid];
				current2++;		
			}
		}
		current0++;
	}
	if (current1<=(sid+eid)/2){
		for (i=current1; i<=(sid+eid)/2; i++){
			A[current0] = tempA[i-sid];
			B[current0] = tempB[i-sid];
			C[current0] = tempC[i-sid];
			current0++;
		}
	}
	if (current2<=eid){
		for (i=current2; i<=eid; i++){
			A[current0] = tempA[i-sid];
			B[current0] = tempB[i-sid];
			C[current0] = tempC[i-sid];
			current0++;
		}
	}

	free(tempA);
	free(tempB);
	free(tempC);
}

void init(void) {
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


/******************************************************************************
Process a keyboard action.  In particular, exit the program when an
"escape" is pressed in the window.
******************************************************************************/

void keyboard(unsigned char key, int x, int y) {
	int i;

  /* set escape key to exit */
  switch (key) {
    case 27:
			poly->finalize();  // finalize_everything
      exit(0);
      break;

		case '0':
			display_mode = 0;
			display();
			break;

		case '1':
			display_mode = 1; // Polygon ID coloring mode
			display();
			break;

		case '2':
			display_mode = 2; // Barycentric coordinates map
			display();
			break;

		case '3':
			display_mode = 3; // Polygon normal coloring
			display();
			break;

		case '4':
			display_mode = 4; // Vertex normal coloring
			display();
			break;

		case '5':
			display_mode = 5; // 3D checkerboard coloring
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

		case 'x':
			switch(ACSIZE){
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
			fprintf(stderr, "ACSIZE=%d\n", ACSIZE);
			display();
			break;

		case '|':
			this_file = fopen("rotmat.txt", "w");
			for (i=0; i<4; i++) 
				fprintf(this_file, "%f %f %f %f\n", rotmat[i][0], rotmat[i][1], rotmat[i][2], rotmat[i][3]);
			fclose(this_file);
			break;

		case '^':
			this_file = fopen("rotmat.txt", "r");
			for (i=0; i<4; i++) 
				fscanf(this_file, "%f %f %f %f ", (&rotmat[i][0]), (&rotmat[i][1]), (&rotmat[i][2]), (&rotmat[i][3]));
			fclose(this_file);
			display();
			break;

		case '+':
			checker_L *= 1.1f; // Increase checker size
			printf("Checker L: %f\n", checker_L);
			display();
			break;

		case '-':
			checker_L /= 1.1f; // Decrease checker size
			printf("Checker L: %f\n", checker_L);
			display();
			break;

		case 'w':  // Toggle watercolor mode
			wcolor_mode = !wcolor_mode;
			printf("Watercolor mode: %s\n", watercolor_mode ? "ON" : "OFF");
			glutPostRedisplay();
			return;

		case 'c':
			color_intensity += 0.1f;
			printf("Color intensity: %f\n", color_intensity);
			display();
			break;

		case 'C':
			color_intensity -= 0.1f;
			printf("Color intensity: %f\n", color_intensity);
			display();
			break;

		case 'e':
			edge_darkening += 0.1f;
			printf("Edge darkening: %f\n", edge_darkening);
			display();
			break;

		case 'E':
			edge_darkening -= 0.1f;
			printf("Edge darkening: %f\n", edge_darkening);
			display();
			break;

		case 'n':
			noise_scale += 0.01f;
			printf("Noise scale: %f\n", noise_scale);
			display();
			break;

		case 'N':
			noise_scale -= 0.01f;
			printf("Noise scale: %f\n", noise_scale);
			display();
			break;

		case 'g':  // Increase granulation
			granulation_scale += 0.1f;
			printf("Granulation scale: %f\n", granulation_scale);
			display();
			break;

		case 'G':  // Decrease granulation
			granulation_scale = std::max(0.0f, granulation_scale - 0.1f);
			printf("Granulation scale: %f\n", granulation_scale);
			display();
			break;

    case 'v':  // Increase viscosity
      VISCOSITY += 0.01f;
      printf("Viscosity: %f\n", VISCOSITY);
      display();
      break;

    case 'V':  // Decrease viscosity
      VISCOSITY = std::max(0.01f, VISCOSITY - 0.01f);
      printf("Viscosity: %f\n", VISCOSITY);
      display();
      break;

    case 'a':  // Increase evaporation
      EVAPORATION_RATE += 0.01f;
      printf("Evaporation rate: %f\n", EVAPORATION_RATE);
      display();
      break;

    case 'A':  // Decrease evaporation
      EVAPORATION_RATE = std::max(0.01f, EVAPORATION_RATE - 0.01f);
      printf("Evaporation rate: %f\n", EVAPORATION_RATE);
      display();
      break;

    case 'r':  // Reset watercolor parameters
      color_intensity = 1.0f;
      granulation_scale = 0.5f;
      VISCOSITY = 0.05f;
      EVAPORATION_RATE = 0.03f;
      printf("Reset watercolor parameters\n");
      display();
      break;
	}
}

Polyhedron::Polyhedron()
{
	nverts = nedges = ntris = 0;
	max_verts = max_tris = 50;

	vlist = new Vertex *[max_verts];
	tlist = new Triangle *[max_tris];		
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

void display_object()
{
	unsigned int i, j;
	Polyhedron *the_patch = poly;
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glShadeModel(GL_SMOOTH);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHT1);
	for (i=0; i<poly->ntris; i++) {
		Triangle *temp_t=poly->tlist[i];
		glBegin(GL_POLYGON);
		GLfloat mat_diffuse[] = {1.0, 1.0, 1.0, 1.0};
		
		glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
   
		glColor3f(1.0, 1.0, 1.0);
		glNormal3d(temp_t->normal.entry[0], temp_t->normal.entry[1], temp_t->normal.entry[2]);
		for (j=0; j<3; j++) {
			Vertex *temp_v = temp_t->verts[j];
			glVertex3d(temp_v->x, temp_v->y, temp_v->z);
		}
		glEnd();
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

	for (i=0; i<this_poly->ntris; i++) {
		if (mode == GL_SELECT)
      glLoadName(i+1);

		Triangle *temp_t=this_poly->tlist[i];

		switch (display_mode) {
		case 0:
			if (i == this_poly->seed) {
				mat_diffuse[0] = 0.0;
				mat_diffuse[1] = 0.0;
				mat_diffuse[2] = 1.0;
				mat_diffuse[3] = 1.0;
			} else {
				mat_diffuse[0] = 1.0;
				mat_diffuse[1] = 1.0;
				mat_diffuse[2] = 0.0;
				mat_diffuse[3] = 1.0;
			}
			glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
			glBegin(GL_POLYGON);
			for (j=0; j<3; j++) {

				Vertex *temp_v = temp_t->verts[j];
				glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
				if (i==this_poly->seed)
					glColor3f(0.0, 0.0, 1.0);
				else
					glColor3f(1.0, 1.0, 0.0);
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();
			break;

		case 1:  // Polygon ID coloring mode
			glBegin(GL_POLYGON);
			{
				generate_polygon_id_color(i, r, g, b);
				mat_diffuse[0] = r;
				mat_diffuse[1] = g;
				mat_diffuse[2] = b;
				mat_diffuse[3] = 1.0;
				glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
				
				for (j = 0; j < 3; j++) {
					Vertex *temp_v = temp_t->verts[j];
					glNormal3d(temp_t->normal.entry[0], temp_t->normal.entry[1], temp_t->normal.entry[2]);
					glColor3f(r, g, b);
					glVertex3d(temp_v->x, temp_v->y, temp_v->z);
				}
			}
			glEnd();
			break;

		case 2: // Barycentric coordinates coloring
			glDisable(GL_LIGHTING); // disabling lighting because OpenGL uses smooth shading and we do not want a color blend
			glBegin(GL_POLYGON);
			// First vertex: Red
			glColor3f(1.0, 0.0, 0.0);
			glNormal3d(temp_t->normal.entry[0], temp_t->normal.entry[1], temp_t->normal.entry[2]);
			glVertex3d(temp_t->verts[0]->x, temp_t->verts[0]->y, temp_t->verts[0]->z);

			// Second vertex: Green
			glColor3f(0.0, 1.0, 0.0);
			glNormal3d(temp_t->normal.entry[0], temp_t->normal.entry[1], temp_t->normal.entry[2]);
			glVertex3d(temp_t->verts[1]->x, temp_t->verts[1]->y, temp_t->verts[1]->z);

			// Third vertex: Blue
			glColor3f(0.0, 0.0, 1.0);
			glNormal3d(temp_t->normal.entry[0], temp_t->normal.entry[1], temp_t->normal.entry[2]);
			glVertex3d(temp_t->verts[2]->x, temp_t->verts[2]->y, temp_t->verts[2]->z);
			glEnd();
			glEnable(GL_LIGHTING); //reenable lighting after draw
			break;

		case 3: // Polygon normal coloring
			glDisable(GL_LIGHTING);
			glBegin(GL_POLYGON);
			// Map normal from [-1,1] to [0,1] for color
			r = 0.5f * (temp_t->normal.entry[0] + 1.0f);
			g = 0.5f * (temp_t->normal.entry[1] + 1.0f);
			b = 0.5f * (temp_t->normal.entry[2] + 1.0f);
			glColor3f(r, g, b);
			for (j = 0; j < 3; j++) {
				Vertex *temp_v = temp_t->verts[j];
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();
			glEnable(GL_LIGHTING);
			break;

		case 4: // Vertex normal coloring
			glDisable(GL_LIGHTING);
			glBegin(GL_POLYGON);
			for (j = 0; j < 3; j++) {
				Vertex *temp_v = temp_t->verts[j];
				r = 0.5f * (temp_v->normal.entry[0] + 1.0f);
				g = 0.5f * (temp_v->normal.entry[1] + 1.0f);
				b = 0.5f * (temp_v->normal.entry[2] + 1.0f);
				glColor3f(r, g, b);
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();
			glEnable(GL_LIGHTING);
			break;

		case 5: // 3D Checkerboard coloring
			glDisable(GL_LIGHTING);
			glBegin(GL_POLYGON);
			for (j = 0; j < 3; j++) {
				Vertex *temp_v = temp_t->verts[j];
				int nx = static_cast<int>(floor(temp_v->x / checker_L));
				int ny = static_cast<int>(floor(temp_v->y / checker_L));
				int nz = static_cast<int>(floor(temp_v->z / checker_L));
				r = checker_f(nx);
				g = checker_f(ny);
				b = checker_f(nz);
				glColor3f(r, g, b);
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();
			glEnable(GL_LIGHTING);
			break;

		case 6:
			glBegin(GL_POLYGON);
			for (j=0; j<3; j++) {
				Vertex *temp_v = temp_t->verts[j];
				glNormal3d(temp_t->normal.entry[0], temp_t->normal.entry[1], temp_t->normal.entry[2]);
				glColor3f(1.0, 1.0, 1.0);
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();
			break;

		case 10:
			glBegin(GL_POLYGON);
			for (j=0; j<3; j++) {
				mat_diffuse[0] = 1.0;
				mat_diffuse[1] = 0.0;
				mat_diffuse[2] = 0.0;
				mat_diffuse[3] = 1.0;
		
				glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);

				Vertex *temp_v = temp_t->verts[j];
				glNormal3d(temp_t->normal.entry[0], temp_t->normal.entry[1], temp_t->normal.entry[2]);

				glColor3f(1.0, 0.0, 0.0);
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();
			break;
		}
		
	}
}

void display()
{
    if (wcolor_mode) {
        display_watercolor();
        return;
    }
    
    GLint viewport[4];
    int jitter;

    glClearColor (1.0, 1.0, 1.0, 1.0);  // background for rendering color coding and lighting
    glGetIntegerv (GL_VIEWPORT, viewport);
   
    glClear(GL_ACCUM_BUFFER_BIT);
    for (jitter = 0; jitter < ACSIZE; jitter++) {
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        set_view(GL_RENDER, poly);
        glPushMatrix ();
        switch(ACSIZE){
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
        display_shape(GL_RENDER, poly);
        glPopMatrix ();
        glAccum(GL_ACCUM, 1.0/ACSIZE);
    }
    glAccum (GL_RETURN, 1.0);
    glFlush();
    glutSwapBuffers();
    glFinish();
}

void display_watercolor()
{
    // Initialize simulation if needed
    if (!watercolor_sim) {
        watercolor_sim = new WatercolorSimulation(GRID_SIZE);
        initialize_pigments_from_geometry();
    }

    // Update simulation parameters
    watercolor_sim->set_parameters(VISCOSITY, EVAPORATION_RATE);

    // Simulate fluid flow
    watercolor_sim->simulate();

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();
    
    // Set up view to match regular display
    set_view(GL_RENDER, poly);
    glPushMatrix();
    set_scene(GL_RENDER, poly);
    
    // Enable necessary OpenGL features
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glDisable(GL_LIGHTING);
    
    // First pass: Base colors
    glBegin(GL_TRIANGLES);
    for (int i = 0; i < poly->ntris; i++) {
        Triangle* t = poly->tlist[i];
        for (int j = 0; j < 3; j++) {
            Vertex* v = t->verts[j];
            
            // Project vertex to 2D grid coordinates
            float x = (v->x + 1.0f) * 0.5f;
            float y = (v->y + 1.0f) * 0.5f;
            
            // Get pigment from simulation
            int grid_x = std::min(GRID_SIZE-1, std::max(0, (int)(x * GRID_SIZE)));
            int grid_y = std::min(GRID_SIZE-1, std::max(0, (int)(y * GRID_SIZE)));
            int idx = grid_y * GRID_SIZE + grid_x;
            
            // Get raw pigment values
            float r = watercolor_sim->grid[idx].pigment[0];
            float g = watercolor_sim->grid[idx].pigment[1];
            float b = watercolor_sim->grid[idx].pigment[2];
            
            // If no pigment, use a default color
            if (r < 0.01f && g < 0.01f && b < 0.01f) {
                r = 0.7f; g = 0.8f; b = 0.9f;  // Light blue default
            }
            
            // Apply noise effect
            float noise_val = octave_noise(x * noise_scale, y * noise_scale, 0.0f, 4, 0.5f);
            float noise_factor = 1.0f + (noise_val * 0.2f * noise_scale);
            
            // Apply granulation effect
            float gran = octave_noise(x * granulation_scale * 0.25f, 
                                    y * granulation_scale * 0.25f, 
                                    0.0f, 4, 0.5f);
            float gran_factor = 1.0f + (gran * 0.05f * granulation_scale);
            
            // Combine noise and granulation
            float final_factor = noise_factor * gran_factor;
            final_factor = std::min(1.2f, std::max(0.8f, final_factor));
            
            r = std::min(1.0f, r * final_factor);
            g = std::min(1.0f, g * final_factor);
            b = std::min(1.0f, b * final_factor);
            
            // Apply color intensity
            r *= color_intensity;
            g *= color_intensity;
            b *= color_intensity;
            
            glColor4f(r, g, b, 0.8f);
            glVertex3f(v->x, v->y, v->z);
        }
    }
    glEnd();
    
    // Second pass: Add subtle highlights
    glEnable(GL_POINT_SMOOTH);
    glPointSize(1.0);
    glBegin(GL_POINTS);
    for (int i = 0; i < poly->nverts; i++) {
        Vertex* v = poly->vlist[i];
        if (rand() % 100 < 20) {  // 20% chance for highlight
            float x = (v->x + 1.0f) * 0.5f;
            float y = (v->y + 1.0f) * 0.5f;
            float highlight = noise(x * 8.0f, y * 8.0f, 0.0f) * 0.2f;
            glColor4f(1.0f, 1.0f, 1.0f, highlight);
            glVertex3f(v->x, v->y, v->z);
        }
    }
    glEnd();
    
    glPopMatrix();
    glutSwapBuffers();
}

void print_help()
{
    printf("\nWatercolor Controls:\n");
    printf("i/k: Increase/decrease color intensity\n");
    printf("g/h: Increase/decrease granulation\n");
    printf("v/b: Increase/decrease viscosity\n");
    printf("e/d: Increase/decrease evaporation rate\n");
    printf("r: Reset watercolor parameters\n");
    printf("w: Toggle watercolor mode\n");
    printf("\n");
}

void Polyhedron::compute_vert_mean_curvature() {
    min_mean_curvature = 1e10;
    max_mean_curvature = -1e10;
    
    for (int i = 0; i < nverts; i++) {
        Vertex* v = vlist[i];
        icVector3 mean_curvature_vector(0.0);
        
        // Compute mean curvature using the cotangent formula
        for (int j = 0; j < v->ntris; j++) {
            Triangle* t = v->tris[j];
            int v_index = face_to_vertex_ref(t, v);
            Vertex* v1 = t->verts[(v_index + 1) % 3];
            Vertex* v2 = t->verts[(v_index + 2) % 3];
            
            icVector3 e1(v1->x - v->x, v1->y - v->y, v1->z - v->z);
            icVector3 e2(v2->x - v->x, v2->y - v->y, v2->z - v->z);
            
            double cot_alpha = dot(e1, e2) / length(cross(e1, e2));
            mean_curvature_vector += t->normal * cot_alpha;
        }
        
        v->mean_curvature = length(mean_curvature_vector) / 2.0;
        min_mean_curvature = std::min(min_mean_curvature, v->mean_curvature);
        max_mean_curvature = std::max(max_mean_curvature, v->mean_curvature);
    }
}

void Polyhedron::compute_vert_gaussian_curvature() {
    min_gauss_curvature = 1e10;
    max_gauss_curvature = -1e10;
    
    for (int i = 0; i < nverts; i++) {
        Vertex* v = vlist[i];
        double angle_sum = 0.0;
        
        // Compute Gaussian curvature using angle deficit
        for (int j = 0; j < v->ntris; j++) {
            Triangle* t = v->tris[j];
            int v_index = face_to_vertex_ref(t, v);
            Vertex* v1 = t->verts[(v_index + 1) % 3];
            Vertex* v2 = t->verts[(v_index + 2) % 3];
            
            icVector3 e1(v1->x - v->x, v1->y - v->y, v1->z - v->z);
            icVector3 e2(v2->x - v->x, v2->y - v->y, v2->z - v->z);
            
            normalize(e1);
            normalize(e2);
            double angle = acos(dot(e1, e2));
            angle_sum += angle;
        }
        
        v->gauss_curvature = 2 * PI - angle_sum;
        min_gauss_curvature = std::min(min_gauss_curvature, v->gauss_curvature);
        max_gauss_curvature = std::max(max_gauss_curvature, v->gauss_curvature);
    }
}

void Polyhedron::compute_silhouette_edges() {
    silhouette.clear();
    
    // Get view direction (assuming camera is at origin looking down -z)
    icVector3 view_dir(0, 0, -1);
    
    for (int i = 0; i < nedges; i++) {
        Edge* e = elist[i];
        if (e->ntris == 2) {  // Only consider edges shared by two triangles
            Triangle* t1 = e->tris[0];
            Triangle* t2 = e->tris[1];
            
            // Check if edge is a silhouette
            double dot1 = dot(t1->normal, view_dir);
            double dot2 = dot(t2->normal, view_dir);
            
            if (dot1 * dot2 < 0) {  // One face is front-facing, other is back-facing
                LineSegment segment;
                segment.start.set(e->verts[0]->x, e->verts[0]->y, e->verts[0]->z);
                segment.end.set(e->verts[1]->x, e->verts[1]->y, e->verts[1]->z);
                silhouette.push_back(segment);
            }
        }
    }
}

// Pigment initialization based on 3D geometry properties
void initialize_pigments_from_geometry() {
    printf("Curvature range: %f to %f\n", poly->min_mean_curvature, poly->max_mean_curvature);
    
    // First pass: Initialize base colors using vertex attributes
    for (int i = 0; i < poly->nverts; i++) {
        Vertex* v = poly->vlist[i];
        
        // Project vertex to grid
        float x = (v->x + 1.0f) * 0.5f * GRID_SIZE;
        float y = (v->y + 1.0f) * 0.5f * GRID_SIZE;
        
        // Map vertex attributes to colors
        float t = (float)i / (float)poly->nverts;  // Vertex index for base variation
        
        // Use curvature for color intensity
        float curvature_t = (v->mean_curvature - poly->min_mean_curvature) / 
                          (poly->max_mean_curvature - poly->min_mean_curvature);
        
        // Create watercolor palette with multiple color channels
        // Reds: Based on vertex index and curvature
        float r = 0.3 + 0.5 * sin(t * 6.28f) + 0.2 * curvature_t;
        
        // Greens: Based on position and curvature
        float pos_t = (v->x + v->y + v->z) / 3.0f;  // Position-based variation
        float g = 0.2 + 0.4 * sin(pos_t * 3.14f) + 0.3 * (1.0f - curvature_t);
        
        // Blues: Based on vertex index and position
        float b = 0.4 + 0.4 * cos(t * 4.0f) + 0.2 * sin(pos_t * 2.0f);
        
        // Add some randomness for natural variation
        r += 0.1f * (rand() / (float)RAND_MAX - 0.5f);
        g += 0.1f * (rand() / (float)RAND_MAX - 0.5f);
        b += 0.1f * (rand() / (float)RAND_MAX - 0.5f);
        
        // Ensure minimum color intensity
        r = std::max(0.2f, std::min(1.0f, r));
        g = std::max(0.2f, std::min(1.0f, g));
        b = std::max(0.2f, std::min(1.0f, b));
        
        // Vary pigment amount based on curvature
        float amount = 0.4f + 0.4f * curvature_t;
        amount = std::max(0.3f, std::min(0.8f, amount));
        
        watercolor_sim->add_pigment((int)x, (int)y, r, g, b, amount);
        
        printf("Vertex %d: curvature=%.3f, colors=(%.2f,%.2f,%.2f), amount=%.2f\n", 
               i, v->mean_curvature, r, g, b, amount);
    }
    
    // Second pass: Add watercolor effects
    for (int i = 0; i < poly->nverts; i++) {
        if (rand() % 100 < 30) {  // 30% chance for each vertex
            Vertex* v = poly->vlist[i];
            float x = (v->x + 1.0f) * 0.5f * GRID_SIZE;
            float y = (v->y + 1.0f) * 0.5f * GRID_SIZE;
            
            // Create splash effect with colors based on vertex attributes
            float t = (float)i / (float)poly->nverts;
            float curvature_t = (v->mean_curvature - poly->min_mean_curvature) / 
                              (poly->max_mean_curvature - poly->min_mean_curvature);
            
            // Splash colors influenced by vertex attributes
            float splash_r = 0.4 + 0.3 * sin(t * 6.28f) + 0.2 * curvature_t;
            float splash_g = 0.2 + 0.3 * cos(t * 4.0f) + 0.2 * (1.0f - curvature_t);
            float splash_b = 0.3 + 0.3 * sin(t * 3.0f) + 0.2 * curvature_t;
            
            // Add splash with random radius
            float radius = 2.0f + 3.0f * (rand() / (float)RAND_MAX);
            for (float dx = -radius; dx <= radius; dx++) {
                for (float dy = -radius; dy <= radius; dy++) {
                    if (dx*dx + dy*dy <= radius*radius) {
                        int splash_x = (int)(x + dx);
                        int splash_y = (int)(y + dy);
                        if (splash_x >= 0 && splash_x < GRID_SIZE && 
                            splash_y >= 0 && splash_y < GRID_SIZE) {
                            float falloff = 1.0f - (dx*dx + dy*dy) / (radius*radius);
                            watercolor_sim->add_pigment(splash_x, splash_y, 
                                                      splash_r, splash_g, splash_b, 
                                                      0.3f * falloff);
                        }
                    }
                }
            }
        }
    }
}

// Forward declarations
float fade(float t);
float lerp(float t, float a, float b);
float grad(int hash, float x, float y, float z);
float perlin_noise(float x, float y, float z);
float octave_noise(float x, float y, float z, int octaves, float persistence);
float noise(float x, float y, float z);

// Helper functions for noise
float fade(float t) { return t * t * t * (t * (t * 6.0f - 15.0f) + 10.0f); }
float lerp(float t, float a, float b) { return a + t * (b - a); }
float grad(int hash, float x, float y, float z) {
    int h = hash & 15;
    float u = h < 8 ? x : y;
    float v = h < 4 ? y : h == 12 || h == 14 ? x : z;
    return ((h & 1) == 0 ? u : -u) + ((h & 2) == 0 ? v : -v);
}

float perlin_noise(float x, float y, float z) {
    // Permutation table
    static const int p[512] = {
        151,160,137,91,90,15,131,13,201,95,96,53,194,233,7,225,140,36,103,30,69,142,8,99,37,240,21,10,23,
        190,6,148,247,120,234,75,0,26,197,62,94,252,219,203,117,35,11,32,57,177,33,88,237,149,56,87,174,20,125,136,171,168,
        68,175,74,165,71,134,139,48,27,166,77,146,158,231,83,111,229,122,60,211,133,230,220,105,92,41,55,46,245,40,244,
        102,143,54,65,25,63,161,1,216,80,73,209,76,132,187,208,89,18,169,200,196,135,130,116,188,159,86,164,100,109,198,173,186,
        3,64,52,217,226,250,124,123,5,202,38,147,118,126,255,82,85,212,207,206,59,227,47,16,58,17,182,189,28,42,223,183,170,213,
        119,248,152,2,44,154,163,70,221,153,101,155,167,43,172,9,129,22,39,253,19,98,108,110,79,113,224,232,178,185,112,104,218,246,97,228,
        251,34,242,193,238,210,144,12,191,179,162,241,81,51,145,235,249,14,239,107,49,192,214,31,181,199,106,157,184,84,204,176,115,121,50,45,127,4,150,254,
        138,236,205,93,222,114,67,29,24,72,243,141,128,195,78,66,215,61,156,180,
        // Duplicate the permutation table
        151,160,137,91,90,15,131,13,201,95,96,53,194,233,7,225,140,36,103,30,69,142,8,99,37,240,21,10,23,
        190,6,148,247,120,234,75,0,26,197,62,94,252,219,203,117,35,11,32,57,177,33,88,237,149,56,87,174,20,125,136,171,168,
        68,175,74,165,71,134,139,48,27,166,77,146,158,231,83,111,229,122,60,211,133,230,220,105,92,41,55,46,245,40,244,
        102,143,54,65,25,63,161,1,216,80,73,209,76,132,187,208,89,18,169,200,196,135,130,116,188,159,86,164,100,109,198,173,186,
        3,64,52,217,226,250,124,123,5,202,38,147,118,126,255,82,85,212,207,206,59,227,47,16,58,17,182,189,28,42,223,183,170,213,
        119,248,152,2,44,154,163,70,221,153,101,155,167,43,172,9,129,22,39,253,19,98,108,110,79,113,224,232,178,185,112,104,218,246,97,228,
        251,34,242,193,238,210,144,12,191,179,162,241,81,51,145,235,249,14,239,107,49,192,214,31,181,199,106,157,184,84,204,176,115,121,50,45,127,4,150,254,
        138,236,205,93,222,114,67,29,24,72,243,141,128,195,78,66,215,61,156,180
    };

    // Find unit cube that contains the point
    int xi = static_cast<int>(floor(x)) & 255;
    int yi = static_cast<int>(floor(y)) & 255;
    int zi = static_cast<int>(floor(z)) & 255;

    // Find relative position in cube
    x -= floor(x);
    y -= floor(y);
    z -= floor(z);

    // Compute fade curves
    float u = fade(x);
    float v = fade(y);
    float w = fade(z);

    // Hash coordinates of the 8 cube corners
    int A = p[xi]+yi;
    int AA = p[A]+zi;
    int AB = p[A+1]+zi;
    int B = p[xi+1]+yi;
    int BA = p[B]+zi;
    int BB = p[B+1]+zi;

    // Add blended results from 8 corners of cube
    return lerp(w, lerp(v, lerp(u, grad(p[AA], x, y, z),
                                   grad(p[BA], x-1, y, z)),
                           lerp(u, grad(p[AB], x, y-1, z),
                                   grad(p[BB], x-1, y-1, z))),
                   lerp(v, lerp(u, grad(p[AA+1], x, y, z-1),
                                   grad(p[BA+1], x-1, y, z-1)),
                           lerp(u, grad(p[AB+1], x, y-1, z-1),
                                   grad(p[BB+1], x-1, y-1, z-1))));
}

float octave_noise(float x, float y, float z, int octaves, float persistence) {
    float total = 0;
    float frequency = 1;
    float amplitude = 1;
    float maxValue = 0;
    
    for(int i = 0; i < octaves; i++) {
        total += perlin_noise(x * frequency, y * frequency, z * frequency) * amplitude;
        maxValue += amplitude;
        amplitude *= persistence;
        frequency *= 2;
    }
    
    return total / maxValue;
}

float noise(float x, float y, float z) {
    return perlin_noise(x, y, z);
}

void Polyhedron::average_normals()
{
    int i, j;

    for (i=0; i<nverts; i++) {
        vlist[i]->normal = icVector3(0.0);
        for (j=0; j<vlist[i]->ntris; j++) 
            vlist[i]->normal += vlist[i]->tris[j]->normal;
        normalize(vlist[i]->normal);
    }
}

