/*

Data structures for learnply

Eugene Zhang 2005

*/

#ifndef __LEARNPLY_H__
#define __LEARNPLY_H__


#include "ply.h"
#include "icVector.H"
#include <vector>

const double EPS = 1.0e-6;
const double PI=3.1415926535898;

/* forward declarations */
class Triangle;
class Corner;

class Vertex {
public:
  double x,y,z;
  int index;

  int ntris;
  Triangle **tris;
  int max_tris;

	icVector3 normal;
  void *other_props;

  //added for corners 
  int ncorners;
  Corner **corners;

public:
  Vertex(double xx, double yy, double zz) 
  {
	x = xx; y = yy; z = zz; 
	ncorners = 0;
	corners = NULL;

	}
};

class Edge {
public:
  int index;
  Vertex *verts[2];
  Corner *corners[2];
  int ntris;
  Triangle **tris;
	double length;
};

class Triangle {
public:
  int index;
  int nverts;
  Vertex *verts[3];
  Edge *edges[3];

  Corner *corners[3];

	double angle[3];
	float area;

	icVector3 normal;
  void *other_props;

};

class Corner {
public:
    int index; 

	Edge *e;
	Vertex *v;
	Triangle *t;

	Corner *o, *n, *p;

	Corner()
	{
		index = 0;
		e = NULL;
		v = NULL;
		t = NULL;
		o = NULL;
		n = NULL;
		p = NULL;
	}
};

struct CornerTableEntry
{
	int c;
	int v_min;
	int v_max;
	int o;

};

class VertexList {

public:

  int num_verts;
  int max_verts;
  Vertex **verts;

  VertexList (int max) {
    max_verts = max;
    verts = new Vertex *[max_verts];
    num_verts = 0;
  }
	void finalize() {
		free(verts);
		free(this);
	}
	void append(Vertex *v)
	{
		int i;

		/* first make sure there is enough room for new vertex */

		if (num_verts >= max_verts) {
			max_verts += 10;
			Vertex **tlist = new Vertex *[max_verts];
			for (i = 0; i < num_verts; i++)
				tlist[i] = verts[i];
			delete (verts);
			verts = tlist;
		}

		/* add new vertex to list */

		verts[num_verts] = v;
		num_verts++;
	}

};

class TriangleList {
public:

  int num_tris;
  int max_tris;
  Triangle **tris;

  TriangleList (int max) {
    max_tris = max;
    tris = new Triangle *[max_tris];
    num_tris = 0;
  }
	void append(Triangle *t)
	{
		int i;

		/* first make sure there is enough room for new triangle */

		if (num_tris >= max_tris) {
			max_tris += 10;
			Triangle **tlist = new Triangle *[max_tris];
			for (i = 0; i < num_tris; i++)
				tlist[i] = tris[i];
			delete (tris);
			tris = tlist;
		}

		/* add new triangle to list */

		tris[num_tris] = t;
		num_tris++;
	}
};

class EdgeList {

public:

  int num_edges;
  int max_edges;
  Edge **edges;

  EdgeList (int max) {
    max_edges = max;
    edges = new Edge *[max_edges];
    num_edges = 0;
  }
	void append(Edge *e)
	{
		int i;

		/* first make sure there is enough room for new edge */

		if (num_edges >= max_edges) {
			max_edges += 10;
			Edge **tlist = new Edge *[max_edges];
			for (i = 0; i < num_edges; i++)
				tlist[i] = edges[i];
			delete (edges);
			edges = tlist;
		}

		/* add new edge to list */

		edges[num_edges] = e;
		num_edges++;
	}
};

class Polyhedron {
public:

	int index;

  Triangle **tlist;  /* list of triangles */
  int ntris;
  int max_tris;

  Vertex **vlist;    /* list of vertices */
  int nverts;
  int max_verts;

  Edge **elist;      /* list of edges */
  int nedges;
  int max_edges;

  Corner **clist;
  int ncorners;
  int max_corners;
  CornerTableEntry *ctable;

	icVector3 center;
	double radius;
	double area;

	int seed;

  PlyOtherProp *vert_other,*face_other;

  std::vector<Corner> corners;

	void average_normals();

  void create_edge(Vertex *, Vertex *);
  void create_edges();
  int face_to_vertex_ref(Triangle *, Vertex *);
  void order_vertex_to_tri_ptrs(Vertex *);
  void vertex_to_tri_ptrs();
  Triangle *find_common_edge(Triangle *, Vertex *, Vertex *);
  Triangle *other_triangle(Edge *, Triangle *);
	void calc_bounding_sphere();
	void calc_face_normals_and_area();
	void calc_edge_length();

	void create_corners();


	Polyhedron();
  Polyhedron(FILE *);
  void write_file(FILE *);

  void create_pointers();

	// initialization and finalization
	void initialize();
	void finalize();

	//3a
	void compute_euler_characteristic();
};


#endif /* __LEARNPLY_H__ */

