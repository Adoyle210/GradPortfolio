/*

Data structures for learnply

Eugene Zhang, 2005

*/


#define _USE_MATH_DEFINES
#include <cassert>
#include "polyhedron.h"
#include "learnply_io.h"
#include <stdlib.h>
#include <math.h>
#include <algorithm>

static PlyFile *in_ply;



;//////////////////////////////////////////////////////////////////////////////
;// Constructors
;//////////////////////////////////////////////////////////////////////////////

// read in a polyhedron from a file
Polyhedron::Polyhedron(FILE *file)
{
	int i,j;
	int elem_count;
	char *elem_name;

	/*** Read in the original PLY object ***/
	in_ply = read_ply (file);

	for (i = 0; i < in_ply->num_elem_types; i++) 
	{
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
	for (i = 0; i < ntris; i++) 
	{
		tlist[i]->verts[0] = vlist[(intptr_t) tlist[i]->verts[0]];
		tlist[i]->verts[1] = vlist[(intptr_t) tlist[i]->verts[1]];
		tlist[i]->verts[2] = vlist[(intptr_t) tlist[i]->verts[2]];
	}

	/* get rid of triangles that use the same vertex more than once */
	for (i = ntris-1; i >= 0; i--) 
	{
		Triangle *tri = tlist[i];
		Vertex *v0 = tri->verts[0];
		Vertex *v1 = tri->verts[1];
		Vertex *v2 = tri->verts[2];

		if (v0 == v1 || v1 == v2 || v2 == v0) 
		{
			free (tlist[i]);
			ntris--;
			tlist[i] = tlist[ntris];
		}
	}
}

// empty constructor
Polyhedron::Polyhedron()
{
	nverts = nedges = ntris = 0;
	max_verts = max_tris = 50;

	vlist = new Vertex * [max_verts];
	tlist = new Triangle * [max_tris];
}



;//////////////////////////////////////////////////////////////////////////////
;// Utility Functions
;//////////////////////////////////////////////////////////////////////////////

// write out a polyhedron to a file
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

	/* 2 */
	for (i = 0; i < ncorners; ++i)
		delete clist[i];
	delete[] clist;
	delete[] ctable;

	free(tlist);
	free(elist);
	free(vlist);
	if (!vert_other)
		free(vert_other);
	if (!face_other)
		free(face_other);

	silhouette.clear();
	major_hatches.clear();
	minor_hatches.clear();
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

  /* 2 */
  e->corners[0] = e->corners[1] = NULL;

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

	for (i=0; i<nverts; i++)
	{
		if (i==0)  {
			min.set(vlist[i]->x, vlist[i]->y, vlist[i]->z);
			max.set(vlist[i]->x, vlist[i]->y, vlist[i]->z);
		}
		else 
		{
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

int Polyhedron::calc_face_normals_and_area()
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
	if (signedvolume < 0)
		return 0; // return orientation
	else 
	{
		for (i=0; i<ntris; i++) // flip normals
			tlist[i]->normal *= -1.0;
		return 1; // return orientation
	}
}



;//////////////////////////////////////////////////////////////////////////////
;// Project 1 functions
;//////////////////////////////////////////////////////////////////////////////

/* 3a */
void Polyhedron::compute_euler_characteristic()
{
	int x = nverts - nedges + ntris;
	printf("Euler characteristic: %d\n", x);
}

/* 3b */
static double compute_angle(Corner *c)
{
	icVector3 v1(c->v->x, c->v->y, c->v->z);
	icVector3 v2(c->n->v->x, c->n->v->y, c->n->v->z);
	icVector3 v3(c->p->v->x, c->p->v->y, c->p->v->z);

	icVector3 v21 = v2 - v1;
	normalize(v21);

	icVector3 v31 = v3 - v1;
	normalize(v31);

	return acos(dot(v21, v31));
}

/* 3b */
void Polyhedron::compute_gaussian_curvature_angle_deficit()
{
	K_min = INFINITY, K_max = -INFINITY;

	double K_disc = 0.0;
	for (int i = 0; i < nverts; ++i)
	{
		double ang_tot = 0.0;

		Vertex *v = vlist[i];
		for (int j = 0; j < v->ncorners; ++j)
			ang_tot += v->corners[j]->angle;

		v->K = 2.0 * M_PI - ang_tot;
		K_disc += v->K;

		K_min = std::min(K_min, v->K);
		K_max = std::max(K_max, v->K);
	}

	printf("Discrete Gaussian Curvature (total angle deficit): %g\n", K_disc);
	printf("Minimum/Maximum Vertex Curvature: %g/%g\n", K_min, K_max);
}

/* 3c */
void Polyhedron::compute_gaussian_curvature_valence_deficit()
{
	min_valence_deficit = 6;

	int total_deficit = 0;
	for (int i = 0; i < nverts; ++i)
	{
		vlist[i]->valence_deficit = 6 - vlist[i]->ntris;
		//vlist[i]->K = M_PI / 3.0 * vlist[i]->valence_deficit;
		//K_disc += vlist[i]->K;
		total_deficit += vlist[i]->valence_deficit;

		min_valence_deficit = std::min(min_valence_deficit, vlist[i]->valence_deficit);
	}

	printf("Total Valence Deficit: %d\n", total_deficit);
	printf("Minimum Valence Deficit: %d\n", min_valence_deficit);
}

/* 3d */
void Polyhedron::compute_handles()
{
	int x = nverts - nedges + ntris;
	int nhandles = (2 - x) / 2;
	printf("Handle count: %d\n", nhandles);
}

/* 2 */
void Polyhedron::create_corners()
{
	int ci;

	/* always 3 corners per triangle */
	max_corners = ncorners = 3 * ntris;
	clist = new Corner *[max_corners];
	for (int i = 0; i < ncorners; ++i)
	{
		clist[i] = new Corner;
		clist[i]->index = i;
	}

	/* setup corners */
	ci = 0;
	for (int i = 0; i < ntris; ++i)
	{
		Triangle *t = tlist[i];
		for (int j = 0; j < 3; ++j)
		{
			Corner *c = clist[ci + j];

			/* assign edge/vertex/triangle */
			c->v = t->verts[j];
			c->e = t->edges[(j+1)%3]; // see edge creation for why this works
			c->t = t;
			t->corners[j] = c;

			/* compute next/previous */
			c->n = clist[ci + (j + 1) % 3];
			c->p = clist[ci + (j + 2) % 3];

			/* store corner in vertex */
			c->v->ncorners++;
			Corner **corners = (Corner **)realloc(c->v->corners, sizeof(Corner *) * c->v->ncorners);
			if (!corners)
				abort();
			corners[c->v->ncorners - 1] = c;
			c->v->corners = corners;
		}

		ci += 3;
	}

	/* construct corner table */
	ctable = new CornerTableEntry[ncorners];
	for (int i = 0; i < ncorners; ++i)
	{
		ctable[i].c = i;
		ctable[i].v_min = std::min(clist[i]->p->v->index, clist[i]->n->v->index);
		ctable[i].v_max = std::max(clist[i]->p->v->index, clist[i]->n->v->index);
		ctable[i].o = 0;
	}

	/* sort corner table */
	std::sort(ctable, ctable + ncorners, [](CornerTableEntry &a, CornerTableEntry &b) {
		/* returns an entry is "less than" another */
		if (a.v_min != b.v_min)
			return a.v_min < b.v_min;
		if (a.v_max != b.v_max)
			return a.v_max < b.v_max;
		return a.c < b.c;
	});

	/* link opposite corners */
	for (int i = 0; i < ncorners; i += 2)
	{
		CornerTableEntry &c = ctable[i];
		CornerTableEntry &o = ctable[i+1];

		/* link in corner table */
		c.o = o.c;
		o.o = c.c;

		/* create link in corner structure */
		clist[c.c]->o = clist[c.o];
		clist[o.c]->o = clist[o.o];

		/* track each corner in its edge */
		assert(clist[c.c]->e == clist[o.c]->e);
		clist[c.c]->e->corners[0] = clist[c.c];
		clist[o.c]->e->corners[1] = clist[c.o];
	}

	/* calculate corner angles */
	for (int i = 0; i < ncorners; ++i)
	{
		Corner *c = clist[i];
		c->angle = compute_angle(c);
	}

	/* print corner table */
	/*for (int i = 0; i < ncorners; ++i)
	{
		CornerTableEntry &c = ctable[i];
		printf("%5d %5d %5d %5d\n", c.c, c.v_min, c.v_max, c.o);
	}*/
}

/* 1c */
void Polyhedron::average_normals()
{
	int i, j;

	for (i = 0; i < nverts; i++) {
		vlist[i]->normal = icVector3(0.0);
		for (j = 0; j < vlist[i]->ntris; j++)
			vlist[i]->normal += vlist[i]->tris[j]->normal;
		normalize(vlist[i]->normal);
	}
}


;//////////////////////////////////////////////////////////////////////////////
;// Project 2 functions
;//////////////////////////////////////////////////////////////////////////////

/* Problem 1 */

void Polyhedron::compute_silhouette_edges(const icMatrix3x3& view, const icVector3& translate)
{
	
}

void Polyhedron::compute_silhouette_faces(const icMatrix3x3& view, const icVector3& translate)
{
	
}


/* Problem 2 */

void Polyhedron::compute_vert_voronoi_areas()
{
	
}

void Polyhedron::compute_vert_mean_curvature()
{
	
}

void Polyhedron::compute_vert_gaussian_curvature()
{

}

void Polyhedron::compute_vert_curvature_tensor()
{
	
}

void Polyhedron::update_vert_global_tensors()
{
	
}

void Polyhedron::smooth_vert_curvature_tensors(int weight_scheme, double step_size, int iterations)
{
	
}

void Polyhedron::compute_vert_principal_curvatures()
{
	
}

void Polyhedron::compute_face_principal_curvatures()
{

}

void Polyhedron::build_curvature_hatch_lines(int principal_direction, int transition_scheme)
{

}

void Polyhedron::hatch_line_step(int principal_direction, int transition_scheme, int iterations,
	Edge* current_edge, Triangle* current_face, icVector3& current_loc, icVector3& current_dir)
{

}