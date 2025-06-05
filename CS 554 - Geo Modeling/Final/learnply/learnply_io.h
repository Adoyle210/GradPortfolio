/*

Data structure for I/O of polygonal models

Eugene Zhang 2005

*/

#ifndef __LEARNPLY_IO_H__
#define __LEARNPLY_IO_H__

#include "ply.h"

typedef struct Vertex_io
{
  float x, y, z;
  int angle_deficit; // Changed back to int32 to match PLY file
  void *other_props; /* other properties */
} Vertex_io;

typedef struct Face_io
{
  unsigned char nverts; /* number of vertex indices in list */
  int *verts;           /* vertex index list */
  int ancestor;         // Add this property
  void *other_props;    /* other properties */
} Face_io;

/* TODO: casting literals to (char *) to surpress warnings, should restructure program to fix this */

char *elem_names[] = {/* list of the kinds of elements in the user's object */
                      (char *)"vertex", (char *)"face"};

PlyProperty vert_props[] = {
    /* list of property information for a vertex */
    {(char *)"x", Float32, Float32, offsetof(Vertex_io, x), 0, 0, 0, 0},
    {(char *)"y", Float32, Float32, offsetof(Vertex_io, y), 0, 0, 0, 0},
    {(char *)"z", Float32, Float32, offsetof(Vertex_io, z), 0, 0, 0, 0},
    {(char *)"angle_deficit", Int32, Int32, offsetof(Vertex_io, angle_deficit), 0, 0, 0, 0}};

PlyProperty face_props[] = {
    /* list of property information for a face */
    {(char *)"vertex_indices", Int32, Int32, offsetof(Face_io, verts),
     1, Uint8, Uint8, offsetof(Face_io, nverts)},
    {(char *)"ancestor", Int32, Int32, offsetof(Face_io, ancestor), 0, 0, 0, 0},
};

#endif /* __LEARNPLY_IO_H__ */
