#ifndef NOISE_H
#define NOISE_H

// Forward declarations for noise functions
float fade(float t);
float lerp(float t, float a, float b);
float grad(int hash, float x, float y, float z);
float perlin_noise(float x, float y, float z);
float octave_noise(float x, float y, float z, int octaves, float persistence);
float noise(float x, float y, float z);

#endif // NOISE_H 