# Watercolor Effect Attempt

Curvature-driven watercolor rendering experiment built on top of a PLY mesh viewer.

## Summary
This final project computes mesh curvature, maps that information into a stylized watercolor palette, and layers in noise and granulation to simulate a non-photorealistic painted look.

## Tech Stack
- C++
- OpenGL / GLUT
- PLY mesh data
- Curvature analysis and watercolor simulation

## Features
- Mean and Gaussian curvature analysis
- Curvature-based pigment placement
- Watercolor simulation with flow, evaporation, and granulation
- Noise-driven surface variation
- Keyboard toggles for watercolor controls

## Controls
- `w` toggles watercolor mode
- `g` / `G` adjust granulation
- `n` / `N` adjust noise scale
- `r` resets watercolor parameters

## Build

### Setup
Install CMake 3.10 or newer and add it to your system `PATH`. Create a `build` directory and work from inside it.

#### Windows
```bat
md build
cd build
```

#### Mac/Linux
```sh
mkdir build
cd build
```

### Generate Build Files
Run this again any time you add or remove files.

#### Windows
```bat
cmake -A Win32 ..
```

#### Mac/Linux
```sh
cmake ..
```

### Build and Run
#### Windows (Visual Studio)
- Set the `learnply` project as the startup project before running.

#### Mac/Linux
```sh
make
./learnply
```
