# Mesh Visualization Modes

PLY mesh viewer with multiple coloring and interaction modes.

## Summary
This assignment loads a bunny mesh from `.ply` data and visualizes it with several display modes, including ID coloring, barycentric coordinates, normal-based shading, and a 3D checkerboard effect.

## Tech Stack
- C++
- OpenGL / GLUT
- PLY mesh data

## Features
- Triangle picking and seed selection
- Polygon ID coloring with generated colors
- Barycentric coordinate coloring
- Polygon and vertex normal coloring
- 3D checkerboard coloring with adjustable cell size

## Controls
- Number keys switch display modes
- `x` toggles anti-aliasing
- `+` and `-` adjust checker size
- Mouse drag rotates the model

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
