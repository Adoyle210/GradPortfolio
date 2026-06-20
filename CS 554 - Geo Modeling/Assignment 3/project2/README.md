# Curvature and Hatch Lines

PLY mesh analysis project focused on curvature computation and stylized rendering.

## Summary
This assignment loads a PLY mesh and computes geometric properties such as Gaussian curvature, mean curvature, curvature tensors, and principal curvatures, then uses those results to render curvature-based hatch lines.

## Tech Stack
- C++
- OpenGL / GLUT
- PLY mesh data
- Curvature analysis routines

## Features
- Mesh loading and triangle reconstruction
- Gaussian and mean curvature computation
- Curvature tensor smoothing and principal curvature extraction
- Curvature-driven hatch line generation
- Silhouette and tensor visualization modes

## Controls
- Number keys switch display modes
- `s` cycles silhouette modes
- `t` cycles tensor display modes
- `k` toggles rendering mode
- `w` toggles wireframe

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
