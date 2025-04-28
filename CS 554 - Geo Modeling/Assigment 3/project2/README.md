# learnply

## Building

### Setup

The project uses CMake to generate a build system. Install a version of CMake of at least version 3.10 or higher and add it to your system `PATH`. CMake requires a build directory, from which it will work out of. Create one called `build` and change the working directory to be in it. You only need to do this step once.

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

### Build System Generation

Any time you add or remove files, you must do this step again. Make sure you are working from within the `build` directory.

#### Windows

```bat
cmake -A Win32 ..
```

**Note:** You must include the `-A Win32` switch in the CMake command, as the executable must be a 32-bit executable.

#### Mac/Linux

```sh
cmake ..
```

### Building

To build, use the build system CMake has generated. For Windows, the default is the latest Visual Studio installation and for Linux/Mac, the default is Unix Makefiles.

#### Windows (Visual Studio)

When building and running on Visual Studio, be sure to select the **learnply** project as the **startup project** through the context menu by right clicking on it in the solution explorer.

#### Mac/Linux

From the build directory, use `make` to build the program and `./learnply` to run the program.
