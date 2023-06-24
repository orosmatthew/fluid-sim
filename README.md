# Fluid Sim

## 3D Simulation

![](https://raw.githubusercontent.com/orosmatthew/fluid-sim/master/media/fluid-sim-3d.gif)

## 2D Simulation

![](https://raw.githubusercontent.com/orosmatthew/fluid-sim/master/media/fluid-sim-2d.gif)

This simulation was created as a final project for my Computational Physics class at Baldwin Wallace University.

This project contains both a 2D and 3D fluid simulation. All 2D related code is in the `src/` directory while all 3D
related code is in the `3d/` directory.

The 2D simulation uses [Raylib](https://www.raylib.com/) for visualization while the 3D simulation uses a custom Vulkan
library called Mini Vulkan Engine (MVE) written from scratch. Thus, the 3D simulation requires a Vulkan capable GPU to
run.

A web export of the simulation is also provided of the 2D simulation without multithreading enabled.

Windows binaries and the web export can be found under `Releases`.

## References

* Fluid Simulation for
  Dummies: [https://mikeash.com/pyblog/fluid-simulation-for-dummies.html](https://mikeash.com/pyblog/fluid-simulation-for-dummies.html)

* The Coding Train - Fluid Simulation: https://thecodingtrain.com/challenges/132-fluid-simulation

* Raymarching Clouds (Part 1): https://davidpeicho.github.io/blog/cloud-raymarching-walkthrough-part1/

* Raymarching Clouds (Part 2): https://davidpeicho.github.io/blog/cloud-raymarching-walkthrough-part2/