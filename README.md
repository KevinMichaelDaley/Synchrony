# Synchrony
Fast, parallel complex network and nonlinear dynamics toolbox in C++ with GPU support.

Developed in conjunction with the [Belykh lab](https://math.gsu.edu/ibelykh/belykh_lab.html "Our gsu.edu lab homepage"), this project aims to be a common, modular framework for simulating complex networks and nonlinear dynamical systems, with a special focus on applications our lab is interested in (pedestrian and suspension bridge modeling, dynamical neuroscience, stochastic switching networks and multilayer networks, etc.)

Synchrony aims to be extremely fast, and to that end supports CUDA via clang's nvidia backend.  We are also planning a separate python/cython mode for researchers who are more comfortable using dynamic languages.

Currently this is a work in heavy progress; however, numerical integration is already fairly stable and well-tested. 

Required Dependencies: g++ or clang++ with c++17 support and GNU make (or equivalent)
Optional Dependencies (for GPU mode): recent clang++ (tested with v6.0) with CUDA and c++17 support, cuda toolkit 8.0 (install the [version](https://developer.nvidia.com/cuda-80-ga2-download-archive "developer.nvidia.com") from the NVIDIA website.  You may have to edit the makefiles if you are on windows or installing to anywhere other than /usr/local; we are working on fixing this).

