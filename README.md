# SunlightHNC

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://hsugawa8651.github.io/SunlightHNC.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://hsugawa8651.github.io/SunlightHNC.jl/dev/)
[![Build Status](https://github.com/hsugawa8651/SunlightHNC.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/hsugawa8651/SunlightHNC.jl/actions/workflows/CI.yml?query=branch%3Amain)


This is a Julia port of SunlightHNC 

General purpose hypernetted chain (HNC) integral equation code

Copyright (C) 2025 Hiroharu Sugawara, Tokyo Metropolitan University <hsugawa@tmu.ac.jp>

Based on SunlightHNC original Version 1.13
https://github.com/patrickbwarren/SunlightHNC/releases/tag/v1.13

SunlightDPD is based on an original code copyright © 2007 Lucian Anton, with modifications copyright © 2008, 2009 Andrey Vlasov, and additional modifications copyright © 2009-2018 Unilever UK Central Resources Ltd (Registered in London number 29140; Registered Office: Unilever House, 100 Victoria Embankment, London EC4Y 0DY, UK). Later modifications copyright © 2020-2025 Patrick B Warren (STFC).

# Status of this Julia port
* This package is written in Julia language only.
* This package refers only to registered Julia packages.
* Matrix calculations depend on proven codes based on LAPACK and FFTW.
* Create as many HNC objects as you like, allowing parallel computations.
* No build process required. Install this package, then you are ready to go.
* Most of commpand line tools located in the original `example` directory have been ported to Julia.
* Julia runs fast on various environments. Develop codes in your PC and run them on high-performance computers.
* No GPU support yet.

## LICENSE
* GPL v2
