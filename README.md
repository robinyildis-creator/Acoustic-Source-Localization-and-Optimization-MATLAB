# Acoustic Source Localization and Optimization (MATLAB)

This project studies **sound waves in a closed room** using numerical methods and solves both forward and inverse problems based on the Helmholtz equation.

The goal is to:

* **identify the location and strength of a sound source** from boundary measurements
* **optimize source placement** to minimize sound in a specific region

---

## Overview

The sound field is modeled using the Helmholtz equation, where the source is assumed to be localized and isotropic.

The project is divided into two main parts:

* **Inverse problem:**
  Estimate the source position and amplitude from measurements of the normal derivative on the room boundary using least squares and Gauss–Newton methods.

* **Optimization problem:**
  Find the optimal placement of a sound source (e.g. a TV) to minimize sound levels in a designated area of the room.

---

## Methods

* Numerical integration (2D trapezoidal rule)
* Gauss–Newton and Levenberg–Marquardt optimization
* Fourier-based integral formulation
* Golden section search
* Nonlinear optimization (`fminsearch`)
* Noise robustness analysis

---

## Files

* `Part-1.m` – Inverse problem (source reconstruction)
* `Part-2.m` – Optimal source placement
* `ProjektCtmat2025.pdf` – Project description

---

## Key Features

* Reconstruction of source location from boundary data
* Comparison of optimization methods (GN vs LM) 
* Robustness analysis under measurement noise
* 2D optimization of sound source placement
* Visualization of acoustic fields in complex room geometry 

---

## Tech

* MATLAB
* Numerical PDE methods
* Optimization & inverse problems

---

## Description

MATLAB project on acoustic wave modeling, inverse source localization, and optimal placement using numerical methods and optimization techniques.
