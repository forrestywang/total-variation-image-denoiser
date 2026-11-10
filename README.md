# Total Variation Image Denoiser

An end-to-end MATLAB implementation of **Total Variation image denoising**. This project discretizes the governing PDE using a 5-point finite difference scheme on a structured 2D grid and solves the linearized system using a custom **Preconditioned Conjugate Gradient (PCG)** solver accelerated by **Incomplete Cholesky factorization**.

For the complete mathematical background, refer to the following [report](docs/theory.pdf).

## Repository Structure

```
total-variation-image-denoising/
├── src/
│   ├── form_matrix.m           % Assembles sparse SPD matrix A(u) from finite difference stencils
│   ├── pcg_solver.m            % Custom Preconditioned Conjugate Gradient solver
│   └── set_image.m             % Generates target image geometry and synthetic Gaussian noise
├── docs/
│   ├── theory.pdf              % Compiled LaTeX technical report detailing stencils and derivations
│   └── theory.tex              % LaTeX source code
├── results/
│   └── TV_Denoiser.pdf         % Exported MATLAB Live Script technical report
│   └── denoised_results.pdf    % Vector graphics export comparing original, noisy, and restored images
├── TV_Denoiser.mlx             % Main MATLAB Live Script runner & benchmark suite
├── LICENSE                     % MIT License
└── README.md

```

## Benchmarks & Performance

Evaluated on a 512 × 512 **grid** (N = 262,144 total unknowns) with a noise-to-signal ratio of 1.00 across K = 8 outer lagged diffusivity iterations:

* **Convergence:** Total iterations bounded under **1,000 cumulative inner PCG steps** across all 8 outer steps.
* **Runtime:** Complete restoration executes in 25 **seconds**.
* **Edge Preservation:** Successfully suppresses heavy Gaussian noise while preserving crisp boundary contours without structural bleeding.

## Getting Started

### Prerequisites

* MATLAB R2021a or newer (Signal Processing / Image Processing Toolbox recommended).

### Running the Code

1. Clone the repository:

   ```
   git clone https://github.com/your-username/total-variation-image-denoising.git
   cd total-variation-image-denoising
   ```

2. Open `TV_Denoiser.mlx` in MATLAB and run all sections to execute the solver, display the iteration convergence table, and export visual results.

## License

Distributed under the [MIT License](LICENSE).
