# Pysemitip: A Python Modernization of the SEMITIP STM Simulation Software

## Overview

The primary goal of the Pysemitip project is to translate the original Fortran Semitip program, developed by R.M. Feenstra's group at Carnegie Mellon University (CMU), into a modern, maintainable, and extensible Python library. Semitip is a physics simulation tool designed for simulating Scanning Tunneling Microscopy (STM) measurement results.

The motivation behind this modernization effort is to facilitate further development, enhance maintainability, and enable seamless integration with the rich ecosystem of modern scientific Python tools and libraries.

## Project Status

This project is currently in the **initial analysis and structured translation phase**.

## Core Methodology

Our development is guided by the **"Analyze First, Translate Second"** philosophy. This means that before any Fortran code is translated into Python, a thorough understanding of its physical purpose and algorithmic structure is paramount.

This understanding is achieved by:
1.  Consulting the original physics papers and technical manuals located in the `/docs/Fortran-semitip/` directory.
2.  Analyzing the structure and interdependencies within the original Fortran source code found in `/src/fortran/`.

This approach is crucial to avoid the pitfalls of superficial, line-by-line translation, ensuring accuracy and preserving the integrity of the underlying physical models.

## Repository Structure

A brief overview of the key directories within this repository:

-   `src/`: The main source code directory.
    -   `fortran/`: Contains the complete original Fortran source code for Semitip, including the main programs `MultInt` and `MultPlane`.
    -   `core/`: Houses existing Python utility modules. Currently, these modules provide functionalities for reading and converting `fort.9` input files to a more human-readable YAML format.
    -   `p_spacecharge/`: This directory contains Pascal code, documentation, and related files from another project. It might hold relevant information or algorithms for understanding and implementing space charge calculations within Semitip.
-   `docs/`: The primary directory for all documentation.
    -   `Fortran-semitip/`: A critical collection of PDF documents, including published papers and technical manuals. These files detail the physical models, algorithms, and theoretical underpinnings of the Semitip software. They serve as the ground truth for the physics.
    -   `modules/`: Intended for documentation related to the new Python modules as they are developed.
-   `environment.yml`: Conda environment file for reproducing the Python development environment.
-   `requirements.txt`: Pip requirements file for setting up the Python development environment.

## Development Setup

You can set up the Python development environment using either Conda or pip.

**Using Conda (recommended):**

1.  Create the Conda environment from the `environment.yml` file:
    ```bash
    conda env create -f environment.yml
    ```
2.  Activate the environment:
    ```bash
    conda activate pysemitip
    ```

**Using pip:**

1.  Install the required packages using `requirements.txt`:
    ```bash
    pip install -r requirements.txt
    ```

## Development Workflow & Strategy

The project will proceed through the following distinct phases:

### Phase 1: Analysis & Documentation (Current Phase)

-   **Objective:** Develop a deep understanding of the original Fortran code and its underlying physics.
-   **Activities:**
    -   Systematically review each Fortran subroutine and function.
    -   Consult the provided PDF documents in `docs/Fortran-semitip/` to understand the physical purpose and mathematical formulation of each routine.
    -   Document findings, including subroutine functionalities, input/output parameters, dependencies, and connections to COMMON blocks.
    -   Map out the overall program flow and data structures.

### Phase 2: Bottom-Up Translation & Unit Testing

-   **Objective:** Translate individual, low-dependency Fortran subroutines into Python functions with verifiable correctness.
-   **Activities:**
    -   Start with utility and mathematical subroutines that have minimal dependencies.
    -   Translate these into equivalent Python functions.
    -   For each translated Python function, create a corresponding unit test. This test will validate the Python function's output against the output of the original Fortran routine using identical inputs.

### Phase 3: Integration

-   **Objective:** Assemble the validated Python functions into larger modules that replicate the functionality of the main Fortran programs.
-   **Activities:**
    -   Begin by reconstructing the `MultInt` program's functionality.
    -   Integrate the translated Python functions, ensuring correct data flow and interaction between them.
    -   Develop higher-level tests to verify the integrated modules.

### Phase 4: Refinement & Modernization

-   **Objective:** Refactor the direct translation into more idiomatic ("Pythonic") code and develop a user-friendly API.
-   **Activities:**
    -   Optimize Python code for clarity, performance, and maintainability.
    -   Leverage libraries like NumPy and SciPy more effectively for array operations and scientific computations.
    -   Design and implement a higher-level API for setting up and running STM simulations, abstracting away the lower-level details.
    -   Develop comprehensive examples and tutorials.

## How to Contribute

Contributions to Pysemitip are welcome. To ensure the quality and integrity of the project, please adhere to the following guidelines:

-   **Follow the Core Methodology:** Prioritize understanding the physics and original code structure before translation.
-   **Structured Translation:** Contribute by translating individual subroutines and providing corresponding unit tests.
-   **Documentation:** All new code, architectural decisions, and significant findings must be clearly documented.
-   **Code Quality:** Ensure your code is clean, well-commented, and adheres to Python best practices (e.g., PEP 8).
-   **Testing:** All contributions must be accompanied by relevant tests.

By following these principles, we aim to create a robust, accurate, and valuable Python library for the STM research community.