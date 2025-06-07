import numpy as np
import matplotlib.pyplot as plt

from src.physics.solvers import generate_prolate_spheroidal_grid
from src.physics.core.poisson import solve_poisson_hyperbolic

if __name__ == "__main__":
    eta_grid, nu_grid = generate_prolate_spheroidal_grid(
        tip_radius=1.0, num_eta=40, num_nu=40, f=1.0
    )

    phi = solve_poisson_hyperbolic(
        eta_grid, nu_grid, phi_tip=1.0, phi_sample=0.0, omega=1.6, max_iter=5000
    )

    f = 1.0
    r = f * np.sqrt((1 - nu_grid ** 2) * (eta_grid ** 2 - 1))
    z = f * eta_grid * nu_grid

    plt.figure(figsize=(6, 5))
    cs = plt.contourf(r, z, phi, levels=20)
    plt.colorbar(cs, label="Potential (V)")
    plt.xlabel("r (nm)")
    plt.ylabel("z (nm)")
    plt.title("Laplace solution in prolate spheroidal coords")
    plt.axis("equal")
    plt.tight_layout()
    plt.show()

