import numpy as np
import matplotlib.pyplot as plt

from src.physics.solvers import generate_prolate_spheroidal_grid

if __name__ == "__main__":
    eta_grid, nu_grid = generate_prolate_spheroidal_grid(
        tip_radius=1.0, num_eta=40, num_nu=40, f=1.0
    )

    # Convert to cylindrical coordinates for visualization
    f = 1.0
    r = f * np.sqrt((1 - nu_grid ** 2) * (eta_grid ** 2 - 1))
    z = f * eta_grid * nu_grid

    plt.figure(figsize=(5, 5))
    plt.scatter(r, z, s=5)
    plt.xlabel("r (nm)")
    plt.ylabel("z (nm)")
    plt.title("Prolate spheroidal grid")
    plt.axis("equal")
    plt.tight_layout()
    plt.show()

