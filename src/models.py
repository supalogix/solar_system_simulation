r"""
Planet class implementation based on the following references:

    1. Bate, R. R., Mueller, D. D., & White, J. E. (1971).
    Fundamentals of Astrodynamics.
    Dover Publications.

    2. Vallado, D. A. (2013).
    Fundamentals of Astrodynamics and Applications (4th ed.).
    Microcosm Press / Springer.
"""


import numpy as np
from constants import G, solar_mass

class Planet:
    def __init__(self, name, type, a, e, i, Omega, omega, nu, mass):
        """
        Create a Planet with its orbital elements and physical mass.

        Args:
            name (str): Identifier for the planet.
            type (str): Type of the planet.
            a (float): Semi-major axis of the orbit (meters).
            e (float): Orbital eccentricity (0=circle, <1=ellipse).
            i (float): Inclination relative to the reference plane (radians).
            Omega (float): Longitude of the ascending node (radians).
            omega (float): Argument of periapsis (radians).
            nu (float): True anomaly at epoch (radians), initial angular position.
            mass (float): Mass of the planet (kilograms).
        """

        self.name  = name
        self.type = type
        self.a     = a  
        self.e     = e  
        self.i     = np.radians(i)
        self.Omega = np.radians(Omega)
        self.omega = np.radians(omega)
        self.nu    = np.radians(nu)
        self.mass  = mass  

    def generate_parameters(self, t): 
        # Compute the mean motion
        n = np.sqrt(G.value * solar_mass.value / self.a**3)

        # Recover the initial eccentric anomaly from the stored true anomaly
        E0 = 2*np.arctan2(
            np.sqrt(1 - self.e) * np.sin(self.nu / 2),
            np.sqrt(1 + self.e) * np.cos(self.nu / 2),
        )

        # Compute the initial mean anomaly 
        M0 = E0 - self.e * np.sin(E0)

        # Compute the current mean anomaly
        M = M0 + n * t

        # Solve Kepler’s equation for the eccentric anomaly
        E = M
        for _ in range(20):
            E = M + self.e * np.sin(E)

        # Compute the true anomaly from the eccentric anomaly
        nu = 2*np.arctan2(
            np.sqrt(1 + self.e) * np.sin(E / 2),
            np.sqrt(1 - self.e) * np.cos(E / 2)
        )

        return n, E0, M0, M, E, nu


    def position(self, t):
        """
        Return the 3D position (m) at time t, correctly accounting for
        initial true anomaly nu, argument of periapsis, node, and inclination.
        """

        n, E0, M0, M, E, nu = self.generate_parameters(t)


        # Compute the orbital radius from the eccentric anomaly
        r = self.a * (1 - self.e * np.cos(E))

        # Coordinates in the orbital plane (before any 3D rotations)
        x_orb = r * np.cos(nu)
        y_orb = r * np.sin(nu)
        z_orb = 0.0

        # Rotate the orbital‐plane coordinates by argument of periapsis ω about the z‐axis
        x1 = x_orb * np.cos(self.omega) - y_orb * np.sin(self.omega)
        y1 = x_orb * np.sin(self.omega) + y_orb * np.cos(self.omega)
        z1 = 0.0

        # Rotate by inclination about the x‐axis
        x2 = x1
        y2 = y1 * np.cos(self.i) - z1 * np.sin(self.i)
        z2 = y1 * np.sin(self.i) + z1 * np.cos(self.i)

        # Rotate by longitude of ascending node Ω about the z‐axis
        x = x2 * np.cos(self.Omega) - y2 * np.sin(self.Omega)
        y = x2 * np.sin(self.Omega) + y2 * np.cos(self.Omega)
        z = z2

        return np.array([x, y, z])

    def velocity(self, t):
        """
        Return the 3D velocity vector (m/s) at time t, correctly accounting for
        the initial true anomaly nu, argument of periapsis, inclination, and node.
        """

        n, E0, M0, M, E, nu = self.generate_parameters(t)

        # Compute the orbital radius and its time derivative
        r = self.a * (1 - self.e * np.cos(E))
        dE_dt = n / (1 - self.e * np.cos(E))
        r_dot = self.a * self.e * np.sin(E) * dE_dt

        # Compute the rate of change of true anomaly nu:
        dnu_dt = (np.sqrt(1 - self.e**2) * dE_dt) / (1 - self.e * np.cos(E))

        # Compute the in‐plane velocity components in the orbital frame
        vx_orb = r_dot * np.cos(nu) - r * dnu_dt * np.sin(nu)
        vy_orb = r_dot * np.sin(nu) + r * dnu_dt * np.cos(nu)
        vz_orb = 0.0

        # Rotate the in‐plane velocity by the argument of periapsis about the z‐axis
        v1x = vx_orb * np.cos(self.omega) - vy_orb * np.sin(self.omega)
        v1y = vx_orb * np.sin(self.omega) + vy_orb * np.cos(self.omega)
        v1z = 0.0

        # Rotate the velocity by the inclination about the x‐axis
        v2x = v1x
        v2y = v1y * np.cos(self.i) - v1z * np.sin(self.i)
        v2z = v1y * np.sin(self.i) + v1z * np.cos(self.i)


        # Rotate the velocity by the longitude of ascending node about the z‐axis
        vx =  v2x * np.cos(self.Omega) - v2y * np.sin(self.Omega)
        vy =  v2x * np.sin(self.Omega) + v2y * np.cos(self.Omega)
        vz =  v2z

        return np.array([vx, vy, vz])


    def momentum(self, t):
        """
        Compute the planet’s linear momentum at a given time.

        Args:
            t (float): Time in seconds at which to evaluate the velocity.

        Returns:
            numpy.ndarray: Momentum vector (kg * m/s), given by mass * velocity.
        """

        # Linear momentum: mass × velocity gives a vector that captures how much
        # “inertia” the planet has in its motion. We use this momentum to compute
        # its speed (and thus kinematic time dilation) without recalculating velocity
        # components separately.
        return self.mass * self.velocity(t)
