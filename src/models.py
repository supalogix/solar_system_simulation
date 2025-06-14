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
from helpers import *
from constants import G, solar_mass
from operator import itemgetter

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
        n = compute_mean_motion(G.value, solar_mass.value, self.a)

        # Recover the initial eccentric anomaly from the stored true anomaly
        E0 = recover_initial_eccentric_anomaly(self.e, self.nu)

        # Compute the initial mean anomaly 
        M0 = compute_initial_mean_anomaly(E0, self.e)

        # Compute the current mean anomaly
        M = compute_current_mean_anomaly(M0, n, t)

        # Solve Kepler’s equation for the eccentric anomaly
        E = solve_keplers_equation_for_the_eccentric_anomaly(M, self.e)

        # Compute the true anomaly from the eccentric anomaly
        nu = compute_true_anomaly(self.e, E)

        # Compute the orbital radius from the eccentric anomaly
        r = compute_orbital_radius(
            a = self.a, 
            e = self.e, 
            E = E
        )

        return {
            'n': n, 
            'E0': E0, 
            'M0': M0, 
            'M': M, 
            'E': E, 
            'nu': nu, 
            'r': r
        }


    def position(self, t):
        """
        Return the 3D position (m) at time t, correctly accounting for
        initial true anomaly nu, argument of periapsis, node, and inclination.
        """

        #nu, r = self.generate_parameters(
        #    t=t
        #)
        nu, r = itemgetter('nu', 'r')(self.generate_parameters(t=t))

        # Coordinates in the orbital plane (before any 3D rotations)
        x_orb, y_orb, z_orb = compute_coordinates_in_orbital_plane(
            r=r, 
            nu=nu
        )

        # Rotate the orbital‐plane coordinates by argument of periapsis about the z‐axis
        x1, y1, z1 = rotate_orbital_plane_coordinates(
            x_orb=x_orb,
            y_orb=y_orb,
            z_orb=z_orb,
            omega=self.omega
        )

        # Rotate by inclination about the x‐axis
        x2, y2, z2 = rotate_by_x_axis_inclination(
            x1=x1,
            y1=y1,
            z1=z1,
            i=self.i
        )

        # Rotate by longitude of ascending node about the z‐axis
        x, y, z = rotate_by_z_axis_longitude_ascending_node(
            x2=x2,
            y2=y2,
            z2=z2,
            Omega=self.Omega
        )

        return np.array([x, y, z])

    def velocity(self, t):
        """
        Return the 3D velocity vector (m/s) at time t, correctly accounting for
        the initial true anomaly nu, argument of periapsis, inclination, and node.
        """

        n, nu, E, r = itemgetter('n', 'nu', 'E', 'r')(self.generate_parameters(t=t))

        # Compute time derivative
        dE_dt, r_dot = compute_time_derivative(
            n=n,
            a=self.a,
            e=self.e,
            E=E
        )

        # Compute the rate of change of true anomaly:
        dnu_dt = compute_true_anomaly_rate_of_change(
            e=self.e,
            E=E,
            dE_dt=dE_dt
        )

        # Compute the in‐plane velocity components in the orbital frame
        vx_orb, vy_orb, vz_orb = compute_orbital_frame_in_plane_velocity_components(
            r=r,
            r_dot=r_dot,
            nu=nu,
            dnu_dt=dnu_dt
        )

        # Rotate the in‐plane velocity by the argument of periapsis about the z‐axis
        v1x, v1y, v1z = rotate_in_plane_velocity(
            vx_orb=vx_orb,
            vy_orb=vy_orb,
            vz_orb=vz_orb,
            omega=self.omega
        )

        # Rotate the velocity by the inclination about the x‐axis
        v2x, v2y, v2z = rotate_x_axis_inclination_velocity(
            v1x=v1x,
            v1y=v1y,
            v1z=v1z,
            i=self.i
        )

        # Rotate the velocity by the longitude of ascending node about the z‐axis
        vx, vy, vz = rotate_velocity_by_longitude(
            v2x=v2x,
            v2y=v2y,
            v2z=v2z,
            Omega=self.Omega
        )

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
