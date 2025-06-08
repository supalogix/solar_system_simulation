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

def compute_mean_motion(G, solar_mass, a):
    n = np.sqrt(G * solar_mass / a**3)
    return n

def recover_initial_eccentric_anomaly(e, nu):
    E0 = 2*np.arctan2(
        np.sqrt(1 - e) * np.sin(nu / 2),
        np.sqrt(1 + e) * np.cos(nu / 2),
    )

    return E0

def compute_initial_mean_anomaly(E0, e):
    M0 = E0 - e * np.sin(E0)

    return M0

def compute_current_mean_anomaly(M0, n, t):
    M = M0 + n * t

    return M

def solve_keplers_equation_for_the_eccentric_anomaly(M, e):
    E = M
    for _ in range(20):
        E = M + e * np.sin(E)

    return E

def compute_true_anomaly(e, E):
    nu = 2*np.arctan2(
        np.sqrt(1 + e) * np.sin(E / 2),
        np.sqrt(1 - e) * np.cos(E / 2)
    )

    return nu

def compute_orbital_radius(a, e, E):
    r = a * (1 - e * np.cos(E))

    return r

def compute_coordinates_in_orbital_plane(r, nu):
    x_orb = r * np.cos(nu)
    y_orb = r * np.sin(nu)
    z_orb = 0.0

    return x_orb, y_orb, z_orb

def rotate_orbital_plane_coordinates(x_orb, y_orb, z_orb, omega):
    x1 = x_orb * np.cos(omega) - y_orb * np.sin(omega)
    y1 = x_orb * np.sin(omega) + y_orb * np.cos(omega)
    z1 = 0.0

    return x1, y1, z1

def rotate_by_x_axis_inclination(x1, y1, z1, i):
    x2 = x1
    y2 = y1 * np.cos(i) - z1 * np.sin(i)
    z2 = y1 * np.sin(i) + z1 * np.cos(i)

    return x2, y2, z2

def rotate_by_z_axis_longitude_ascending_node(x2, y2, z2, Omega):
    x = x2 * np.cos(Omega) - y2 * np.sin(Omega)
    y = x2 * np.sin(Omega) + y2 * np.cos(Omega)
    z = z2

    return x, y, z
    
def compute_time_derivative(n, a, e, E):
    dE_dt = n / (1 - e * np.cos(E))
    r_dot = a * e * np.sin(E) * dE_dt

    return dE_dt, r_dot

def compute_true_anomaly_rate_of_change(e, E, dE_dt):
    dnu_dt = (np.sqrt(1 - e**2) * dE_dt) / (1 - e * np.cos(E))

    return dnu_dt

def compute_orbital_frame_in_plane_velocity_components(r, r_dot, nu, dnu_dt):
    vx_orb = r_dot * np.cos(nu) - r * dnu_dt * np.sin(nu)
    vy_orb = r_dot * np.sin(nu) + r * dnu_dt * np.cos(nu)
    vz_orb = 0.0

    return vx_orb, vy_orb, vz_orb

def rotate_in_plane_velocity(vx_orb, vy_orb, vz_orb, omega):
    v1x = vx_orb * np.cos(omega) - vy_orb * np.sin(omega)
    v1y = vx_orb * np.sin(omega) + vy_orb * np.cos(omega)
    v1z = 0.0

    return v1x, v1y, v1z

def rotate_x_axis_inclination_velocity(v1x, v1y, v1z, i):
    v2x = v1x
    v2y = v1y * np.cos(i) - v1z * np.sin(i)
    v2z = v1y * np.sin(i) + v1z * np.cos(i)

    return v2x, v2y, v2z

def rotate_velocity_by_longitude(v2x, v2y, v2z, Omega):
    vx =  v2x * np.cos(Omega) - v2y * np.sin(Omega)
    vy =  v2x * np.sin(Omega) + v2y * np.cos(Omega)
    vz =  v2z

    return vx, vy, vz


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
        #n = np.sqrt(G.value * solar_mass.value / self.a**3)
        n = compute_mean_motion(G.value, solar_mass.value, self.a)

        # Recover the initial eccentric anomaly from the stored true anomaly
        # E0 = 2*np.arctan2(
        #     np.sqrt(1 - self.e) * np.sin(self.nu / 2),
        #     np.sqrt(1 + self.e) * np.cos(self.nu / 2),
        # )
        E0 = recover_initial_eccentric_anomaly(self.e, self.nu)

        # Compute the initial mean anomaly 
        #M0 = E0 - self.e * np.sin(E0)
        M0 = compute_initial_mean_anomaly(E0, self.e)

        # Compute the current mean anomaly
        #M = M0 + n * t
        M = compute_current_mean_anomaly(M0, n, t)

        # Solve Kepler’s equation for the eccentric anomaly
        #E = M
        #for _ in range(20):
        #    E = M + self.e * np.sin(E)
        E = solve_keplers_equation_for_the_eccentric_anomaly(M, self.e)

        # Compute the true anomaly from the eccentric anomaly
        #nu = 2*np.arctan2(
        #    np.sqrt(1 + self.e) * np.sin(E / 2),
        #    np.sqrt(1 - self.e) * np.cos(E / 2)
        #)
        nu = compute_true_anomaly(self.e, E)

        # Compute the orbital radius from the eccentric anomaly
        #r = self.a * (1 - self.e * np.cos(E))
        r = compute_orbital_radius(
            a = self.a, 
            e = self.e, 
            E = E
        )

        return n, E0, M0, M, E, nu, r


    def position(self, t):
        """
        Return the 3D position (m) at time t, correctly accounting for
        initial true anomaly nu, argument of periapsis, node, and inclination.
        """

        n, E0, M0, M, E, nu, r = self.generate_parameters(
            t=t
        )

        # Compute the orbital radius from the eccentric anomaly
        #r = self.a * (1 - self.e * np.cos(E))
        #r = compute_orbital_radius(
        #    a = self.a, 
        #    e = self.e, 
        #    E = E
        #)

        # Coordinates in the orbital plane (before any 3D rotations)
        #x_orb = r * np.cos(nu)
        #y_orb = r * np.sin(nu)
        #z_orb = 0.0
        x_orb, y_orb, z_orb = compute_coordinates_in_orbital_plane(
            r=r, 
            nu=nu
        )

        # Rotate the orbital‐plane coordinates by argument of periapsis about the z‐axis
        #x1 = x_orb * np.cos(self.omega) - y_orb * np.sin(self.omega)
        #y1 = x_orb * np.sin(self.omega) + y_orb * np.cos(self.omega)
        #z1 = 0.0
        x1, y1, z1 = rotate_orbital_plane_coordinates(
            x_orb=x_orb,
            y_orb=y_orb,
            z_orb=z_orb,
            omega=self.omega
        )

        # Rotate by inclination about the x‐axis
        #x2 = x1
        #y2 = y1 * np.cos(self.i) - z1 * np.sin(self.i)
        #z2 = y1 * np.sin(self.i) + z1 * np.cos(self.i)
        x2, y2, z2 = rotate_by_x_axis_inclination(
            x1=x1,
            y1=y1,
            z1=z1,
            i=self.i
        )

        # Rotate by longitude of ascending node about the z‐axis
        #x = x2 * np.cos(self.Omega) - y2 * np.sin(self.Omega)
        #y = x2 * np.sin(self.Omega) + y2 * np.cos(self.Omega)
        #z = z2
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

        n, E0, M0, M, E, nu, r = self.generate_parameters(t)

        # Compute the orbital radius and its time derivative
        #r = self.a * (1 - self.e * np.cos(E))
        #dE_dt = n / (1 - self.e * np.cos(E))
        #r_dot = self.a * self.e * np.sin(E) * dE_dt

        dE_dt, r_dot = compute_time_derivative(
            n=n,
            a=self.a,
            e=self.e,
            E=E
        )

        # Compute the rate of change of true anomaly:
        #dnu_dt = (np.sqrt(1 - self.e**2) * dE_dt) / (1 - self.e * np.cos(E))
        dnu_dt = compute_true_anomaly_rate_of_change(
            e=self.e,
            E=E,
            dE_dt=dE_dt
        )

        # Compute the in‐plane velocity components in the orbital frame
        #vx_orb = r_dot * np.cos(nu) - r * dnu_dt * np.sin(nu)
        #vy_orb = r_dot * np.sin(nu) + r * dnu_dt * np.cos(nu)
        #vz_orb = 0.0
        vx_orb, vy_orb, vz_orb = compute_orbital_frame_in_plane_velocity_components(
            r=r,
            r_dot=r_dot,
            nu=nu,
            dnu_dt=dnu_dt
        )

        # Rotate the in‐plane velocity by the argument of periapsis about the z‐axis
        #v1x = vx_orb * np.cos(self.omega) - vy_orb * np.sin(self.omega)
        #v1y = vx_orb * np.sin(self.omega) + vy_orb * np.cos(self.omega)
        #v1z = 0.0
        v1x, v1y, v1z = rotate_in_plane_velocity(
            vx_orb=vx_orb,
            vy_orb=vy_orb,
            vz_orb=vz_orb,
            omega=self.omega
        )

        # Rotate the velocity by the inclination about the x‐axis
        #v2x = v1x
        #v2y = v1y * np.cos(self.i) - v1z * np.sin(self.i)
        #v2z = v1y * np.sin(self.i) + v1z * np.cos(self.i)
        v2x, v2y, v2z = rotate_x_axis_inclination_velocity(
            v1x=v1x,
            v1y=v1y,
            v1z=v1z,
            i=self.i
        )


        # Rotate the velocity by the longitude of ascending node about the z‐axis
        #vx =  v2x * np.cos(self.Omega) - v2y * np.sin(self.Omega)
        #vy =  v2x * np.sin(self.Omega) + v2y * np.cos(self.Omega)
        #vz =  v2z
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
