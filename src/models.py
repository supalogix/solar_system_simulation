import numpy as np
from constants import G, solar_mass

# Planet class definition with velocity and momentum
class Planet:
    def __init__(self, name, a, e, i, Omega, omega, nu, mass):
        """
        Create a Planet with its orbital elements and physical mass.

        Args:
            name (str): Identifier for the planet.
            a (float): Semi-major axis of the orbit (meters).
            e (float): Orbital eccentricity (0=circle, <1=ellipse).
            i (float): Inclination relative to the reference plane (radians).
            Omega (float): Longitude of the ascending node (radians).
            omega (float): Argument of periapsis (radians).
            nu (float): True anomaly at epoch (radians), initial angular position.
            mass (float): Mass of the planet (kilograms).
        """

        self.name  = name
        self.a     = a  
        self.e     = e  
        self.i     = np.radians(i)
        self.Omega = np.radians(Omega)
        self.omega = np.radians(omega)
        self.nu    = np.radians(nu)
        self.mass  = mass  

    def position(self, t):
        """
        Return the planet’s instantaneous position vector at a given time.

        Args:
            t (float): Time in seconds at which to evaluate the orbital motion.

        Returns:
            numpy.ndarray: 3D position vector (m), computed from the orbit’s
                        semi-major axis, eccentric anomaly, and true anomaly.
        """

        # Kepler’s equation for eccentric anomaly E:
        # We start with the "mean anomaly" M, which grows uniformly in time
        # (as if the planet moved in a circle). But real orbits are ellipses,
        # so we need the "eccentric anomaly" E, which is an auxiliary angle that 
        # lets us compute the true position on the ellipse. Kepler’s equation,
        #     M = E - e * sin(E),
        # is transcendental (i.e. can’t be solved in closed form), so we use
        # simple iteration: initialize E = M and repeatedly replace
        #     M + e * sin(E) -> E
        # until it converges. The resulting E then feeds into the formulas
        # for the planet’s x/y coordinates on its orbit.
        n = np.sqrt(G.value * solar_mass.value / (self.a**3))
        M = n * t
        E = M
        for _ in range(10):
            E = M + self.e * np.sin(E)

        # True anomaly \nu: the actual angular position of the planet along its elliptical path,
        # converted from the auxiliary eccentric anomaly so we know exactly where it is in orbit.
        nu = 2 * np.arctan2(
            np.sqrt(1 + self.e) * np.sin(E / 2),
            np.sqrt(1 - self.e) * np.cos(E / 2)
        )

        # Orbital radius r: the planet’s distance from the central body at this point in its orbit,
        # determined by the ellipse’s size (semi-major axis) and shape (eccentricity) for the current anomaly.
        r = self.a * (1 - self.e * np.cos(E))

        x = r * np.cos(nu)
        y = r * np.sin(nu)
        return np.array([x, y, 0])

    def velocity(self, t):
        """
        Return the planet’s instantaneous velocity vector at a given time.

        Args:
            t (float): Time in seconds at which to evaluate the orbital motion.

        Returns:
            numpy.ndarray: 3D velocity vector (m/s), combining radial and tangential
                        components computed from the orbital elements.
        """

        # Compute mean motion (n): the planet’s average angular speed if it circled
        # at constant rate in a perfect circle the size of its orbit. This comes from
        # Kepler’s third law relating orbital period to the semi-major axis.
        n = np.sqrt(G.value * solar_mass.value / (self.a**3))

        # Mean anomaly (M): the "clock angle" that ticks forward uniformly in time,
        # as if the planet moved at that average rate in a circle.
        M = n * t

        # Eccentric anomaly (E) iteration: we start with the mean anomaly and
        # repeatedly correct it by adding the orbital eccentricity’s sine term.
        # This simple loop converges on the auxiliary angle that accounts for
        # the ellipse’s shape, bridging the gap between uniform circular motion
        # and the true elliptical path.
        E = M
        for _ in range(10):
            E = M + self.e * np.sin(E)

        # Orbital radius r: the planet’s current distance from the focus, 
        # given by how the ellipse "stretches" at anomaly E
        r = self.a * (1 - self.e * np.cos(E))

        # Rate of change of eccentric anomaly dE_dt: 
        # converts the uniform mean motion into the actual speed at which E advances 
        # along the ellipse
        dE_dt = n / (1 - self.e * np.cos(E))

        # Radial velocity r_dot: 
        # how fast the planet moves closer to or farther from the focus, 
        # coming from the ellipse’s shape (eccentricity) and the speed of E
        r_dot = self.a * self.e * np.sin(E) * dE_dt

        # True anomaly: the planet’s actual angular position along its ellipse,
        # converted from the auxiliary angle E so we know exactly where it sits in its orbit.
        nu = 2 * np.arctan2(
            np.sqrt(1 + self.e) * np.sin(E / 2),
            np.sqrt(1 - self.e) * np.cos(E / 2)
        )

        # Angular rate dnu_dt: how fast that true anomaly is changing at this moment,
        # combining the orbit’s eccentricity and the speed of E to give the instantaneous
        # angular velocity around the focus.
        dnu_dt = (np.sqrt(1 - self.e**2) * dE_dt) / (1 - self.e * np.cos(E))

        # Radial velocity: how fast the planet moves closer to or farther from the focus
        v_rad = r_dot

        # Tangential velocity: how fast the planet sweeps around the focus perpendicular to the radius
        v_tan = r * dnu_dt

        # Decompose into x/y components in the orbital plane:
        #   - vx uses the current angle to project radial motion along x and subtract tangential motion
        #   - vy projects radial motion along y and adds tangential motion
        vx = v_rad * np.cos(nu) - v_tan * np.sin(nu)
        vy = v_rad * np.sin(nu) + v_tan * np.cos(nu)

        # Return the full 3D velocity vector (zero out-of-plane component)
        return np.array([vx, vy, 0])


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



class Satellite:
    def __init__(self, name, distance_from_sun, mass=1.0):
        """
        Initialize a Satellite with its name, orbital distance, and mass.

        Args:
            name (str): Identifier for the satellite.
            distance_from_sun (float): Distance from the Sun (meters).
            mass (float, optional): Mass of the satellite (kilograms). Defaults to 1.0.
        """

        self.name               = name
        self.distance_from_sun  = distance_from_sun
        self.mass               = mass

    def position(self, angle):
        """
        Return the satellite's position vector at a given orbital angle.

        Args:
            angle (float): Angle in radians from the reference direction within the orbital plane.

        Returns:
            numpy.ndarray: 3D position vector [x, y, z] in meters, computed by placing the
                        satellite at its constant distance from the sun and rotating that
                        vector by the specified angle around the sun.
        """

        # Compute the satellite’s position in its orbital plane:
        # - Multiply the constant radius by the cosine of the angle to get the x-coordinate.
        # - Multiply the radius by the sine of the angle to get the y-coordinate.
        # This places the satellite on a circle of radius distance_from_sun at the given angle.
        x = self.distance_from_sun * np.cos(angle)
        y = self.distance_from_sun * np.sin(angle)
        return np.array([x, y, 0])


    def velocity(self, angle):
        """
        Return the satellite’s instantaneous velocity vector at a given orbital angle.

        Args:
            angle (float): Angle in radians from the reference direction within the orbital plane.

        Returns:
            numpy.ndarray: 3D velocity vector [vx, vy, vz] in m/s, representing the
                        tangential speed around the Sun at the satellite’s constant
                        orbital radius, perpendicular to the radial direction.
        """

        # Compute the satellite’s angular speed for a circular orbit:
        # gravity provides the centripetal force
        omega = np.sqrt(G.value * solar_mass.value / (self.distance_from_sun**3))

        # Convert that angular speed into linear velocity at this angle:
        # - vx and vy are the tangential components of the orbital velocity,
        #   perpendicular to the radius vector.
        # The minus sign on vx ensures the velocity is orthogonal and follows
        # the right-hand rule around the Sun.
        vx = -self.distance_from_sun * omega * np.sin(angle)
        vy =  self.distance_from_sun * omega * np.cos(angle)

        # Return the 3D velocity vector (zero out-of-plane component)
        return np.array([vx, vy, 0])

    def momentum(self, angle):
        """
        Compute the satellite’s linear momentum vector at a given orbital angle.

        Args:
            angle (float): Angle in radians from the reference direction within the orbital plane.

        Returns:
            numpy.ndarray: Momentum vector (kg * m/s), given by mass * velocity(angle).
        """

        return self.mass * self.velocity(angle)
