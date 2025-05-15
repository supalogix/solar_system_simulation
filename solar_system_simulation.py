import numpy as np
import matplotlib.pyplot as plt
from einsteinpy.geodesic import Geodesic
import astropy.units as u
import mpmath as mp
import matplotlib.ticker as ticker


# Set high precision (50 decimal places)
mp.dps = 50  

# Constants
G = 6.674e-11 * u.m**3 / (u.kg * u.s**2)  # Gravitational constant
c = 3.0e8 * u.m / u.s  # Speed of light
afs_conversion = 1e18  # Convert seconds to attoseconds
solar_mass = 1.989e30 * u.kg  # Solar mass



# Planet class definition
class Planet:
    def __init__(self, name, a, e, i, Omega, omega, nu, mass):
        self.name = name
        self.a = a  # Semi-major axis in meters
        self.e = e  # Eccentricity
        self.i = np.radians(i)  # Inclination in radians
        self.Omega = np.radians(Omega)  # Longitude of ascending node in radians
        self.omega = np.radians(omega)  # Argument of periapsis in radians
        self.nu = np.radians(nu)  # True anomaly in radians
        self.mass = mass  # Mass in kg

    def position(self, t):
        n = np.sqrt(G.value * solar_mass.value / (self.a**3))  # Mean motion
        M = n * t  # Mean anomaly in radians
        E = M  # First approximation
        for _ in range(10):
            E = M + self.e * np.sin(E)
        nu = 2 * np.arctan2(np.sqrt(1 + self.e) * np.sin(E / 2),
                             np.sqrt(1 - self.e) * np.cos(E / 2))
        r = self.a * (1 - self.e * np.cos(E))
        x_orb = r * np.cos(nu)
        y_orb = r * np.sin(nu)
        return np.array([x_orb, y_orb, 0])



# Satellite class definition
class Satellite:
    def __init__(self, name, distance_from_sun):
        self.name = name
        self.distance_from_sun = distance_from_sun

    def position(self, angle):
        x = self.distance_from_sun * np.cos(angle)
        y = self.distance_from_sun * np.sin(angle)
        return np.array([x, y, 0])



# Function to calculate time dilation considering all planetary contributions
def calculate_time_dilation(planets, position):
    G_mp = mp.mpf(G.to_value(u.m**3 / (u.kg * u.s**2)))  # High-precision G
    c_mp = mp.mpf(c.to_value(u.m / u.s))  # High-precision c

    total_time_dilation = mp.mpf(0)  # Initialize total time dilation

    for planet in planets:
        mass_kg = mp.mpf(planet.mass)
        planet_pos = planet.position(0)  # Get planet's position at t=0
        distance = np.linalg.norm(position - planet_pos)

        if distance == 0:
            continue  # Avoid singularities

        # Compute Schwarzschild time dilation for each planet
        time_dilation = (G_mp * mass_kg) / (distance * c_mp**2)
        total_time_dilation += time_dilation

    return float(total_time_dilation * afs_conversion)  # Convert to attoseconds




# Example parameters for planets
planets = [
    Planet("Mercury", 5.791e10, 0.2056, 7.0, 48.3, 29.1, 0, 3.285e23),
    Planet("Venus", 1.082e11, 0.0067, 3.4, 76.7, 54.9, 0, 4.867e24),
    Planet("Earth", 1.496e11, 0.0167, 0.0, 0.0, 102.9, 0, 5.972e24),
    Planet("Mars", 2.279e11, 0.0935, 1.9, 49.6, 286.5, 0, 6.417e23),
    Planet("Jupiter", 7.785e11, 0.0484, 1.3, 100.5, 273.6, 0, 1.898e27),
    Planet("Saturn", 1.429e12, 0.0565, 2.5, 113.6, 339.6, 0, 5.683e26),
    Planet("Uranus", 2.871e12, 0.0463, 0.8, 74.0, 96.6, 0, 8.681e25),
    Planet("Neptune", 4.495e12, 0.0086, 1.8, 131.8, 276.5, 0, 1.024e26)
]





jupiter = planets[4]
l4_satellite = Satellite("L4 Satellite", jupiter.a)
l5_satellite = Satellite("L5 Satellite", jupiter.a)
l1_satellite = Satellite("L1 Satellite", jupiter.a * 0.99)
l2_satellite = Satellite("L2 Satellite", jupiter.a * 1.01)
l3_satellite = Satellite("L3 Satellite", jupiter.a)



days_in_year = 365 * 24
times = np.arange(days_in_year)
colors = ["red", "green", "purple", "orange", "blue"]



# Dictionary to store proper time offsets (in attoseconds)
proper_time_offset_data = {key: [] for key in ["L1", "L2", "L3", "L4", "L5"]}

# Iterate over time steps
for day in times:
    t = day * 86400  # Convert days to seconds
    jupiter_pos = jupiter.position(t)  # Get Jupiter's position
    angle = np.arctan2(jupiter_pos[1], jupiter_pos[0])  # Compute current orbital angle

    # Compute positions for Lagrange points
    l1_pos = l1_satellite.position(angle)
    l2_pos = l2_satellite.position(angle)
    l3_pos = l3_satellite.position(angle + np.pi)
    l4_pos = l4_satellite.position(angle + 2 * np.pi / 3)
    l5_pos = l5_satellite.position(angle - 2 * np.pi / 3)

    # Store positions in a dictionary
    positions = {"L1": l1_pos, "L2": l2_pos, "L3": l3_pos, "L4": l4_pos, "L5": l5_pos}

    # Compute time dilation at each Lagrange point
    time_dilations = {key: calculate_time_dilation(planets, pos) for key, pos in positions.items()}

    # Select L4 as the reference time dilation
    reference_time_dilation = time_dilations["L4"]

    # Compute proper time offset for each point relative to L4
    for key in time_dilations:
        offset = (time_dilations[key] - reference_time_dilation) * afs_conversion  # Convert to attoseconds
        proper_time_offset_data[key].append(offset)









# Plot time dilation graph
plt.figure(figsize=(10, 5))
for key, color in zip(proper_time_offset_data.keys(), colors):
    plt.plot(times, proper_time_offset_data[key], label=f"{key} Time Dilation", color=color)

plt.xlabel("Days")
plt.ylabel("Time Dilation (attoseconds)")
plt.title("Time Dilation at Lagrange Points")
plt.legend()
plt.grid(True)
plt.show()


#Convert attoseconds to millions of attoseconds
proper_time_offsets_millions = {key: [value / 1e6 for value in values] for key, values in proper_time_offset_data.items()}

for key, color in zip(proper_time_offset_data.keys(), colors):
    plt.figure(figsize=(10, 4))
    plt.plot(times, proper_time_offsets_millions[key], label=key, color=color)

    plt.xlabel("Days")
    plt.ylabel("Proper Time Offset (millions of attoseconds)")  # Make it clear
    plt.title(f"Proper Time Offset: {key} Over a Year")
    plt.legend()
    plt.grid(True)

    # Format the y-axis to show clear million markers
    plt.gca().yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{x:.0f}M"))

    plt.show()
