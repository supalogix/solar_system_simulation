import astropy.units as u
import mpmath as mp

# Universal gravitational constant (G): determines the strength of gravity
# in Newton’s law of gravitation (attracting masses over a distance).
G = 6.674e-11 * u.m**3 / (u.kg * u.s**2)

# Speed of light in vacuum (c): the ultimate speed limit for any signal
# or object; appears in relativity formulas for time dilation.
c = 3.0e8 * u.m / u.s

# Attosecond conversion factor: to turn seconds into 10^{-18} seconds,
# making tiny relativistic offsets easier to work with.
afs_conversion = 1e18

# Mass of the Sun: central mass in the solar system, used for computing
# orbital motions and gravitational fields of planets.
solar_mass  = 1.989e30 * u.kg

# Convert constants to high-precision mpmath floats:
# We use mpf to preserve many significant digits, since time dilation
# effects are extremely small and require high precision summation.
G_mp = mp.mpf(G.to_value(u.m**3 / (u.kg * u.s**2)))  # High-precision G
c_mp = mp.mpf(c.to_value(u.m / u.s))                 # High-precision c
