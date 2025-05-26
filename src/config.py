# Rather than using abstract units, we anchor each planet’s orbital elements and mass to Earth’s 
# baseline so it’s easy to see how each world differs:
#  - semi-major axis (a) is given in meters but described as a multiple of 1 AU (Earth’s orbit)
#  - eccentricity (e) measures how stretched the orbit is (Earth’s e≈0.0167 is nearly circular)
#  - inclination (i) is the tilt of the orbital plane relative to Earth’s (i=0 by definition)
#  - mass is quoted in kilograms but noted as a multiple of Earth’s mass
# Choosing Earth as the reference makes these comparisons intuitive: we know what 1 AU,
# e≈0.017, i=0° and 1 Earth-mass feel like, so we can immediately grasp how Mercury’s
# 0.39 AU or Jupiter’s 318× Earth-mass stand out.
planets_data = [
    {
        "name": "Mercury",
        "a": 5.791e10,
        "e": 0.2056,
        "i": 7.0,
        "Omega": 48.3,
        "omega": 29.1,
        "nu": 0.0,
        "mass": 3.285e23
    },
    {
        "name": "Venus",
        "a": 1.082e11,
        "e": 0.0067,
        "i": 3.4,
        "Omega": 76.7,
        "omega": 54.9,
        "nu": 0.0,
        "mass": 4.867e24
    },
    {
        "name": "Earth",
        "a": 1.496e11,
        "e": 0.0167,
        "i": 0.0,
        "Omega": 0.0,
        "omega": 102.9,
        "nu": 0.0,
        "mass": 5.972e24
    },
    {
        "name": "Mars",
        "a": 2.279e11,
        "e": 0.0935,
        "i": 1.9,
        "Omega": 49.6,
        "omega": 286.5,
        "nu": 0.0,
        "mass": 6.417e23
    },
    {
        "name": "Jupiter",
        "a": 7.785e11,
        "e": 0.0484,
        "i": 1.3,
        "Omega": 100.5,
        "omega": 273.6,
        "nu": 0.0,
        "mass": 1.898e27
    },
    {
        "name": "Saturn",
        "a": 1.429e12,
        "e": 0.0565,
        "i": 2.5,
        "Omega": 113.6,
        "omega": 339.6,
        "nu": 0.0,
        "mass": 5.683e26
    },
    {
        "name": "Uranus",
        "a": 2.871e12,
        "e": 0.0463,
        "i": 0.8,
        "Omega": 74.0,
        "omega": 96.6,
        "nu": 0.0,
        "mass": 8.681e25
    },
    {
        "name": "Neptune",
        "a": 4.495e12,
        "e": 0.0086,
        "i": 1.8,
        "Omega": 131.8,
        "omega": 276.5,
        "nu": 0.0,
        "mass": 1.024e26
    }
]

colors = ["red", "green", "purple", "orange", "blue"]

#planets = [
#    Planet(
#        name="Mercury",
#        a=5.791e10,    # semi-major axis: 0.39 AU (~0.39 x Earth's 1 AU)
#        e=0.2056,      # eccentricity: 0.2056 (~12.3 x Earth's e=0.0167)
#        i=7.0,         # inclination: 7.0 deg vs Earth's 0 deg
#        Omega=48.3,    # longitude of ascending node: 48.3 deg vs 0 deg
#        omega=29.1,    # argument of periapsis: 29.1 deg vs Earth's 102.9 deg
#        nu=0,          # true anomaly at epoch: 0 deg
#        mass=3.285e23  # mass: 0.055 x Earth mass
#    ),

#    Planet(
#        name="Venus",
#        a=1.082e11,    # semi-major axis: 0.72 AU (~0.72 x Earth's 1 AU)
#        e=0.0067,      # eccentricity: 0.0067 (~0.4 x Earth's e)
#        i=3.4,         # inclination: 3.4 deg vs Earth's 0 deg
#        Omega=76.7,    # longitude of ascending node: 76.7 deg vs 0 deg
#        omega=54.9,    # argument of periapsis: 54.9 deg vs Earth's 102.9 deg
#        nu=0,          # true anomaly at epoch: 0 deg
#        mass=4.867e24  # mass: 0.82 x Earth mass
#    ),

#    Planet(
#        name="Earth",
#        a=1.496e11,    # semi-major axis: 1.00 AU (baseline)
#        e=0.0167,      # eccentricity: 0.0167 (baseline)
#        i=0.0,         # inclination: 0 deg (baseline)
#        Omega=0.0,     # longitude of ascending node: 0 deg (baseline)
#        omega=102.9,   # argument of periapsis: 102.9 deg
#        nu=0,          # true anomaly at epoch: 0 deg
#        mass=5.972e24  # mass: 1.00 x Earth mass (baseline)
#    ),

#    Planet(
#        name="Mars",
#        a=2.279e11,    # semi-major axis: 1.52 AU (~1.52 x Earth's 1 AU)
#        e=0.0935,      # eccentricity: 0.0935 (~5.6 x Earth's e)
#        i=1.9,         # inclination: 1.9 deg vs 0 deg
#        Omega=49.6,    # longitude of ascending node: 49.6 deg vs 0 deg
#        omega=286.5,   # argument of periapsis: 286.5 deg
#        nu=0,          # true anomaly at epoch: 0 deg
#        mass=6.417e23  # mass: 0.11 x Earth mass
#    ),

#    Planet(
#        name="Jupiter",
#        a=7.785e11,    # semi-major axis: 5.20 AU (~5.20 x Earth's 1 AU)
#        e=0.0484,      # eccentricity: 0.0484 (~2.9 x Earth's e)
#        i=1.3,         # inclination: 1.3 deg vs 0 deg
#        Omega=100.5,   # longitude of ascending node: 100.5 deg vs 0 deg
#        omega=273.6,   # argument of periapsis: 273.6 deg
#        nu=0,          # true anomaly at epoch: 0 deg
#        mass=1.898e27  # mass: 318 x Earth mass
#    ),

#    Planet(
#        name="Saturn",
#        a=1.429e12,    # semi-major axis: 9.58 AU (~9.58 x Earth's 1 AU)
#        e=0.0565,      # eccentricity: 0.0565 (~3.4 x Earth's e)
#        i=2.5,         # inclination: 2.5 deg vs 0 deg
#        Omega=113.6,   # longitude of ascending node: 113.6 deg vs 0 deg
#        omega=339.6,   # argument of periapsis: 339.6 deg
#        nu=0,          # true anomaly at epoch: 0 deg
#        mass=5.683e26  # mass: 95 x Earth mass
#    ),

#    Planet(
#        name="Uranus",
#        a=2.871e12,    # semi-major axis: 19.19 AU (~19.19 x Earth's 1 AU)
#        e=0.0463,      # eccentricity: 0.0463 (~2.8 x Earth's e)
#        i=0.8,         # inclination: 0.8 deg vs 0 deg
#        Omega=74.0,    # longitude of ascending node: 74.0 deg vs 0 deg
#        omega=96.6,    # argument of periapsis: 96.6 deg
#        nu=0,          # true anomaly at epoch: 0 deg
#        mass=8.681e25  # mass: 14 x Earth mass
#    ),

#    Planet(
#        name="Neptune",
#        a=4.495e12,    # semi-major axis: 30.07 AU (~30.07 x Earth's 1 AU)
#        e=0.0086,      # eccentricity: 0.0086 (~0.5 x Earth's e)
#        i=1.8,         # inclination: 1.8 deg vs 0 deg
#        Omega=131.8,   # longitude of ascending node: 131.8 deg vs 0 deg
#        omega=276.5,   # argument of periapsis: 276.5 deg
#        nu=0,          # true anomaly at epoch: 0 deg
#        mass=1.024e26  # mass: 17 x Earth mass
#    ),
#]