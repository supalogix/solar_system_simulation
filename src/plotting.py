import numpy as np
import mpmath as mp
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from constants import G_mp, c_mp, afs_conversion, G, solar_mass


def plot_all_offsets(times, data, colors):
    """
    Plot the proper-time offsets for each satellite over the given time series.

    Args:
        times (array-like of float): Sequence of time points (e.g., days).
        data (dict of str -> sequence of float): Mapping from satellite name to its
            offset values at each time point.
        colors (list of str): Color strings for each satellite curve, in the same
            order as data.keys().

    Returns:
        None
    """

    plt.figure(figsize=(10, 5))
    for key, col in zip(data.keys(), colors):
        plt.plot(times, data[key], label=f"{key} Time Dilation", color=col)
    plt.xlabel("Days")
    plt.ylabel("Time Dilation (as)")
    plt.title("Time Dilation at Lagrange Points")
    plt.legend(); plt.grid(True)
    plt.show()


def plot_individual(times, data, colors):
    """
    Plot each satellite’s proper-time offset in a separate subplot.

    Args:
        times (array-like of float): Sequence of time points (e.g., days).
        data (dict of str -> sequence of float): Mapping from satellite name to its
            offset values at each time point.
        colors (list of str): Color strings for each satellite plot, matching data.keys().

    Returns:
        None
    """

    for key, col in zip(data.keys(), colors):
        fig, ax = plt.subplots(figsize=(10,4))
        ax.plot(times, [v/1e6 for v in data[key]], label=key, color=col)
        ax.set(xlabel="Days",
               ylabel="Proper Time Offset (millions of as)",
               title=f"Proper Time Offset: {key} Over a Year")
        ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{x:.0f}M"))
        ax.legend(); ax.grid(True)
        plt.show()

def plot_satellite_contributions_over_orbit(planets, satellite, steps=360):
    """
    Plot each planet's time-dilation contribution at a single satellite as it orbits the sun.

    Parameters:
        planets   : list of Planet objects
        satellite : a Satellite object with .position(angle) and .name
        steps     : number of points around one full orbit (default 360)

    Returns: 
        None
    """

    # Generate 'steps' evenly spaced angles from 0 to 2*\pi radians,
    # so we can pick out points uniformly around the full circular orbit.
    angles = np.linspace(0, 2 * np.pi, steps)

    # Prepare storage for each planet's contributions
    contributions = {p.name: [] for p in planets}

    # Angular speed for a circular orbit:
    # In a perfectly circular path, the Sun’s gravity provides the centripetal force,
    # so we have
    #   $$\omega_{\text{sat}} = \sqrt{\frac{G\,M}{r^3}}$$
    # where \(M\) is the Sun’s mass and \(r\) is the orbital radius.
    r = satellite.distance_from_sun
    omega_sat = np.sqrt(G.value * solar_mass.value / (r**3))


    for angle in angles:
        # Compute the time at which the satellite reaches this orbital angle:
        #   $$t = \frac{\theta}{\omega_{\mathrm{sat}}}$$
        # where θ is the current angle (radians) and ω_sat is the angular speed (rad/s).
        t = float(angle / omega_sat)

        # Get the satellite’s 3D position vector at angle θ along its orbit.
        pos = satellite.position(angle)

        # Loop over each planet and convert its mass to a high-precision mpf:
        # time dilation effects are vanishingly small, so using arbitrary-precision
        # arithmetic ensures we don’t lose those tiny contributions to rounding.
        for planet in planets:

            # High-precision mass
            m_mp = mp.mpf(planet.mass)

            # Find how far the planet is from the satellite right now:
            planet_pos = planet.position(t)
            dist = np.linalg.norm(pos - planet_pos)
            if dist == 0:
                contributions[planet.name].append(0.0)
                continue

            # Gravitational time dilation:
            # A mass "dips" spacetime, so clocks nearer the mass run slower.
            # In the weak-field limit, the tiny fractional slowdown per second is
            #   $$\frac{G\,m}{r\,c^2}$$
            # where r is the planet–satellite distance. We use that approximation here.
            td_grav = (G_mp * m_mp) / (dist * c_mp**2)

            # Kinematic slowdown:
            # Start from momentum $\mathbf{p} = m\,\mathbf{v}$,
            # compute $v^2 = \mathbf{p}\cdot\mathbf{p} \,/\, m^2$,
            # and apply the weak‐speed approximation
            #   $\displaystyle \Delta t - \tau \approx \tfrac12\,\frac{v^2}{c^2}\,$.
            p_vec = planet.momentum(t)
            p_sq  = mp.mpf(np.dot(p_vec, p_vec))
            v_sq  = p_sq / (m_mp**2)
            td_kin = v_sq / (2 * c_mp**2)

            # Total contribution:
            # Combine the gravitational term $\displaystyle \frac{G\,m}{r\,c^2}$
            # with the kinematic term $\tfrac12\,\frac{v^2}{c^2}$,
            # then convert the sum into attoseconds.
            total_td = td_grav + td_kin
            contributions[planet.name].append(float(total_td * afs_conversion))

    # Plot each planet's contribution over the orbit
    plt.figure(figsize=(10, 6))
    for name, vals in contributions.items():
        plt.plot(angles, vals, label=name)

    plt.xlabel('Orbit Angle (rad)')
    plt.ylabel('Time Dilation (attoseconds)')
    plt.title(f"Planet Contributions to Time Dilation of {satellite.name}")
    plt.legend(); plt.grid(True)
    plt.tight_layout()
    plt.show()