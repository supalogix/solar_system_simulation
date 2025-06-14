import numpy as np

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