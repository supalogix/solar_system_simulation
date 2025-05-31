from typing import List, Tuple, Dict
import numpy as np

def generate_lagrange_orbital_params( planet, m_sun):

    a0: float = planet['a']
    mu: float = planet['mass'] / (m_sun + planet['mass'])

    # copy through all the non‐nu elements
    common: Dict[str, float] = {
        "e":     planet['e'],
        "i":     planet['i'],
        "Omega": planet['Omega'],
        "omega": planet['omega'],
    }

    # Define the five Lagrange points by:
    # 1) identifier of the point
    # 2) shift in semi-major axis (meters) relative to the planet’s orbit
    # 3) true anomaly offset 
    specs: List[Tuple[str, float, float]] = [
        ("L1", -a0 * (mu/3)**(1/3),      0.0      ),
        ("L2",  a0 * (mu/3)**(1/3),      0.0      ),
        ("L3",  a0 * (7*mu/12),          np.pi    ),
        ("L4",  0.0,                     np.pi/3  ),
        ("L5",  0.0,                    -np.pi/3  ),
    ]

    L_params: List[Dict[str, float]] = []
    for name, da, dnu in specs:
        aL = a0 + da

        # wrap true anomaly into [0, 2*pi)
        nuL = (planet['nu'] + dnu) % (2 * np.pi)

        L_params.append({
            "name":  name,         # str
            "type":  "SATELLITE",  # str
            "a":      aL,          # float (m)
            **common,              # e, i, Omega, omega 
            "nu":     nuL,         # float (rad)
            "mass":   100.0,       # float (kg)
        })

    return L_params

