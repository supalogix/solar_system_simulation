# IO Documentation

This README describes how to use the `planets_data` array as the input to our solar‐system simulation program to output time dilation data between lagrange points for the duration of the simulation. 

It explains the structure of each planet’s data, the meaning of every field, and how the program consumes this input to construct `Planet` objects for simulation, and the structure of the output.

---

## Program Input Overview

The simulation expects a Python list of dictionaries called `planets_data`. Each dictionary represents one solar‐system body (in this case, a planet) and provides its orbital elements plus physical mass. The program will:

1. Read `planets_data` at startup.
2. Compute the Sun’s barycentric semi‐major axis from the combined planetary masses and semi‐major axes.
3. Insert the Sun (as a special “planet” entry) into the data.
4. Instantiate a `Planet` object for each entry (Sun plus all planets).
5. Proceed with orbital propagation (position, velocity, momentum, etc.) for planets and any additional satellites/Lagrange points.

### Data Structure

Each entry in `planets_data` can be viewed as JSON‐style dictionary with the following keys:

| Key     | Type   | Units / Format            | Description                                                                            |
| ------- | ------ | ------------------------- | -------------------------------------------------------------------------------------- |
| `name`  | string | —                         | Identifier for the body (e.g. `"Mercury"`, `"Earth"`).                                 |
| `type`  | string | —                         | Body type, currently always `"PLANET"`. (Future use: `"SATELLITE"`, etc.)              |
| `a`     | float  | meters                    | Semi‐major axis of the orbit (distance from central body to orbit’s center).           |
| `e`     | float  | dimensionless (0 ≤ e < 1) | Orbital eccentricity. `e = 0` means circular; `0 < e < 1` means elliptical.            |
| `i`     | float  | degrees                   | Inclination relative to the reference plane (e.g., ecliptic).                          |
| `Omega` | float  | degrees                   | Longitude of ascending node (angle from reference direction to ascending node).        |
| `omega` | float  | degrees                   | Argument of periapsis (angle from ascending node to periapsis).                        |
| `nu`    | float  | degrees                   | True anomaly at epoch (initial angular position measured from periapsis).              |
| `mass`  | float  | kilograms                 | Physical mass of the body (used to compute barycenter and gravitational interactions). |

**Note:** All angular quantities are specified in **degrees** in the input. Internally, the program converts them to **radians** when constructing `Planet` objects.


#### Example (Mercury)

```jsonc
{
  "name": "Mercury",
  "type": "PLANET",
  "a": 5.791e10,           // 0.387 AU in meters
  "e": 0.2056,             // Mercury’s eccentricity
  "i": 7.0,                // 7 degree inclination
  "Omega": 48.3,           // 48.3 degree node longitude
  "omega": 29.1,           // 29.1 degree argument of periapsis
  "nu": 0.0,               // starting at true anomaly 0°
  "mass": 3.285e23         // Mercury’s mass in kg
}
```

---

## Program Output Overview

The output contains the following top‐level fields:

* **time** (`number`): Simulation time in seconds.

* **planets** (`object`): Dictionary mapping planet names (e.g. “Sun”, “Earth”, etc.) to their state at time `time`.

  * Each planet has:

    * **position** (`object`): Cartesian coordinates in meters (`px`, `py`, `pz`).
    * **velocity** (`object`): Cartesian velocity components in m/s (`vx`, `vy`, `vz`).
    * **momentum** (`object`): Cartesian momentum components in kg·m/s (`mx`, `my`, `mz`).

* **satellites** (`object`): Dictionary mapping Trojan/Lagrange‐point names (e.g. “L1”, “L2”, etc.) to their state.

  * Each satellite entry has:

    * **position** (`object`): Cartesian coordinates in meters (`px`, `py`, `pz`).
    * **velocity** (`object`): Cartesian velocity in m/s (`vx`, `vy`, `vz`).
    * **momentum** (`object`): Cartesian momentum in kg·m/s (`mx`, `my`, `mz`).
    * **time\_dilation\_per\_planet** (`object`): Dictionary of time‐dilation factors (in seconds) contributed by each planet. Keys are planet names; values are numeric.
    * **time\_dilation\_total** (`number`): Sum of all per‐planet time‐dilations for this satellite (in seconds).

* **day** (`integer`): Current simulation day count (e.g. integer day index).

* **time\_dilations** (`object`): Nested dictionary that records relative time‐dilation differences among satellites. Top‐level keys are satellite names (e.g. “L4”, “L5”). Each sub‐dictionary maps every satellite name to a signed numeric difference in seconds.

### Example Output

```jsonc
{
  "time": 0,
  "planets": {
    "Sun": {
      "position": { "px": 1509076438.2524517, "py": 0.0, "pz": 0.0 },
      "velocity": { "vx": 0.0, "vy": 296588.8865114812, "vz": 0.0 },
      "momentum": { "mx": 0.0, "my": 5.897670008280804e+35, "mz": 0.0 }
    },
    "Mercury": {
      "position": { "px": 10159911459.266905, "py": 44784847603.54687, "pz": 2726610714.156952 },
      "velocity": { "vx": -57274.20423034981, "vy": 12610.868087461728, "vz": 6280.7053262140225 },
      "momentum": { "mx": -1.881457608966991e+28, "my": 4.1426701667311773e+27, "mz": 2.0632116996613064e+27 }
    },
    "Venus": {
      "position": { "px": -94400428462.4021, "py": 51096231964.440384, "pz": 6156351101.110938 },
      "velocity": { "vx": -16819.463747900652, "vy": -30971.27241364595, "vz": 549.1599615501833 },
      "momentum": { "mx": -8.186033006103248e+28, "my": -1.5073718283721485e+29, "mz": 2.6727615328647425e+27 }
    },
    "Earth": {
      "position": { "px": -141763036119.5926, "py": 43612004483.333786, "pz": 0.0 },
      "velocity": { "vx": -9245.138120324325, "vy": -28586.42727665787, "vz": 0.0 },
      "momentum": { "mx": -5.521196485457688e+28, "my": -1.707181436962008e+29, "mz": 0.0 }
    },
    "Mars": {
      "position": { "px": 91551443362.14935, "py": 206514092191.48663, "pz": 2127277244.8948302 },
      "velocity": { "vx": -21233.8934403583, "vy": 11884.562557296389, "vz": 791.9504601557471 },
      "momentum": { "mx": -1.3625789420677923e+28, "my": 7.626323793017093e+27, "mz": 5.08194610281943e+26 }
    },
    "Jupiter": {
      "position": { "px": -553791978100.0188, "py": 571604681782.7864, "pz": 9992984829.955616 },
      "velocity": { "vx": -9539.776488025473, "vy": -8483.782082842274, "vz": 247.9474077158403 },
      "momentum": { "mx": -1.8106495774272346e+31, "my": -1.6102218393234635e+31, "mz": 4.706041798446649e+29 }
    },
    "Saturn": {
      "position": { "px": 84734859227.70465, "py": -1507183936470.914, "pz": 22954817750.32677 },
      "velocity": { "vx": 9086.511063296744, "vy": 505.1781797564645, "vz": -372.3749117873166 },
      "momentum": { "mx": 5.16386423727154e+30, "my": 2.8709275955559878e+29, "mz": -2.11620662368732e+29 }
    },
    "Uranus": {
      "position": { "px": 2556377694907.7095, "py": 1512005473946.877, "pz": -28493762056.515335 },
      "velocity": { "vx": -3516.092903563602, "vy": 5548.05985189801, "vz": 68.54895159433241 },
      "momentum": { "mx": -3.0523202495835628e+29, "my": 4.816270757432663e+29, "mz": 5.950734487903996e+27 }
    },
    "Neptune": {
      "position": { "px": 4381457273329.3457, "py": -908434514502.6641, "pz": -83618070296.43863 },
      "velocity": { "vx": 1065.5523173844413, "vy": 5351.268013149529, "vz": -137.05426942553254 },
      "momentum": { "mx": 1.0911255730016677e+29, "my": 5.479698445465117e+29, "mz": -1.4034357189174532e+28 }
    }
  },
  "satellites": {
    "L1": {
      "position": {
        "px": 667440670899.889,
        "py": 175355813906.02612,
        "pz": -15617905101.177258
      },
      "velocity": {
        "vx": -3600.715845599989,
        "vy": 13735.067848989218,
        "vz": 23.541995051746635
      },
      "momentum": {
        "mx": -360071.5845599989,
        "my": 1373506.7848989218,
        "mz": 2354.1995051746635
      },
      "time_dilation_per_planet": {
        "Sun": 490835035326.70087,
        "Mercury": 19326754383.012787,
        "Venus": 6902313853.488105,
        "Earth": 5014763196.803622,
        "Mars": 3293046621.0684485,
        "Jupiter": 906892531.9195603,
        "Saturn": 461117463.7291687,
        "Uranus": 239742132.0615704,
        "Neptune": 165521040.43982935
      },
      "time_dilation_total": 527145186549.22394
    },
    "L2": {
      "position": {
        "px": 765205084197.1158,
        "py": 201041330255.61597,
        "pz": -17905562110.585636
      },
      "velocity": {
        "vx": -3362.8404291571273,
        "vy": 12827.6830053777,
        "vz": 21.98673156610581
      },
      "momentum": {
        "mx": -336284.04291571275,
        "my": 1282768.30053777,
        "mz": 2198.673156610581
      },
      "time_dilation_per_planet": {
        "Sun": 490561020202.0841,
        "Mercury": 19326754335.456257,
        "Venus": 6902313314.311848,
        "Earth": 5014762606.13526,
        "Mars": 3293046502.4130464,
        "Jupiter": 906823602.514624,
        "Saturn": 461109973.53922796,
        "Uranus": 239743314.7536977,
        "Neptune": 165521488.21403912
      },
      "time_dilation_total": 526871095339.4221
    },
    "L3": {
      "position": {
        "px": 705391491927.355,
        "py": 227331573929.55014,
        "pz": -16679654313.196138
      },
      "velocity": {
        "vx": -4162.317834393062,
        "vy": 13053.499277659448,
        "vz": 38.89174726003649
      },
      "momentum": {
        "mx": -416231.7834393062,
        "my": 1305349.9277659447,
        "mz": 3889.174726003649
      },
      "time_dilation_per_planet": {
        "Sun": 490687293351.0584,
        "Mercury": 19326754358.410007,
        "Venus": 6902313584.905869,
        "Earth": 5014762904.032745,
        "Mars": 3293046570.7710238,
        "Jupiter": 906874478.0885895,
        "Saturn": 461109548.43899083,
        "Uranus": 239742884.13841873,
        "Neptune": 165521149.6270339
      },
      "time_dilation_total": 526997418829.47107
    },
    "L4": {
      "position": {
        "px": 712775331036.1454,
        "py": 201265527491.95685,
        "pz": -16736615006.760975
      },
      "velocity": {
        "vx": -3706.1623344030145,
        "vy": 13195.334123399998,
        "vz": 28.126911261921737
      },
      "momentum": {
        "mx": -370616.23344030144,
        "my": 1319533.4123399998,
        "mz": 2812.6911261921737
      },
      "time_dilation_per_planet": {
        "Sun": 490688606623.66895,
        "Mercury": 19326754357.924034,
        "Venus": 6902313573.92141,
        "Earth": 5014762891.613674,
        "Mars": 3293046561.985054,
        "Jupiter": 906862881.4993763,
        "Saturn": 461112309.530122,
        "Uranus": 239742771.2112363,
        "Neptune": 165521225.1985109
      },
      "time_dilation_total": 526998723196.55237
    },
    "L5": {
      "position": {
        "px": 696323162102.9601,
        "py": 252861045220.71075,
        "pz": -16582888275.128662
      },
      "velocity": {
        "vx": -4615.431051867221,
        "vy": 12902.297574148104,
        "vz": 49.62743693975065
      },
      "momentum": {
        "mx": -461543.1051867221,
        "my": 1290229.7574148104,
        "mz": 4962.743693975065
      },
      "time_dilation_per_planet": {
        "Sun": 490688072942.17523,
        "Mercury": 19326754359.251247,
        "Venus": 6902313600.666928,
        "Earth": 5014762921.97501,
        "Mars": 3293046580.530011,
        "Jupiter": 906887248.9013387,
        "Saturn": 461106964.785198,
        "Uranus": 239742972.58735296,
        "Neptune": 165521066.5930591
      },
      "time_dilation_total": 526998208657.4654
    }
  },
  "day": 0,
  "time_dilations": {
    "L4": {
      "L1": 146463352.67156982,
      "L2": -127627857.13024902,
      "L3": -1304367.0812988281,
      "L4": 0.0,
      "L5": -514539.08697509766
    },
    "L5": {
      "L1": 146977891.75854492,
      "L2": -127113318.04327393,
      "L3": -789827.9943237305,
      "L4": 514539.08697509766,
      "L5": 0.0
    }
  }
}
```

