
# Exercise 2: orbiting asteroid simulation

**Due:** September 8, before class.

**Collaboration:** you may work in teams of two.  Your team may consult online
resources, others in the class, and the instructor and TA, but your work must
be your own.  Copying other students' code will be considered an honor code
violation.

## Step 1: Get logged into `coc-ice.pace.gatech.edu`

Follow the steps in the [orientation slides](http://docs.pace.gatech.edu/training/img/PACE_orientation_Sep2021.pdf).

## Step 2: Getting the code

This project is hosted in a repository on <https://gitlab.gatech.edu/cse6230fa21/>.

Your team should clone the repository into a local copy that you will use to develop.

## About the code

The file `elements.txt` contains the [orbital elements](https://en.wikipedia.org/wiki/Orbital_elements)
of 585,962 asteroids from <https://ssd.jpl.nasa.gov/?sb_elem>:

```
$ head elements.txt
 Num   Name              Epoch      a          e        i         w        Node        M         H     G   Ref
------ ----------------- ----- ---------- ---------- --------- --------- --------- ----------- ----- ----- ----------
     1 Ceres             59396  2.7656553 0.07839202  10.58820  73.73827  80.26764 247.5499723  3.53  0.12 JPL 48
     2 Pallas            59396  2.7738148 0.22975842  34.89778 310.43686 172.92018 229.2297176  4.22  0.11 JPL 46
     3 Juno              59396  2.6681427 0.25696659  12.99149 247.99924 169.85223 215.0926965  5.28  0.32 JPL 121
     4 Vesta             59396  2.3616594 0.08835130   7.14154 151.01560 103.80606 311.6920605  3.31  0.32 JPL 36
     5 Astraea           59396  2.5739090 0.19062063   5.36759 358.65915 141.57082 112.3804292  6.99  0.15 JPL 121
     6 Hebe              59396  2.4254273 0.20324026  14.73993 239.60615 138.64186 294.2096859  5.65  0.24 JPL 98
     7 Iris              59396  2.3861991 0.22940210   5.51808 145.31202 259.52652 353.1265064  5.60  0.15 JPL 117
     8 Flora             59396  2.2017060 0.15596202   5.88901 285.53203 110.87603  74.7425269  6.54  0.28 JPL 126
```

The program `./simulation PRECISION N_STEPS [N_ASTEROIDS]` simulates the orbits of `N_ASTEROIDS` for 1 year (if `N_ASTEROIDS` is larger that 585,962, some are duplicated).  It does so using `N_STEPS` time steps, each of length `1.0 / N_STEPS`.

Each asteroid is represented by the six components of its state vector:

- its position `x`, `y`, `z`
- its velocity `u`, `v`, `w`

A time step first updates each asteroids velocity from the effect of the Sun's gravity, then
uses the velocity to change the asteroids position.

The `PRECISION` is the type of floating point number used for the time step calculations: `float` (32-bit) or `double` (64-bit).

The ideal simulation should be fast and accurate, but these goals are in opposition:

- `float` arithmetic is faster than `double` arithmetic, but is less accurate
- fewer time steps makes the simulation faster, but bigger errors are made each time step

## Step 3: Parallelize the code with OpenMP

Use OpenMP to parallelize the simulation.  At a minimum, you should parallelize:

- The functions `elems_to_state_arrays()` and `compute_accuracy()` in `simulation.cc`.
- The functions `timestep()` and `position_timestep()` in `timestep.cc`.

You are welcome to make additional parallelization changes.

Things you are not allowed to change:

- You cannot change the math in the timestep.
- All asteroids must be updated at the same time: you cannot simulate them individually or in chunks.

## Step 4: Figure 1, asteroid time steps per second
