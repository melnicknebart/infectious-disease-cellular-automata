# Infectious Disease Cellular Automata Simulator

A Python cellular-automata project for exploring how an infectious disease spreads across a spatial population.

## What this project demonstrates

- Object-oriented Python with `Cell` and `DiseaseMap`
- Stochastic infection, recovery, and mortality
- Reproducible random seeds
- Numerical approximation of a mortality distribution
- Synchronous cellular-automata updates
- NumPy + Matplotlib visualization
- CSV input and history export
- Command-line parameters with `argparse`

## Repository structure

```text
infectious-disease-cellular-automata/
├── README.md
├── simulator.py
├── generate_demo_map.py
├── requirements.txt
├── .gitignore
├── data/
│   └── README.md
├── examples/
└── outputs/
```

## Install

```bash
python -m venv .venv
```

Activate the environment, then:

```bash
pip install -r requirements.txt
```

## Quick start

Create the synthetic public demo map:

```bash
python generate_demo_map.py
```

Run the default simulation:

```bash
python simulator.py
```

The program prints state counts and saves:

- `outputs/final_state.png`
- `outputs/history.csv`

## Try different parameters

```bash
python simulator.py --virality 0.9 --recovery-time 4 --death-mean 5 --death-sd 1.2 --steps 35 --seed 7
```

```bash
python simulator.py --virality 0.25 --recovery-time 2 --death-mean 4 --death-sd 1 --steps 35 --seed 7
```

Using the same seed makes parameter comparisons easier to reproduce.

## Model logic

Each active cell is one of three states:

- `S`: susceptible
- `I`: infected
- `R`: resistant/dead

At each time step:

1. Newly infected cells wait before becoming infectious.
2. Infected cells first check for recovery.
3. Remaining infected cells evaluate a time-dependent mortality probability.
4. Surviving infected cells may infect north/east/south/west susceptible neighbors.
5. State changes are applied synchronously.

## Data

The original project used a course-supplied 150 × 150 New York City coordinate map.

That dataset is not redistributed here. `generate_demo_map.py` instead creates a synthetic 150 × 150 grid so the public portfolio project can run independently.

If you have permission to use the original map locally, place it at:

```text
data/nyc_map.csv
```

and run:

```bash
python simulator.py --map data/nyc_map.csv
```

## Improvements over the original course version

- command-line controls
- reproducible random seeds
- synchronous updates
- input validation
- exported simulation history
- saved figures
- cleaner separation of simulation configuration and model logic
- docstrings and clearer naming
- publicly shareable synthetic demo data

## Future improvements

- outbreak animation
- S/I/R population curves
- automatic parameter sweeps
- heterogeneous population density
- distance-weighted transmission
- unit tests
