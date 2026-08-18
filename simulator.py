import argparse
import csv
import math
import random
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def normpdf(x, mean, sd):
    if sd <= 0:
        raise ValueError("Standard deviation must be greater than 0.")
    var = float(sd) ** 2
    denom = (2 * math.pi * var) ** 0.5
    num = math.exp(-(float(x) - float(mean)) ** 2 / (2 * var))
    return num / denom


def pdeath(x, mean, sd, step=0.01):
    """Approximate death probability during infection time step x."""
    start = x - 0.5
    end = x + 0.5
    integral = 0.0

    while start <= end:
        integral += step * (
            normpdf(start, mean, sd) + normpdf(start + step, mean, sd)
        ) / 2
        start += step

    return min(max(integral, 0.0), 1.0)


@dataclass
class SimulationConfig:
    virality: float = 0.7
    recovery_time: int = 2
    death_mean: float = 3.0
    death_sd: float = 1.0

    def validate(self):
        if not 0 <= self.virality <= 1:
            raise ValueError("virality must be between 0 and 1")
        if self.recovery_time < 1:
            raise ValueError("recovery_time must be at least 1")
        if self.death_sd <= 0:
            raise ValueError("death_sd must be greater than 0")


@dataclass
class Cell:
    x: int
    y: int
    state: str = "S"
    infection_time: int = 0

    def infect(self):
        self.state = "I"
        self.infection_time = 0


class DiseaseMap:
    def __init__(self, width=150, height=150, config=None, rng=None):
        self.width = width
        self.height = height
        self.config = config or SimulationConfig()
        self.config.validate()
        self.rng = rng or random.Random()
        self.cells = {}

    def add_cell(self, cell):
        if not (0 <= cell.x < self.height and 0 <= cell.y < self.width):
            raise ValueError(f"Cell {(cell.x, cell.y)} is outside the map.")
        self.cells[(cell.x, cell.y)] = cell

    def adjacent_cells(self, x, y):
        positions = [
            (x, y - 1),
            (x + 1, y),
            (x, y + 1),
            (x - 1, y),
        ]
        return [self.cells[p] for p in positions if p in self.cells]

    def infect_starting_cell(self, x, y):
        if (x, y) not in self.cells:
            raise ValueError(f"Starting cell {(x, y)} is not on the map.")
        self.cells[(x, y)].infect()

    def time_step(self):
        """Advance the simulation synchronously by one time step."""
        infected = [c for c in self.cells.values() if c.state == "I"]
        susceptible = {
            coord for coord, c in self.cells.items() if c.state == "S"
        }

        recoveries = set()
        deaths = set()
        survivors = set()
        new_infections = set()

        for cell in infected:
            coord = (cell.x, cell.y)

            if cell.infection_time < 1:
                survivors.add(coord)
                continue

            # Required model order: recover -> die -> infect neighbors.
            if cell.infection_time >= self.config.recovery_time:
                recoveries.add(coord)
                continue

            death_probability = pdeath(
                cell.infection_time,
                self.config.death_mean,
                self.config.death_sd,
            )

            if self.rng.random() <= death_probability:
                deaths.add(coord)
                continue

            survivors.add(coord)

            for neighbor in self.adjacent_cells(cell.x, cell.y):
                ncoord = (neighbor.x, neighbor.y)
                if (
                    ncoord in susceptible
                    and self.rng.random() <= self.config.virality
                ):
                    new_infections.add(ncoord)

        for coord in recoveries:
            self.cells[coord].state = "S"
            self.cells[coord].infection_time = 0

        for coord in deaths:
            self.cells[coord].state = "R"
            self.cells[coord].infection_time = 0

        for coord in survivors:
            self.cells[coord].infection_time += 1

        for coord in new_infections:
            self.cells[coord].infect()

    def count_states(self):
        s = sum(c.state == "S" for c in self.cells.values())
        i = sum(c.state == "I" for c in self.cells.values())
        r = sum(c.state == "R" for c in self.cells.values())
        return s, i, r

    def to_image(self):
        image = np.zeros((self.height, self.width, 3), dtype=float)
        colors = {
            "S": [0.0, 1.0, 0.0],
            "I": [1.0, 0.0, 0.0],
            "R": [0.5, 0.5, 0.5],
        }

        for cell in self.cells.values():
            image[cell.x, cell.y] = colors[cell.state]

        return image

    def display(self, title="Disease Simulation", save_path=None, show=True):
        plt.figure(figsize=(7, 7))
        plt.imshow(self.to_image())
        plt.title(title)
        plt.axis("off")
        plt.tight_layout()

        if save_path:
            save_path = Path(save_path)
            save_path.parent.mkdir(parents=True, exist_ok=True)
            plt.savefig(save_path, dpi=200, bbox_inches="tight")

        if show:
            plt.show()
        else:
            plt.close()


def read_map(filename, config, rng, width=150, height=150):
    path = Path(filename)
    if not path.exists():
        raise FileNotFoundError(
            f"{path} not found. Run `python generate_demo_map.py` first "
            "or provide your own permitted map CSV."
        )

    disease_map = DiseaseMap(
        width=width,
        height=height,
        config=config,
        rng=rng,
    )

    with path.open("r", encoding="utf-8-sig", newline="") as f:
        reader = csv.reader(f)
        for line_number, row in enumerate(reader, start=1):
            if not row:
                continue
            if len(row) < 2:
                raise ValueError(f"Expected x,y on line {line_number}: {row}")

            x = int(row[0].strip())
            y = int(row[1].strip())
            disease_map.add_cell(Cell(x, y))

    return disease_map


def write_history(history, filename):
    path = Path(filename)
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["step", "susceptible", "infected", "resistant"])
        writer.writerows(history)


def run_simulation(args):
    config = SimulationConfig(
        virality=args.virality,
        recovery_time=args.recovery_time,
        death_mean=args.death_mean,
        death_sd=args.death_sd,
    )
    config.validate()

    rng = random.Random(args.seed)
    disease_map = read_map(
        args.map,
        config=config,
        rng=rng,
        width=args.width,
        height=args.height,
    )
    disease_map.infect_starting_cell(args.start_x, args.start_y)

    history = []
    s, i, r = disease_map.count_states()
    history.append((0, s, i, r))
    print(f"Step 0 | S: {s} | I: {i} | R: {r}")

    for step in range(1, args.steps + 1):
        disease_map.time_step()
        s, i, r = disease_map.count_states()
        history.append((step, s, i, r))
        print(f"Step {step} | S: {s} | I: {i} | R: {r}")

        if i == 0:
            print("No infected cells remain; ending early.")
            break

    write_history(history, args.history)

    disease_map.display(
        title=(
            f"Final State | virality={config.virality}, "
            f"recovery={config.recovery_time}, seed={args.seed}"
        ),
        save_path=args.output,
        show=not args.no_show,
    )


def build_parser():
    parser = argparse.ArgumentParser(
        description="Cellular-automata infectious-disease simulator."
    )

    parser.add_argument("--map", default="data/demo_map.csv")
    parser.add_argument("--steps", type=int, default=20)
    parser.add_argument("--start-x", type=int, default=75)
    parser.add_argument("--start-y", type=int, default=75)
    parser.add_argument("--width", type=int, default=150)
    parser.add_argument("--height", type=int, default=150)

    parser.add_argument("--virality", type=float, default=0.7)
    parser.add_argument("--recovery-time", type=int, default=2)
    parser.add_argument("--death-mean", type=float, default=3.0)
    parser.add_argument("--death-sd", type=float, default=1.0)
    parser.add_argument("--seed", type=int, default=42)

    parser.add_argument("--output", default="outputs/final_state.png")
    parser.add_argument("--history", default="outputs/history.csv")
    parser.add_argument("--no-show", action="store_true")

    return parser


if __name__ == "__main__":
    run_simulation(build_parser().parse_args())
