from dataclasses import dataclass, field
from math import ceil

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from execution_comparison import ComparisonGraph, time_formatter_ms

colours = ["cadetblue", "lightcoral", "darkseagreen", "thistle"]

font = {"weight": "normal", "size": 15}

mpl.rc("font", **font)


@dataclass
class ComparisonLineGraph:
    data: dict[str, tuple[str, list[float]]]
    title: str
    x_as: str
    y_as: str
    datasets: list[str] = (
        field(
            default_factory=lambda: [
                "Human-Prot",
                "Swiss-Prot zonder missed cleavage",
                "Swiss-Prot met missed cleavage",
                "SIHUMI 03",
                "SIHUMI 05",
                "SIHUMI 07",
                "SIHUMI 08",
                "SIHUMI 11",
                "SIHUMI 14",
            ],
        ),
    )


def create_relative_comparison(
    data: ComparisonLineGraph, output_name: str | None = None
):
    fig, ax = plt.subplots(layout="constrained")

    for i, (key, (linestyle, values)) in enumerate(data.data.items()):
        ax.plot(
            data.datasets,
            values,
            color=colours[i],
            marker="o",
            linestyle=linestyle,
            label=key,
        )

    # flatten all the values into a 1d array and take max
    max_val = max([x for row in data.data.values() for x in row[1]])

    ax.set_ylim([0, max_val * 1.05])
    ax.set_xlabel(data.x_as)
    ax.set_ylabel(data.y_as)
    # ax.set_title(data.title)
    ax.set_xticks(np.arange(min(data.datasets), max(data.datasets) + 1, 1.0))
    ax.grid(color="lightgray", linestyle="-", linewidth=1)
    ax.legend(loc="lower right")
    ax.margins(0.1, 0.05)
    height = 3 if len(data.datasets) == 2 else 6
    plt.gcf().set_size_inches(10, height)

    if output_name is not None:
        plt.savefig(output_name)
    plt.show()


if __name__ == "__main__":
    manual_threads = [
        28953.89794921875,
        15479.1708984375,
        10773.0380859375,
        7757.5771484375,
        6313.8115234375,
        5045.465576171875,
        4561.252685546875,
        4009.536865234375,
        3693.316650390625,
        3283.839599609375,
        2954.5703125,
        2681.176025390625,
    ]

    rayon_threads = [
        25481.59521484375,
        13808.319580078125,
        9103.64990234375,
        6982.69580078125,
        5301.76708984375,
        4554.90966796875,
        3817.50537109375,
        3413.670654296875,
        3035.99169921875,
        2764.484375,
        2505.06689453125,
        2333.380615234375,
    ]

    # make the speedup relative
    manual_threads_start = manual_threads[0]
    for i, val in enumerate(manual_threads):
        manual_threads[i] = manual_threads_start / val

    rayon_threads_start = rayon_threads[0]
    for i, val in enumerate(rayon_threads):
        rayon_threads[i] = rayon_threads_start / val

    data = [
        ComparisonLineGraph(
            {
                "Manueel aanmaken van threads": ("solid", manual_threads),
                "Threading m.b.v. Rayon ": ("solid", rayon_threads),
                "Perfecte schaling": (
                    "dashed",
                    [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
                ),
            },
            "Versnelling van de uitvoeringstijd en opzichte van het aantal threads",
            "Aantal threads",
            "Relatieve versnelling",
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
        )
    ]

    for graph in data:
        create_relative_comparison(graph)
