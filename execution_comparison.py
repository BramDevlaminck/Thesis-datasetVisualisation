from dataclasses import dataclass, field
from math import ceil

import matplotlib.pyplot as plt
import numpy as np

from time_output_parser import parse_and_aggregate_time_output_file

"""File that generates grouped bar charts to compare the execution of different implementations"""

colours = ["cadetblue", "blanchedalmond"]


@dataclass
class ComparisonGraph:
    data: dict[str, list[float]]
    title: str
    x_as: str
    y_as: str
    datasets: list[str] = field(
        default_factory=lambda: [
            "Human-Prot",
            "Swiss-Prot zonder missed cleavage",
            "Swiss-Prot met missed cleavage",
            "SIHUMI S03",
            "SIHUMI S05",
            "SIHUMI S07",
            "SIHUMI S08",
            "SIHUMI S11",
            "SIHUMI S14",
        ]
    )


def create_speed_comparison(data: ComparisonGraph, output_name: str | None = None):
    fig, ax = plt.subplots(layout="constrained")
    y_pos = np.arange(len(data.datasets))

    width = (
        1 / 3
    )  # the width of the bars (divide by 1 bigger than the number of algorithms we are comparing!)
    multiplier = 0
    for i, (key, values) in enumerate(data.data.items()):
        offset = width * multiplier
        bars = ax.barh(
            y_pos + offset,
            width=values,
            height=width,
            label=key,
            color=colours[i],
            edgecolor="black",
            alpha=0.7,
        )
        ax.bar_label(bars, padding=3, labels=["{:.2f}".format(val) for val in values])
        multiplier += 1

    ax.set_yticks(y_pos + width / 2, labels=data.datasets)
    ax.invert_yaxis()  # labels read top-to-bottom
    ax.set_xlabel(data.x_as)
    ax.set_ylabel(data.y_as)
    ax.set_title(data.title)
    ax.set_xlim(right=ceil(max(np.array(list(data.data.values())).flatten()) * 1.14))

    ax.legend()
    ax.margins(0.1, 0.05)
    height = 4 if len(data.datasets) == 2 else 8
    plt.gcf().set_size_inches(12, height)

    if output_name is not None:
        plt.savefig(output_name)
    plt.show()


if __name__ == "__main__":
    # the files that contain the output from the time command about building the tree
    cpp_tree_build_time_files = [
        "/Users/brdvlami/Documents/Ugent/MA2/Thesis/BenchmarkResults/CPP_ukkonen_CCB/TreeBuild/immunopeptidomics/1_output_ukkonen_cpp_28children.txt",
        "/Users/brdvlami/Documents/Ugent/MA2/Thesis/BenchmarkResults/CPP_ukkonen_CCB/TreeBuild/swissprot/output_cpp_28children_swissprot.txt",
    ]
    rust_tree_build_time_files = [
        "/Users/brdvlami/Documents/Ugent/MA2/Thesis/BenchmarkResults/Rust_ukkonen/TreeBuild/output_rust_immunopeptidomics_build.txt",
        "/Users/brdvlami/Documents/Ugent/MA2/Thesis/BenchmarkResults/Rust_ukkonen/TreeBuild/output_rust_swissprot_build.txt",
    ]

    # parse the files to extract the data
    cpp_execution_data = [
        parse_and_aggregate_time_output_file(file) for file in cpp_tree_build_time_files
    ]
    rust_execution_data = [
        parse_and_aggregate_time_output_file(file)
        for file in rust_tree_build_time_files
    ]

    all_data = [
        ComparisonGraph(  # until match found
            {
                "C++": [
                    204.374,
                    147.355,
                    134.187,
                    17.3605,
                    16.7635,
                    17.8944,
                    17.6533,
                    17.249,
                    17.3359,
                ],
                "Rust": [
                    219.67120659179687,
                    152.98342421875,
                    140.635694921875,
                    17.862667626953126,
                    17.937821044921876,
                    17.9280654296875,
                    18.086377294921874,
                    17.827000048828126,
                    18.37514921875,
                ],
            },
            "Tijd in milliseconden om een match voor alle peptiden te zoeken",
            "Tijd in milliseconden",
            "Zoekbestand",
        ),
        ComparisonGraph(  # withs subtree
            {
                "C++": [
                    1.45269e06,
                    1.57777e07,
                    1.54936e07,
                    3.67146e06,
                    3.67822e06,
                    3.56324e06,
                    3.60447e06,
                    3.63128e06,
                    3.65424e06,
                ],
                "Rust": [
                    861712.7821858724,
                    287.7971191406,
                    210.38191731770834,
                    41.4705810546875,
                    37.624462890625,
                    51.1537353515625,
                    50.9006103515625,
                    50.213818359375,
                    35.260107421875,
                ],
            },
            "Tijd in milliseconden om met doorzoeken van subboom alle peptiden te zoeken",
            "Tijd in milliseconden",
            "Zoekbestand",
        ),
        ComparisonGraph(
            {
                "C++": [val.execution_time_seconds for val in cpp_execution_data],
                "Rust": [val.execution_time_seconds for val in rust_execution_data],
            },
            "Tijd in seconden voor het opbouwen van de suffixboom via Ukkonen",
            "Tijd in seconden",
            "Proteïne databank",
            ["Human-Prot", "Swiss-Prot"],
        ),
        ComparisonGraph(
            {
                "C++": [val.max_mem_size * 1e-6 for val in cpp_execution_data],
                "Rust": [val.max_mem_size * 1e-6 for val in rust_execution_data],
            },
            "Maximale gebruikte hoeveelheid geheugen in GB voor het opbouwen van de suffixboom via Ukkonen",
            "Geheugengebruik in GB",
            "Proteïne databank",
            ["Human-Prot", "Swiss-Prot"],
        ),
    ]

    for data in all_data:
        create_speed_comparison(data)
