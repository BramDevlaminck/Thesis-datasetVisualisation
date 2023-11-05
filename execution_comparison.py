from dataclasses import dataclass, field

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
            "Swiss-Prot without missed cleavage",
            "Swiss-Prot with missed cleavage",
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

    ax.legend()
    ax.margins(0.1, 0.05)
    plt.gcf().set_size_inches(12, 4)

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
                "C++": [209.522, 170.814, 179.346],
                "Rust": [226.07584635416666, 175.210205078125, 151.49943033854166],
            },
            "Gemiddelde Tijd in ms om een match voor alle Peptiden te zoeken",
            "Tijd in ms",
            "Zoekbestand",
        ),
        ComparisonGraph(
            {
                "C++": [
                    1.45269e06,
                    1.57777e07,
                    0,
                ],  # TODO: laatste waarde is een placeholder, wordt atm uitgerekend
                "Rust": [861712.7821858724, 287.7971191406, 210.38191731770834],
            },
            "Gemiddelde Tijd in ms om met doorzoeken van subboom alle Peptiden te zoeken",
            "Tijd in ms",
            "Zoekbestand",
        ),
        ComparisonGraph(
            {
                "C++": [val.execution_time_seconds for val in cpp_execution_data],
                "Rust": [val.execution_time_seconds for val in rust_execution_data],
            },
            "Gemiddelde Tijd in seconden voor het opbouwen van de suffixboom via Ukkonen",
            "Tijd in seconden",
            "proteinen databank",
            ["Human-Prot", "Swiss-Prot"],
        ),
        ComparisonGraph(
            {
                "C++": [val.max_mem_size * 1e-6 for val in cpp_execution_data],
                "Rust": [val.max_mem_size * 1e-6 for val in rust_execution_data],
            },
            "Maximale gebruikte hoeveelheid geheugen in GB voor het opbouwen van de suffixboom via Ukkonen",
            "Gebruikt geheugen in GB",
            "proteinen databank",
            ["Human-Prot", "Swiss-Prot"],
        ),
    ]

    for data in all_data:
        create_speed_comparison(data)
