from dataclasses import dataclass, field
from math import ceil
from typing import Callable

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from time_output_parser import parse_and_aggregate_time_output_file

"""File that generates grouped bar charts to compare the execution of different implementations"""

colours = ["cadetblue", "blanchedalmond"]

font = {"weight": "normal", "size": 15}

mpl.rc("font", **font)


def time_formatter_ms(time_in_ms: float) -> str:
    """
    aid function to transform a time in milliseconds into a string

    use >= 1.2 since we don't want to display 61 minutes as 1 hour,...
    also show minutes when we are calculating in hours, otherwise only show the largest "unit"
    """
    hour_in_ms = 3600000
    minute_in_ms = 60000
    second_in_ms = 1000
    hours = round(time_in_ms / hour_in_ms, 2)
    if hours >= 1.2:
        hours = int(hours)
        rest_min = int(round((time_in_ms - hours * hour_in_ms) / minute_in_ms, 0))
        return f"{hours} h, {rest_min} min"
    minutes = round(time_in_ms / minute_in_ms, 2)
    if minutes >= 1.2:
        minutes = int(round(minutes, 0))
        return f"{minutes} min"
    seconds = round(time_in_ms / second_in_ms, 2)
    if seconds >= 1.2:
        return f"{seconds} s"
    time_in_ms = round(time_in_ms, 2)
    return f"{time_in_ms} ms"


def time_formatter_sec(time_in_sec: float) -> str:
    """aid function to transform a time in seconds into a string"""
    hour_in_ms = 3600
    minute_in_ms = 60
    second_in_ms = 1
    hours = round(time_in_sec / hour_in_ms, 2)
    if hours >= 1:
        return f"{hours} h"
    minutes = round(time_in_sec / minute_in_ms, 2)
    if minutes >= 1:
        return f"{minutes} min"
    seconds = round(time_in_sec / second_in_ms, 2)
    if seconds >= 1:
        return f"{seconds} s"


def memory_formatter_gb(memory: float) -> str:
    """aid function to return string but with 'GB' added to it"""
    return f"{round(memory, 2)} GB"


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
            "SIHUMI 03",
            "SIHUMI 05",
            "SIHUMI 07",
            "SIHUMI 08",
            "SIHUMI 11",
            "SIHUMI 14",
        ]
    )
    label_formatter: Callable[[float], str] | None = None


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

        if data.label_formatter is not None:
            labels = [data.label_formatter(val) for val in values]
        else:
            labels = ["{:.2f}".format(val) for val in values]
        ax.bar_label(bars, padding=3, labels=labels)
        multiplier += 1

    ax.set_yticks(y_pos + width / 2, labels=data.datasets)
    ax.invert_yaxis()  # labels read top-to-bottom
    ax.set_xlabel(data.x_as)
    # ax.set_ylabel(data.y_as)
    ax.set_title(data.title)
    ax.set_xlim(right=ceil(max(np.array(list(data.data.values())).flatten()) * 1.14))

    # ax.legend(loc="lower right")
    ax.legend()
    ax.margins(0.1, 0.05)
    height = 4 if len(data.datasets) == 2 else 8
    plt.gcf().set_size_inches(15.2, height)

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

    rust_array_build_time_files = [
        "/Users/brdvlami/Documents/Ugent/MA2/Thesis/BenchmarkResults/Rust_suffix_array/BuildOnly/output_suffix_array_rust_immunopeptidomics_build.txt",
        "/Users/brdvlami/Documents/Ugent/MA2/Thesis/BenchmarkResults/Rust_suffix_array/BuildOnly/output_suffix_array_rust_swissprot_build.txt",
    ]

    # parse the files to extract the data
    cpp_tree_execution_data = [
        parse_and_aggregate_time_output_file(file) for file in cpp_tree_build_time_files
    ]
    rust_tree_execution_data = [
        parse_and_aggregate_time_output_file(file)
        for file in rust_tree_build_time_files
    ]

    rust_array_execution_data = [
        parse_and_aggregate_time_output_file(file)
        for file in rust_array_build_time_files
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
            "Tijd (ms) om een match voor alle peptiden te zoeken",
            "Tijd (ms)",
            "Peptidebestand",
            label_formatter=time_formatter_ms,
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
            "Tijd (ms) om met doorzoeken van subboom alle peptiden te zoeken",
            "Tijd (ms)",
            "Peptidebestand",
            label_formatter=time_formatter_ms,
        ),
        ComparisonGraph(
            {
                "C++": [val.execution_time_seconds for val in cpp_tree_execution_data],
                "Rust": [
                    val.execution_time_seconds for val in rust_tree_execution_data
                ],
            },
            "Tijd (s) voor het opbouwen van de suffixboom via Ukkonen",
            "Tijd (s)",
            "Proteïnedatabank",
            ["Human-Prot", "Swiss-Prot"],
            label_formatter=time_formatter_sec,
        ),
        ComparisonGraph(
            {
                "C++": [val.max_mem_size * 1e-6 for val in cpp_tree_execution_data],
                "Rust": [val.max_mem_size * 1e-6 for val in rust_tree_execution_data],
            },
            "Maximale gebruikte hoeveelheid geheugen (GB) voor het opbouwen van de suffixboom via Ukkonen",
            "Geheugengebruik (GB)",
            "Proteïnedatabank",
            ["Human-Prot", "Swiss-Prot"],
            label_formatter=memory_formatter_gb,
        ),
        ComparisonGraph(  # until match found comparison between suffix tree and suffix array
            {
                "Suffixboom": [
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
                "Suffix array": [
                    459.671376953125,
                    522.13283203125,
                    383.06416748046877,
                    92.2581396484375,
                    88.7735302734375,
                    91.50265380859375,
                    92.354150390625,
                    95.3837451171875,
                    90.72038330078125,
                ],
            },
            "Tijd (ms) om een match voor alle peptiden te zoeken",
            "Tijd (ms)",
            "Peptidebestand",
            label_formatter=time_formatter_ms,
        ),
        ComparisonGraph(  # withs subtree
            {
                "Suffixboom": [
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
                "Suffix array": [
                    73284.92388427735,
                    625.5600708007812,
                    461.5416577148437,
                    99.98901611328125,
                    96.99223388671875,
                    108.72378662109375,
                    109.3398388671875,
                    109.3398388671875,
                    100.088798828125,
                ],
            },
            "Tijd (ms) om all voorkomens te vinden waarmee een peptide matcht",
            "Tijd (ms)",
            "Peptidebestand",
            label_formatter=time_formatter_ms,
        ),
        ComparisonGraph(
            {
                "Suffixboom": [
                    val.execution_time_seconds for val in rust_tree_execution_data
                ],
                "Suffix array": [
                    val.execution_time_seconds for val in rust_array_execution_data
                ],
            },
            "Tijd (s) voor het opbouwen van de indexstructuur in Rust",
            "Tijd (s)",
            "Proteïnedatabank",
            ["Human-Prot", "Swiss-Prot"],
            label_formatter=time_formatter_sec,
        ),
        ComparisonGraph(
            {
                "Suffixboom": [
                    val.max_mem_size * 1e-6 for val in rust_tree_execution_data
                ],
                "Suffix array": [
                    val.max_mem_size * 1e-6 for val in rust_array_execution_data
                ],
            },
            "Maximale gebruikte hoeveelheid geheugen (GB) voor het opbouwen van de indexstructuur in Rust",
            "Geheugengebruik (GB)",
            "Proteïnedatabank",
            ["Human-Prot", "Swiss-Prot"],
            label_formatter=memory_formatter_gb,
        ),
        # grafiek die de R-index en suffix array vergelijkt tijdens het opbouwen, ZONDER mapping van tekst naar proteïne of taxon informatie
        ComparisonGraph(
            {
                "R-index": [11.648724, 19.594848, 33.090360, 87.313960, 110.846824],
                "Suffix array": [4.167052, 8.158876, 16.130820, 32.423844, 40.800236],
            },
            "Maximale gebruikte hoeveelheid geheugen (GB) voor het opbouwen van de indexstructuur",
            "Geheugengebruik (GB)",
            "Proteïnedatabank",
            ["0.5%", "1%", "2%", "4%", "5%"],
            label_formatter=memory_formatter_gb,
        ),
        # grafiek die de R-index en suffix array vergelijkt tijdens het opbouwen, ZONDER mapping van tekst naar proteïne of taxon informatie
        ComparisonGraph(
            {
                "R-index": [
                    1.954170576,
                    3.320975480,
                    5.624046578,
                    15.798992742,
                    20.087319696,
                ],
                "Suffix array": [
                    3.837553812,
                    7.516684503,
                    14.863689414,
                    29.87913753,
                    37.598805144,
                ],
            },
            "Indexgrootte na opbouwen in GB",
            "Geheugengebruik (GB)",
            "Proteïnedatabank",
            ["0.5%", "1%", "2%", "4%", "5%"],
            label_formatter=memory_formatter_gb,
        ),
    ]

    for data in all_data:
        create_speed_comparison(data)
