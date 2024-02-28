from dataclasses import dataclass, field
from math import ceil, floor
from typing import Callable

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from time_output_parser import parse_and_aggregate_time_output_file

"""File that generates grouped bar charts to compare the execution of different implementations"""

colours = ["cadetblue", "blanchedalmond", "lightcoral", "darkseagreen", "thistle"]

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
        minutes = int(floor(minutes))
        res = f"{minutes} min"
        seconds = int(round((time_in_ms - minutes * minute_in_ms) / second_in_ms, 0))
        if minutes < 10:
            res += f", {seconds} s"
        return res
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

    width = 1 / (
        len(data.data) + 1
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

    ax.set_yticks(y_pos + width * len(data.data) / 2 - width / 2, labels=data.datasets)
    ax.invert_yaxis()  # labels read top-to-bottom
    ax.set_xlabel(data.x_as)
    # ax.set_ylabel(data.y_as)
    # ax.set_title(data.title)
    ax.set_xlim(right=ceil(max(np.array(list(data.data.values())).flatten()) * 1.14))

    ax.legend(loc="lower right")
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

    rust_array_libdivsufsort_build_time_files = [
        "/Users/brdvlami/Documents/Ugent/MA2/Thesis/BenchmarkResults/Rust_suffix_array/BuildOnly_libdivsufsort/output_suffix_array_rust_immunopeptidomics_build.txt",
        "/Users/brdvlami/Documents/Ugent/MA2/Thesis/BenchmarkResults/Rust_suffix_array/BuildOnly_libdivsufsort/output_suffix_array_rust_swissprot_build.txt",
    ]

    rust_array_libsais_build_time_files = [
        "/Users/brdvlami/Documents/Ugent/MA2/Thesis/BenchmarkResults/Rust_suffix_array/BuildOnly_libsais/output_suffix_array_rust_immunopeptidomics_build.txt",
        "/Users/brdvlami/Documents/Ugent/MA2/Thesis/BenchmarkResults/Rust_suffix_array/BuildOnly_libsais/output_suffix_array_rust_swissprot_build.txt",
    ]

    # parse the files to extract the data
    cpp_tree_execution_data = [
        parse_and_aggregate_time_output_file(file) for file in cpp_tree_build_time_files
    ]
    rust_tree_execution_data = [
        parse_and_aggregate_time_output_file(file)
        for file in rust_tree_build_time_files
    ]

    rust_array_libdivsufsort_execution_data = [
        parse_and_aggregate_time_output_file(file)
        for file in rust_array_libdivsufsort_build_time_files
    ]

    rust_array_libsais_execution_data = [
        parse_and_aggregate_time_output_file(file)
        for file in rust_array_libsais_build_time_files
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
                    val.execution_time_seconds
                    for val in rust_array_libdivsufsort_execution_data
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
                    val.max_mem_size * 1e-6
                    for val in rust_array_libdivsufsort_execution_data
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
        # results for swissprot + no missed cleavage file zoektijden + DENSE mapping van suffix naar protein
        ComparisonGraph(
            {
                "Match": [
                    396.708173828125,
                    637.1782958984375,
                    814.3014501953126,
                    1071.4788232421874,
                    2390.039921875,
                ],
                "All occurrences": [
                    503.968720703125,
                    900.22390625,
                    1826.64822265625,
                    9463.048955078126,
                    133851.62409179687,
                ],
                "Find LCA": [
                    3524.2559912109373,
                    3948.636240234375,
                    4941.030161132812,
                    12628.776108398437,
                    135133.98046875,
                ],
            },
            "Zoektijd Swiss-Prot zonder missed cleavages (dense mapping van suffix naar proteïne)",
            "Tijd (ms)",
            "Sample factor",
            ["1", "2", "3", "4", "5"],
            label_formatter=time_formatter_ms,
        ),
        ComparisonGraph(  # swissprot + no missed cleavage file zoektijden
            {
                "Dense": [
                    3524.2559912109373,
                    3948.636240234375,
                    4941.030161132812,
                    12628.776108398437,
                    135133.98046875,
                ],
                "Sparse": [
                    3680.764990234375,
                    4154.237607421875,
                    5239.153872070313,
                    13705.790971679688,
                    144266.00616210938,
                ],
            },
            "Zoektijd naar LCA (taxon id) voor dense vs sparse mapping voor verschillende sample factors",
            "Tijd (ms)",
            "Sample factor",
            ["1", "2", "3", "4", "5"],
            label_formatter=time_formatter_ms,
        ),
        ComparisonGraph(
            {
                "Manueel aanmaken van threads": [
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
                ],
                "Threading m.b.v. Rayon ": [
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
                ],
            },
            "Tijd om de LCA* van het Swiss-Prot peptidebestand zonder missed cleavage te zoeken in 5% van UniProt",
            "Tijd (ms)",
            "Aantal threads",
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
            label_formatter=time_formatter_ms,
        ),
        ComparisonGraph(
            {
                "Suffix array (libsais)": [
                    val.execution_time_seconds
                    for val in rust_array_libsais_execution_data
                ],
                "Suffixboom (Rust)": [
                    val.execution_time_seconds for val in rust_tree_execution_data
                ],
            },
            "Tijd (s) voor het opbouwen van de indexstructuur",
            "Tijd (s)",
            "Proteïnedatabank",
            ["Human-Prot", "Swiss-Prot"],
            label_formatter=time_formatter_sec,
        ),
        ComparisonGraph(
            {
                "Suffix array (libsais)": [
                    val.max_mem_size * 1e-6 for val in rust_array_libsais_execution_data
                ],
                "Suffixboom (Rust)": [
                    val.max_mem_size * 1e-6 for val in rust_tree_execution_data
                ],
            },
            "Maximale gebruikte hoeveelheid geheugen (GB) voor het opbouwen van de indexstructuur",
            "Geheugengebruik (GB)",
            "Proteïnedatabank",
            ["Human-Prot", "Swiss-Prot"],
            label_formatter=memory_formatter_gb,
        ),
        ComparisonGraph(
            {
                "Suffix array (libsais)": [
                    val.execution_time_seconds
                    for val in rust_array_libsais_execution_data
                ],
                "Suffix array (libdivsufsort)": [
                    val.execution_time_seconds
                    for val in rust_array_libdivsufsort_execution_data
                ],
            },
            "Tijd (s) voor het opbouwen van de indexstructuur",
            "Tijd (s)",
            "Proteïnedatabank",
            ["Human-Prot", "Swiss-Prot"],
            label_formatter=time_formatter_sec,
        ),
        ComparisonGraph(
            {
                "Suffix array (libsais)": [
                    val.max_mem_size * 1e-6 for val in rust_array_libsais_execution_data
                ],
                "Suffix array (libdivsufsort)": [
                    val.max_mem_size * 1e-6
                    for val in rust_array_libdivsufsort_execution_data
                ],
            },
            "Maximale gebruikte hoeveelheid geheugen (GB) voor het opbouwen van de indexstructuur",
            "Geheugengebruik (GB)",
            "Proteïnedatabank",
            ["Human-Prot", "Swiss-Prot"],
            label_formatter=memory_formatter_gb,
        ),
        ComparisonGraph(  # suffixboom (met lca voorberekend) vs suffixarray (met lca*), met cutoff op 10000 matches voor de SA
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
                "Suffix array 1 thread (niet sparse)": [
                    1750.2314453125,
                    924.6883544921875,
                    727.555029296875,
                    159.974169921875,
                    159.2865478515625,
                    178.650390625,
                    183.7022705078125,
                    184.5429931640625,
                    153.8831787109375,
                ],
                "Suffix array 12 threads (niet sparse)": [
                    417.0977783203125,
                    185.4984375,
                    173.9829833984375,
                    40.562841796875,
                    44.768310546875,
                    42.59638671875,
                    43.5115966796875,
                    42.554638671875,
                    40.7529541015625,
                ],
            },
            "Tijd (ms) om het geaggregeerde taxon ID te vinden voor elke peptide met cutoff op 10000 matches",
            "Tijd (ms)",
            "Peptidebestand",
            label_formatter=time_formatter_ms,
        ),
        ComparisonGraph(
            # suffixboom (met lca voorberekend) vs suffixarray (met lca*), zonder cutoff voor de SA
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
                "Suffix array 1 thread (niet sparse)": [
                    572101.4310058594,
                    947.8179931640625,
                    729.418505859375,
                    160.544873046875,
                    162.6793212890625,
                    178.7932373046875,
                    184.179833984375,
                    187.9348876953125,
                    151.7920654296875,
                ],
                "Suffix array 12 threads (niet sparse)": [
                    163968.31882324218,
                    193.1651611328125,
                    188.2863037109375,
                    40.9485595703125,
                    42.1366455078125,
                    42.0172119140625,
                    42.5451904296875,
                    42.07587890625,
                    39.64833984375,
                ],
            },
            "Tijd (ms) om het geaggregeerde taxon ID te vinden voor elke peptide zonder cutoff",
            "Tijd (ms)",
            "Peptidebestand",
            label_formatter=time_formatter_ms,
        ),
        ComparisonGraph(
            # evolutie van zoektijd naar taxonid (1 tread)
            {
                "sparseness factor 1": [910.463134765625, 717.8666015625],
                "sparseness factor 2": [1409.4352783203126, 1054.492822265625],
                "sparseness factor 3": [2553.0102783203124, 1568.374853515625],
                "sparseness factor 4": [11808.8814453125, 4321.583813476563],
                "sparseness factor 5": [155939.5229003906, 48842.2990234375],
            },
            "Tijd (ms) om het geaggregeerde taxon ID te vinden voor elke peptide zonder cutoff",
            "Tijd (ms)",
            "Peptidebestand",
            ["Swiss-Prot zonder missed cleavages", "Swiss-Prot met missed cleavages"],
            label_formatter=time_formatter_ms,
        ),
    ]

    for data in all_data:
        create_speed_comparison(data)
