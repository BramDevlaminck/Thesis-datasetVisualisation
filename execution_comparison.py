from dataclasses import dataclass, field
from math import ceil, floor
from typing import Callable

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from time_output_parser import parse_and_aggregate_time_output_file

"""File that generates grouped bar charts to compare the execution of different implementations"""

colours = ["cadetblue", "blanchedalmond", "lightcoral", "darkseagreen", "thistle"]

font = {"weight": "normal", "size": 21}

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
        hours = int(floor(hours))
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
    hour_in_sec = 3600
    minute_in_sec = 60
    second_in_sec = 1
    hours = round(time_in_sec / hour_in_sec, 2)
    if hours >= 1:
        hours = int(floor(hours))
        rest_min = int(round((time_in_sec - hours * hour_in_sec) / minute_in_sec, 0))
        return f"{hours} h, {rest_min} min"
    minutes = round(time_in_sec / minute_in_sec, 2)
    if minutes >= 1:
        minutes = int(floor(minutes))
        res = f"{minutes} min"
        seconds = int(round((time_in_sec - minutes * minute_in_sec) / second_in_sec, 0))
        # if minutes < 10:
        res += f", {seconds} s"
        return res
    seconds = round(time_in_sec / second_in_sec, 2)
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
            "Swiss-Prot tryptisch",
            "Swiss-Prot missed cleavage",
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
    ax.set_title(data.title)
    ax.set_xlim(right=ceil(max(np.array(list(data.data.values())).flatten()) * 1.10))

    ax.legend(loc="lower right")
    # ax.legend()
    ax.margins(0.1, 0.05)
    height = 4.5 if len(data.datasets) <= 2 else 15
    plt.gcf().set_size_inches(20, height)

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
        ComparisonGraph(
            # suffixboom (met lca voorberekend) vs suffixarray (met lca*), met cutoff op 10000 matches voor de SA
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
                "Suffix array 1 thread (k = 1)": [
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
                "Suffix array 12 threads (k = 1)": [
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
            "Tijd (ms) voor zoeken + taxonomische analyse\nvoor elke peptide met cutoff op 10000 matches",
            "Tijd (ms)",
            "Peptidebestand",
            label_formatter=time_formatter_ms,
        ),
        ComparisonGraph(
            # suffixboom (met lca voorberekend) vs suffixarray (met lca*), zonder cutoff voor de SA, de SA heeft altijd sparseness factor 1
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
                "Suffix array 1 thread": [
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
                "Suffix array 12 threads": [
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
                "k = 1": [910.463134765625, 717.8666015625],
                "k = 2": [1409.4352783203126, 1054.492822265625],
                "k = 3": [2553.0102783203124, 1568.374853515625],
                "k = 4": [11808.8814453125, 4321.583813476563],
                "k = 5": [155939.5229003906, 48842.2990234375],
            },
            "Tijd (ms) om het geaggregeerde taxon ID te vinden voor elke peptide zonder cutoff",
            "Tijd (ms)",
            "Peptidebestand",
            ["Swiss-Prot tryptisch", "Swiss-Prot missed cleavages"],
            label_formatter=time_formatter_ms,
        ),
        ComparisonGraph(
            {
                "Dense": [
                    52261.39501953125,
                    379.025146484375,
                    291.041015625,
                    61.3359375,
                    60.760009765625,
                    65.18994140625,
                    68.41796875,
                    66.320068359375,
                    66.740966796875,
                ],
                "Sparse": [
                    182869.44775390625,
                    331.035888671875,
                    270.232177734375,
                    72.0390625,
                    64.677978515625,
                    65.0830078125,
                    75.007080078125,
                    74.037109375,
                    66.64599609375,
                ],
            },
            "Tijd (ms) om alle voorkomens te vinden gebruik makende van een dense of sparse mapping van suffix naar proteïne",
            "Tijd (ms)",
            "Peptidebestand",
            label_formatter=time_formatter_ms,
        ),
        ComparisonGraph(
            {
                "Dense": [
                    val * 1e-9
                    for val in [
                        1208909824,
                        3357868032,
                        3221602304,
                        3326017536,
                        3307388928,
                        3249684480,
                        3317972992,
                        3317022720,
                        3260497920,
                    ]
                ],
                "Sparse": [
                    val * 1e-9
                    for val in [
                        1035321344,
                        2491056128,
                        2501902336,
                        2567421952,
                        2594766848,
                        2501951488,
                        2498560000,
                        2498691072,
                        2591326208,
                    ]
                ],
            },
            "Geheugen (GB) om alle voorkomens te vinden gebruik makende van een dense of sparse mapping van suffix naar proteïne",
            "Geheugen (GB)",
            "Peptidebestand",
            label_formatter=memory_formatter_gb,
        ),
        ComparisonGraph(
            {
                "Dense": [
                    val * 1e-9
                    for val in [
                        1208909824,
                        3294756864,
                    ]
                ],
                "Sparse": [
                    val * 1e-9
                    for val in [
                        1035321344,
                        2530709504,
                    ]
                ],
            },
            "Geheugen (GB) gebruik tijdens het zoeken in de Human-Prot of Swiss-Prot databank voor de sparse en dense mapping",
            "Geheugen (GB)",
            "Peptidebestand",
            ["Human-prot", "Swiss-Prot"],
            label_formatter=memory_formatter_gb,
        ),
        ComparisonGraph(  # dit gaat over swissprot
            {
                "SA grootte (zonder tekst)": [
                    1.652189544,
                    0.826094776,
                    0.550729848,
                    0.413047392,
                    0.330437912,
                ],
                "Index grootte (SA + tekst)": [
                    val + 0.206523693
                    for val in [
                        1.652189544,
                        0.826094776,
                        0.550729848,
                        0.413047392,
                        0.330437912,
                    ]
                ],
            },
            "Grootte van de indexstructuur (GB) voor verschillende sparseness factoren",
            "Grootte (GB)",
            "Sparseness factor",
            [
                "sparseness factor 1",
                "sparseness factor 2",
                "sparseness factor 3",
                "sparseness factor 4",
                "sparseness factor 5",
            ],
            label_formatter=memory_formatter_gb,
        ),
        ComparisonGraph(  # zoektijd in volledige Uniprot database, dit is in de index waar we I en L apart houden, en dan via een boom proberen zoeken naar waar ze gelijk zijn
            {
                "Zoeken I = L": [
                    58.5436666666667,
                    85.698,
                    23.9133333333333,
                    1.9925,
                    1.7813333333,
                    3.046,
                    5.801,
                    4.353,
                    1.796,
                ],
                "Zoeken I ≠ L": [
                    68.9836666667,
                    103.06666666667,
                    24.294,
                    4.8796666666667,
                    4.2696666666667,
                    6.283,
                    9.1523333333333,
                    8.2953333333333,
                    4.4906666666667,
                ],
            },
            "Tijd (in s) om peptides te zoeken in UniProtKB",
            "Tijd (s)",
            "",
            label_formatter=time_formatter_sec,
        ),
        ComparisonGraph(  # zoektijd in volledige Uniprot database, dit is door I en L in de index gelijk te stellen, en als ze dit niet zijn achteraf te filteren
            {
                "Zoeken I = L": [
                    54.257,
                    72.1886666667,
                    17.221,
                    2.014,
                    1.661,
                    2.9173333333333,
                    4.781,
                    3.809,
                    1.7293333333333,
                ],
                "Zoeken I ≠ L": [
                    50.577,
                    58.9803333333333,
                    15.9306666666667,
                    3.118,
                    2.8753333333333,
                    4.088,
                    6.0046666666667,
                    5.1613333333333,
                    2.824,
                ],
            },
            "Tijd (in s) om peptides te zoeken in UniProtKB",
            "Tijd (s)",
            "",
            label_formatter=time_formatter_sec,
        ),
        ComparisonGraph(
            # grafiek voor de les computationele biologie, oude unipept op swissprot vs suffixtree
            {
                "Huidige Unipept index": [132, 1537],
                "Suffixboom": [0.15298, 0.14064],
            },
            "",
            "Tijd (s)",
            "",
            ["zonder missed cleavages", "met missed cleavages"],
            label_formatter=time_formatter_sec,
        ),
        ComparisonGraph(
            # voor de les computationele biologie, oude unipept op swissprot vs suffixtree
            {
                "Huidige Unipept index": [6.7],
                "Suffixboom": [85],
            },
            "",
            "Geheugengebruik (GB)",
            "Proteïnedatabank",
            ["Swiss-Prot"],
            label_formatter=memory_formatter_gb,
        ),
        ComparisonGraph(
            # voor de les computationele biologie, huidige uniprot vs suffixboom vs suffix array
            {
                "Huidige Unipept index": [6.7],
                "Suffixboom": [85],
                "Suffix array": [2.64],
            },
            "",
            "Geheugengebruik (GB)",
            "Proteïnedatabank",
            ["Swiss-Prot"],
            label_formatter=memory_formatter_gb,
        ),
        ComparisonGraph(
            # voor de les computationele biologie, huidige uniprot vs suffixboom vs suffix array
            {
                "Suffixboom": [152.98, 140.64],
                "Suffix array": [947.82, 729.42],
            },
            "",
            "Tijd (ms)",
            "Proteïnedatabank",
            ["zonder missed cleavages", "met missed cleavages"],
            label_formatter=time_formatter_ms,
        ),
        ComparisonGraph(
            # voor de les computationele biologie, huidige uniprot vs suffixboom vs suffix array
            {
                "Suffix array": [59, 16],
            },
            "",
            "Tijd (s)",
            "Proteïnedatabank",
            ["zonder missed cleavages", "met missed cleavages"],
            label_formatter=time_formatter_sec,
        ),
        ComparisonGraph(  # dit gaat over uniprotKB
            {
                "SA grootte (zonder tekst)": [705.62, 352.81, 235.21, 176.40, 141.124],
                "Index grootte (SA + tekst)": [
                    val + 87.055 for val in [705.62, 352.81, 235.21, 176.40, 141.124]
                ],
            },
            "Grootte van de indexstructuur (GB) voor verschillende sparseness factoren",
            "Grootte (GB)",
            "Sparseness factor",
            [
                "k = 1",
                "k = 2",
                "k = 3",
                "k = 4",
                "k = 5",
            ],
            label_formatter=memory_formatter_gb,
        ),
        ComparisonGraph(  # dit gaat over 5%, 10%,... van uniprot
            {
                "Suffix array": [
                    4.167052,
                    8.158876,
                    16.130820,
                    32.423844,
                    40.800236,
                ],
                "R-index": [
                    11.648724,
                    19.594848,
                    33.090360,
                    87.313960,
                    110.846824,
                ],
            },
            "Maximaal geheugengebruik tijdens het opbouwen van de suffix array en R-index voor verschillende databanken",
            "Grootte (GB)",
            "Fractie van UniProtKB",
            [
                "0.5%",
                "1%",
                "2%",
                "4%",
                "5%",
            ],
            label_formatter=memory_formatter_gb,
        ),
        ComparisonGraph(  # dit gaat over 5%, 10%,... van uniprot
            {
                "Suffix array (SA + tekst)": [
                    3.837553812,
                    7.516684503,
                    14.863689414,
                    29.87913753,
                    37.598805144,
                ],
                "R-index": [
                    1.954170576,
                    3.320975480,
                    5.624046578,
                    15.798992742,
                    20.087319696,
                ],
            },
            "Grootte van de resulterende suffix array en R-index voor verschillende databanken",
            "Grootte (GB)",
            "Fractie van UniProtKB",
            [
                "0.5%",
                "1%",
                "2%",
                "4%",
                "5%",
            ],
            label_formatter=memory_formatter_gb,
        ),
        ComparisonGraph(  # font size 21, width 10, height 15, x_limit factor 1.21
            {
                "Unipept 6.x (I ≠ L)": [
                    3.2986666667,
                    3.1353333333,
                    5.1593333333,
                    8.08333333333,
                    7.5296666667,
                    4.05066666667,
                ],
                "Unipept 6.x (I = L)": [
                    8.05,
                    6.9283333333,
                    9.979,
                    13.7416666666667,
                    11.9933333333333,
                    6.716,
                ],
                "Unipept 5.x (I ≠ L)": [230, 185, 435, 216, 610, 101],
                "Unipept 5.x (I = L)": [
                    202,
                    143,
                    464,
                    330,
                    845,
                    84,
                ],
            },
            "Tijd (s) nodig voor het verwerken van een peptidebestand\nmet missed cleavage handling",
            "Time (s)",
            "Proteïnedatabank",
            [
                "SIHUMI 03",
                "SIHUMI 05",
                "SIHUMI 07",
                "SIHUMI 08",
                "SIHUMI 11",
                "SIHUMI 14",
            ],
            label_formatter=time_formatter_sec,
        ),
        ComparisonGraph(  # font size 21, width 10, height 15, x_limit factor 1.21
            {
                "Unipept 6.x": [145.66, 115.81],
                "Unipept 5.x": [93, 110],
                # "New Unipept (I ≠ L)": [
                #     145.66,
                # ],
                # "New Unipept (I = L)": [
                #     115.81,
                # ],
                # "Old Unipept (I ≠ L)": [93],
                # "Old Unipept (I = L)": [
                #     110
                # ],
            },
            "Tijd (s) nodig voor het verwerken van 100 000 tryptische peptiden",
            "Time (s)",
            "Proteïnedatabank",
            ["I ≠ L", "I = L"],
            label_formatter=time_formatter_sec,
        ),
        ComparisonGraph(  # search only Unipept 6.x
            {
                "I ≠ L": [
                    117.936,
                    24.500,
                    3.1743333333333,
                    2.4996666666667,
                    5.005,
                    8.3986666666667,
                    6.139,
                    2.5673333333333,
                ],
                "I = L": [
                    95.4246666667,
                    31.3756666666667,
                    4.562,
                    4.2593333333333,
                    6.7183333333333,
                    9.9463333333333,
                    6.319,
                    3.4423333333333,
                ],
            },
            "Tijd (s) nodig voor het verwerken van een peptidebestand met missed cleavage handling",
            "Time (s)",
            "Proteïnedatabank",
            [
                "Swiss-Prot tryptisch",
                "Swiss-Prot missed cleavages",
                "SIHUMI 03",
                "SIHUMI 05",
                "SIHUMI 07",
                "SIHUMI 08",
                "SIHUMI 11",
                "SIHUMI 14",
            ],
            label_formatter=time_formatter_sec,
        ),
    ]

    libsais_build_time_uniprot = [
        31883639.763671875,
        35007841.790771484,
        34330708.56616211,
        27523439.407958984,
        29120545.077392578,
        27644291.99584961,
    ]

    libdivsufsort_build_time_uniprot = [
        17697356.158203125,
        16287824.045898438,
        18641389.829833984,
        18943567.30859375,
        18042967.510253906,
        19415305.08618164,
    ]

    all_data.extend(
        [
            ComparisonGraph(  # graph with execution time for building uniprot
                {
                    "Libdivsufsort": [
                        sum(libdivsufsort_build_time_uniprot)
                        * 0.001
                        / len(libdivsufsort_build_time_uniprot)
                    ],
                    "Libsais": [
                        sum(libsais_build_time_uniprot)
                        * 0.001
                        / len(libsais_build_time_uniprot)
                    ],
                },
                "Tijd (s) om de suffixarray op te bouwen voor UniProtKB",
                "Tijd (s)",
                "",
                [""],
                label_formatter=time_formatter_sec,
            ),
            ComparisonGraph(
                {
                    "Libdivsufsort": [734.47],
                    "Libsais": [731.10],
                },
                "Geheugengebruik (GB) om de suffixarray op te bouwen voor UniProtKB",
                "Geheugengebruik (GB)",
                "",
                [""],
                label_formatter=memory_formatter_gb,
            ),
        ]
    )

    for data in all_data:
        create_speed_comparison(data)
