from dataclasses import dataclass
from math import floor

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import ticker


def round_up(num_to_round: float, multiple_to_round_to: float) -> float:
    """Aid function to round up to a certain multiple"""
    if multiple_to_round_to == 0:
        return num_to_round

    remainder = num_to_round % multiple_to_round_to
    if remainder == 0:
        return num_to_round

    return round(num_to_round + multiple_to_round_to - remainder)


def calculate_character_statistics(proteins: list[str]) -> dict[str, int]:
    """Creates a list of all the individual characters of the proteins"""
    alphabet = "abcdefghijklmnopqrstuvwxyz".upper()
    statistiscs: dict[str, int] = {letter: 0 for letter in alphabet}
    for line in proteins:
        for char in line:
            statistiscs[char.upper()] += 1

    return statistiscs


def calculate_length_statistics(proteins: list[str]) -> list[int]:
    """Creates a list containing the lengths off all the proteins"""
    statistiscs: list[int] = []
    for line in proteins:
        statistiscs.append(len(line))

    return statistiscs


def make_alphabet_complete(
    data: tuple[list[str], list[int]]
) -> tuple[list[str], list[int]]:
    """
    Data is a tuple of (char, count) in practice. These chars are already sorted alphabetically, BUT not all the characters are always there!

    This function adds the missing characters from the alphabet at the right place with count 0 as value
    """
    alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"

    labels, counts = data
    for i, letter in enumerate(alphabet):
        if i >= len(labels) or labels[i] != letter:
            labels.insert(i, letter)
            counts.insert(i, 0)

    return labels, counts


def create_barh_acid_code_occurrences(
    proteins: list[str], short_name: str, output_file: str | None = None
):
    """Creates a horizontal barchart of the distribution of the letters in the proteins"""

    statistics = calculate_character_statistics(proteins)
    data_tuples = [(k, v) for k, v in statistics.items()]
    data_tuples.sort(key=lambda x: x[0])
    labels, counts = list(zip(*data_tuples))
    # counts for char in statistics how many times it occurred, the char will be places in labels (sorted alphabetically) and the count in counts
    # labels, counts = np.unique(statistics, return_counts=True)
    # some letters from the alphabet will not occur in this dataset, add them at the right place with count 0
    # labels, counts = make_alphabet_complete((list(labels), list(counts)))

    create_barh(
        labels,
        counts,
        "Aminozuur",
        "Aantal voorkomens",
        f"Aantal voorkomens per \naminozuur voor {short_name}",
        output_file,
    )


def create_barchart_vertical(statistics: list[any]):
    labels, counts = np.unique(statistics, return_counts=True)
    labels, counts = make_alphabet_complete((list(labels), list(counts)))
    ticks = range(len(counts))
    fig, ax = plt.subplots()

    bar_container = ax.bar(
        ticks, counts, color="cadetblue", edgecolor="black", alpha=0.7, align="center"
    )
    ax.bar_label(bar_container, fmt="{:,.0f}", rotation="vertical", padding=3)
    plt.xticks(ticks, labels)
    # add enough y-label hight to make sure the labels above the bar don't overlap with the title
    ax.set_ylim(0, max(counts) * 1.2)
    # Add labels and a title
    plt.xlabel("Amino Acid Code")
    plt.ylabel("Number of Occurrences")
    plt.title("Number of Occurrences per\nAmino Acid code for the swissprot dataset")

    # Add a grid
    plt.grid(True, linestyle="--", alpha=0.5, axis="y")
    plt.gcf().set_size_inches(11.69, 8.27)

    # Show the plot
    # plt.show()
    plt.savefig("Graphs/character_occurrences_swissprot_var1_protein_database.jpg")


def create_barh_protein_length(
    proteins: list[str],
    bin_size: int,
    short_name: str,
    max_allowed_protein_length: float,
    sequence_type: str,
    output_file: str | None = None,
):
    """From the given proteins calculate their length and create a horizontal barchart with the data"""
    statistics = calculate_length_statistics(proteins)
    max_bin = int(
        min(
            max_allowed_protein_length - bin_size,
            (max(statistics) // bin_size) * bin_size,
        )
    )
    frequency_dict: dict[int, int] = {
        key: 0 for key in range(0, max_bin + bin_size, bin_size)
    }

    for val in statistics:
        if val < max_allowed_protein_length:
            frequency_dict[(val // bin_size) * bin_size] += 1

    labels, counts = list(zip(*frequency_dict.items()))
    labels = [f"{label}-{label + bin_size - 1}" for label in labels]

    create_barh(
        labels,
        counts,
        f"{sequence_type}lengte",
        "Aantal voorkomens",
        f"Verdeling van de {sequence_type.lower()}lengte voor {short_name}",
        output_file,
    )


def create_barh(
    labels: list[any],
    counts: list[int],
    ylabel: str,
    x_label: str,
    title: str,
    file_name: str | None = None,
):
    """
    creates a horizontal barchart with the provided labels counts and axis labels

    if file_name is not none the barchart is stored to a file with the given name
    """
    reversed(labels)
    ticks = range(len(counts), 0, -1)
    fig, ax = plt.subplots()

    bar_container = ax.barh(
        ticks, counts, color="cadetblue", edgecolor="black", alpha=0.7, align="center"
    )
    ax.bar_label(bar_container, fmt="{:,.0f}", rotation="horizontal", padding=3)
    plt.yticks(ticks, labels)
    # add enough y-label height to make sure the labels above the bar don't overlap with the title
    ax.set_xlim(0, max(counts) * 1.15)
    # Add labels and a title
    plt.ylabel(ylabel)
    plt.xlabel(x_label)
    # plt.title(title)

    # Add a grid and set image size
    plt.grid(True, linestyle="--", alpha=0.5, axis="x")

    # edge case for graph with only 1 bar, reduce the height to prevent a really thiccc bar
    image_height = 4 if len(labels) == 1 else 8.27

    plt.gcf().set_size_inches(11.69, image_height)
    plt.tight_layout()

    if file_name is not None:
        plt.savefig(file_name)
    else:
        plt.show()


def format_float_to(value: float, precision: int) -> str:
    return "{:.{precision}f}".format(value, precision=precision)


def decrease_end_bound(value: float, precision: int) -> float:
    return value - 1 / 10**precision


def protein_search_distribution(
    searchResults: list[tuple[bool, int, float]],
    bin_size: float,
    short_name: str,
    programming_language: str,
):
    found_list, proteins_size_list, search_time_list = list(zip(*searchResults))

    max_bin = (max(search_time_list) // bin_size) * bin_size
    frequency_dict: dict[float, int] = {
        key: 0 for key in np.arange(0, max_bin + bin_size, bin_size)
    }

    for val in search_time_list:
        frequency_dict[(val // bin_size) * bin_size] += 1

    precision = (
        0 if bin_size >= 1 else len(str(bin_size).split(".")[-1])
    )  # calculate number of digits after the . in the bin size
    labels, counts = list(zip(*frequency_dict.items()))
    labels = [
        f"{format_float_to(label, precision)}-{format_float_to(decrease_end_bound(label + bin_size, precision), precision)}"
        for label in labels
    ]

    create_barh(
        labels,
        counts,
        "Search time in ms",
        "Number of Occurrences",
        f"Distribution of search time {short_name} ({programming_language})",
    )


def create_search_time_distribution(searchResults: list[tuple[bool, int, float]]):
    detail_range = list(range(5, 11))
    min_detailed_value = detail_range[0]
    max_detailed_value = detail_range[-1]
    list_num_occurrences = [0 for _ in range(len(detail_range) + 1)]
    list_total_time_for_group = [0 for _ in range(len(detail_range) + 1)]
    for found, prot_size, search_time in searchResults:
        # determine the group
        index = -1
        if prot_size <= max_detailed_value:
            index = prot_size - min_detailed_value
        # add value to the appropriate group
        list_num_occurrences[index] += 1
        list_total_time_for_group[index] += search_time

    for index, occ in enumerate(list_num_occurrences):
        list_total_time_for_group[index] /= occ

    labels = [str(val) for val in detail_range]
    labels.append(f"{labels[-1]}+")
    create_barh(
        labels=labels,
        counts=list_total_time_for_group,
        x_label="Time (ms)",
        ylabel="peptide length",
        title="test",
    )

    create_barh(
        labels=labels,
        counts=list_num_occurrences,
        x_label="Number of occurrences",
        ylabel="peptide length",
        title="test",
    )


@dataclass
class ProteinLengthGraphSettings:
    max_allowed_length: float
    bin_size: int


@dataclass
class GraphInfoForfile:
    """Basic data class that keeps the data that is needed for creating all the needed graphs for a file together"""

    short_name: str
    proteinLengthGraphSettings: list[ProteinLengthGraphSettings]
    input_file_location: str
    sequence_type: str


@dataclass
class ConfigurationForSearchTimeGraphs:
    bin_size: float
    input_file_location: str
    short_name: str
    language: str


def heatmap(
    data,
    row_labels,
    col_labels,
    short_name: str,
    programming_language: str,
    ax=None,
    cbar_kw=None,
    cbarlabel="",
    **kwargs,
):
    if ax is None:
        ax = plt.gca()

    if cbar_kw is None:
        cbar_kw = {}

    # Plot the heatmap
    im = ax.imshow(data, **kwargs)

    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
    cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")

    # Show all ticks and label them with the respective list entries.
    ax.set_xticks(np.arange(data.shape[1]), labels=col_labels)
    ax.set_yticks(np.arange(data.shape[0]), labels=row_labels)
    ax.set_ylabel("Peptide length", rotation=90)
    ax.set_xlabel("Search time in ms")
    ax.set_title(
        f"Heatmap of the needed search time \n {short_name} \n in combination with the peptide length ({programming_language})"
    )

    # Let the horizontal axes labeling appear on top.
    # ax.tick_params(top=True, bottom=False,
    #                labeltop=True, labelbottom=False)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")

    # Turn spines off and create white grid.
    ax.spines[:].set_visible(False)

    ax.set_xticks(np.arange(data.shape[1] + 1) - 0.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0] + 1) - 0.5, minor=True)
    ax.grid(which="minor", color="w", linestyle="-", linewidth=3)
    ax.tick_params(which="minor", bottom=False, left=False)

    return im, cbar


def annotate_heatmap(
    im,
    data=None,
    valfmt="{x:.2f}",
    textcolors=("black", "white"),
    threshold=None,
    **textkw,
):
    if not isinstance(data, (list, np.ndarray)):
        data = im.get_array()

    # Normalize the threshold to the images color range.
    if threshold is not None:
        threshold = im.norm(threshold)
    else:
        threshold = im.norm(data.max()) / 2.0

    # Set default alignment to center, but allow it to be
    # overwritten by textkw.
    kw = dict(horizontalalignment="center", verticalalignment="center")
    kw.update(textkw)

    # Get the formatter in case a string is supplied
    if isinstance(valfmt, str):
        valfmt = ticker.StrMethodFormatter(valfmt)

    # Loop over the data and create a `Text` for each "pixel".
    # Change the text's color depending on the data.
    texts = []
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            kw.update(color=textcolors[int(im.norm(data[i, j]) > threshold)])
            text = im.axes.text(j, i, valfmt(data[i, j], None), **kw)
            texts.append(text)

    return texts


def create_heatmap_of_search_execution(
    search_results: list[tuple[bool, int, float]],
    short_name: str,
    programming_language: str,
):
    fig, ax = plt.subplots()

    num_x_bins = 5
    num_y_bins = 5
    max_protein_length = max(
        map(lambda l: l[1], search_results)
    )  # to be just a bit bigger than the actual max value, is useful during mapping to the bins
    max_protein_length += max_protein_length * 0.001
    max_protein_length = round_up(max_protein_length, 10)
    max_search_time = max(map(lambda l: l[2], search_results))
    precision = 0
    # edge case where every sequence is found in 1 ms, we only want 1 bin in the time dimension in that case
    if max_search_time >= 1:
        max_search_time += max_search_time * 0.001
        max_search_time = round_up(max_search_time, 50)
    else:
        num_digits = 0
        max_copy = max_search_time
        while max_copy < 1:
            max_copy *= 10
            num_digits += 1
        factor = 5
        if max_copy >= 5:
            factor = 10

        max_search_time = factor * 1 / 10**num_digits
        precision = num_digits + 1

    x_labels = []  # searchTimes
    y_labels = []  # protein length
    binned_data = []

    # create empty binned_data matrix and fill in the axis labels
    for i in range(num_x_bins):
        x_labels.append(
            f"{format_float_to(i * max_search_time / num_x_bins, precision)}-{format_float_to(decrease_end_bound(((i + 1) * max_search_time / num_x_bins), precision), precision)}"
        )
    for i in range(num_y_bins):
        binned_data.append([0 for _ in range(num_x_bins)])
        y_labels.append(
            f"{i * max_protein_length // num_y_bins}-{((i + 1) * max_protein_length // num_y_bins) - 1}"
        )

    # fill in the binned_data matrix
    for _, size, search_time in search_results:
        y = floor(search_time / max_search_time * num_x_bins)
        x = floor(size / max_protein_length * num_y_bins)
        binned_data[x][y] += 1

    im, cbar = heatmap(
        np.array(binned_data),
        y_labels,
        x_labels,
        short_name,
        programming_language,
        ax=ax,
        cmap="Blues",
        cbarlabel="number of occurrences",
        cbar_kw={
            "shrink": 0.75
        },  # add shrink to make vertical height of the bar a bit smaller
    )
    texts = annotate_heatmap(im, valfmt="{x:,.0f}")

    plt.gcf().set_size_inches(8.27, 8.27)
    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    inputForGraphs = [
        GraphInfoForfile(
            "de Human-Prot databank",
            [
                ProteinLengthGraphSettings(float("inf"), 1000),
                ProteinLengthGraphSettings(1000, 50),
            ],
            "/Users/brdvlami/Documents/Ugent/MA2/Thesis/Dataset/BenchmarkData/humanprot/protein_database.tsv",
            "Proteïne",
        ),
        GraphInfoForfile(
            "de Swiss-Prot databank",
            [
                ProteinLengthGraphSettings(float("inf"), 1000),
                ProteinLengthGraphSettings(1000, 50),
            ],
            "/Users/brdvlami/Documents/Ugent/MA2/Thesis/Dataset/BenchmarkData/swissprot_var1/protein_database.tsv",
            "Proteïne",
        ),
        GraphInfoForfile(
            "het Human-Prot peptidebestand",
            [ProteinLengthGraphSettings(float("inf"), 5)],
            "/Users/brdvlami/Documents/Ugent/MA2/Thesis/Dataset/BenchmarkData/humanprot/search_file.tsv",
            "Peptide",
        ),
        GraphInfoForfile(
            "het Swiss-Prot peptidebestand met missed cleavage",
            [ProteinLengthGraphSettings(float("inf"), 5)],
            "/Users/brdvlami/Documents/Ugent/MA2/Thesis/Dataset/BenchmarkData/swissprot_var2/search_file_mch.tsv",
            "Peptide",
        ),
        GraphInfoForfile(
            "het Swiss-Prot peptidebestand zonder missed cleavage",
            [ProteinLengthGraphSettings(float("inf"), 5)],
            "/Users/brdvlami/Documents/Ugent/MA2/Thesis/Dataset/BenchmarkData/swissprot_var1/search_file_no_mch.tsv",
            "Peptide",
        ),
        GraphInfoForfile(
            "het SIHUMI SO3 peptidebestand",
            [ProteinLengthGraphSettings(float("inf"), 5)],
            "/Users/brdvlami/Documents/Ugent/MA2/Thesis/Dataset/BenchmarkData/SIHUMI/S03.txt",
            "Peptide",
        ),
        GraphInfoForfile(
            "het SIHUMI SO5 peptidebestand",
            [ProteinLengthGraphSettings(float("inf"), 5)],
            "/Users/brdvlami/Documents/Ugent/MA2/Thesis/Dataset/BenchmarkData/SIHUMI/S05.txt",
            "Peptide",
        ),
        GraphInfoForfile(
            "het SIHUMI SO7 peptidebestand",
            [ProteinLengthGraphSettings(float("inf"), 5)],
            "/Users/brdvlami/Documents/Ugent/MA2/Thesis/Dataset/BenchmarkData/SIHUMI/S07.txt",
            "Peptide",
        ),
        GraphInfoForfile(
            "het SIHUMI SO8 peptidebestand",
            [ProteinLengthGraphSettings(float("inf"), 5)],
            "/Users/brdvlami/Documents/Ugent/MA2/Thesis/Dataset/BenchmarkData/SIHUMI/S08.txt",
            "Peptide",
        ),
        GraphInfoForfile(
            "het SIHUMI S11 peptidebestand",
            [ProteinLengthGraphSettings(float("inf"), 5)],
            "/Users/brdvlami/Documents/Ugent/MA2/Thesis/Dataset/BenchmarkData/SIHUMI/S11.txt",
            "Peptide",
        ),
        GraphInfoForfile(
            "het SIHUMI S14 peptidebestand",
            [ProteinLengthGraphSettings(float("inf"), 5)],
            "/Users/brdvlami/Documents/Ugent/MA2/Thesis/Dataset/BenchmarkData/SIHUMI/S14.txt",
            "Peptide",
        ),
    ]

    # for graphInput in inputForGraphs:
    #     with open(graphInput.input_file_location) as fp:
    #         data = []
    #         for line in fp:
    #             data.append(line.split("\t")[-1].rstrip("\n"))
    #         print(graphInput.input_file_location)
    #         create_barh_acid_code_occurrences(data, graphInput.short_name)
    #         # do this in a loop since we can provide mutliple max-value and bin_size combinations to provide a better view
    #         for proteinLengthGraphSetting in graphInput.proteinLengthGraphSettings:
    #             create_barh_protein_length(
    #                 data,
    #                 proteinLengthGraphSetting.bin_size,
    #                 graphInput.short_name,
    #                 proteinLengthGraphSetting.max_allowed_length,
    #                 graphInput.sequence_type,
    #             )

    search_time_graph_configurations = [
        ConfigurationForSearchTimeGraphs(
            0.005,
            "/Users/brdvlami/Documents/Ugent/MA2/Thesis/BenchmarkResults/CPP_ukkonen_CCB/SearchTime/Hit/search_cpp_onlyMatch_immunopeptidomics_avg10.txt",
            "until match for the Human-Prot search file",
            "C++",
        ),
        ConfigurationForSearchTimeGraphs(
            0.005,
            "/Users/brdvlami/Documents/Ugent/MA2/Thesis/BenchmarkResults/CPP_ukkonen_CCB/SearchTime/Hit/search_cpp_onlyMatch_swissprot_var1_avg10.txt",
            "until match for the Swiss-Prot search file without missed cleavage",
            "C++",
        ),
        ConfigurationForSearchTimeGraphs(
            0.0005,
            "/Users/brdvlami/Documents/Ugent/MA2/Thesis/BenchmarkResults/CPP_ukkonen_CCB/SearchTime/Hit/search_cpp_onlyMatch_swissprot_var2_avg10.txt",
            "until match for the Swiss-Prot search file with missed cleavage",
            "C++",
        ),
        ConfigurationForSearchTimeGraphs(
            100,
            "/Users/brdvlami/Documents/Ugent/MA2/Thesis/BenchmarkResults/CPP_ukkonen_CCB/SearchTime/HitAndSearch/search_cpp_matchAndSearchTree_immunopeptidomics_avg10.txt",
            "with traversal of subtree for the Human-Prot search file",
            "C++",
        ),
        ConfigurationForSearchTimeGraphs(
            50,
            "/Users/brdvlami/Documents/Ugent/MA2/Thesis/BenchmarkResults/CPP_ukkonen_CCB/SearchTime/HitAndSearch/search_cpp_matchAndSearchTree_swissprot_var1_avg10.txt",
            "with traversal of subtree for the Swiss-Prot search file without missed cleavage",
            "C++",
        ),
        ConfigurationForSearchTimeGraphs(
            50,
            "/Users/brdvlami/Documents/Ugent/MA2/Thesis/BenchmarkResults/CPP_ukkonen_CCB/SearchTime/HitAndSearch/search_cpp_matchAndSearchTree_swissprot_var2_avg10.txt",
            "with traversal of subtree for the Swiss-Prot search file with missed cleavage",
            "C++",
        ),
        ConfigurationForSearchTimeGraphs(
            0.005,
            "/Users/brdvlami/Documents/Ugent/MA2/Thesis/BenchmarkResults/Rust_ukkonen/SearchTime/Hit/search_rust_onlyMatch_immunopeptidomics_avg10.txt",
            "until match for the Human-Prot search file",
            "Rust",
        ),
        ConfigurationForSearchTimeGraphs(
            0.0005,
            "/Users/brdvlami/Documents/Ugent/MA2/Thesis/BenchmarkResults/Rust_ukkonen/SearchTime/Hit/search_rust_onlyMatch_swissprot_var1_avg10.txt",
            "until match for the Swiss-Prot search file without missed cleavage",
            "Rust",
        ),
        ConfigurationForSearchTimeGraphs(
            0.0005,
            "/Users/brdvlami/Documents/Ugent/MA2/Thesis/BenchmarkResults/Rust_ukkonen/SearchTime/Hit/search_rust_onlyMatch_swissprot_var2_avg10.txt",
            "until match for the Swiss-Prot search file with missed cleavage",
            "Rust",
        ),
        ConfigurationForSearchTimeGraphs(
            100,
            "/Users/brdvlami/Documents/Ugent/MA2/Thesis/BenchmarkResults/Rust_ukkonen/SearchTime/HitAndSearch/search_rust_matchAndSearchTree_immunopeptidomics_avg10.txt",
            "with traversal of subtree for the Human-Prot search file",
            "Rust",
        ),
        ConfigurationForSearchTimeGraphs(
            0.05,
            "/Users/brdvlami/Documents/Ugent/MA2/Thesis/BenchmarkResults/Rust_ukkonen/SearchTime/HitAndSearch/search_rust_matchAndSearchTree_swissprot_var1_avg10.txt",
            "with traversal of subtree for the Swiss-Prot search file without missed cleavage",
            "Rust",
        ),
        ConfigurationForSearchTimeGraphs(
            0.05,
            "/Users/brdvlami/Documents/Ugent/MA2/Thesis/BenchmarkResults/Rust_ukkonen/SearchTime/HitAndSearch/search_rust_matchAndSearchTree_swissprot_var2_avg10.txt",
            "with traversal of subtree for the Swiss-Prot search file with missed cleavage",
            "Rust",
        ),
    ]

    # for configuration in search_time_graph_configurations:
    #     with open(configuration.input_file_location) as fp:
    #         data = []
    #         for line in fp:
    #             found, protein_size, search_time = line.split(";")
    #             # only use the time information if there was a match.
    #             # If we don't match the information is not representative since the mismatch could be after e.g. 1 character => fast, but also nothing that is useful
    #             if found == "1":
    #                 search_time = search_time.rstrip("\n")
    #                 data.append((bool(found), int(protein_size), float(search_time)))
    #         protein_search_distribution(
    #             data,
    #             configuration.bin_size,
    #             configuration.short_name,
    #             configuration.language,
    #         )
    #         create_heatmap_of_search_execution(
    #             data, configuration.short_name, configuration.language
    #         )

    sparse_files = [
        "/Users/brdvlami/Documents/Ugent/MA2/Thesis/BenchmarkResults/Rust_suffix_array/SparseSuffixArrays/swissprot_var1_sample_factor1_avg10.txt",
        "/Users/brdvlami/Documents/Ugent/MA2/Thesis/BenchmarkResults/Rust_suffix_array/SparseSuffixArrays/swissprot_var1_sample_factor2_avg10.txt",
        "/Users/brdvlami/Documents/Ugent/MA2/Thesis/BenchmarkResults/Rust_suffix_array/SparseSuffixArrays/swissprot_var1_sample_factor3_avg10.txt",
        "/Users/brdvlami/Documents/Ugent/MA2/Thesis/BenchmarkResults/Rust_suffix_array/SparseSuffixArrays/swissprot_var1_sample_factor4_avg10.txt",
    ]

    for file in sparse_files:
        with open(file) as fp:
            data = []
            for line in fp:
                found, protein_size, search_time = line.split(";")
                # only use the time information if there was a match.
                # If we don't match the information is not representative since the mismatch could be after e.g. 1 character => fast, but also nothing that is useful
                if found == "1":
                    search_time = search_time.rstrip("\n")
                    data.append((bool(found), int(protein_size), float(search_time)))
            create_search_time_distribution(data)
