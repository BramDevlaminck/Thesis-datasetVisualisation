from dataclasses import dataclass
from math import floor

import matplotlib.pyplot as plt
from matplotlib import ticker
import numpy as np


def round_up(num_to_round: int, multiple_to_round_to: int) -> int:
    """Aid function to round up to a certain multiple"""
    if multiple_to_round_to == 0:
        return num_to_round

    remainder = num_to_round % multiple_to_round_to
    if remainder == 0:
        return num_to_round

    return round(num_to_round + multiple_to_round_to - remainder)


def calculate_character_statistics(proteins: list[str]) -> list[str]:
    """Creates a list of all the individual characters of the proteins"""
    statistiscs: list[str] = []
    for line in proteins:
        statistiscs += [*line]

    return statistiscs


def calculate_length_statistics(proteins: list[str]) -> list[int]:
    """Creates a list containing the lengths off all the proteins"""
    statistiscs: list[int] = []
    for line in proteins:
        statistiscs.append(len(line))

    return statistiscs


def make_alphabet_complete(data: tuple[list[str], list[int]]) -> tuple[list[str], list[int]]:
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


def create_barh_acid_code_occurrences(proteins: list[str], short_name: str, output_file: str | None = None):
    """Creates a horizontal barchart of the distribution of the letters in the proteins"""

    statistics = calculate_character_statistics(proteins)
    # counts for char in statistics how many times it occurred, the char will be places in labels (sorted alphabetically) and the count in counts
    labels, counts = np.unique(statistics, return_counts=True)
    # some letters from the alphabet will not occur in this dataset, add them at the right place with count 0
    labels, counts = make_alphabet_complete((list(labels), list(counts)))

    create_barh(
        labels,
        counts,
        'Amino Acid Code',
        'Number of Occurrences',
        f'Number of Occurrences per\nAmino Acid code for the {short_name}',
        output_file
    )


def create_barchart_vertical(statistics: list[any]):
    labels, counts = np.unique(statistics, return_counts=True)
    labels, counts = make_alphabet_complete((list(labels), list(counts)))
    ticks = range(len(counts))
    fig, ax = plt.subplots()

    bar_container = ax.bar(ticks, counts, color='cadetblue', edgecolor='black', alpha=0.7, align="center")
    ax.bar_label(bar_container, fmt='{:,.0f}', rotation='vertical', padding=3)
    plt.xticks(ticks, labels)
    # add enough y-label hight to make sure the labels above the bar don't overlap with the title
    ax.set_ylim(0, max(counts) * 1.2)
    # Add labels and a title
    plt.xlabel('Amino Acid Code')
    plt.ylabel('Number of Occurrences')
    plt.title('Number of Occurrences per\nAmino Acid code for the swissprot_var1 dataset')

    # Add a grid
    plt.grid(True, linestyle='--', alpha=0.5, axis='y')
    plt.gcf().set_size_inches(11.69, 8.27)

    # Show the plot
    # plt.show()
    plt.savefig("Graphs/character_occurrences_swissprot_var1_protein_database.jpg")


def create_barh_protein_length(proteins: list[str], bin_size: int, short_name: str, max_allowed_protein_length: float,
                               output_file: str | None = None):
    """From the given proteins calculate their length and create a horizontal barchart with the data"""
    statistics = calculate_length_statistics(proteins)
    max_bin = int(min(max_allowed_protein_length - bin_size, (max(statistics) // bin_size) * bin_size))
    frequency_dict: dict[int, int] = {key: 0 for key in range(0, max_bin + bin_size, bin_size)}

    for val in statistics:
        if val < max_allowed_protein_length:
            frequency_dict[(val // bin_size) * bin_size] += 1

    labels, counts = list(zip(*frequency_dict.items()))
    labels = [f"{label}-{label + bin_size}" for label in labels]

    create_barh(
        labels,
        counts,
        "Protein Length",
        "Number of Occurrences",
        f"Distribution of protein length for the {short_name}",
        output_file
    )


def create_barh(labels: list[any], counts: list[int], ylabel: str, x_label: str, title: str,
                file_name: str | None = None):
    """
    creates a horizontal barchart with the provided labels counts and axis labels

    if file_name is not none the barchart is stored to a file with the given name
    """
    reversed(labels)
    ticks = range(len(counts), 0, -1)
    fig, ax = plt.subplots()

    bar_container = ax.barh(ticks, counts, color='cadetblue', edgecolor='black', alpha=0.7, align="center")
    ax.bar_label(bar_container, fmt='{:,.0f}', rotation='horizontal', padding=3)
    plt.yticks(ticks, labels)
    # add enough y-label height to make sure the labels above the bar don't overlap with the title
    ax.set_xlim(0, max(counts) * 1.15)
    # Add labels and a title
    plt.ylabel(ylabel)
    plt.xlabel(x_label)
    plt.title(title)

    # Add a grid and set image size
    plt.grid(True, linestyle='--', alpha=0.5, axis='x')
    plt.gcf().set_size_inches(11.69, 8.27)

    if file_name is not None:
        plt.savefig(file_name)
    else:
        plt.show()


def protein_search_distribution(searchResults: list[tuple[bool, int, float]], bin_size: int, short_name: str):
    found_list, proteins_size_list, search_time_list = list(zip(*searchResults))

    max_bin = int((max(search_time_list) // bin_size) * bin_size)
    frequency_dict: dict[int, int] = {key: 0 for key in range(0, max_bin + bin_size, bin_size)}

    for val in search_time_list:
        frequency_dict[(val // bin_size) * bin_size] += 1

    labels, counts = list(zip(*frequency_dict.items()))
    labels = [f"{label}-{label + bin_size}" for label in labels]

    create_barh(
        labels,
        counts,
        "Search time in ms",
        "Number of Occurrences",
        f"Distribution of search time for the {short_name}"
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


@dataclass
class ConfigurationForSearchTimeGraphs:
    bin_size: int
    input_file_location: str
    short_name: str


def heatmap(data, row_labels, col_labels, short_name: str, ax=None,
            cbar_kw=None, cbarlabel="", **kwargs):
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
    ax.set_ylabel("Protein length", rotation=90)
    ax.set_xlabel("Search time in ms")
    ax.set_title(f"Heatmap of the needed search time in combination with the protein lengths\n for {short_name}")

    # Let the horizontal axes labeling appear on top.
    # ax.tick_params(top=True, bottom=False,
    #                labeltop=True, labelbottom=False)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
             rotation_mode="anchor")

    # Turn spines off and create white grid.
    ax.spines[:].set_visible(False)

    ax.set_xticks(np.arange(data.shape[1] + 1) - .5, minor=True)
    ax.set_yticks(np.arange(data.shape[0] + 1) - .5, minor=True)
    ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
    ax.tick_params(which="minor", bottom=False, left=False)

    return im, cbar


def annotate_heatmap(im, data=None, valfmt="{x:.2f}",
                     textcolors=("black", "white"),
                     threshold=None, **textkw):
    if not isinstance(data, (list, np.ndarray)):
        data = im.get_array()

    # Normalize the threshold to the images color range.
    if threshold is not None:
        threshold = im.norm(threshold)
    else:
        threshold = im.norm(data.max()) / 2.

    # Set default alignment to center, but allow it to be
    # overwritten by textkw.
    kw = dict(horizontalalignment="center",
              verticalalignment="center")
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


def create_heatmap_of_search_execution(search_results: list[tuple[bool, int, float]], short_name: str):
    fig, ax = plt.subplots()

    num_x_bins = 5
    num_y_bins = 5
    max_protein_length = max(map(lambda l: l[1],
                                 search_results))  # to be just a bit bigger than the actual max value, is useful during mapping to the bins
    max_protein_length += max_protein_length * 0.001
    max_protein_length = round_up(max_protein_length, 10)
    max_search_time = max(map(lambda l: l[2], search_results))
    # edge case where every sequence is found in 1 ms, we only want 1 bin in the time dimension in that case
    if max_search_time >= 1:
        max_search_time += max_search_time * 0.001
        max_search_time = round_up(max_search_time, 50)
    else:
        num_x_bins = 1
        max_search_time = 1

    x_labels = []  # searchTimes
    y_labels = []  # protein length
    binned_data = []

    # create empty binned_data matrix and fill in the axis labels
    for i in range(num_x_bins):
        x_labels.append(f"{i * max_search_time // num_x_bins}-{(i + 1) * max_search_time // num_x_bins}")
    for i in range(num_y_bins):
        binned_data.append([0 for _ in range(num_x_bins)])
        y_labels.append(f"{i * max_protein_length // num_y_bins}-{(i + 1) * max_protein_length // num_y_bins}")

    # fill in the binned_data matrix
    for (_, size, search_time) in search_results:
        y = floor(search_time / max_search_time * num_x_bins)
        x = floor(size / max_protein_length * num_y_bins)
        binned_data[x][y] += 1

    im, cbar = heatmap(
        np.array(binned_data),
        y_labels, x_labels,
        short_name,
        ax=ax,
        cmap="Blues",
        cbarlabel="number of occurrences",
        cbar_kw={"shrink": 0.9} # add shrink to make vertical height of the bar a bit smaller
    )
    texts = annotate_heatmap(im, valfmt="{x:,.0f}")

    plt.gcf().set_size_inches(8.27, 7.27)
    fig.tight_layout()
    plt.show()


if __name__ == '__main__':
    inputForGraphs = [
        GraphInfoForfile(
            "immunopeptidomics database",
            [
                ProteinLengthGraphSettings(float('inf'), 1000),
                ProteinLengthGraphSettings(1000, 50)
            ],
            "/Users/brdvlami/Documents/Ugent/MA2/Thesis/Dataset/BenchmarkData/immunopeptidomics/protein_database.tsv"),
        GraphInfoForfile(
            "swissprot database",
            [
                ProteinLengthGraphSettings(float('inf'), 1000),
                ProteinLengthGraphSettings(1000, 50)
            ],
            "/Users/brdvlami/Documents/Ugent/MA2/Thesis/Dataset/BenchmarkData/swissprot_var1/protein_database.tsv"),
        GraphInfoForfile(
            "immunopeptidomics searchfile",
            [ProteinLengthGraphSettings(float('inf'), 5)],
            "/Users/brdvlami/Documents/Ugent/MA2/Thesis/Dataset/BenchmarkData/immunopeptidomics/search_file.tsv"),
        GraphInfoForfile(
            "swissprot searchfile without missed cleavage (var 2)",
            [ProteinLengthGraphSettings(float('inf'), 5)],
            "/Users/brdvlami/Documents/Ugent/MA2/Thesis/Dataset/BenchmarkData/swissprot_var2/search_file_mch.tsv"),
        GraphInfoForfile(
            "swissprot searchfile with missed cleavage (var 1)",
            [ProteinLengthGraphSettings(float('inf'), 5)],
            "/Users/brdvlami/Documents/Ugent/MA2/Thesis/Dataset/BenchmarkData/swissprot_var1/search_file_no_mch.tsv")
    ]

    for graphInput in inputForGraphs:
        with open(graphInput.input_file_location) as fp:
            data = []
            for line in fp:
                data.append(line.split("\t")[-1].rstrip('\n'))

            create_barh_acid_code_occurrences(data, graphInput.short_name)
            # do this in a loop since we can provide mutliple max-value and bin_size combinations to provide a better view
            for proteinLengthGraphSetting in graphInput.proteinLengthGraphSettings:
                create_barh_protein_length(data, proteinLengthGraphSetting.bin_size, graphInput.short_name,
                                           proteinLengthGraphSetting.max_allowed_length)

    search_time_graph_configurations = [
        ConfigurationForSearchTimeGraphs(
            1,
             "/Users/brdvlami/Documents/Ugent/MA2/Thesis/BenchmarkResults/CPP_ukkonen_CCB/SearchTime/Hit/search_cpp_onlyMatch_immunopeptidomics_avg10.txt",
             "immunopeptidomics search file for only checking if a match exists",
        ),
        ConfigurationForSearchTimeGraphs(
            1,
            "/Users/brdvlami/Documents/Ugent/MA2/Thesis/BenchmarkResults/CPP_ukkonen_CCB/SearchTime/Hit/search_cpp_onlyMatch_swissprot_var1_avg10.txt",
            "swissprot with missed cleavage (var1) for only checking if a match exists",
        ),
        ConfigurationForSearchTimeGraphs(
            1,
            "/Users/brdvlami/Documents/Ugent/MA2/Thesis/BenchmarkResults/CPP_ukkonen_CCB/SearchTime/Hit/search_cpp_onlyMatch_swissprot_var2_avg10.txt",
            "swissprot without missed cleavage (var2) for only checking if a match exists",
        ),
        ConfigurationForSearchTimeGraphs(
            100,
            "/Users/brdvlami/Documents/Ugent/MA2/Thesis/BenchmarkResults/CPP_ukkonen_CCB/SearchTime/HitAndSearch/search_cpp_matchAndSearchTree_immunopeptidomics_avg10.txt",
            "immunopeptidomics search file with traversal of the whole subtree after match",
        ),
        ConfigurationForSearchTimeGraphs(
            50,
            "/Users/brdvlami/Documents/Ugent/MA2/Thesis/BenchmarkResults/CPP_ukkonen_CCB/SearchTime/HitAndSearch/search_cpp_matchAndSearchTree_swissprot_var1_avg10.txt",
            "swissprot with missed cleavage (var1) with traversal of the whole subtree after match",
        )
    ]

    for configuration in search_time_graph_configurations:
        with open(configuration.input_file_location) as fp:
            data = []
            for line in fp:
                found, protein_size, search_time = line.split(";")
                search_time = search_time.rstrip("\n")
                data.append((bool(found), int(protein_size), float(search_time)))
            protein_search_distribution(data, configuration.bin_size, configuration.short_name)
            create_heatmap_of_search_execution(data, configuration.short_name)
