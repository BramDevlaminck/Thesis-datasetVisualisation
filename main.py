import matplotlib.pyplot as plt
import numpy as np


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


def create_barh_acid_code_occurrences(proteins: list[str]):
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
        'Number of Occurrences per\nAmino Acid code for the swissprot_var2 search file with missed cleavage'
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


def create_barh_protein_length(proteins: list[str], bin_size: int):
    """From the given proteins calculate their length and create a horizontal barchart with the data"""
    statistics = calculate_length_statistics(proteins)
    max_bin = (max(statistics) // bin_size) * bin_size
    frequency_dict: dict[int, int] = {key: 0 for key in range(0, max_bin + bin_size, bin_size)}

    for val in statistics:
        frequency_dict[(val // bin_size) * bin_size] += 1

    labels, counts = list(zip(*frequency_dict.items()))
    labels = [f"{label}-{label + bin_size}" for label in labels]

    create_barh(
        labels,
        counts,
        "Protein Length",
        "Number of Occurrences",
        "Distribution of protein length for the swissprot_var2 search file with missed cleavage"
    )


def create_barh(labels: list[any], counts: list[int], ylabel: str, x_label: str, title: str,
                file_name: str | None = None):
    """
    creates a horizontal barchart with the provided labels counts and axis labels

    if filen_name is not none the barchart is stored to a file with the given name
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


if __name__ == '__main__':
    filename = "/Users/brdvlami/Documents/Ugent/MA2/Thesis/Dataset/BenchmarkData/swissprot_var2/search_file_mch.tsv"
    with open(filename) as fp:
        data = []
        for line in fp:
            data.append(line.split("\t")[-1].rstrip('\n'))

        create_barh_acid_code_occurrences(data)
        create_barh_protein_length(data, 5)
