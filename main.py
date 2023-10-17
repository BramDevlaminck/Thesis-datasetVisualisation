import matplotlib.pyplot as plt
import numpy as np


def calculate_character_statistics(data: list[str]) -> list[str]:
    statistiscs: list[str] = []
    for line in data:
        statistiscs += [*line]

    return statistiscs


def calculate_length_statistics(data: list[str]) -> dict[int, int]:
    statistiscs: dict[int, int] = dict()
    for line in data:
        size = len(line)
        current_count = statistiscs.get(size, 0)
        statistiscs[size] = current_count + 1

    return statistiscs


def make_alphabet_complete(data: tuple[list[str], list[int]]) -> tuple[list[str], list[int]]:
    alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"

    labels, counts = data
    for i, letter in enumerate(alphabet):
        if i >= len(labels) or labels[i] != letter:
            labels.insert(i, letter)
            counts.insert(i, 0)

    return labels, counts


def create_histogram_horizontal(statistics: list[any]):
    labels, counts = np.unique(statistics, return_counts=True)
    labels, counts = make_alphabet_complete((list(labels), list(counts)))
    reversed(labels)
    ticks = range(len(counts), 0, -1)
    fig, ax = plt.subplots()

    bar_container = ax.barh(ticks, counts, color='cadetblue', edgecolor='black', alpha=0.7, align="center")
    ax.bar_label(bar_container, fmt='{:,.0f}', rotation='horizontal', padding=3)
    plt.yticks(ticks, labels)
    # add enough y-label hight to make sure the labels above the bar don't overlap with the title
    ax.set_xlim(0, max(counts) * 1.15)
    # Add labels and a title
    plt.ylabel('Amino Acid Code')
    plt.xlabel('Number of Occurrences')
    plt.title('Number of Occurrences per\nAmino Acid code for the swissprot_var2 search file with missed cleavage')

    # Add a grid
    plt.grid(True, linestyle='--', alpha=0.5, axis='x')
    plt.gcf().set_size_inches(11.69, 8.27)

    # Show the plot
    # plt.show()
    plt.savefig("Graphs/character_occurrences_swissprot_var2_search_file_mch_horizontal.jpg")


def create_histogram_of_dict(statistics: list[any]):
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


if __name__ == '__main__':
    filename = "/Users/brdvlami/Documents/Ugent/MA2/Thesis/Dataset/BenchmarkData/swissprot_var2/search_file_mch.tsv"
    with open(filename) as fp:
        data = []
        for line in fp:
            data.append(line.split("\t")[-1].rstrip('\n'))

        statistics = calculate_character_statistics(data)
        create_histogram_horizontal(statistics)
