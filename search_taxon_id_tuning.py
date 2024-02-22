from math import ceil

import matplotlib.pyplot as plt
import numpy as np


def create_2d_heatmap(data: list[tuple[int, int, float]]):
    max_shown_prot_size = 54
    max_matches = max(v[1] for v in data)
    occ_bucket_size = int(ceil(max_matches / 100))  # we want 100 buckets
    image = [[0 for _ in range(100)] for _ in range(50)]

    # TODO: count for each interval in the square how many times it occurs and show it as an image
    for size, occ, _ in data:
        if size <= max_shown_prot_size:
            size_bucket = size - 5  # - 5 since the min value is 5
            occ_bucket = occ // occ_bucket_size
            image[size_bucket][occ_bucket] += 1
    print(max_matches)

    fig, ax = plt.subplots()

    # set 0 values to NaN to be displayed as white in the image
    for i, row in enumerate(image):
        for j, val in enumerate(row):
            if val == 0:
                image[i][j] = np.NAN

    hist = ax.imshow(
        image,
        interpolation="none",
        cmap="cividis",
    )
    x_tick_labels = [0, 500000, 1000000, 1500000, 2000000]
    x_ticks = [label / max_matches * 100 for label in x_tick_labels]
    ax.set_xticks(
        x_ticks, labels=[f"{label:,}".replace(",", ".") for label in x_tick_labels]
    )
    y_ticks = [0, 10, 20, 30, 40]
    ax.set_yticks(y_ticks, labels=[str(v + 5) for v in y_ticks])
    fig.colorbar(hist)
    ax.set_ylabel("Peptide length")
    ax.set_xlabel("Number of matching proteins")
    ax.set_title(
        "Distribution of number of matches in Uniprot\nfor the Swiss-Prot peptide file\n without missed cleavages"
    )
    plt.show()


def create_2d_heatmap_searchtime(data: list[tuple[int, int, float]]):
    data_filtered = []
    for length, id, time in data:
        if time is not None:
            data_filtered.append((length, id, time))

    max_shown_prot_size = 54
    max_time = max(v[2] for v in data_filtered) * 1000000
    time_bucket_size = int(ceil(max_time / 100))  # we want 100 buckets
    image = [[0 for _ in range(100)] for _ in range(50)]

    for size, _, time in data_filtered:
        if size <= max_shown_prot_size:
            size_bucket = size - 5  # - 5 since the min value is 5
            time_bucket = int((time * 1000000) // time_bucket_size)
            image[size_bucket][time_bucket] += 1
    print(max_time)
    fig, ax = plt.subplots()

    # set 0 values to NaN to be displayed as white in the image
    for i, row in enumerate(image):
        for j, val in enumerate(row):
            if val == 0:
                image[i][j] = np.NAN

    hist = ax.imshow(
        image,
        interpolation="none",
        cmap="cividis",
    )
    x_tick_labels = [0, 50000, 100000, 150000, 200000]
    x_ticks = [label / (max_time / 1000000) * 100 for label in x_tick_labels]
    ax.set_xticks(
        x_ticks, labels=[f"{label:,}".replace(",", ".") for label in x_tick_labels]
    )
    y_ticks = [0, 10, 20, 30, 40]
    ax.set_yticks(y_ticks, labels=[str(v + 5) for v in y_ticks])
    fig.colorbar(hist)
    ax.set_ylabel("Peptide length")
    ax.set_xlabel("search time (ms)")
    ax.set_title(
        "Distribution of number of search time in Uniprot\nfor the Swiss-Prot peptide file\n without missed cleavages"
    )
    plt.show()


def plot_occurrences_distribution(data: list[tuple[int, int, int, float]]):
    filter_cond = lambda v: v[3] is not None and v[1] < 10000

    filtered_data = filter(filter_cond, data)
    filtered_data_copy_iter = filter(
        filter_cond, data
    )  # copy needed since filtered_data is an iterator and we need it twice
    occ_distribution = [val[1] for val in filtered_data]
    binwidth = 1000
    plt.hist(
        occ_distribution,
        bins=range(min(occ_distribution), max(occ_distribution) + binwidth, binwidth),
        log=True,
        label="All entries",
        color="cadetblue",
        # edgecolor="black",
        alpha=0.7,
    )
    distribution_root = []
    for _, num_matches, taxon_id, _ in filtered_data_copy_iter:
        if taxon_id == 1:
            distribution_root.append(num_matches)
            print("found!")
    plt.hist(
        distribution_root,
        bins=range(min(occ_distribution), max(occ_distribution) + binwidth, binwidth),
        log=True,
        label="Entries with root as LCA (aggregation method: LCA)",
        color="blanchedalmond",
        # edgecolor="black",
        # alpha=0.7,
    )
    plt.title("Logarithmic distribution of number of matching proteins")
    plt.xlabel("Number of matching proteins")
    plt.ylabel("Number of occurrences")
    plt.legend()
    plt.show()


if __name__ == "__main__":
    data_occ_count = []
    with open(
        "/Users/brdvlami/Documents/Ugent/MA2/Thesis/BenchmarkResults/Rust_suffix_array/OccurrencesPerLengthUniprot/uniprot_searchtime_all_occ_count.txt"
    ) as fp:
        for line in fp:
            protein_size, num_matches, search_time = line.rstrip("\n").split(";")
            protein_size = int(protein_size)
            num_matches = int(num_matches)
            search_time = float(search_time)
            data_occ_count.append((protein_size, num_matches, search_time))

        create_2d_heatmap(data_occ_count)

    data_taxon_id = []
    with open(
        "/Users/brdvlami/Documents/Ugent/MA2/Thesis/BenchmarkResults/Rust_suffix_array/TaxonIdUniprot/uniprot_searchtime_taxon_id.txt"
    ) as fp:
        for line in fp:
            protein_size, taxon_id, search_time = line.rstrip("\n").split(";")
            protein_size = int(protein_size)
            if taxon_id == "" or taxon_id == "/":
                taxon_id = None
            else:
                taxon_id = int(taxon_id)
            if search_time == "":
                search_time = None
            else:
                search_time = float(search_time)

            data_taxon_id.append((protein_size, taxon_id, search_time))
        create_2d_heatmap_searchtime(data_taxon_id)

    combined_data = []
    for (prot_size1, num_matches, search_time_occ), (
        prot_size2,
        taxon_id,
        search_time_taxonid,
    ) in zip(data_occ_count, data_taxon_id):
        assert prot_size1 == prot_size2
        combined_data.append((prot_size1, num_matches, taxon_id, search_time_taxonid))

    plot_occurrences_distribution(combined_data)

    print("------------------")
    filtered = []
    total_time = 0
    for size, matches, tax_id, search_time_taxonid in combined_data:
        if search_time_taxonid is not None and search_time_taxonid > 1000:
            filtered.append((size, matches, tax_id, search_time_taxonid))
            total_time += search_time_taxonid

    filtered.sort(key=lambda x: x[3], reverse=True)
    for line in filtered:
        print(line)

    print(total_time / 1000)

    print(len(filtered))

    print("------------------")
    filtered = []
    total_time = 0
    for size, matches, tax_id, search_time_taxonid in combined_data:
        if matches is not None and matches > 4000 and search_time_taxonid is not None:
            filtered.append((size, matches, tax_id, search_time_taxonid))
            total_time += search_time_taxonid

    filtered.sort(key=lambda x: x[1], reverse=True)
    for line in filtered:
        print(line)

    print(total_time / 1000)

    print(len(filtered))
