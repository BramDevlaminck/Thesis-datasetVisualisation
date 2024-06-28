import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import plotly.express as px
import squarify

font = {"weight": "normal", "size": 15}

mpl.rc("font", **font)

colours = [
    "cadetblue",
    "blanchedalmond",
    "RosyBrown",
    "darkseagreen",
    "SlateGray",
    "lightcoral",
    "thistle",
]


def create_pie_chart(data, labels):
    total = sum(data)
    scaled_data = [val / total * 100 for val in data]

    fig, ax = plt.subplots()
    ax.pie(scaled_data, labels=labels, colors=colours)
    plt.show()


def create_treemap():
    # unneeded info: protein to text index: 4 GB
    # Create a data frame with fake data
    df = pd.DataFrame(
        {
            "nb_people": [235.21, 88.2, 0.6, 2, 2, 10.55, 8.4],
            "group": [
                "Sparse suffix array (I = L)\n(sparseness factor 3)\n(235 GB, 68%)",
                "Tekst (I ≠ L)\n(Proteïne sequenties)\n(88 GB, 25.4%)",
                "NCBI taxonomy aggregator (0.6 GB, 0.2%)",
                "Taxonomische annotaties (2 GB, 0.6%)",
                "Suffix → proteïne (2 GB, 0.6%)",
                "Functionele\nannotaties\n(10.5 GB, 3%)",
                "UniProt\naccession\nnumbers\n(8.4 GB, 2.4%)",
            ],
        }
    )

    num_labels_in_legend = 2
    labels = list(df["group"])
    labels_on_plot = labels[:2] + ["" for _ in range(3)] + labels[5:]
    # plot it
    ax = squarify.plot(
        sizes=df["nb_people"],
        color=colours,
        label=labels_on_plot,
        alpha=0.8,
    )
    plt.axis("off")
    plt.gcf().set_size_inches(15.2, 8)

    plt.annotate(
        "Suffix → proteïne (2 GB, 0.6%)",
        xy=(70.8, 95.3),
        fontsize=15,
        xytext=(45.5, 94.5),
        arrowprops=dict(facecolor="black", width=0.8, headwidth=8),
        color="black",
    )
    plt.annotate(
        "Taxonomische annotaties (2 GB, 0.6%)",
        xy=(70.8, 86.3),
        fontsize=15,
        xytext=(39.65, 85.5),
        arrowprops=dict(facecolor="black", width=0.8, headwidth=8),
        color="black",
    )
    plt.annotate(
        "NCBI taxonomy aggregator (0.6 GB, 0.2%)",
        xy=(70.8, 80.3),
        fontsize=15,
        xytext=(36.9, 79.5),
        arrowprops=dict(facecolor="black", width=0.8, headwidth=8),
        color="black",
    )
    plt.text(
        0,
        101,
        "Totale index (346.5 GB) - UniProtKB 2024_01",
        fontsize=15,
        color="black",
    )
    plt.tight_layout()

    plt.show()


def create_treemap_abstract_format():
    font = {"weight": "normal", "size": 21}

    mpl.rc("font", **font)

    # unneeded info: protein to text index: 4 GB
    # Create a data frame with fake data
    df = pd.DataFrame(
        {
            "nb_people": [235.21, 88.2, 0.6, 2, 2, 10.55, 8.4],
            "group": [
                "Sparse suffix array (I = L)\n(sparseness factor 3)\n(235 GB, 68%)",
                "Text (I ≠ L)\n(Protein sequences)\n(88 GB, 25.4%)",
                "NCBI taxonomy aggregator (0.6 GB, 0.2%)",
                "Taxonomic annotations (2 GB, 0.6%)",
                "Suffix → protein (2 GB, 0.6%)",
                "Functional\nannotations\n(10.5 GB, 3%)",
                "UniProt\naccession\nnumbers\n(8.4 GB,\n2.4%)",
            ],
        }
    )

    num_labels_in_legend = 2
    labels = list(df["group"])
    labels_on_plot = labels[:2] + ["" for _ in range(3)] + labels[5:]
    # plot it
    ax = squarify.plot(
        sizes=df["nb_people"],
        color=colours,
        label=labels_on_plot,
        alpha=0.8,
    )
    plt.axis("off")
    plt.gcf().set_size_inches(15.2, 10)

    plt.annotate(
        "Suffix → protein (2 GB, 0.6%)",
        xy=(70.8, 95.3),
        fontsize=23,
        xytext=(34.6, 94.4),
        arrowprops=dict(facecolor="black", width=0.8, headwidth=8),
        color="black",
    )
    plt.annotate(
        "Taxonomic annotations (2 GB, 0.6%)",
        xy=(70.8, 86.3),
        fontsize=23,
        xytext=(26.7, 85.4),
        arrowprops=dict(facecolor="black", width=0.8, headwidth=8),
        color="black",
    )
    plt.annotate(
        "NCBI taxonomy aggregator (0.6 GB, 0.2%)",
        xy=(70.8, 80.3),
        fontsize=23,
        xytext=(20.1, 79.4),
        arrowprops=dict(facecolor="black", width=0.8, headwidth=8),
        color="black",
    )
    plt.text(
        0,
        102,
        "Total index (346.5 GB) - UniProtKB 2024_01",
        fontsize=23,
        color="black",
    )
    plt.tight_layout()

    plt.show()


def plotly_treemap():
    # Define data
    values = [235, 87, 0, 11, 1]
    labels = [
        "Suffix array (235 GB)",
        "Text (87 GB)",
        "Annotations (11 GB)",
        "FA (11 GB)",
        "suffix -> prot (1 GB)",
    ]

    # Define the parent label and value
    parent_label = "full index (335 GB)"

    # Create dataframe
    parents = [parent_label, parent_label, parent_label, labels[2], parent_label]

    # Create treemap
    fig = px.treemap(
        names=labels,
        parents=parents,
        values=values,
        color=labels,
        color_discrete_sequence=colours,
    )
    fig.update_traces(root_color="lightgrey")
    fig.update_layout(hovermode=False)

    fig.show()
    # img_bytes = fig.to_image(format="png", width=1000, height=800, scale=10)
    # fp = io.BytesIO(img_bytes)
    # with fp:
    #     i = mpimg.imread(fp, format="png")
    # plt.axis("off")
    # plt.imshow(i, interpolation="nearest")
    # plt.show()


if __name__ == "__main__":
    # create_pie_chart(
    #     [235, 87, 11, 1], ["SA", "Text", "functionele annotaties", "mapping"]
    # )
    create_treemap()
    create_treemap_abstract_format()
    # plotly_treemap()
