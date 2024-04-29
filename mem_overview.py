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
                "Sparse suffix array\n(sparseness factor 3)\n(235 GB)",
                "Tekst\n(Proteïne sequenties)\n(88 GB)",
                "NCBI taxonomy aggregator (0.6 GB)",
                "Taxonomische annotaties (2 GB)",
                "Suffix → proteïne (2 GB)",
                "Functionele\nannotaties\n(10.5 GB)",
                "UniProt\naccession\nnumbers\n(8.4 GB)",
            ],
        }
    )

    num_labels_in_legend = 2
    labels = list(df["group"])
    labels_on_plot = labels[:2] + ["" for _ in range(3)] + labels[5:]
    print(labels_on_plot)
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
        "Suffix → proteïne (2 GB)",
        xy=(70.8, 95.3),
        fontsize=15,
        xytext=(49.5, 94.5),
        arrowprops=dict(facecolor="black", width=0.8, headwidth=8),
        color="black",
    )
    plt.annotate(
        "Taxonomische annotaties (2 GB)",
        xy=(70.8, 86.3),
        fontsize=15,
        xytext=(43.65, 85.5),
        arrowprops=dict(facecolor="black", width=0.8, headwidth=8),
        color="black",
    )
    plt.annotate(
        "NCBI taxonomy aggregator (0.6 GB)",
        xy=(70.8, 80.3),
        fontsize=15,
        xytext=(40.9, 79.5),
        arrowprops=dict(facecolor="black", width=0.8, headwidth=8),
        color="black",
    )
    plt.text(0, 101, "Totale index (346.5 GB)", fontsize=15, color="black")
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
    # plotly_treemap()
