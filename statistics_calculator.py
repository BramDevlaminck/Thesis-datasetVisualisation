from dataclasses import dataclass


@dataclass
class TSVRow:
    taxon_id: int
    sequence: str


def read_tsv_file(filename) -> list[TSVRow]:
    rows = []
    with open(filename) as fp:
        for row in fp:
            (_, _, taxon_id, _, _, sequence) = row.rstrip().split("\t")
            rows.append(TSVRow(int(taxon_id), sequence))

    return rows


def get_sequence_statistic(rows: list[TSVRow], func) -> int:
    return func([len(row.sequence) for row in rows])


def avg(data: list[float]) -> float:
    return sum(data) / len(data)


def med(data: list[float]) -> float:
    in_order = sorted(data)
    if len(in_order) % 2 == 1:
        return in_order[len(in_order) // 2]
    else:
        index_1 = len(in_order) // 2
        index_2 = index_1 + 1
        return (in_order[index_1] + in_order[index_2]) / 2


if __name__ == "__main__":

    files = [
        "/Users/brdvlami/Documents/Ugent/MA2/Thesis/Dataset/BenchmarkData/immunopeptidomics/protein_database.tsv",
        "/Users/brdvlami/Documents/Ugent/MA2/Thesis/Dataset/BenchmarkData/swissprot_var1/protein_database.tsv"
    ]
    for file in files:
        rows = read_tsv_file(file)

        minimum = get_sequence_statistic(rows, min)
        maximum = get_sequence_statistic(rows, max)
        average = get_sequence_statistic(rows, avg)
        median = get_sequence_statistic(rows, med)
        print(f"file: {file.split('/')[-2]}")
        print(f"minimum sequence length: {minimum}")
        print(f"maximum sequence length: {maximum}")
        print(f"average sequence length: {average}")
        print(f"median sequence length: {median}")
        print()
