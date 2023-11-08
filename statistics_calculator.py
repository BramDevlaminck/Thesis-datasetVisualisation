from dataclasses import dataclass
from enum import Enum

"""
File that prints some basic information about the proteins in the given files

This information contains the min, max, avg, median length, total length and number of entries
"""


class FileType(str, Enum):
    DATABASE = ("database",)
    SEARCHFILE = "searchfile"


@dataclass
class TSVRow:
    taxon_id: int
    sequence: str


def read_tsv_file(filename: str, filetype: FileType) -> list[str]:
    rows = []
    with open(filename) as fp:
        for row in fp:
            if filetype == FileType.DATABASE:
                (_, _, taxon_id, _, _, sequence) = row.rstrip().split("\t")
            else:
                sequence = row.rstrip()
            rows.append(sequence)

    return rows


def get_sequence_statistic(rows: list[str], func) -> int:
    return func([len(row) for row in rows])


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


def output_statistics(rows: list[str]):
    minimum = get_sequence_statistic(rows, min)
    maximum = get_sequence_statistic(rows, max)
    average = get_sequence_statistic(rows, avg)
    median = get_sequence_statistic(rows, med)
    print(f"file: {file.split('/')[-2:]}")
    print(f"total length: {len(''.join(rows))}")
    print(f"number of sequences: {len(rows)}")
    print(f"minimum sequence length: {minimum}")
    print(f"maximum sequence length: {maximum}")
    print(f"average sequence length: {average}")
    print(f"median sequence length: {median}")
    print()


if __name__ == "__main__":
    database_files = [
        "/Users/brdvlami/Documents/Ugent/MA2/Thesis/Dataset/BenchmarkData/immunopeptidomics/protein_database.tsv",
        "/Users/brdvlami/Documents/Ugent/MA2/Thesis/Dataset/BenchmarkData/swissprot_var1/protein_database.tsv",
    ]

    search_files = [
        "/Users/brdvlami/Documents/Ugent/MA2/Thesis/Dataset/BenchmarkData/immunopeptidomics/search_file.tsv",
        "/Users/brdvlami/Documents/Ugent/MA2/Thesis/Dataset/BenchmarkData/swissprot_var1/search_file_no_mch.tsv",
        "/Users/brdvlami/Documents/Ugent/MA2/Thesis/Dataset/BenchmarkData/swissprot_var2/search_file_mch.tsv",
        "/Users/brdvlami/Documents/Ugent/MA2/Thesis/Dataset/BenchmarkData/SIHUMI/S03.txt",
        "/Users/brdvlami/Documents/Ugent/MA2/Thesis/Dataset/BenchmarkData/SIHUMI/S05.txt",
        "/Users/brdvlami/Documents/Ugent/MA2/Thesis/Dataset/BenchmarkData/SIHUMI/S07.txt",
        "/Users/brdvlami/Documents/Ugent/MA2/Thesis/Dataset/BenchmarkData/SIHUMI/S08.txt",
        "/Users/brdvlami/Documents/Ugent/MA2/Thesis/Dataset/BenchmarkData/SIHUMI/S11.txt",
        "/Users/brdvlami/Documents/Ugent/MA2/Thesis/Dataset/BenchmarkData/SIHUMI/S14.txt",
    ]
    for file in database_files:
        rows = read_tsv_file(file, FileType.DATABASE)

        output_statistics(rows)

    for file in search_files:
        rows = read_tsv_file(file, FileType.SEARCHFILE)

        output_statistics(rows)
