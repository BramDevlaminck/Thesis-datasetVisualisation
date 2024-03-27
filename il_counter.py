import json

import regex as re

from statistics_calculator import FileType, read_tsv_file


def count_successive_il(proteins: list[str]) -> dict[int, int]:
    result: dict[int, int] = dict()
    for protein in proteins:
        for match in re.findall("[IL]+", protein, overlapped=True):
            key = len(match)
            for i in range(1, key + 1):
                current_count = result.get(i, 0)
                result[i] = current_count + 1

    return result


# returns a list which contains the number of I's or L's for each protein
def count_il_occurrences(proteins: list[str]) -> list[int]:
    result = []
    for protein in proteins:
        count = 0
        for char in protein:
            if char == "I" or char == "L":
                count += 1
        result.append(count)

    return result


if __name__ == "__main__":
    proteins = read_tsv_file(
        "/mnt/hdd/uniprot/uniprotKB_protein_database.tsv", FileType.DATABASE
    )

    counts = count_successive_il(proteins)

    with open("il_multiples.json", "w") as fp:
        json.dump(counts, fp)

    il_count_per_protein = count_il_occurrences(proteins)
    with open("il_occurrences_per_sequence.json", "w") as fp:
        json.dump(il_count_per_protein, fp)
    print()
