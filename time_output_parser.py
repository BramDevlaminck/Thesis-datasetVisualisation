from dataclasses import dataclass

"""
File That contains function to extract the extended output from the /usr/bin/time -v command

The (useful) data is then stored in a TimeOutput object for each iteration
"""


@dataclass
class TimeOutput:
    command: str
    execution_time_seconds: float
    max_mem_size: float


def get_value_from_line(line: str) -> str:
    """
    The time -v command has as output always some keys, then a ':' and then the actual data

    This function simply extracts the data from the line
    """
    return "".join(line.split(":")[1:]).strip()


def progress_iterator(iterator: iter, count: int):
    """Progresses the given iterator `count` times"""
    for _ in range(count):
        next(iterator)


def parse_time_output_file(filename: str) -> list[TimeOutput]:
    """Parses the file"""
    data = []
    with open(filename) as fp:
        while (command_line := next(fp, None)) is not None:
            command = get_value_from_line(command_line).strip('"')
            user_time = float(get_value_from_line(next(fp)))
            system_time = float(get_value_from_line(next(fp)))
            progress_iterator(fp, 6)
            max_mem_size = float(get_value_from_line(next(fp)))
            progress_iterator(fp, 13)
            data.append(TimeOutput(command, user_time + system_time, max_mem_size))

    return data


def aggregate_time_output(data: list[TimeOutput]) -> TimeOutput:
    """Takes the average value of each field of the timeOutput objects and returns this average object"""
    command = data[0].command
    execution_time = sum(map(lambda ex: ex.execution_time_seconds, data)) / len(data)
    memory = sum(map(lambda ex: ex.max_mem_size, data)) / len(data)

    return TimeOutput(command, execution_time, memory)


def parse_and_aggregate_time_output_file(filename: str) -> TimeOutput:
    return aggregate_time_output(parse_time_output_file(filename))
