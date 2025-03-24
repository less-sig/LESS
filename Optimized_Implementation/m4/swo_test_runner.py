#!/usr/bin/env python3

from swo_parser import swo_parser_main
from subprocess import Popen, PIPE
from os import path, pipe
from sys import argv, stdout
from re import match
from queue import Queue, Empty
from threading import Thread


def enqueue_output(out, queue):
    for line in out:
        queue.put(line)


if __name__ == "__main__":
    platformio_package_dir = argv[1]
    project_dir = path.abspath(argv[2])
    openocd_path = path.join(platformio_package_dir, "tool-openocd", "bin", "openocd")
    server_args = ["-s", path.join(platformio_package_dir, "tool-openocd", "scripts")]
    server_args.extend(["-f", "board/stm32f4discovery.cfg"])
    swo_trace_clkin_freq = 24000000
    swo_trace_freq = 2000000
    server_args.extend(
        [
            "-c",
            "init; tpiu config internal - uart false %s %s; itm ports on"
            % (swo_trace_clkin_freq, swo_trace_freq),
        ]
    )
    server_args.insert(0, openocd_path)
    print(
        f"Starting OpenOCD with SWO Trace clock-in frequency %s, SWO trace frequency %s. Invocation:"
        % (swo_trace_clkin_freq, swo_trace_freq),
        flush=True,
    )
    print(server_args)
    # start openocd and client processes in parallel
    openocd_process = Popen(server_args)

    # start client process in parallel
    (swo_parser_r_fd, swo_parser_w_fd) = pipe()
    swo_parser_r = open(swo_parser_r_fd, "r")
    swo_parser_w = open(swo_parser_w_fd, "w")

    swo_parser_t = Thread(target=swo_parser_main, args=(swo_parser_w,))
    swo_parser_t.daemon = True
    swo_parser_t.start()

    # https://stackoverflow.com/a/4896288/523079
    q = Queue()
    enqueue_output_t = Thread(target=enqueue_output, args=(swo_parser_r, q))
    enqueue_output_t.daemon = True
    enqueue_output_t.start()

    timeout = int(argv[3])
    tests_done = False

    while True:
        try:
            line = q.get(timeout=timeout)
        except Empty:
            break

        if match(r".*[0-9]* Tests [0-9]* Failures [0-9]* Ignored.*", line):
            tests_done = True
            stdout.write(line)
            stdout.flush()
            continue
        else:
            stdout.write(line)
            stdout.flush()

        if tests_done:
            break

    openocd_process.terminate()
