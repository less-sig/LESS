#!/usr/bin/env python3
""" """
import os
import os.path
import re

__all__ = ["parse", "parse_file", "ParseError"]

# Precompiled regex used to parse comments.
_COMMENT_RE = re.compile("\s*(#|$)")

# Precompiled regexes used to parse header fields.
_FIELD_DESC_RE = re.compile("desc:\s(?P<data>.*)$")
_FIELD_CMD_RE = re.compile("cmd:\s(?P<data>.*)$")
_FIELD_TIME_UNIT_RE = re.compile("time_unit:\s(?P<data>ms|B|i)$")

# Precompiled regexes used to parse snaphot fields.
_FIELD_SNAPSHOT_RE = re.compile("snapshot=(?P<data>\d+)")
_FIELD_TIME_RE = re.compile("time=(?P<data>\d+)")
_FIELD_MEM_HEAP_RE = re.compile("mem_heap_B=(?P<data>\d+)")
_FIELD_MEM_EXTRA_RE = re.compile("mem_heap_extra_B=(?P<data>\d+)")
_FIELD_MEM_STACK_RE = re.compile("mem_stacks_B=(?P<data>\d+)")
_FIELD_HEAP_TREE_RE = re.compile("heap_tree=(?P<data>\w+)")

# Precompiled regex to parse heap entries. Matches three things:
#   - the number of children,
#   - the number of bytes,
#   - and the details section.
_HEAP_ENTRY_RE = re.compile("""
    \s*n                    # skip zero or more spaces, then 'n'
    (?P<num_children>\d+)   # match number of children, 1 or more digits
    :\s                     # skip ':' and one space
    (?P<num_bytes>\d+)      # match the number of bytes, 1 or more digits
    \s                      # skip one space
    (?P<details>.*)         # match the details
""", re.VERBOSE)

# Precompiled regex to check if the details section is below threshold.
_HEAP_BELOW_THRESHOLD_RE = re.compile(r"""in.*places?.*""")

# Precompiled regex to parse the details section of entries above threshold.
# This should match four things:
#   - the hexadecimal address,
#   - the function name,
#   - the file name or binary path, i.e. file.cpp or usr/local/bin/foo.so,
#   - and a line number if present.
# Last two parts are optional to handle entries without a file name or binary
# path.
_HEAP_DETAILS_RE = re.compile(r"""
    (?P<address>[a-fA-F0-9x]+)  # match the hexadecimal address
    :\s                         # skip ': '
    (?P<function>.+?)           # match the function's name, non-greedy
    (?:                         # don't capture fname/line group
        \s
        \(
        (?:in\s)?               # skip 'in ' if present
        (?P<fname>[^:]+)        # match the file name
        :?                      # skip ':', if present
        (?P<line>\d+)?          # match the line number, if present
        \)
    )?                          # fname/line group is optional
    $                           # should have reached the EOL
""", re.VERBOSE)


class ParseContext:
    """
    A simple context for parsing. Dumbed down version of fileinput.
    """
    def __init__(self, fd):
        self._fd = fd
        self._line = 0

    def line(self):
        return self._line

    def readline(self):
        self._line += 1
        return self._fd.readline()

    def filename(self):
        return os.path.abspath(self._fd.name)


class ParseError(Exception):
    """
    Error raised when a parsing error is encountered.
    """
    def __init__(self, msg, ctx):
        self.msg = msg
        self.line = ctx.line()
        self.filename = ctx.filename()

    def __str__(self):
        return " ".join([str(self.msg), 'at line', str(self.line), 'in',
                        str(self.filename)])


def parse_file(filepath):
    """
    Convenience function taking a file path instead of a file descriptor.
    """
    with open(filepath) as fd:
        return parse(fd)


def parse(fd):
    """
    Parse an already opened massif output file.
    """
    mdata = {}
    ctx = ParseContext(fd)
    _parse_header(ctx, mdata)
    _parse_snapshots(ctx, mdata)
    return mdata


def _match_unconditional(ctx, regex, string):
    """
    Unconditionaly match a regular expression against a string, i.e. if there
    is no match we raise a ParseError.
    """
    match = regex.match(string)
    if match is None:
        raise ParseError("".join(["can't match '", string, "' against '",
                         regex.pattern, "'"]), ctx)
    return match


def _get_next_line(ctx, may_reach_eof=False):
    """
    Read another line from ctx. If may_reach_eof is False, reaching EOF will
    be considered as an error.
    """
    line = ctx.readline()  # Returns an empty string on EOF.

    if len(line) == 0:
        if may_reach_eof is False:
            raise ParseError("unexpected EOF", ctx)
        else:
            return None
    else:
        return line.strip("\n")


def _get_next_field(ctx, field_regex, may_reach_eof=False):
    """
    Read the next data field. The field_regex arg is a regular expression that
    will be used to match the field. Data will be extracted from the match
    object by calling m.group('data'). If may_reach_eof is False, reaching EOF
    will be considered as an error.
    """
    line = _get_next_line(ctx, may_reach_eof)
    while line is not None:
        if _COMMENT_RE.match(line):
            line = _get_next_line(ctx, may_reach_eof)
        else:
            match = _match_unconditional(ctx, field_regex, line)
            return match.group("data")

    return None


def _parse_header(ctx, mdata):
    mdata["desc"] = _get_next_field(ctx, _FIELD_DESC_RE)
    mdata["cmd"] = _get_next_field(ctx, _FIELD_CMD_RE)
    mdata["time_unit"] = _get_next_field(ctx, _FIELD_TIME_UNIT_RE)


def _parse_snapshots(ctx, mdata):
    index = 0
    snapshots = []
    detailed_snapshot_indices = []
    peak_snapshot_index = None

    snapshot = _parse_snapshot(ctx)

    while snapshot is not None:
        if snapshot["is_detailed"]:
            detailed_snapshot_indices.append(index)
        if snapshot["is_peak"]:
            peak_snapshot_index = index
        snapshots.append(snapshot["data"])
        snapshot = _parse_snapshot(ctx)
        index += 1

    mdata["snapshots"] = snapshots
    mdata["detailed_snapshot_indices"] = detailed_snapshot_indices

    if peak_snapshot_index is not None:
        mdata["peak_snapshot_index"] = peak_snapshot_index


def _parse_snapshot(ctx):
    """
    Parse another snapshot, appending it to the mdata["snapshots"] list. On
    EOF, False will be returned.
    """
    snapshot_id = _get_next_field(ctx, _FIELD_SNAPSHOT_RE, may_reach_eof=True)

    if snapshot_id is None:
        return None

    snapshot_id = int(snapshot_id)
    time = int(_get_next_field(ctx, _FIELD_TIME_RE))
    mem_heap = int(_get_next_field(ctx, _FIELD_MEM_HEAP_RE))
    mem_heap_extra = int(_get_next_field(ctx, _FIELD_MEM_EXTRA_RE))
    mem_stacks = int(_get_next_field(ctx, _FIELD_MEM_STACK_RE))
    heap_tree_field = _get_next_field(ctx, _FIELD_HEAP_TREE_RE)

    heap_tree = None
    is_detailed = False
    is_peak = False

    if heap_tree_field != "empty":
        is_detailed = True
        if heap_tree_field == "peak":
            is_peak = True
        heap_tree = _parse_heap_tree(ctx)

    return {
        "is_detailed": is_detailed,
        "is_peak": is_peak,
        "data": {
            "id": snapshot_id,
            "time": time,
            "mem_heap": mem_heap,
            "mem_heap_extra": mem_heap_extra,
            "mem_stack": mem_stacks,
            "heap_tree": heap_tree
        }
    }


def _parse_heap_tree(ctx):
    """
    Parse a heap tree.
    """
    line = _get_next_line(ctx)

    entry_match = _match_unconditional(ctx, _HEAP_ENTRY_RE, line)
    details_group = entry_match.group("details")

    details = None
    details_match = _HEAP_DETAILS_RE.match(details_group)

    if details_match:
        # The 'line' field could be None if the binary/library wasn't compiled
        # with debug info. To avoid errors on this condition, we need to make
        # sure that the 'line' field is not None before trying to convert it to
        # an integer.
        linum = details_match.group(4)
        if linum is not None:
            linum = int(linum)

        details = {
            "address": details_match.group("address"),
            "function": details_match.group("function"),
            "file": details_match.group("fname"),
            "line": linum
        }

    children = []
    for i in range(0, int(entry_match.group("num_children"))):
        children.append(_parse_heap_tree(ctx))

    heap_node = {}
    heap_node["nbytes"] = int(entry_match.group("num_bytes"))
    heap_node["children"] = children
    heap_node["details"] = details

    return heap_node


def massif_analyse(massif_filename):
    data = parse_file(massif_filename)
    
    #header = 'time,stack,heap,heap_extra'
    max_heap = 0
    max_stack = 0

    for snapshot in data['snapshots']:
        time = snapshot['time']
        stack = snapshot['mem_stack']
        heap = snapshot['mem_heap']
        heap_extra = snapshot['mem_heap_extra']
        #print(time, stack, heap, heap_extra)
        if heap > max_heap:
            max_heap = heap

        if stack > max_stack:
            max_stack = stack
    
    # print("max heap:", max_heap, "max_stack", max_stack)
    return max_heap, max_stack



# Execute a command
def command(cmd):
    print(cmd)
    os.system(cmd)

values = {
        252 : {"192": {"keygen": 0, "sign": 0, "verify": 0},
             "68": {"keygen": 0, "sign": 0, "verify": 0},
             "45": {"keygen": 0, "sign": 0, "verify": 0},},
        400 : {"220": {"keygen": 0, "sign": 0, "verify": 0},
             "102": {"keygen": 0, "sign": 0, "verify": 0},},
        548 : {"345": {"keygen": 0, "sign": 0, "verify": 0},
             "137": {"keygen": 0, "sign": 0, "verify": 0},}
}

# create all the massif files:
# command("./bench.sh")

file = "build/LESS_memory_1_BALANCED_sign.massif"
rootdir = "./build"
for subdir, dirs, files in os.walk(rootdir):
    for file in files:
        filepath = subdir + os.sep + file
        if filepath.endswith(".massif"):
            # this is stupid: but one version is called `SHORT_SIG`
            try:
                _, _, cat, kind, method = file.split("_")
            except:
                _, _, cat, kind, kind_2, method = file.split("_")
                kind = kind + "_" + kind_2

            method = method.split(".")[0]
            cat = int(cat)
            max_heap, max_stack = massif_analyse(filepath)
            #print(cat, kind, method, max_stack)

            values[cat][kind][method] = max_stack

print(values)
print("Everything's done!\n")
