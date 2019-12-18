from collections import namedtuple

MemInfoEntry = namedtuple('MemInfoEntry', ['value', 'unit'])
meminfo = {}
with open('/proc/meminfo') as file:
    for line in file:
        key, value, *unit = line.strip().split()
        meminfo[key.rstrip(':')] = MemInfoEntry(value, unit)
memtotal = round((int(meminfo['MemTotal'].value) / 10000000) - 2)