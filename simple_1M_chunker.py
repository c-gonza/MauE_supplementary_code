#!/usr/bin/env python
# coding: utf-8
import re
import sys
import time
entries_seen = 0
start = time.time()
lines = []
count = 0
for line in sys.stdin.buffer:
    #print(line)
    if line.startswith(b"OX") or line.startswith(b"ID") or line.startswith(b"DR   InterPro") or line.startswith(b"OS") or line.startswith(b"OC") or line.startswith(b"AC"):
        lines.append(line)
    if line.startswith(b"//\n"):
        lines.append(line)
        entries_seen += 1
        if entries_seen >= 1000000:
            #print(lines)
            with open(r"PATH\1M_entry_chunks_260\line_chunk_{}.txt".format(count),"wb") as file:
                file.writelines(lines)
                file.close()
            lines.clear
            entries_seen = 0
            lines = []
            end = time.time()
            print(end - start)
            start = None
            end = None
            start = time.time()
            count+=1
    if count > 260:
        break
#print(lines)
#print(sys.getsizeof(lines)/1000000)

del lines
del start
del end
del entries_seen
sys.exit()