#!/usr/bin/python
import getopt
import sys
import math
import copy
import os
import socket
import subprocess
import re
import types

class ScalFMMConfig(object):
    num_threads    = 1
    num_nodes      = 1
    algorithm      = "implicit"
    model          = "cube"
    num_particules = 10000
    height         = 4
    bloc_size      = 100
    order          = 5

    def show(self):
        print ("=== Simulation parameters ===")
        print ("Number of nodes: " + str(self.num_nodes))
        print ("Number of threads: " + str(self.num_threads))
        print ("Model: " + str(self.model))
        print ("Number of particules: " + str(self.num_particules))
        print ("Height: " + str(self.height))
        print ("Bloc size: " + str(self.bloc_size))
        print ("Order: " + str(self.order))

    def gen_header(self):
        columns = [
            "model",
            "algo",
            "nnode",
            "nthreads",
            "npart",
            "height",
            "bsize",
            "global_time",
            "runtime_time",
            "task_time",
            "idle_time",
            "scheduling_time",
            "communication_time",
            "rmem",
        ]
        header = ""
        for i in range(len(columns)):
            if not i == 0:
                header += ","
            header += "\"" + columns[i] + "\""
        header += "\n"
        return header


    def gen_record(self, global_time, runtime_time, task_time, idle_time, scheduling_time, rmem):
        columns = [
            self.model,
            self.algorithm,
            self.num_nodes,
            self.num_threads,
            self.num_particules,
            self.height,
            self.bloc_size,
            global_time,
            runtime_time,
            task_time,
            idle_time,
            scheduling_time,
            0.0,
            rmem,
        ]
        record = ""
        for i in range(len(columns)):
            if not i == 0:
                record += ","
            if (type(columns[i]) is bool or
                type(columns[i]) == str):
                record += "\""
            record += str(columns[i])
            if (type(columns[i]) == bool or
                type(columns[i]) == str):
                record += "\""
        record += "\n"
        return record

def get_times_from_trace_file(filename):
    cmd = "starpu_trace_state_stats.py " + filename
    proc = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    if not proc.returncode == 0:
        sys.exit("FATAL: Failed to parse trace.rec!")
        return proc.returncode
    task_time = 0.0
    idle_time = 0.0
    runtime_time = 0.0
    scheduling_time = 0.0
    for line in stdout.decode().splitlines():
        arr = line.replace("\"", "").split(",")
        if arr[0] == "Name":
            continue
        if len(arr) >= 4:
            if arr[2] == "Runtime":
                if arr[0] == "Scheduling":
                    scheduling_time = float(arr[3])
                else:
                    runtime_time = float(arr[3])
            elif arr[2] == "Task":
                task_time += float(arr[3])
            elif arr[2] == "Other":
                idle_time = float(arr[3])
            # sys.exit("Invalid time!")
    return runtime_time, task_time, idle_time, scheduling_time

def main():
    output_trace_file=""
    trace_filename="trace.rec"
    output_filename="loutre.db"

    long_opts = ["help",
                 "trace-file=",
                 "output-trace-file=",
                 "output-file="]

    opts, args = getopt.getopt(sys.argv[1:], "ht:i:o:", long_opts)
    for o, a in opts:
        if o in ("-h", "--help"):
            # usage()
            print("No help")
            sys.exit()
        elif o in ("-t", "--trace-file"):
            trace_filename = str(a)
        elif o in ("-i", "--output-trace-file"):
            output_trace_file = str(a)
        elif o in ("-o", "--output-file"):
            output_filename = str(a)
        else:
            assert False, "unhandled option"

    config=ScalFMMConfig()
    rmem = 0
    global_time = 0.0
    runtime_time = 0.0
    task_time = 0.0
    idle_time = 0.0
    scheduling_time = 0.0

    if (os.path.isfile(output_filename)): #Time in milli
        output_file = open(output_filename, "a")
    else:
        output_file = open(output_filename, "w")
        output_file.write(config.gen_header())

    with open(output_trace_file, "r") as ins:
        for line in ins:
            if re.search("Average", line):
                a = re.findall("[-+]?\d*\.\d+|\d+", line)
                if len(a) == 1:
                    global_time = a[0]
            elif re.search("Total particles", line):
                a = re.findall("[-+]?\d*\.\d+|\d+", line)
                if len(a) == 1:
                    config.num_particules = int(a[0])
            elif re.search("Group size", line):
                a = re.findall("[-+]?\d*\.\d+|\d+", line)
                if len(a) == 1:
                    config.bloc_size = int(a[0])
            elif re.search("Nb node", line):
                a = re.findall("[-+]?\d*\.\d+|\d+", line)
                if len(a) == 1:
                    config.num_nodes = int(a[0])
            elif re.search("Tree height", line):
                a = re.findall("[-+]?\d*\.\d+|\d+", line)
                if len(a) == 1:
                    config.height = int(a[0])
            elif re.search("Nb thread", line):
                a = re.findall("[-+]?\d*\.\d+|\d+", line)
                if len(a) == 1:
                    config.num_threads = int(a[0])
            elif re.search("Model", line):
                config.model = line[line.index(":")+1:].strip()
            elif re.search("Algorithm", line):
                config.algorithm = line[line.index(":")+1:].strip()

    if (os.path.isfile(trace_filename)): #Time in milli
        runtime_time, task_time, idle_time, scheduling_time = get_times_from_trace_file(trace_filename)
    else:
        print("File doesn't exist " + trace_filename)

    # Write a record to the output file.
    output_file.write(config.gen_record(float(global_time),
                      float(runtime_time),
                      float(task_time),
                      float(idle_time),
                      float(scheduling_time),
                      int(rmem)))

main()
