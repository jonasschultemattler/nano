import subprocess
import psutil
import time
import re
import os

import matplotlib.pyplot as plt
import numpy as np


def track_memory_and_runtime(executable_path, argument):
    try:
        start_time = time.time()
        process = subprocess.Popen([executable_path, argument], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        pid = process.pid

        # Monitor memory usage
        max_memory = 0
        ps_process = psutil.Process(pid)
        while process.poll() is None:  # While the process is running
            try:
                memory_info = ps_process.memory_info()
                max_memory = max(max_memory, memory_info.rss)  # Track peak memory usage
            except psutil.NoSuchProcess:
                break  # Process has terminated
            time.sleep(0.1)  # Poll every 100ms

        # Wait for the process to finish and capture output
        stdout, stderr = process.communicate()
        end_time = time.time()

        # Calculate runtime
        runtime = end_time - start_time

        # Convert memory usage to MB
        max_memory_mb = max_memory / (1024 * 1024)

        print(f"Runtime: {runtime:.2f} seconds")
        print(f"Peak Memory Usage: {max_memory_mb:.2f} MB")
        if stdout:
            output = stdout.decode()
            print("Program Output:")
            print(output)

            result = int(re.findall(r"Distinct kmers:\s*(\d+)", output)[0])

        if stderr:
            print("Program Errors:")
            print(stderr.decode())

    except FileNotFoundError:
        print(f"Error: The executable '{executable_path}' was not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

    return runtime, max_memory_mb, result



def compare_performance(results, algorithms, file_sizes):
    file_sizes = np.array(file_sizes)/(1024*1024)
    sorted_indices = np.argsort(file_sizes)
    results = results[:,sorted_indices,:]
    file_sizes = file_sizes[sorted_indices]

    fig, (ax0, ax1) = plt.subplots(1,2)
    for i in range(results.shape[1]):
        ax0.plot(file_sizes, results[i,:,0], label=algorithms[i])
        ax0.set_title("Runtime")
        ax0.set_ylabel("Time [s]")
        ax0.set_xlabel("Filesize [Mb]")
    ax0.legend()
    for i in range(results.shape[1]):
        ax1.plot(file_sizes, results[i,:,1], label=algorithms[i])
        ax1.set_title("Memory")
        ax1.set_ylabel("Size [Mb]")
        ax1.set_xlabel("Filesize [Mb]")
    ax1.legend()
    plt.show()


def compare_quality(results, algorithms):
    # todo: violinplot error per dataset
    # todo: hist2d plot truth againt estimate
    fig, ax = plt.subplots()
    ax.errorbar(np.arange(results.shape[1])-0.1, results[1,:,2], np.abs(results[0,:,2] - results[1,:,2]), fmt='o', linewidth=2, capsize=6, label="FM")
    ax.plot(np.arange(results.shape[1]), results[0,:,2], 'kx', label="Ground truth")
    ax.errorbar(np.arange(results.shape[1])+0.1, results[2,:,2], np.abs(results[0,:,2] - results[2,:,2]), fmt='o', linewidth=2, capsize=6, label="HLL")
    ax.legend()
    plt.show()


if __name__ == "__main__":
    files = ["data/ecoli1_k31_ust.fa.gz", "data/ecoli2_k31_ust.fa.gz", "data/ecoli4_k31_ust.fa.gz"]

    naive_results, fm_results, hll_results = [], [], []
    for file in files:
        naive_results.append(track_memory_and_runtime("build/source/naivecounting", file))
        fm_results.append(track_memory_and_runtime("build/source/flajoletmartin", file))
        hll_results.append(track_memory_and_runtime("build/source/hyperloglog", file))

    file_sizes = [os.path.getsize(file) for file in files]
    results = np.array([naive_results, fm_results, hll_results])
    algorithms = ["Naive", "Flajolet Martin", "HLL"]
    compare_performance(results, algorithms, file_sizes)
    compare_quality(results, algorithms)



