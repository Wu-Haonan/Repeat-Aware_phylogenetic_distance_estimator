import argparse
import os
import csv
import subprocess
from itertools import combinations, product
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import numpy as np
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceMatrix
from Bio import Phylo
from sourmash import load_file_as_signatures
from scipy.optimize import brentq
import shutil

def get_script_dir():
    """Get the directory where this script is located"""
    return os.path.dirname(os.path.abspath(__file__))

def get_cpp_executable():
    """Get the absolute path to the cpp histogram executable"""
    script_dir = get_script_dir()
    cpp_path = os.path.join(script_dir, "cpp", "kmc_histogram")
    
    if not os.path.isfile(cpp_path):
        raise FileNotFoundError(
            f"kmc_histogram executable not found at: {cpp_path}\n"
            f"Please ensure the cpp/kmc_histogram is compiled in the same directory as this script."
        )
    
    if not os.access(cpp_path, os.X_OK):
        raise PermissionError(
            f"kmc_histogram executable found but not executable: {cpp_path}\n"
            f"Please check file permissions."
        )
    
    return cpp_path

def parse_args():
    parser = argparse.ArgumentParser(description="Repeat-Aware Mutation Distance Estimator")
    parser.add_argument("--config", required=True, help="CSV with columns: label,filepath")
    parser.add_argument("--mode", choices=["half", "all"], default="half", help="Whether to compute one-direction or both")
    parser.add_argument("--k", type=int, default=31, help="k-mer size")
    parser.add_argument("--theta", type=float, default=0.01, help="FracMinHash theta")
    parser.add_argument("--cores", type=int, default=4, help="Number of CPU cores to use (processes)")
    parser.add_argument("--kmc-threads", type=int, default=1, help="Threads per KMC task (default: 1)")
    parser.add_argument("--outdir", default="results", help="Output directory")
    parser.add_argument("--cleanup", action="store_true", help="Clean intermediate files")
    return parser.parse_args()

def read_config(config_path):
    labels, paths = [], []
    with open(config_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            labels.append(row["label"])
            paths.append(row["filepath"])
    return labels, paths

def run_kmc(label, fasta_path, k, middle_dir, kmc_threads):
    kmc_dir = os.path.join(middle_dir, f"kmc_{label}_k{k}")
    Path(kmc_dir).mkdir(parents=True, exist_ok=True)
    out_prefix = os.path.join(kmc_dir, f"{label}_k{k}")
    tmp_dir = os.path.join(kmc_dir, "tmp")
    Path(tmp_dir).mkdir(parents=True, exist_ok=True)
    subprocess.run("ulimit -n 2048", shell=True, check=True)
    if not os.path.exists(out_prefix + ".kmc_pre"):
        cmd = f"kmc -k{k} -ci1 -t{kmc_threads} -fm {fasta_path} {out_prefix} {tmp_dir}"
        subprocess.run(cmd, shell=True, check=True)
    return out_prefix

def run_kmc_histogram(label, out_prefix, k, middle_dir, cpp_exec_path):
    hist_out = os.path.join(middle_dir, f"{label}_k{k}_hist.csv")
    if not os.path.exists(hist_out):
        # Use the provided cpp executable path
        subprocess.run([cpp_exec_path, out_prefix, hist_out], check=True)
    return hist_out

def run_sourmash_sketch(label, fasta_path, k, theta, middle_dir):
    sketch_path = os.path.join(middle_dir, f"{label}_k{k}_sig.sig")
    if not os.path.exists(sketch_path):
        scaled = int(1 / theta)
        subprocess.run(["sourmash", "sketch", "dna", "-p", f"k={k},scaled={scaled}", "-o", sketch_path, fasta_path], check=True)
    return sketch_path

def compute_intersection(sig1_file, sig2_file):
    try:
        s1 = next(load_file_as_signatures(sig1_file))
        s2 = next(load_file_as_signatures(sig2_file))
        intersection_size = s1.minhash.count_common(s2.minhash)
        print(f"[INFO] Sketch intersection : {intersection_size}")
        return intersection_size
    except StopIteration:
        print("[ERROR] One of the files contains no signatures")
        return None
    except Exception as e:
        print(f"[ERROR] {e}")
        return None

def solve_histogram_equation(hist_dict, I, total_kmers):
    if I == 0:
        return 1
    else:
        def poly(x):
            return (total_kmers - I) - sum(count * x**i for i, count in hist_dict.items())
        try:
            root = brentq(poly, 0, 1, xtol=1e-15)
        except Exception as e:
            print(f"[ERROR] Newton solver failed: {e}")
            return None
        return root

def parse_histogram_csv(path):
    hist = {}
    with open(path, 'r') as f:
        reader = csv.reader(f)
        next(reader)  # skip header
        
        for row in reader:
            if len(row) >= 2:
                count = int(row[0])
                num_kmers = int(row[1])
                hist[count] = num_kmers
    return hist

def compute_distance(labelA, sigA, histA_path, labelB, sigB, histB_path, k, theta, middle_dir):
    out_prefix = os.path.join(middle_dir, f"{labelA}___{labelB}")
    out_file = out_prefix + ".csv"
    if os.path.exists(out_file):
        return out_file

    I = compute_intersection(sigA, sigB)
    if I is None:
        return None

    I /= theta
    histA = parse_histogram_csv(histA_path)
    total_kmers = sum(histA.values())

    print(f"[INFO] Total number of kmers of {labelA}: {total_kmers}")
    print(f"[INFO] Estimated intersection I: {I}")

    if I >= total_kmers:
        r = 0.0
    else:
        q = solve_histogram_equation(histA, I, total_kmers)
        r = 1 - (1 - q) ** (1/k)

    print(f"[RESULT] {labelA} , {labelB} Estimated r : {r}")
    with open(out_file, 'w') as f:
        f.write(f"label1,label2,distance\n{labelA},{labelB},{r}\n")
    return out_file

def cleanup_middle_files(middle_dir):
    """Clean up all intermediate files in the middle directory"""
    if os.path.exists(middle_dir):
        print(f"[INFO] Cleaning up intermediate files in {middle_dir}")
        shutil.rmtree(middle_dir)

def build_nj_tree(labels, distances, out_file):
    """
    Build NJ tree, using DistanceMatrix API
    """
    n = len(labels)
    
    print(f"[INFO] Building NJ tree for {n} species")
    
    full_matrix = np.zeros((n, n))
    
    for i in range(n):
        for j in range(n):
            if i == j:
                full_matrix[i][j] = 0.0  # set diagonal as 0
            else:
                label_i, label_j = labels[i], labels[j]
                dist = distances.get((label_i, label_j), distances.get((label_j, label_i), 0.0))
                full_matrix[i][j] = dist
    
    dm = DistanceMatrix(names=labels)
    for i in range(n):
        for j in range(i):
            dm[labels[i], labels[j]] = full_matrix[i][j]
    
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dm)
    
    Phylo.write(tree, out_file, "newick")
    print(f"[INFO] NJ tree successfully written to {out_file}")

def run_kmc_histogram_wrapper(args):
    """Wrapper function for multiprocessing"""
    return run_kmc_histogram(*args)

def main():
    args = parse_args()
    
    # Get cpp executable path automatically
    cpp_exec_path = get_cpp_executable()
    print(f"[INFO] Using cpp executable: {cpp_exec_path}")
    
    # Create output directory for final results
    Path(args.outdir).mkdir(parents=True, exist_ok=True)
    
    # Create middle directory for intermediate files
    middle_dir = "./middle_results"
    Path(middle_dir).mkdir(parents=True, exist_ok=True)

    labels, fasta_paths = read_config(args.config)

    print("[Step 1.1] Running KMC and Sourmash Sketch in parallel...")
    kmc_outs = [None] * len(labels)
    sigs = [None] * len(labels)

    from concurrent.futures import as_completed
    with ProcessPoolExecutor(max_workers=args.cores) as executor:
        futures = []
        for i, (label, path) in enumerate(zip(labels, fasta_paths)):
            futures.append((i, 'kmc', executor.submit(run_kmc, label, path, args.k, middle_dir, args.kmc_threads)))
            futures.append((i, 'sour', executor.submit(run_sourmash_sketch, label, path, args.k, args.theta, middle_dir)))

        for i, task_type, fut in futures:
            result = fut.result()
            if task_type == 'kmc':
                kmc_outs[i] = result
            else:
                sigs[i] = result

    print("[Step 1.2] Parsing KMC Histograms...")
    with ProcessPoolExecutor(max_workers=args.cores) as executor:
        # Pass cpp_exec_path to each histogram computation
        arg_list = [(label, out_prefix, args.k, middle_dir, cpp_exec_path) for label, out_prefix in zip(labels, kmc_outs)]
        hists = list(executor.map(run_kmc_histogram_wrapper, arg_list))

    print("[Step 2] Computing pairwise distances...")
    label_map = {l: (sig, hist) for l, sig, hist in zip(labels, sigs, hists)}
    if args.mode == "half":
        pairs = combinations(labels, 2)
    else:
        pairs = [(A, B) for A in labels for B in labels if A != B]

    distances = {}
    with ProcessPoolExecutor(max_workers=args.cores) as executor:
        futures = []
        for A, B in pairs:
            futures.append(executor.submit(compute_distance, A, label_map[A][0], label_map[A][1],
                                                        B, label_map[B][0], label_map[B][1],
                                                        args.k, args.theta, middle_dir))
        for f in futures:
            path = f.result()
            if path:
                with open(path) as pf:
                    next(pf)  # skip header
                    row = pf.readline().strip().split(',')
                    A, B, val = row[0], row[1], float(row[2])
                    distances[(A, B)] = val

    print("[Step 3] Building NJ Tree...")
    
    # Save pairwise distances to final output directory
    dist_file = os.path.join(args.outdir, "pairwise_distances.csv")
    with open(dist_file, 'w') as f:
        f.write("label1,label2,distance\n")
        for (A, B), val in distances.items():
            f.write(f"{A},{B},{val}\n")
    print(f"Pairwise distances written to {dist_file}")

    # Save NJ tree to final output directory
    tree_path = os.path.join(args.outdir, "nj_tree.nwk")
    unique_distances = {}
    for (A, B), val in distances.items():
        key = tuple(sorted([A, B]))
        if key not in unique_distances:
            unique_distances[key] = val
    
    tree_distances = {(A, B): val for (A, B), val in unique_distances.items()}
    build_nj_tree(labels, tree_distances, tree_path)

    # Clean up intermediate files if requested
    if args.cleanup:
        cleanup_middle_files(middle_dir)
    else:
        print(f"[INFO] Intermediate files preserved in {middle_dir}")

if __name__ == "__main__":
    main()