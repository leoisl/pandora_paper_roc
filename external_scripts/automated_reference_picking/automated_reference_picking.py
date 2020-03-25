import argparse
import subprocess
import collections


def run_command(command):
    subprocess.check_call(command, shell=True)

def get_args():
    parser = argparse.ArgumentParser(description='Pick references from a set of references and assemblies to run single-reference callers.')
    parser.add_argument('--refs', type=str, help='Path to a file where each line is a path to a ref', required=True)
    parser.add_argument('--assemblies', type=str, help='Path to a file where each line is a path to an assembly', required=True)
    parser.add_argument('--nb_best_ref_per_assembly', type=int, help='How many refs to output per assembly (n bests are output)',
                        required=True)
    parser.add_argument('--nb_best_refs_for_all_assemblies', type=int,
                        help='How many refs to output for all assemblies (n bests to all assemblies are output)',
                        required=True)
    parser.add_argument('--output', type=str,
                        help='Where to save the best refs found',
                        required=True)
    args = parser.parse_args()
    return args


def choose_best_refs(sorted_mash_pb_fh, nb_best_ref_per_assembly, nb_best_refs_for_all_assemblies):
    assembly_to_best_ref = collections.defaultdict(list)
    ref_to_total_mash_distance = collections.defaultdict(int)

    for line in sorted_mash_pb_fh:
        # read input
        line_split = line.strip().split()
        ref, assembly, mash_distance, p_value, matching_hashes = line_split
        mash_distance = float(mash_distance)

        assembly_to_best_ref[assembly].append(
            {"ref": ref, "mash_distance": mash_distance})
        ref_to_total_mash_distance[ref] += mash_distance

    # sort assembly_to_best_ref according to mash distance
    for _, list_of_refs in assembly_to_best_ref.items():
        list_of_refs.sort(key=lambda x: x["mash_distance"])
    with open("assembly_to_best_ref.debug", "w") as assembly_to_best_ref_filehandler:
        print(assembly_to_best_ref, file=assembly_to_best_ref_filehandler)


    ref_to_total_mash_distance_sorted_by_best_overall = list(ref_to_total_mash_distance.items())
    ref_to_total_mash_distance_sorted_by_best_overall.sort(key=lambda x: x[1])
    with open("ref_to_total_mash_distance_sorted_by_best_overall.debug", "w") as ref_to_total_mash_distance_sorted_by_best_overall_filehandler:
        print(ref_to_total_mash_distance_sorted_by_best_overall, file=ref_to_total_mash_distance_sorted_by_best_overall_filehandler)

    # get best_refs
    best_refs_for_each_sample = set()
    for assembly, list_of_refs in assembly_to_best_ref.items():
        for ref_mash_distance in list_of_refs[:nb_best_ref_per_assembly]:
            print(f"Best ref for {assembly} is {ref_mash_distance['ref']}")
            best_refs_for_each_sample.add(ref_mash_distance["ref"])

    best_refs_for_all_samples = set()
    for ref_mash_distance in ref_to_total_mash_distance_sorted_by_best_overall:
        if len(best_refs_for_all_samples) >= nb_best_refs_for_all_assemblies: break
        print(f"Best global ref: {ref_mash_distance[0]}")
        best_refs_for_all_samples.add(ref_mash_distance[0])

    best_refs = best_refs_for_each_sample.union(best_refs_for_all_samples)

    print("Best refs:")
    print("\n".join(best_refs))

    print(f"Number of samples: {len(assembly_to_best_ref)}")
    print(f"Number of best refs for each sample: {len(best_refs_for_each_sample)}")
    print(f"Number of collisions: {len(assembly_to_best_ref) - len(best_refs_for_each_sample)}")
    print(f"Number of global best refs: {len(best_refs_for_all_samples)}")
    print(f"Number of collisions between local and global best refs: {len(best_refs_for_each_sample.intersection(best_refs_for_all_samples))}")
    print(f"Number of best refs: {len(best_refs)}")


    return best_refs


def main():
    args = get_args()
    refs = args.refs
    assemblies = args.assemblies
    nb_best_ref_per_assembly = args.nb_best_ref_per_assembly
    nb_best_refs_for_all_assemblies = args.nb_best_refs_for_all_assemblies
    output = args.output

    print("Checking if mash is working...")
    run_command("mash -h")
    print("mash works!")

    print("Running pairwise mash dists...")
    run_command(f"mash sketch -l {refs}")
    run_command(f"mash sketch -l {assemblies}")
    run_command(f"mash dist {refs}.msh {assemblies}.msh > mash_pb.out")
    run_command(f"cat mash_pb.out | sort -k3 -n > sorted_mash_pb.out")
    print("Done!")

    print("Choosing best refs...")
    with open("sorted_mash_pb.out") as sorted_mash_pb_fh:
        best_refs = choose_best_refs(sorted_mash_pb_fh, nb_best_ref_per_assembly, nb_best_refs_for_all_assemblies)

    with open(output, "w") as output_filehandler:
        print("\n".join(best_refs), file=output_filehandler)
    print("All done!")


if __name__ == "__main__":
    main()