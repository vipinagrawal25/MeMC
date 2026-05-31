"""
makesolidpts.py  --  generate solid_index.h5 for semisolid simulations

Usage:
    python makesolidpts.py <input.h5> <num_solid_points> [--bdry_type 1] [--nghst 12] [--seed 42]

Reads mesh connectivity from input.h5, excludes frame particles and particles
in direct contact with the frame, then randomly selects num_solid_points from
the remaining bulk and writes their indices to solid_index.h5 in the same folder.
"""

import numpy as np
import h5py, sys, os, argparse

def get_nframe(N, bdry_type):
    nf1 = int(np.sqrt(N))
    if bdry_type == 0:
        return 2 * nf1
    elif bdry_type == 1:
        return 4 * nf1
    else:
        return 0

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input',       help='path to input.h5')
    parser.add_argument('num_solid',   type=int, help='number of solid points')
    parser.add_argument('--bdry_type', type=int, default=1)
    parser.add_argument('--nghst',     type=int, default=12)
    parser.add_argument('--seed',      type=int, default=42)
    args = parser.parse_args()

    with h5py.File(args.input, 'r') as f:
        pos      = f['pos'][()]
        node_nbr = f['node_nbr'][()]

    N      = len(pos) // 3
    nframe = get_nframe(N, args.bdry_type)
    print(f"N = {N},  nframe = {nframe},  bulk = {N - nframe}")

    # exclude frame particles and their direct neighbours
    excluded = set(range(nframe))
    for i in range(nframe):
        for k in range(args.nghst):
            nbr = node_nbr[i * args.nghst + k]
            if nbr >= 0:
                excluded.add(nbr)

    candidates = [i for i in range(nframe, N) if i not in excluded]
    print(f"Candidates (bulk minus frame-neighbours): {len(candidates)}")

    if args.num_solid > len(candidates):
        print(f"ERROR: num_solid ({args.num_solid}) > available candidates ({len(candidates)})")
        sys.exit(1)

    rng = np.random.default_rng(args.seed)
    solid_indices = np.sort(
        rng.choice(candidates, size=args.num_solid, replace=False)
    ).astype(np.int32)

    outfile = os.path.join(os.path.dirname(os.path.abspath(args.input)), 'solid_index.h5')
    with h5py.File(outfile, 'w') as f:
        f.create_dataset('solid_idx', data=solid_indices)

    print(f"Wrote {args.num_solid} solid indices -> {outfile}")

if __name__ == '__main__':
    main()
