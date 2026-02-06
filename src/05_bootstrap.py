#!/usr/bin/env python3
"""Bootstrap CI for instability rate."""

import csv
import argparse
import numpy as np

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--comparison', required=True)
    parser.add_argument('--n-bootstrap', type=int, default=10000)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    
    differs = []
    with open(args.comparison, 'r') as f:
        for row in csv.DictReader(f, delimiter='\t'):
            differs.append(int(row['differs']))
    
    differs = np.array(differs)
    n = len(differs)
    point_est = differs.mean()
    
    print(f"Computing {args.n_bootstrap:,} bootstrap resamples...")
    boot_rates = []
    for _ in range(args.n_bootstrap):
        sample = np.random.choice(differs, size=n, replace=True)
        boot_rates.append(sample.mean())
    
    ci_low = np.percentile(boot_rates, 2.5)
    ci_high = np.percentile(boot_rates, 97.5)
    
    with open(args.output, 'w') as f:
        f.write(f"Instability Rate Bootstrap Analysis\n")
        f.write(f"====================================\n")
        f.write(f"N variants: {n:,}\n")
        f.write(f"Point estimate: {point_est:.4f} ({point_est*100:.2f}%)\n")
        f.write(f"95% CI: [{ci_low:.4f}, {ci_high:.4f}]\n")
        f.write(f"        [{ci_low*100:.2f}%, {ci_high*100:.2f}%]\n")
        f.write(f"\nStop rule check:\n")
        if ci_low > 0.05:
            f.write(f"  CI lower bound > 5%: EFFECT IS NON-TRIVIAL\n")
        elif ci_high < 0.05:
            f.write(f"  CI upper bound < 5%: EFFECT IS TRIVIAL\n")
        else:
            f.write(f"  CI spans 5% threshold: INCONCLUSIVE\n")
    
    print(f"\nInstability: {point_est:.2%} [{ci_low:.2%}, {ci_high:.2%}]")
    print(f"Results: {args.output}")

if __name__ == "__main__":
    main()
