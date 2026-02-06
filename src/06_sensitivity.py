#!/usr/bin/env python3
"""Window size sensitivity analysis."""

import csv
import argparse
import subprocess
import os

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--variants', required=True)
    parser.add_argument('--tss-v29', required=True)
    parser.add_argument('--tss-v45', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    
    windows = [250, 500, 1000]
    results = []
    
    for w in windows:
        print(f"Testing window ±{w} kb...")
        tmp_out = f"/tmp/map_{w}kb.tsv"
        subprocess.run([
            'python', 'src/04_map_variants.py',
            '--variants', args.variants,
            '--tss-v29', args.tss_v29,
            '--tss-v45', args.tss_v45,
            '--window-kb', str(w),
            '--output', tmp_out
        ], capture_output=True)
        
        n, n_diff = 0, 0
        with open(tmp_out, 'r') as f:
            for row in csv.DictReader(f, delimiter='\t'):
                n += 1
                n_diff += int(row['differs'])
        
        rate = n_diff / n if n > 0 else 0
        results.append({'window_kb': w, 'n': n, 'n_differ': n_diff, 'rate': rate})
        os.remove(tmp_out)
    
    with open(args.output, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['window_kb','n','n_differ','rate'], delimiter='\t')
        writer.writeheader()
        writer.writerows(results)
    
    print(f"\nWindow Sensitivity:")
    for r in results:
        print(f"  ±{r['window_kb']} kb: {r['rate']:.2%}")
    
    rates = [r['rate'] for r in results]
    print(f"\nRange: {(max(rates)-min(rates))*100:.2f} percentage points")

if __name__ == "__main__":
    main()
