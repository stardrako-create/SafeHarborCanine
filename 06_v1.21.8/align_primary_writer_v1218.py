from pathlib import Path
P=Path(r'D:\Biblioteca\Bioquímica 1º ano 2025-2026\2º Semestre\Vasco M. Barreto\Safe Harbor Canino');p=P/'scripts/build_mother_track_v2.py';s=p.read_text(encoding='utf-8')
old='''        starts = list(range(0, n_bins * bin_size, bin_size))[:n_bins]
        ends = [min(st + bin_size, size) for st in starts]
        chroms_col = [chrom] * n_bins
        mean_bw.addEntries(chroms_col, starts, ends=ends, values=[float(v) for v in fold])'''
new='''        starts = np.arange(0, size, bin_size, dtype=np.int64)
        ends = np.minimum(starts + bin_size, size)
        mask = buffered_evidence[chrom]
        if mask.any():
            mean_bw.addEntries([chrom] * int(mask.sum()), starts[mask].tolist(),
                               ends=ends[mask].tolist(), values=fold[mask].astype(float).tolist())'''
assert old in s;s=s.replace(old,new);p.write_text(s,encoding='utf-8');print('Primary mean writer now preserves no-evidence gaps, matching validated rebuild.')
