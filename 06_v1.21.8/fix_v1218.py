from pathlib import Path
P=Path(r'D:\Biblioteca\Bioquímica 1º ano 2025-2026\2º Semestre\Vasco M. Barreto\Safe Harbor Canino');S=P/'scripts'
p=S/'score_ship_candidates_v2.py';s=p.read_text(encoding='utf-8');old='''        if c["final_score"] is None:
            missing.append("score_component")''';new='''        # A partial weighted score is descriptive, not evidence that every
        # requested component was measured. Keep each absent component explicit.
        for component, weight in weights.items():
            if weight > 0 and c.get(f"score_{component}") is None:
                missing.append(f"score_component:{component}")
        if c["final_score"] is None and not any(x.startswith("score_component:") for x in missing):
            missing.append("score_component")''';assert old in s;s=s.replace(old,new);p.write_text(s,encoding='utf-8')
p=S/'build_mother_track_v2.py';s=p.read_text(encoding='utf-8');s=s.replace('''        return np.zeros(n_bins, dtype=np.float32)''','''        raise RuntimeError(f"{sample_name}: missing chromosome {chrom}; cannot substitute measured zero")''',1)
start=s.index('            except RuntimeError:\n                print(f"  [warn]');end=s.index('    arr = np.nan_to_num',start)
s=s[:start]+'''            except RuntimeError as exc:
                raise RuntimeError(
                    f"{sample_name}: unreadable {chrom}:{chunk_start}-{chunk_end}; "
                    "rebuild required, zero substitution forbidden"
                ) from exc
'''+s[end:]
start=s.index('            except RuntimeError as e:\n');end=s.index('\n        raw_mat',start)
s=s[:start]+'''            except RuntimeError as e:
                raise RuntimeError(f"ATAC build stopped for {s} on {chrom}") from e
'''+s[end:]
s=s.replace('''        # Localize instead: retry chunk-by-chunk, zero only the specific
        # chunks that individually fail.''','''        # Retry chunk-by-chunk, but abort if a block remains unreadable.
        # I/O failure must never masquerade as measured zero.''')
p.write_text(s,encoding='utf-8');print('Patched missing-component status and fail-closed ATAC reads.')
