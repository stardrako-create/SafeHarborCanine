from pathlib import Path
s=Path('regional_fragments_v1217.py').read_text(encoding='utf-8')
s=s.replace("O=Path(r'D:","OLD=Path(r'D:",1).replace(";C=O/'range_cache';C.mkdir(exist_ok=True)",";O=OLD.parent/'06_v1.22.0';C=OLD/'range_cache';C.mkdir(exist_ok=True)",1)
s=s.replace("partial=O/", "partial=OLD/",1)
out=Path(r'D:\Biblioteca\Bioquímica 1º ano 2025-2026\2º Semestre\Vasco M. Barreto\Safe Harbor Canino\06_v1.22.0\regional_fragments_v1220.py');out.write_text(s,encoding='utf-8')
