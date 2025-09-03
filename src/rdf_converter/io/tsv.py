import csv
from typing import Iterable, Dict, List
from pathlib import Path

def write_tsv(path: str, rows: Iterable[Dict]):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    rows = list(rows)
    if not rows:
        # create empty file with no header
        Path(path).write_text("")
        return
    headers: List[str] = list(rows[0].keys())
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=headers, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)
