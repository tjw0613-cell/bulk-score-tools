"""Verify the public archive's file inventory and SHA-256 checksums."""
from pathlib import Path
import argparse
import csv
import hashlib
import json

def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parent.parent)
    parser.add_argument('--strict', action='store_true', help='Also reject unlisted files; omit after generating new outputs.')
    args = parser.parse_args()
    root = args.root.resolve()
    manifest = root / 'FILE_MANIFEST_SHA256.csv'
    rows = list(csv.DictReader(manifest.open(encoding='utf-8-sig')))
    seen, errors = set(), []
    for row in rows:
        name = row['path']
        if name in seen:
            errors.append({'path': name, 'error': 'duplicate_entry'})
            continue
        seen.add(name)
        path = (root / name).resolve()
        if not path.is_relative_to(root) or not path.is_file():
            errors.append({'path': name, 'error': 'missing_or_outside_root'})
        elif path.stat().st_size != int(row['bytes']) or digest(path) != row['sha256']:
            errors.append({'path': name, 'error': 'size_or_hash_mismatch'})
    if args.strict:
        extras = sorted(p.relative_to(root).as_posix() for p in root.rglob('*')
                        if p.is_file() and p != manifest and p.relative_to(root).as_posix() not in seen)
        errors.extend({'path': x, 'error': 'unlisted_file'} for x in extras)
    print(json.dumps({'status': 'PASS' if not errors else 'FAIL', 'files_checked': len(rows), 'errors': errors}, indent=2))
    return 1 if errors else 0

if __name__ == '__main__':
    raise SystemExit(main())
