#!/usr/bin/env python3
"""
patch_ghost_dlr.py  — fixes two bugs in astro_ghost DLR.py
Run once, in the environment where astro_ghost is installed:  python patch_ghost_dlr.py
"""
import sys
import shutil
from pathlib import Path

try:
    import astro_ghost
    dlr_path = Path(astro_ghost.__file__).parent / "DLR.py"
except ImportError:
    print("ERROR: astro_ghost not importable. Activate the environment where astro_ghost is installed.")
    sys.exit(1)

backup = dlr_path.with_suffix(".py.bak")
if not backup.exists():
    shutil.copy2(dlr_path, backup)
    print(f"Backed up original to {backup}")
else:
    shutil.copy2(backup, dlr_path)
    print("Restored from backup before patching")

content = dlr_path.read_text()
original = content
n = 0

# -----------------------------------------------------------------------
# Patch 1 — line ~213
# float(dist/DLR) crashes when SkyMapper returns multiple rows
# -----------------------------------------------------------------------
OLD1 = "    R = float(dist/DLR)\n"
NEW1 = "    R = float(np.atleast_1d(dist/DLR).ravel()[0])\n"

if OLD1 in content:
    content = content.replace(OLD1, NEW1, 1)
    n += 1
    print("Patch 1 applied")
elif NEW1 in content:
    print("Patch 1 already applied")
else:
    print("WARNING: Patch 1 target not found around line 213")

# -----------------------------------------------------------------------
# Patch 2 — line 304
# Exact indentation is 32 spaces (8 levels of 4-space indent)
# -----------------------------------------------------------------------
OLD2 = '                                elong = host_df[band + "_elong"].values[0]\n'
NEW2 = (
    '                                _elong_col = band + "_elong"\n'
    '                                elong = host_df[_elong_col].values[0] if _elong_col in host_df.columns else 1.0\n'
)

if OLD2 in content:
    content = content.replace(OLD2, NEW2, 1)
    n += 1
    print("Patch 2 applied")
elif "_elong_col" in content:
    print("Patch 2 already applied")
else:
    print("WARNING: Patch 2 target not found around line 304")

if content != original:
    dlr_path.write_text(content)
    print(f"\nWrote {n} patch(es) to {dlr_path}")
    try:
        import ast
        ast.parse(content)
        print("Syntax check PASSED.")
    except SyntaxError as e:
        print(f"Syntax check FAILED: {e}")
        print("Restoring backup...")
        shutil.copy2(backup, dlr_path)
        sys.exit(1)
else:
    print("\nNo changes written.")

print("Done.")