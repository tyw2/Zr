#!/usr/bin/env python3


import os
import shutil
from pathlib import Path
from zr_validation.settings import *

def pair_block(pot_path: str, table_points: int = 1000) -> str:
    return f"""pair_style  hybrid/overlay mlip load_from={pot_path} table linear {table_points}
pair_coeff  * * mlip
pair_coeff  1 1 table {TABLE_POTENTIAL} {TABLE_NAME}"""

def write_file(path: Path, content: str, executable: bool = False):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")
    if executable:
        path.chmod(0o755)

def write_file_if_missing(path: Path, content: str, executable: bool = False):

    if path.exists():
        return
    write_file(path, content, executable)

def make_symlink(src: Path, dst: Path):
    if dst.exists() or dst.is_symlink():
        return
    dst.parent.mkdir(parents=True, exist_ok=True)
    os.symlink(src, dst)

def copy_dir_replace(src: Path, dst: Path):


    if not src.exists():
        print(f"    [warn] {src} not found, replica copy skipped.")
        return
    if dst.exists():
        shutil.rmtree(dst)
    shutil.copytree(src, dst)

