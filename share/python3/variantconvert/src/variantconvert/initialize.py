# -*- coding: utf-8 -*-
"""
Use example for testing
variantconvert init
"""

import glob
import logging as log
import os
import variantconvert
import shutil

from os.path import join as osj


def main_init(args):
    if args.dir != None:
        local_dir = args.dir
    else:
        local_dir = variantconvert.__local_config__

    log.info(f"Copying default configs to: {local_dir}")

    if not os.path.isdir(local_dir):
        os.makedirs(local_dir)

    for genome in glob.glob(osj(variantconvert.__default_config__, "**")):
        basename = os.path.basename(genome)
        if basename not in ("__init__.py", "__pycache__"):
            shutil.copytree(
                genome, osj(local_dir, basename), ignore=shutil.ignore_patterns("__init__.py")
            )
    log.info("Done.")
