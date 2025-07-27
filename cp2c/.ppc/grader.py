#!/usr/bin/env python3

from ppcgrader.cli import cli
import ppccp

if __name__ == "__main__":
    cli(
        ppccp.Config(code='cp2c',
                     gpu=False,
                     openmp=False,
                     single_precision=False,
                     vectorize=False))
