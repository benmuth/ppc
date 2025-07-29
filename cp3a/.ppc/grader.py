#!/usr/bin/env python3

from ppcgrader.cli import cli
import ppccp

if __name__ == "__main__":
    cli(
        ppccp.Config(code='cp3a',
                     gpu=False,
                     openmp=True,
                     single_precision=False,
                     vectorize=True))
