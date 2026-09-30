"""The desktop launcher, started as `python -m adcp_launcher`.
"""
import sys

from .entry import main

if __name__ == "__main__":
    sys.exit(main())
