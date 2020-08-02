from pathlib import Path
import os

HEADPATH = Path(__file__).parent

def mfp(path):
    """Make relative filepath"""
    return os.path.join(HEADPATH, path)