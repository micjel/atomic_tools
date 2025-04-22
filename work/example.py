import os
import sys

wdir = os.path.dirname(os.path.abspath(__file__))
root = os.path.abspath(os.path.join(wdir, "../"))
sys.path.append(os.path.join(root, "src"))

import atomic_tools

He1 = atomic_tools.Atom(2, 2)