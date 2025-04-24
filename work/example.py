import os
import sys

wdir = os.path.dirname(os.path.abspath(__file__))
root = os.path.abspath(os.path.join(wdir, "../"))
sys.path.append(os.path.join(root, "src"))

from atom import Atom, Snt

He = Atom(2, 2)
Fe_p1 = Atom(26, 25)