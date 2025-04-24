from atom import Atom, Snt

class Snt_scripts(list):
    def __init__(self, *args):
        super().__init__(args)

    def add_Snt(self, obj):
        if isinstance(obj,Snt):
            self.append(Snt)
        else:
            raise TypeError('Object must be of type Snt')