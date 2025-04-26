import sys
import os

wdir = os.path.dirname(os.path.abspath(__file__))
root = os.path.abspath(os.path.join(wdir, "../"))
sys.path.append(os.path.join(root, "submodules"))

from mylib_python.Nucl import kshell_scripts
import reference

class Atom:
    def __init__(self, Z, Ne):
        """
        Atomic system with proton number Z and N_e electrons
        Z (int) - atomic number (number of protons)
        Ne (int) - number of electrons
        """
        self.element = reference.ELEM[Z]
        self.Z = Z
        self.Ne = Ne
    
    @classmethod
    def from_symbol(cls, symbol, Ne):
        try:
            Z = reference.ELEM.index(symbol.title())
            Ne = Z
            return cls(Z, Ne)
        except ValueError:
            print(f'Element {symbol} not found.')
            return None

    def get_ref(self):
        """
        Get target-space "reference" ARG from element. This is used to generate
        the snt file name for the atomic system. 
        """
        try:
            return reference.element[self.Ne] + str(self.Ne)
        except:
            print(f'Element {self.element} not found.')
            return None
    
    def get_elem(self):
        """
        Get element corresponding to the atom with Z protons. 
        """
        return reference.ELEM[self.Z]

    def get_Z_from_elem(self, ref):
        """
        Get the atomic number from the reference string. 
        """
        try:
            Z = reference.ELEM.index(ref)
            return Z
        except ValueError:
            print(f'Element {ref} not found.')
            return None

    def get_ion_ref(self, q):
        """
        q - charge in units of e < 0
        """
        Z_ion = int(self.get_Z()) - q
        if Z_ion < 0:
            raise ValueError('Ionization charge too large.')
        return reference.element[Z_ion] + str(Z_ion)

    def is_closed(self):
        return (int(self.get_Z()) in [2, 10, 18, 36, 54, 86])

    def get_core(self):
        """
        Return the closed-shell core electrons in the atom.
        """
        for i in [2, 10, 18, 36, 54, 86]:
            if self.Z >= i:
                return i

    def get_valence_occupation(self):
        """
        Get the valence occupation in terms of occupied single-orbital subshells. 

        Return list [c, n_s, n_p, n_d, n_f] where n_s, n_p, n_d, and n_f are the number
        of occupied states in s, p, d, and f shells respectively, and c is the number 
        of electrons in the closed shell core. 
        """
        core = self.get_core()
        Nval = self.Ne - core
        n_s = n_p = n_d = n_f = 0

        if Nval > 0:
            n_s = min(2, Nval)
            Nval -= n_s

        if Nval > 0 and self.Ne >= 56:
            n_f = min(14, Nval)
            Nval -= n_f

        if Nval > 0 and self.Ne >= 20:
            n_d = min(10, Nval)
            Nval -= n_d

        if Nval > 0:
            n_p = min(6, Nval)

        return [core, n_s, n_p, n_d, n_f]

    def get_valence_subshell_l(self):
        """
        Return l for the highest occupied subshell in the atom
        """
        [c, n_s, n_p, n_d, n_f] = self.get_valence_occupation()
        if n_f > 0:
            return 3
        elif n_p > 0:
            return 1
        elif n_d > 0:
            return 2
        elif n_s > 0:
            return 0
        else: return 0
    
    def get_full_shell_occupation(self, l):
        filled = [2, 6, 10, 14]
        return filled[l]
    
    def get_gs_Jp_term(self):
        """
        Return the total angular momentum (2J) and parity (+/-) of the ground state
        based on the shell model in the form used for KSHELL input. 
        """
        occs = self.get_valence_occupation()
        highest_l = self.get_valence_subshell_l()
        print(highest_l)
        filling = self.get_full_shell_occupation(highest_l)
        subshell_occ = occs[highest_l + 1]

        # filling parameters for loop
        atomic_ml = highest_l
        occ_remaining = subshell_occ
        spin = 1
        L = 0
        S = 0

        while occ_remaining > 0:
            occ_remaining -= 1
            L += atomic_ml
            S += spin
            atomic_ml -= 1
            if atomic_ml < -highest_l:
                atomic_ml = highest_l
                spin = -spin

        # calculate J
        if subshell_occ <= filling/2:
            J = abs(S - L)
        else: J = L + S

        # calculate P
        if highest_l % 2 == 0:
            P = '+'
        else:
            if subshell_occ % 2 == 0:
                P = '+'
            else: P = '-'
        term = str(2*J) + '-'
        return term 
    
    def states_list(self):
        """
        Guess states list for valence-space diagonalization
        """
        if self.Ne % 2 == 0: return '0+2,1+2,1-2,2+2,2-2'
        else: return '0.5+1,1.5+1,0.5-1'

class Snt(Atom):
    def __init__(self, Z, Ne, emax, zeta, shell='0hw-shell', hamil='Coulomb', 
    method='magnus', s=500):
        Atom.__init__(self, Z, Ne)
        self.emax = emax
        self.zeta = zeta
        self.shell = shell
        self.hamil = hamil
        self.method = method
        self.s = s
    
    @classmethod
    def from_file(self, fname):
        """
        Using the snt file format {shell}_{hamil}_{args}_{ref}_{element}_{emax}
        _{s}_eta{eta}.snt
        """
        params = fname.split('_')
        print(params)
        print(params[5][1:])
        self.shell = params[0]
        self.hamil = params[1]
        self.method = params[2]
        self.emax = int(params[5][1:])
        self.s = int(params[6][1:])
        self.zeta = float(params[7][3:].replace('.snt', ''))

        ref = params[3]
        Z_elem = params[4]
        
        Atom.__init__(self, self.get_Z_from_elem(Z_elem), reference.ELEM.index(ref))

    def get_snt_fname(self):
        return f'{self.shell}_{self.hamil}_{self.method}_{self.get_ref()}_{self.element}_e{self.emax}_s{self.s}_eta{self.zeta}.snt'
    
    def get_ion_snt_fname(self, q):
        return f'{self.shell}_{self.hamil}_{self.method}_{self.get_ion_ref(q)}_{self.element}_e{self.emax}_s{self.s}_eta{self.zeta}.snt'
    
    def get_summary_fname(self):
        return f'summary_{self.get_ref()}_{self.shell}_{self.hamil}_{self.method}_{self.get_ref()}_{self.element}_e{self.emax}_s{self.s}_eta{self.zeta}.txt'

    def get_ion_summary_fname(self, q):
        return f'''summary_{self.get_ion_ref()}_{self.shell}_{self.hamil}_{self.method}_{self.get_ion_ref(q)}_{self.element}_e{self.emax}_s{self.s}_eta{self.zeta}.txt'''
    
    def diagonalize_valence_space(self, path_to_kshell=''):
        """
        Wrapper for running KSHELL. Option to create a new scratch directory locally
        will be added in the Snt_scripts class. 
        """
        ksh_atom = kshell_scripts(path_to_kshell, self.get_summary_fname(), 
                                  self.get_ref(), self.states_list())
        ksh_atom.run_kshell()
        return None

    def calc_gs(self, path_to_snt='', path_to_sum=''):
        """
        Get the ground state energy of the diagonalized valence-space 
        Hamiltonian in atomic units.
        path_to_snt (str): path to the directory with the snt file with / directory
        separator not including trailing /. 
        path_to_sum (str): path to the directory with the summary file with / 
        directory separator not including trailing /.

        E.g. calc_gs(self, path_to_snt='snt_files', path_to_sum='sum_files')
        """
        with open(f'{path_to_snt}/{self.get_snt_fname()}') as fsnt:
            for _ in range(4):
                next(fsnt)
            E_0 = float(next(fsnt).split()[-1])
        
        if not self.is_closed():
            with open(f'{path_to_sum}/{self.get_summary_fname()}') as fsum:
                for _ in range(5):
                    next(fsum)
                E_1 = float(next(fsum).split()[-3])
                return E_0 + E_1

    def calc_ion_gs(self, path_to_snt='', path_to_sum=''):
        with open(f'{path_to_snt}/{self.get_ion_snt_fname()}') as fsnt:
            for _ in range(4):
                next(fsnt)
            E_0 = float(next(fsnt).split()[-1])
        
        if not self.is_closed():
            with open(f'{path_to_sum}/{self.get_ion_summary_fname()}') as fsum:
                for _ in range(5):
                    next(fsum)
                E_1 = float(next(fsum).split()[-3])
                return E_0 + E_1
            
    def calc_ionization(self, path_to_atom_snt='', path_to_atom_sum='', path_to_ion_snt='',
                        path_to_ion_sum='', eV=True):
        """
        Ionization energy given diagonalized valence spaces. eV=False gives result in 
        atomic units. 
        """
        E_atom = self.calc_gs(path_to_atom_snt, path_to_atom_sum)
        if not self.is_closed():
            E_ion = self.calc_ion_gs(path_to_ion_snt, path_to_ion_sum)
        if eV:
            return (E_ion - E_atom)*reference.au_to_eV
        else:
            return E_ion - E_atom
    
    def get_ex_spectrum(self, path_to_snt='', path_to_sum=''):
        """
        Return excitation energies and terms as list
        """
        return 0
    
He_params = Snt.from_file('0hw-shell_Coulomb_magnus_He2_He_e0_s500_eta0.snt')
print(He_params.get_snt_fname())
print(He_params.get_ion_snt_fname(-1))