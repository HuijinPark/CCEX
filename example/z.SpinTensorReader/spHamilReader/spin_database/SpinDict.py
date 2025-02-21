import copy
import warnings
from collections.abc import Mapping
from collections import UserDict
import numpy as np
from numpy.lib.recfunctions import repack_fields
import sys
from spHamilReader.spin_database.constants import HBAR, ELECTRON_GYRO, HBAR_SI, NUCLEAR_MAGNETON, PI2

class SpinType:
    r"""
    Class which contains properties of each spin type in the bath.

    Args:
        name (str): Name of the bath spin.
        s (float): Total spin of the bath spin.

            Default 0.

        gyro (float): Gyromagnetic ratio in rad * kHz / G.

            Default 0.

        q (float): Quadrupole moment in barn (for s > 1/2).

            Default 0.

        detuning (float): Energy detuning from the zeeman splitting in kHz,
            included as an extra :math:`+\omega \hat S_z` term in the Hamiltonian,
            where :math:`\omega` is the detuning.

            Default 0.


    Attributes:

        name (str): Name of the bath spin.
        s (float): Total spin of the bath spin.
        dim (int): Spin dimensionality = 2s + 1.

        gyro (float): Gyromagnetic ratio in rad/(ms * G).
        q (float): Quadrupole moment in barn (for s > 1/2).
        detuning (float): Energy detuning from the zeeman splitting in kHz.

    """

    def __init__(self, name, s=0., gyro=0., q=0., detuning=0.):

        self.name = name
        self.s = s

        try:
            self.dim = int(2 * s + 1 + 1e-8)

        except TypeError:
            self.dim = (2 * s + 1 + 1e-8).astype(np.int32)

        self.gyro = gyro
        self.q = q
        self.detuning = detuning

    def __eq__(self, obj):
        if not isinstance(obj, SpinType):
            return False

        checks = (self.name == obj.name) & (self.s == obj.s) & (
                self.gyro == obj.gyro) & (self.q == obj.q) & (self.detuning == obj.detuning)

        return checks

    def __repr__(self):
        try:
            base_message = f'{self.name}: ({self.s:.1f}, {self.gyro:.4f}'
        except TypeError:
            base_message = f'{self.name}: ({self.s}, {self.gyro}'

        if np.asarray(self.q).any():
            try:
                m = f', {self.q:.4f}'
            except TypeError:
                m = f', {self.q}'
            base_message += m

        if np.asarray(self.detuning).any():
            try:
                m = f', {self.detuning:.4f}'
            except TypeError:
                m = f', {self.detuning}'
            base_message += m

        base_message += ')'

        return base_message


class SpinDict(UserDict):
    """
    Wrapper class for dictionary tailored for containing properties of the spin types.
    Can take ``np.void`` or ``BathArray`` instances as keys.
    Every entry is instance of the ``SpinType``.

    Each entry of the ``SpinDict`` can be initianlized as follows:

        * As a Tuple containing name (optional), spin, gyromagnetic ratio, quadrupole constant (optional)
          and detuning (optional).
        * As a ``SpinType`` instance.

    Examples:

        >>> types = SpinDict()
        >>> types['1H'] = ('1H', 1 / 2, 26.7519)
        >>> types['2H'] = 1, 4.1066, 0.00286
        >>> types['3H'] = SpinType('3H', 1 / 2, 28.535, 0)
        >>> print(types)
        SpinDict({'1H': (1H, 0.5, 26.7519, 0.0), '2H': (2H, 1, 4.1066, 0.00286), '3H': (3H, 0.5, 28.535, 0)})

    If ``SpinType`` of the given bath spin is not provided, when requested
    ``SpinDict`` will try to find information about the bath spins in the ``common_isotopes``.

    If found, adds an entry to the given ``SpinDict`` instance and
    returns it. Otherwise ``KeyError`` is raised.

    To initiallize several ``SpinType`` entries one can use ``add_types`` method.

    Args:
        *args: Any numbers of arguments which could initialize ``SpinType`` instances.
        **kwargs: Any numbers of keyword arguments which could initialize ``SpinType`` instances.
            For details see ``SpinDict.add_type`` method.

    """

    def __init__(self, *args, **kwargs):
        super(SpinDict, self).__init__()
        self.add_type(*args, **kwargs)

    def __delitem__(self, key):
        try:
            super().__delitem__(key)
        except TypeError:
            if key.shape:
                try:
                    names = key['N']
                except IndexError:
                    names = key
                for k in names:
                    super().__delitem__(k)
                return

            k = key[()]
            try:
                k = k['N']
            except TypeError:
                pass
            super().__delitem__(k)

    def __setitem__(self, key, value):
        try:
            value = _check_key_spintype(key, value)
            super().__setitem__(key, value)

        except (TypeError, ValueError):
            key = np.asarray(key)
            if key.shape:
                try:
                    names = key['N']
                except IndexError:
                    names = key

                for k, v in zip(names, value):
                    v = _check_key_spintype(k, v)
                    super().__setitem__(k, v)
                return

            k = key[()]

            try:
                k = k['N']
            except TypeError:
                pass

            value = _check_key_spintype(k, value)
            super().__setitem__(k, value)

    def __getitem__(self, key):

        try:
            key = key[()]
            return self._super_get_item(key['N'])
            # self._super_get_item(key)

        except (TypeError, IndexError):
            try:
                return self._super_get_item(key)
            except TypeError:

                if key.dtype.names:
                    key = key['N']

                unique_names = np.unique(key)

                if unique_names.size == 1:
                    n = unique_names[0]
                    # ones = np.ones(key.shape, dtype=np.float64)

                    spins = self._super_get_item(n).s  # * ones
                    gyros = self._super_get_item(n).gyro  # * ones
                    quads = self._super_get_item(n).q  # * ones
                    detus = self._super_get_item(n).detuning  # * ones

                else:
                    spins = np.empty(key.shape, dtype=np.float64)
                    gyros = np.empty(key.shape, dtype=np.float64)
                    quads = np.empty(key.shape, dtype=np.float64)
                    detus = np.empty(key.shape, dtype=np.float64)

                    for n in unique_names:
                        spins[key == n] = self._super_get_item(n).s
                        gyros[key == n] = self._super_get_item(n).gyro
                        quads[key == n] = self._super_get_item(n).q
                        detus[key == n] = self._super_get_item(n).detuning

                return SpinType(key, s=spins, gyro=gyros, q=quads, detuning=detus)

                # params = {}
                # unique_names = np.unique(key)
                # first_name = unique_names[0]
                # first_dict = vars(self._super_get_item(first_name))
                #
                # if unique_names.size == 1:
                #     for k in first_dict:
                #         params[k] = np.array([first_dict[k]]*key.size).reshape(key.shape)
                #     return SpinType(**params)
                #
                # for k in first_dict:
                #     params[k] = np.empty(key.shape, dtype=type(first_dict[k]))
                #     params[k][key == first_name] = getattr(self._super_get_item(first_name), k)
                #
                # for n in unique_names[1:]:
                #     for k in params:
                #         params[k][key == n] = getattr(self._super_get_item(n), k)
                # return SpinType(**params)

    # adding two objects
    def __add__(self, obj):
        new_obj = SpinDict()
        keys_1 = list(self.keys())
        keys_2 = list(obj.keys())

        for k in {*keys_1, *keys_2}:

            if (k in keys_1) and (k in keys_2):
                assert obj[k] == self[k], f'Error, type {k} has different properties in provided types'
                new_obj[k] = self[k]

            elif k in keys_1:
                new_obj[k] = self[k]
            else:
                new_obj[k] = obj[k]
        return new_obj

    def __repr__(self):
        message = f"{type(self).__name__}("
        for k in self.data:
            message += f"{self.data[k]}, "
            if len(message) > 75:
                message += '..., '
                break

        message = message[:-2] + ')'

        return message

    def _super_get_item(self, n):
        try:
            return super().__getitem__(n)
        except KeyError:
            if not n in common_isotopes:
                raise KeyError(_spin_not_found_message(n))
            super().__setitem__(n, copy.copy(common_isotopes[n]))
            return super().__getitem__(n)

    def add_type(self, *args, **kwargs):
        """
        Add one or several spin types to the spin dictionary.

        Args:
            *args:
                Any numbers of arguments which could initialize ``SpinType`` instances.
                Accepted arguments:

                    * Tuple containing name, spin, gyromagnetic ratio, quadrupole constant (optional)
                      and detuning (optional).
                    * ``SpinType`` instance.

                Can also initialize one instance of ``SpinType`` if each argument corresponds to
                each positional argument necessary to initiallize.

            **kwargs: Any numbers of keyword arguments which could initialize ``SpinType`` instances.
                Usefull as an alternative for updating the dictionary. for each keyword argument adds an entry
                to the ``SpinDict`` with the same name as keyword.

        Examples:
            >>> types = SpinDict()
            >>> types.add_type('1H', 1 / 2, 26.7519)
            >>> types.add_type(('1H_det', 1 / 2, 26.7519, 10), ('2H', 1, 4.1066, 0.00286),
            >>>                 SpinType('3H', 1 / 2, 28.535, 0), e=(1 / 2, 6.7283, 0))
            >>> print(types)
            SpinDict(1H: (1H, 0.5, 26.7519), 1H_det: (1H_det, 0.5, 26.7519, 10),
            2H: (2H, 1, 4.1066, 0.00286), 3H: (3H, 0.5, 28.535), e: (e, 0.5, 6.7283))

        """
        keys = []
        try:
            for nuc in args:
                if isinstance(nuc, SpinType):
                    key = nuc.name
                    self[key] = nuc
                    keys.append(key)
                elif isinstance(nuc, Mapping):
                    self.update(nuc)
                    keys += list(nuc.keys())
                else:
                    key = nuc[0]
                    self[key] = SpinType(*nuc)
                    keys.append(key)
        except TypeError:
            for k in keys:
                self.pop(k)
            self[args[0]] = SpinType(*args)
        for nuc in kwargs:
            self[nuc] = kwargs[nuc]


def _check_key_spintype(k, v):
    if isinstance(v, SpinType):
        return v

    if v[0] == k:
        v = SpinType(*v)
    else:
        v = SpinType(k, *v)
    return v




_spin_not_found_message = lambda x: 'Spin type for {} was not provided and was not found in common isotopes.'.format(x)
