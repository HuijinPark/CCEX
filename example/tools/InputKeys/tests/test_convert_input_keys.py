"""Regression tests of convert_input_keys.py (review 20261013-180049 C1).

    python3 -m unittest discover -s example/tools/InputKeys/tests      (or python3 -m pytest example/tools/InputKeys/tests)

fixtures/<name>.old.json are tutorial inputs as they were before A47 (git show 87cd7365a^:<path>), <name>.new.json
what the converter made of them in 87cd7365a, the commit that converted the repository's inputs, and <name>.current.json
what it makes of either today (A64: the configuration block; A65: the execution, cluster, output and frame blocks).
"""
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
TOOL = os.path.join(os.path.dirname(HERE), 'convert_input_keys.py')
FIX = os.path.join(HERE, 'fixtures')
sys.path.insert(0, os.path.dirname(HERE))
from convert_input_keys import convert_text, strip_comments   # noqa: E402

NAMES = sorted(f[:-len('.old.json')] for f in os.listdir(FIX) if f.endswith('.old.json'))


def read(path):
    with open(path) as f:
        return f.read()


def comments(text):
    """The comments of a CCEX input (# or ! to the end of the line, outside strings), stripped."""
    out = []
    for line in text.splitlines():
        code = strip_comments(line)
        if code != line:
            out.append(line[len(code):].strip())
    return out


def run(*args):
    return subprocess.run([sys.executable, TOOL, *args], capture_output=True, text=True)


class Fixtures(unittest.TestCase):
    def test_there_are_fixtures(self):
        self.assertGreaterEqual(len(NAMES), 5)

    def test_converted_as_in_87cd7365a(self):
        """The A47 output, and the inputs before A47, become today's form."""
        for name in NAMES:
            for form in ('old', 'new'):
                with self.subTest(name=name, form=form):
                    new, log = convert_text(read(os.path.join(FIX, f'{name}.{form}.json')))
                    self.assertEqual(new, read(os.path.join(FIX, name + '.current.json')))
                    self.assertTrue(log)

    def test_comments_kept(self):
        for name in NAMES:
            with self.subTest(name):
                old = read(os.path.join(FIX, name + '.old.json'))
                new, _ = convert_text(old)
                self.assertTrue(comments(old))
                self.assertEqual(sorted(comments(old)), sorted(comments(new)))

    def test_second_pass_changes_nothing(self):
        for name in NAMES:
            with self.subTest(name):
                new = read(os.path.join(FIX, name + '.current.json'))
                again, log = convert_text(new)
                self.assertEqual(again, new)
                self.assertEqual(log, [])


class SaveMode(unittest.TestCase):
    """The savemode names before 2026-10-13 and the top-level save_nodiv (left unconverted until 2026-10-14, A57); output.savemode
    since A65."""

    def conv(self, text):
        return json.loads(strip_comments(convert_text(text)[0]))['output']

    def test_names(self):
        for old, new in (('normal', 'coherence'), ('info', 'coherence_tagged'), ('all', 'cluster_coherence'),
                         ('allfull', 'cluster_density_matrix')):
            with self.subTest(old=old):
                self.assertEqual(self.conv('{"method": "CCE", "savemode": "%s"}' % old)['savemode'], new)

    def test_save_nodiv_into_the_block(self):
        d = self.conv('{"method": "CCE", "savemode": "normal", "save_nodiv": false}')
        self.assertEqual(d['savemode'], {'mode': 'coherence', 'nodiv': False})
        self.assertNotIn('save_nodiv', json.loads(convert_text('{"savemode": "normal", "save_nodiv": false}')[0]))
        d = self.conv('{"method": "CCE", "save_nodiv": true}')
        self.assertEqual(d['savemode'], {'mode': 'coherence', 'nodiv': True})

    def test_new_names_left_alone(self):
        new, log = convert_text('{"method": "cCCE", "savemode": "coherence"}')
        self.assertEqual(json.loads(strip_comments(new))['output']['savemode'], 'coherence')
        self.assertFalse([l for l in log if l.startswith('savemode "')])      # the value is left alone (the key moves, A65)

    def test_second_pass_changes_nothing(self):
        new, _ = convert_text('{"method": "CCE", "savemode": "all", "save_nodiv": false}')
        again, log = convert_text(new)
        self.assertEqual(again, new)

    def test_save_nodiv_0_1_refused(self):     # CCEX refuses them (readStrictBool before A47, savemode.nodiv after)
        for v in ('0', '1', '"true"'):
            with self.subTest(v), self.assertRaisesRegex(ValueError, 'save_nodiv must be true or false'):
                convert_text('{"method": "CCE", "save_nodiv": %s}' % v)


class SwitchBlocks(unittest.TestCase):
    """gcce_refinement and cluster_contributions: true|false -> {"enabled": ...}; CCEX took only a JSON boolean (A47); in the
    execution and output blocks since A65."""

    def test_boolean_into_the_block(self):
        for name, block in (('gcce_refinement', 'execution'), ('cluster_contributions', 'output')):
            for v in (True, False):
                with self.subTest(name=name, v=v):
                    new, _ = convert_text('{"%s": %s}' % (name, json.dumps(v)))
                    self.assertEqual(json.loads(new), {block: {name: {'enabled': v}}})

    def test_0_1_refused(self):
        for name in ('gcce_refinement', 'cluster_contributions'):
            for v in ('0', '1'):
                with self.subTest(name=name, v=v), self.assertRaisesRegex(ValueError, f'{name} must be true or false'):
                    convert_text('{"%s": %s}' % (name, v))


class OldAndNewKeyTogether(unittest.TestCase):
    """An old key beside its new name was written twice (`"count":1,"count":2`): CCEX reads the first, Python the last."""

    def test_different_values_refused(self):
        for text in ('{"npulse": 2, "pulse": {"count": 1}}',
                     '{"npulse": 2, "pulse": {"Count": 1}}',
                     '{"statefile": "a", "bath_state_file": "b"}',
                     '{"Adaptive": {"policy": "merge"}, "adaptive": {"policy": "hub"}}',
                     '{"Qubit": {"qubit": [{"name": "q", "Spin": 1, "spin": 0.5}]}}',
                     '{"addsubclus": 1, "add_subclusters": 1}'):          # 1 -> true, beside a 1 that CCEX refuses
            with self.subTest(text), self.assertRaisesRegex(ValueError, 'written twice with different values'):
                convert_text(text)

    def test_conflict_names_the_key(self):
        with self.assertRaisesRegex(ValueError, r'pulse\.count .*\(1 and 2\)'):
            convert_text('{"npulse": 2, "pulse": {"count": 1}}')
        with self.assertRaisesRegex(ValueError, r'Qubit\.qubit\[0\]\.Spin'):
            convert_text('{"Qubit": {"qubit": [{"name": "q", "Spin": 1, "spin": 0.5}]}}')

    def test_same_value_kept_once(self):
        new, log = convert_text('{"npulse": 2, "pulse": {"count": 2}}')
        self.assertEqual(json.loads(new), {'pulse': {'count': 2}})
        self.assertEqual(new.count('count'), 1)
        self.assertIn('pulse.count: the same value was written twice, kept once', log)
        new, _ = convert_text('{"Adaptive": {"policy": "merge"}, "adaptive": {"policy": "merge"}}')
        self.assertEqual(json.loads(new), {'cluster': {'adaptive': {'policy': 'merge'}}})

    def test_same_value_keeps_the_comment(self):
        old = '{\n    "npulse": 2,     # Hahn echo\n    "pulse": {\n        "count": 2\n    }\n}\n'
        new, _ = convert_text(old)
        self.assertEqual(json.loads(strip_comments(new)), {'pulse': {'count': 2}})
        self.assertEqual(comments(new), ['# Hahn echo'])


class Configuration(unittest.TestCase):
    """A64 (2026-10-15): gyrofile, bathfile, bath_shift, bath_state_file, nstate and seed in the configuration block, and
    Defect.axis_file and Defect.state_file as its defect_axis_file and defect_state_file; the state and axis files one per
    bath file."""

    def conv(self, text):
        new, log = convert_text(text)
        return json.loads(strip_comments(new)), new, log

    def test_into_the_block(self):
        d, new, _ = self.conv('{\n "method": "cCCE",\n # the bath\n "gyrofile": "g",\n "bathfile": ["b"],   # one file\n'
                              ' "bath_shift": [[0, 0, 1]],\n "bath_state_file": "s_",\n "order": 2,\n "nstate": 3,\n'
                              ' "seed": 7\n}\n')
        self.assertEqual(d, {'method': 'cCCE', 'order': 2,
                             'configuration': {'gyrofile': 'g', 'bathfile': ['b'], 'bath_shift': [[0, 0, 1]],
                                               'bath_state_file': ['s_'], 'nstate': 3, 'seed': 7}})
        self.assertEqual(comments(new), ['# the bath', '# one file'])
        self.assertEqual(convert_text(new), (new, []))

    def test_names_before_a47_too(self):
        d, _, _ = self.conv('{"bathfile": ["b"], "statefile": "s_", "bathadjust": [[0, 0, 0]], "Defect": [{"dfname": "P1"}],'
                            ' "avaaxfile": "ax", "exstatefile": "ex_"}')
        self.assertEqual(d['configuration'], {'bathfile': ['b'], 'bath_shift': [[0, 0, 0]], 'bath_state_file': ['s_'],
                                              'defect_axis_file': ['ax'], 'defect_state_file': ['ex_']})
        self.assertEqual(d['Defect'], {'species': [{'name': 'P1'}]})

    def test_one_name_per_bath_file(self):
        d, _, log = self.conv('{"bathfile": ["a", "b"], "bath_state_file": "s_", "Defect": {"axis_file": "ax", '
                              '"state_file": ["e1_", "e2_"], "species": []}}')
        # with two bath files CCEX did not read a single state or axis file (it drew them at random): removed
        self.assertNotIn('bath_state_file', d['configuration'])
        self.assertNotIn('defect_axis_file', d['configuration'])
        self.assertEqual(d['configuration']['defect_state_file'], ['e1_', 'e2_'])
        self.assertEqual(d['Defect'], {'species': []})
        self.assertTrue([l for l in log if 'bath_state_file removed' in l])

    def test_defect_files_into_the_block(self):
        d, new, log = self.conv('{"configuration": {"bathfile": ["b"]},\n "Defect": {\n  # the axes\n  "axis_file": ["a_{k}"],\n'
                                '  "species": []}}')
        self.assertEqual(d, {'configuration': {'bathfile': ['b'], 'defect_axis_file': ['a_{k}']}, 'Defect': {'species': []}})
        self.assertEqual(comments(new), ['# the axes'])
        self.assertIn('Defect.axis_file -> configuration.defect_axis_file', log)
        self.assertEqual(convert_text(new), (new, []))

    def test_new_form_left_alone(self):
        text = '{"configuration": {"path": "data", "k_range": [1, 3], "bathfile": ["c_{k}"], "seed": 1}}'
        self.assertEqual(convert_text(text), (text, []))


class Blocks(unittest.TestCase):
    """A65 (2026-10-15): execution, cluster (with selection), output, frame; the Defect species keys; the auto_convergence
    block without launcher, its configurations as sampling and its file patterns in the configuration block."""

    def conv(self, text):
        new, log = convert_text(text)
        return json.loads(strip_comments(new)), new, log

    def test_into_the_blocks(self):
        d, new, _ = self.conv('{\n "method": "cCCE",\n "cluster_mode": "pcce",   # partitions\n "backend": "auto",\n'
                              ' "gpu_options": {"min_dim": 64},\n "propagator": "eigen",\n "evolution_update": "exact",\n'
                              ' "profile": true,\n "add_subclusters": false,\n "max_clusters_per_order": [0, 0, 0],\n'
                              ' "pcce": {"partition_size": 2},\n "outfile": "out",\n "savemode": "coherence",\n'
                              ' "cluster_only": false,\n "order": 2\n}\n')
        self.assertEqual(d, {'method': 'cCCE', 'order': 2,
                             'execution': {'backend': 'auto', 'gpu_options': {'min_dim': 64}, 'propagator': 'eigen',
                                           'evolution_update': 'exact', 'profile': True},
                             'cluster': {'mode': 'pcce', 'selection': {'add_subclusters': False, 'max_per_order': [0, 0, 0]},
                                         'pcce': {'partition_size': 2}},
                             'output': {'file': 'out', 'savemode': 'coherence', 'cluster_only': False}})
        self.assertEqual(comments(new), ['# partitions'])
        self.assertEqual(convert_text(new), (new, []))

    def test_probe_memory_into_the_adaptive_block(self):
        d, _, _ = self.conv('{"cluster_mode": "bundle", "adaptive": {"output": "a"}, "probe_memory_gb": 2}')
        self.assertEqual(d['cluster']['adaptive'], {'output': 'a', 'probe_memory_gb': 2})
        d, _, log = self.conv('{"probe_memory_gb": 2}')
        self.assertEqual(d, {})
        self.assertTrue([l for l in log if 'probe_memory_gb removed' in l])

    def test_frame(self):
        on = '{"coordinate_frame_rotation": {"enabled": true, "qubit_axis": [1, 1, 1]}}'
        self.assertEqual(self.conv(on)[0], {'frame': {'rotation': {'enabled': True, 'qubit_axis': [1, 1, 1]}}})
        for off in ('{"coordinate_frame_rotation": {"enabled": false, "qubit_axis": [1, 1, 1]}}',
                    '{"coordinate_frame_rotation": {"qubit_axis": [1, 1, 1]}}',      # enabled was false by default
                    '{"bath_coordinate_rotation": {"enabled": false}}'):
            with self.subTest(off):
                self.assertEqual(self.conv(off)[0], {})
        self.assertEqual(self.conv('{"defect_axis_reference": "qubit_axis"}')[0], {})

    def test_defect_species(self):
        d, _, _ = self.conv('{"Defect": {"species": [{"name": "P1", "n_nuclear_spins": 1, "n_configurations": 4}]}}')
        self.assertEqual(d['Defect']['species'], [{'name': 'P1', 'n_onsite_spins': 1, 'n_axis_states': 4}])

    def test_auto_convergence(self):
        d, _, log = self.conv('{"auto_convergence": {"workdir": "w", "launcher": ["mpirun", "-n", "4"], '
                              '"configurations": {"bathfile": "data/c_{k}", "statefile": "data/s_{k}_", '
                              '"avaaxfile": "data/a_{k}", "exstatefile": "data/e_{k}_", "number": 10}}}')
        self.assertEqual(d, {'auto_convergence': {'workdir': 'w', 'sampling': {'number': 10}},
                             'configuration': {'bathfile': ['data/c_{k}'], 'bath_state_file': ['data/s_{k}_'],
                                               'defect_axis_file': ['data/a_{k}'], 'defect_state_file': ['data/e_{k}_']}})
        self.assertTrue([l for l in log if 'launcher removed' in l])


class CommandLine(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.dir)

    def copy(self, name):
        path = os.path.join(self.dir, name + '.json')
        shutil.copy(os.path.join(FIX, name + '.old.json'), path)
        return path

    def write(self, text):
        path = os.path.join(self.dir, 'input.json')
        with open(path, 'w') as f:
            f.write(text)
        return path

    def test_rewrites_in_place(self):
        path = self.copy(NAMES[0])
        r = run(path)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(read(path), read(os.path.join(FIX, NAMES[0] + '.current.json')))
        r = run(path)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn('already in the new form', r.stdout)

    def test_check_and_stdout_leave_the_file(self):
        for name in NAMES:
            with self.subTest(name):
                path = self.copy(name)
                before, stamp = read(path), os.stat(path).st_mtime_ns
                r = run('--check', path)
                self.assertEqual(r.returncode, 1)                 # it still has old keys
                r = run('--stdout', path)
                self.assertEqual(r.returncode, 0, r.stderr)
                self.assertEqual(r.stdout, read(os.path.join(FIX, name + '.current.json')))
                self.assertEqual(read(path), before)
                self.assertEqual(os.stat(path).st_mtime_ns, stamp)

    def test_check_of_a_new_input(self):
        path = os.path.join(self.dir, 'new.json')
        shutil.copy(os.path.join(FIX, NAMES[0] + '.current.json'), path)
        r = run('--check', path)
        self.assertEqual(r.returncode, 0)
        self.assertIn('already in the new form', r.stdout)

    def test_conflict_leaves_the_file(self):
        text = '{\n    "npulse": 2,\n    "pulse": {"count": 1}\n}\n'
        for option in ((), ('--check',), ('--stdout',)):
            with self.subTest(option):
                path = self.write(text)
                r = run(*option, path)
                self.assertEqual(r.returncode, 2)
                self.assertIn('pulse.count is written twice with different values', r.stderr)
                self.assertEqual(read(path), text)

    def test_unreadable_file_reported_and_the_others_converted(self):
        missing = os.path.join(self.dir, 'missing.json')
        latin = os.path.join(self.dir, 'latin.json')
        with open(latin, 'wb') as f:
            f.write(b'{"npulse": 1  # \xb0\xa1\n}\n')               # a cp949 comment: not UTF-8
        path = self.copy(NAMES[0])
        r = run(missing, latin, path)
        self.assertEqual(r.returncode, 2)
        self.assertIn('missing.json: not converted', r.stderr)
        self.assertIn('latin.json: not converted', r.stderr)
        self.assertNotIn('Traceback', r.stderr)
        self.assertEqual(read(path), read(os.path.join(FIX, NAMES[0] + '.current.json')))

    def test_crlf_kept(self):
        old = read(os.path.join(FIX, NAMES[0] + '.old.json')).replace('\n', '\r\n')
        path = os.path.join(self.dir, 'crlf.json')
        with open(path, 'w', newline='') as f:
            f.write(old)
        r = run(path)
        self.assertEqual(r.returncode, 0, r.stderr)
        with open(path, 'rb') as f:
            data = f.read()
        self.assertEqual(data.count(b'\n'), data.count(b'\r\n'))          # no bare \n, also on the lines the converter added
        self.assertEqual(data.decode().replace('\r\n', '\n'), read(os.path.join(FIX, NAMES[0] + '.current.json')))
        r = run('--check', path)
        self.assertEqual(r.returncode, 0)
        self.assertIn('already in the new form', r.stdout)

    def test_write_is_atomic(self):
        import convert_input_keys as tool
        from unittest import mock
        path = self.copy(NAMES[0])
        os.chmod(path, 0o640)
        before = read(path)
        with mock.patch.object(tool.os, 'replace', side_effect=OSError('disk full')), \
                mock.patch('sys.stdout'), mock.patch('sys.stderr'):
            self.assertEqual(tool.main([path]), 2)
        self.assertEqual(read(path), before)                                 # the original is whole
        self.assertEqual(os.listdir(self.dir), [os.path.basename(path)])     # and no temporary file is left
        with mock.patch('sys.stdout'):
            self.assertEqual(tool.main([path]), 0)
        self.assertEqual(read(path), read(os.path.join(FIX, NAMES[0] + '.current.json')))
        self.assertEqual(os.stat(path).st_mode & 0o777, 0o640)              # the mode of the original
        self.assertEqual(os.listdir(self.dir), [os.path.basename(path)])


if __name__ == '__main__':
    unittest.main()
