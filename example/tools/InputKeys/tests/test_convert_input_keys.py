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
    """The comments of a CCEX input (#, ! or // to the end of the line, outside strings), stripped."""
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
        new, _ = convert_text('{"Adaptive": {"atol": 1}, "adaptive": {"atol": 1}}')
        self.assertEqual(json.loads(new), {'cluster': {'adaptive': {'atol': 1}}})

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
                                               'bath_state_file': ['s_'], 's_range': [1, 3], 'seed': 7}})
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
        self.assertEqual(d['cluster']['adaptive'], {'workdir': 'a', 'probe_memory_gb': 2})
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

    def test_frame_keys_before_a65(self):
        """The flat coordinate_frame_rotation with the A90 keys inside goes to frame.rotation and then on to the frames."""
        d, _, _ = self.conv('{"hf_readmode": 2, "hf_tensorfile": "A", "coordinate_frame_rotation": {"enabled": true, '
                            '"bath_axis": [0, 0, 1], "qubit_axis": [1, 1, 1], "hf_tensor_frame": "bath"}}')
        self.assertEqual(d['frame'], {'rotation': {'enabled': True, 'qubit_axis': [1, 1, 1]}})
        self.assertEqual(d['hamiltonian']['common']['hyperfine']['frame'], 'file')

    def test_defect_species(self):
        d, _, _ = self.conv('{"Defect": {"species": [{"name": "P1", "n_nuclear_spins": 1, "n_configurations": 4}]}}')
        self.assertEqual(d['Defect']['species'], [{'name': 'P1', 'n_onsite_spins': 1, 'n_axis_states': 4}])

    def test_auto_convergence(self):
        d, _, log = self.conv('{"auto_convergence": {"workdir": "w", "launcher": ["mpirun", "-n", "4"], '
                              '"configurations": {"bathfile": "data/c_{k}", "statefile": "data/s_{k}_", '
                              '"avaaxfile": "data/a_{k}", "exstatefile": "data/e_{k}_", "number": 10}}}')
        self.assertEqual(d, {'auto_convergence': {'workdir': 'w', 'sampling': {'configurations': 10}},
                             'configuration': {'bathfile': ['data/c_{k}'], 'bath_state_file': ['data/s_{k}_'],
                                               'defect_axis_file': ['data/a_{k}'], 'defect_state_file': ['data/e_{k}_']}})
        self.assertTrue([l for l in log if 'launcher removed' in l])


class SlashComments(unittest.TestCase):
    """A92 (2026-10-16): // comments (.jsonc) are read like # and !, and stay where they were; // inside a string is text."""

    def conv(self, text):
        new, log = convert_text(text)
        return json.loads(strip_comments(new)), new, log

    def test_old_keys(self):
        text = ('{\n  // the method\n  "method": "cCCE",\n  "cluster_mode": "pcce",   // partitions\n'
                '  "pcce": {"partition_size": 2},\n  // where it goes\n  "outfile": "out",\n  "order": 2   // last\n}\n')
        d, new, log = self.conv(text)
        self.assertEqual(d, {'method': 'cCCE', 'order': 2, 'cluster': {'mode': 'pcce', 'pcce': {'partition_size': 2}},
                             'output': {'file': 'out'}})
        self.assertTrue(log)
        self.assertEqual(sorted(comments(new)), sorted(comments(text)))
        self.assertIn('  // the method\n  "method"', new)
        self.assertIn('"order": 2   // last\n', new)
        self.assertEqual(convert_text(new), (new, []))

    def test_inside_strings(self):
        text = '{\n  "outfile": "http://host//out",   // a // comment\n  "bathfile": ["a//b.txt"]\n}\n'
        d, new, _ = self.conv(text)
        self.assertEqual(d['output']['file'], 'http://host//out')
        self.assertIn('"a//b.txt"', new)
        self.assertEqual(comments(new), ['// a // comment'])

    def test_mixed_with_hash(self):
        text = ('{\n  # hash above\n  "cluster_mode": "pcce",   // slash after\n  // slash above\n'
                '  "outfile": "out",   # hash after\n  ! bang above\n  "order": 2\n}\n')
        d, new, _ = self.conv(text)
        self.assertEqual(d, {'order': 2, 'cluster': {'mode': 'pcce'}, 'output': {'file': 'out'}})
        self.assertEqual(sorted(comments(new)), sorted(comments(text)))
        self.assertEqual(len(comments(text)), 5)

    def test_new_form_left_byte_for_byte(self):
        text = '{\n  // c\n  "method": "cCCE",   // m\n  "output": {"file": "http://x"}   // o\n}\n'
        self.assertEqual(convert_text(text), (text, []))


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

    def test_jsonc_name(self):
        path = os.path.join(self.dir, 'input.jsonc')
        with open(path, 'w') as f:
            f.write('{\n  // old keys\n  "cluster_mode": "pcce",   // partitions\n  "outfile": "out"\n}\n')
        r = run('--check', path)
        self.assertEqual(r.returncode, 1, r.stderr)
        r = run(path)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(read(path), '{\n  "cluster": {\n      // old keys\n      "mode": "pcce"   // partitions\n  },\n'
                                     '  "output": {"file": "out"}\n}\n')     # as with # comments
        r = run(path)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn('already in the new form', r.stdout)
        self.assertEqual(os.listdir(self.dir), ['input.jsonc'])

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


class Frames(unittest.TestCase):
    """A90 (2026-10-16): frame.rotation keeps qubit_axis and enabled; one key "frame" ("file" | "qubit") beside the
    numbers; with the rotation on, the hidden defaults written out ("qubit")."""
    A90 = os.path.join(FIX, 'a90')

    def conv(self, text):
        new, log = convert_text(text)
        return json.loads(strip_comments(new)), new, log

    def test_fixture(self):
        before = read(os.path.join(self.A90, 'rotation_frames.before.json'))
        current = read(os.path.join(self.A90, 'rotation_frames.current.json'))
        new, log = convert_text(before)
        self.assertEqual(new, current)
        self.assertEqual(sorted(comments(before)), sorted(comments(new)))
        self.assertEqual(convert_text(current), (current, []))
        d = json.loads(strip_comments(new))
        self.assertEqual(d['frame'], {'rotation': {'qubit_axis': [1, 1, 1]}})
        self.assertEqual(d['Qubit']['intmap'][0]['frame'], 'qubit')        # A93: Qubit.zfs and Qubit.frame -> intmap
        self.assertEqual(d['hamiltonian']['common']['hyperfine']['frame'], 'file')
        self.assertEqual(d['hamiltonian']['common']['quadrupole']['frame'], 'qubit')
        self.assertEqual([s['frame'] for s in d['Defect']['species']], ['qubit', 'file'])
        self.assertTrue([l for l in log if 'bath_axis removed (it was the default' in l])
        self.assertTrue([l for l in log if 'reference_qubit "q0" removed' in l])

    def test_bath_axis_not_default_noted(self):
        d, _, log = self.conv('{"frame": {"rotation": {"bath_axis": [1, 0, 0], "qubit_axis": [1, 1, 1]}}}')
        self.assertEqual(d, {'frame': {'rotation': {'qubit_axis': [1, 1, 1]}}})
        self.assertTrue([l for l in log if 'bath_axis [1, 0, 0] removed' in l and 'can change (A90)' in l])

    def test_intmap(self):
        one = ('{"Qubit": {"qubit": [{"name": "q0", "xyz": [0, 0, 0]}], "intmap": [{"between": ["q0", "q0"], '
               '"tensor": [1, 0, 0, 0, 1, 0, 0, 0, 1]}, {"between": ["q0", "q0"], "tensor_frame": "bath", '
               '"tensor": [1, 0, 0, 0, 1, 0, 0, 0, 1]}]}, "frame": {"rotation": {"qubit_axis": [1, 1, 1]}}}')
        d, _, _ = self.conv(one)
        self.assertEqual([e['frame'] for e in d['Qubit']['intmap']], ['qubit', 'file'])
        # two qubits: an entry without tensor_frame was refused before, so nothing is guessed for it
        two = ('{"Qubit": {"qubit": [{"name": "a", "xyz": [0, 0, 0]}, {"name": "b", "xyz": [0, 0, 5]}], '
               '"intmap": [{"between": ["a", "a"], "tensor": [1, 0, 0, 0, 1, 0, 0, 0, 1]}]}, '
               '"frame": {"rotation": {"qubit_axis": [1, 1, 1], "reference_qubit": "a", "qubit_position_frame": "bath"}}}')
        d, _, _ = self.conv(two)
        self.assertNotIn('frame', d['Qubit']['intmap'][0])
        self.assertEqual(d['frame'], {'rotation': {'qubit_axis': [1, 1, 1]}})

    def test_rotation_off_writes_no_frame(self):
        text = ('{"Qubit": {"file": "q", "zfs": [1, 0, 0, 0, 1, 0, 0, 0, 1]}, "Defect": {"species": [{"name": "P1", '
                '"coordinate_frame": "qubit"}, {"name": "X"}]}, "frame": {"rotation": {"enabled": false, '
                '"bath_axis": [0, 0, 1], "qubit_axis": [1, 1, 1]}}}')
        d, _, _ = self.conv(text)
        self.assertNotIn('frame', d['Qubit'])
        self.assertEqual(d['Defect']['species'], [{'name': 'P1', 'frame': 'qubit'}, {'name': 'X'}])
        self.assertEqual(d['frame'], {'rotation': {'enabled': False, 'qubit_axis': [1, 1, 1]}})

    def test_tensor_frame_without_its_block(self):
        d, _, log = self.conv('{"frame": {"rotation": {"qubit_axis": [1, 1, 1], "qd_tensor_frame": "bath"}}}')
        self.assertEqual(d, {'frame': {'rotation': {'qubit_axis': [1, 1, 1]}}})
        self.assertTrue([l for l in log if 'qd_tensor_frame removed' in l])

    def test_old_value_in_a_new_key(self):
        d, _, _ = self.conv('{"Defect": {"species": [{"name": "P1", "frame": "bath"}]}}')
        self.assertEqual(d['Defect']['species'][0]['frame'], 'file')


class QubitStatesAndZfs(unittest.TestCase):
    """A93 (2026-10-16): Qubit.qubit[].alphams/betams and Qubit.alphams/betams (Qubit.file) -> Qubit.alpha_state/beta_state
    {"ms": [...]}; Qubit.zfs and Qubit.frame (Qubit.file) -> a Qubit.intmap self entry."""
    A93 = os.path.join(FIX, 'a93')

    def conv(self, text):
        new, log = convert_text(text)
        return json.loads(strip_comments(new)), log

    def test_fixtures(self):
        for name in ('states', 'file_zfs'):
            with self.subTest(name):
                before = read(os.path.join(self.A93, name + '.before.json'))
                current = read(os.path.join(self.A93, name + '.current.json'))
                new, log = convert_text(before)
                self.assertEqual(new, current)
                self.assertTrue(log)
                self.assertEqual(sorted(comments(before)), sorted(comments(new)))
                self.assertEqual(convert_text(current), (current, []))

    def test_per_qubit_with_defaults(self):
        d = json.loads(strip_comments(read(os.path.join(self.A93, 'states.current.json'))))
        self.assertEqual(d['Qubit']['alpha_state'], {'ms': [-1, 1]})      # NV2 wrote no alphams: +1
        self.assertEqual(d['Qubit']['beta_state'], {'ms': [0, 0]})
        self.assertFalse([q for q in d['Qubit']['qubit'] if 'alphams' in q or 'betams' in q])

    def test_one_side_only(self):
        d, _ = self.conv('{"Qubit": {"qubit": [{"name": "a"}, {"name": "b", "betams": -1}]}}')
        self.assertEqual(d['Qubit'], {'qubit': [{'name': 'a'}, {'name': 'b'}], 'beta_state': {'ms': [0, -1]}})

    def test_written_state_wins(self):
        d, log = self.conv('{"Qubit": {"qubit": [{"name": "a", "alphams": 1, "betams": 0}], "alpha_state": {"index": 2}}}')
        self.assertEqual(d['Qubit'], {'qubit': [{'name': 'a'}], 'alpha_state': {'index': 2}, 'beta_state': {'ms': [0]}})
        self.assertTrue([l for l in log if 'alphams removed: Qubit.alpha_state is written' in l and 'per_qubit' in l])
        d, log = self.conv('{"Qubit": {"file": "q", "alphams": 1, "alpha_state": {"ms": [-1]}}}')
        self.assertEqual(d['Qubit'], {'qubit': [{'file': 'q'}], 'alpha_state': {'ms': [-1]}})
        self.assertTrue([l for l in log if 'Qubit.alphams removed' in l])

    def test_file_form(self):
        d, _ = self.conv('{"Qubit": {"file": "q", "alphams": -1, "betams": 0}}')
        self.assertEqual(d['Qubit'], {'qubit': [{'file': 'q'}], 'alpha_state': {'ms': [-1]}, 'beta_state': {'ms': [0]}})

    def test_before_a47(self):
        d, _ = self.conv('{"qubitfile": "q", "qzfs": [1, 0, 0, 0, 1, 0, 0, 0, -2], "qalphams": -1}')
        self.assertEqual(d['Qubit'], {'qubit': [{'file': 'q'}], 'intmap': [{'between': ['q0', 'q0'], 'tensor': [1, 0, 0, 0, 1, 0, 0, 0, -2],
                                                                 'unit': 'kHz'}], 'alpha_state': {'ms': [-1]}})

    def test_zfs(self):
        z = [[-960000, 0, 0], [0, -960000, 0], [0, 0, 1920000]]
        d, _ = self.conv('{"Qubit": {"file": "q", "name": "NV", "zfs": %s}}' % json.dumps(z))
        self.assertEqual(d['Qubit'], {'qubit': [{'file': 'q', 'name': 'NV'}], 'intmap': [{'between': ['NV', 'NV'], 'tensor': z, 'unit': 'kHz'}]})
        # the rotation on: A90 writes "qubit" (Qubit.zfs was never turned), which goes into the entry
        d, _ = self.conv('{"Qubit": {"file": "q", "zfs": %s}, "frame": {"rotation": {"qubit_axis": [1, 1, 1]}}}' % json.dumps(z))
        self.assertEqual(d['Qubit']['intmap'][0]['frame'], 'qubit')
        self.assertNotIn('frame', d['Qubit'])
        d, log = self.conv('{"Qubit": {"file": "q", "frame": "file"}}')
        self.assertEqual(d['Qubit'], {'qubit': [{'file': 'q'}]})
        self.assertTrue([l for l in log if 'Qubit.frame removed' in l])

    def test_zfs_beside_intmap_refused(self):
        with self.assertRaisesRegex(ValueError, 'both "file" and "intmap"'):
            convert_text('{"Qubit": {"file": "q", "zfs": [1, 0, 0, 0, 1, 0, 0, 0, -2], "intmap": []}}')


class InputKeyReview(unittest.TestCase):
    """A94 (2026-10-16, (internal review)): s_range, the seed of a k_range run, the Qubit.file form as an
    entry, pulse.bath.levels, print_trials, Defect zfs.axis and unit case, the hamiltonian sections, cluster.adaptive.workdir,
    the auto_convergence keys."""
    A94 = os.path.join(FIX, 'a94')

    def conv(self, text):
        new, log = convert_text(text)
        return json.loads(strip_comments(new)), log

    def test_fixture(self):
        before = read(os.path.join(self.A94, 'keys.before.jsonc'))
        current = read(os.path.join(self.A94, 'keys.current.jsonc'))
        new, log = convert_text(before)
        self.assertEqual(new, current)
        self.assertEqual(sorted(comments(before)), sorted(comments(new)))
        self.assertEqual(convert_text(current), (current, []))

    def test_nstate(self):
        d, _ = self.conv('{"configuration": {"bathfile": ["b"], "nstate": 5}}')
        self.assertEqual(d['configuration'], {'bathfile': ['b'], 's_range': [1, 5]})
        d, log = self.conv('{"configuration": {"bathfile": ["b"], "nstate": 0}}')
        self.assertEqual(d['configuration'], {'bathfile': ['b']})
        self.assertTrue([l for l in log if 'nstate 0 removed' in l])
        d, _ = self.conv('{"nstate": 2, "bathfile": ["b"]}')                 # before A64: top level
        self.assertEqual(d['configuration'], {'bathfile': ['b'], 's_range': [1, 2]})

    def test_seed_of_a_range(self):
        # draws (sampled states): the old meaning, the same seed for every k, as a list
        d, _ = self.conv('{"configuration": {"k_range": [2, 4], "bathfile": ["c{k}"], "nstate": 1, "seed": 5}}')
        self.assertEqual(d['configuration']['seed'], [5, 5, 5])
        # the ensemble of a bath without Defect species draws nothing: left as it is
        d, _ = self.conv('{"configuration": {"k_range": [2, 4], "bathfile": ["c{k}"], "seed": 5}}')
        self.assertEqual(d['configuration']['seed'], 5)
        # Defect species and pcce draw too
        d, _ = self.conv('{"configuration": {"k_range": [1, 2], "bathfile": ["c{k}"], "seed": 5}, "Defect": {"species": [{"name": "P1"}]}}')
        self.assertEqual(d['configuration']['seed'], [5, 5])
        d, _ = self.conv('{"cluster": {"mode": "pcce"}, "configuration": {"k_range": [1, 2], "bathfile": ["c{k}"], "seed": 5}}')
        self.assertEqual(d['configuration']['seed'], [5, 5])
        # one configuration, or the automatic convergence (its runner reads the seed itself): left as it is
        d, _ = self.conv('{"configuration": {"k_range": [3, 3], "bathfile": ["c{k}"], "nstate": 1, "seed": 5}}')
        self.assertEqual(d['configuration']['seed'], 5)
        d, _ = self.conv('{"auto_convergence": {"workdir": "w"}, "configuration": {"k_range": [1, 2], "bathfile": ["c{k}"], '
                         '"nstate": 1, "seed": 5}}')
        self.assertEqual(d['configuration']['seed'], 5)

    def test_qubit_file_form(self):
        d, _ = self.conv('{"Qubit": {"file": "q", "spin": 0.5, "alpha_state": {"ms": [0.5]}}}')
        self.assertEqual(d['Qubit'], {'qubit': [{'file': 'q', 'spin': 0.5}], 'alpha_state': {'ms': [0.5]}})
        d, _ = self.conv('{"qubitfile": "q", "qspin": 0.5}')                  # before A47
        self.assertEqual(d['Qubit'], {'qubit': [{'file': 'q', 'spin': 0.5}]})

    def test_detuning_list(self):
        d, _ = self.conv('{"Qubit": {"qubit": [{"xyz": [0, 0, 0], "detuning": [2.5]}]}}')
        self.assertEqual(d['Qubit']['qubit'][0]['detuning'], 2.5)

    def test_bath_pulse(self):
        d, _ = self.conv('{"pulse": {"count": 1, "bath": {"count": 0, "spin": 1, "alphams": 1, "betams": 0}}}')
        self.assertEqual(d['pulse']['bath'], {'count': 0, 'levels': [1, 0]})
        d, _ = self.conv('{"pulse": {"count": 1, "bath": {"count": 0, "betams": 0.5}}}')
        self.assertEqual(d['pulse']['bath'], {'count': 0, 'levels': [0.5, 0.5]})

    def test_print_trials(self):
        pcce = '{"cluster": {"mode": "pcce", "pcce": {"partition_size": 2, %s}}}'
        for written, now in (('"print_trials": true', {'print_trials': 32}),
                             ('"print_trials": true, "print_limit": 5', {'print_trials': 5}),
                             ('"print_trials": false, "print_limit": 5', {}), ('"print_limit": 5', {}),
                             ('"print_trials": 7', {'print_trials': 7})):
            with self.subTest(written):
                d, _ = self.conv(pcce % written)
                self.assertEqual(d['cluster']['pcce'], dict(partition_size=2, **now))

    def test_defect(self):
        d, log = self.conv('{"Defect": {"species": [{"name": "NV", "axis": [0, 0, 1], "detuning": [[1, "e", 0]]}]}}')
        self.assertEqual(d['Defect']['species'][0], {'name': 'NV', 'detuning': [[1, 'e', 0]]})
        self.assertTrue([l for l in log if 'only the D form' in l])
        with self.assertRaisesRegex(ValueError, 'zfs.E is 5'):
            convert_text('{"Defect": {"species": [{"name": "NV", "zfs": {"D": 1, "E": 5}}]}}')
        d, _ = self.conv('{"Defect": {"species": [{"name": "P1", "gyros": {"values": [1], "unit": "mhz/t"}}]}}')
        self.assertEqual(d['Defect']['species'][0]['gyros']['unit'], 'MHz/T')

    def test_hamiltonian(self):
        # the sections by method stay (A94's G1 reverted, 2026-10-17)
        text = '{"hamiltonian": {"cCCE": {"knight": true}, "meCCE": {"jump_operators": []}, "gCCE": {}, "noise": {}}}'
        self.assertEqual(convert_text(text), (text, []))

    def test_hamiltonian_sections_of_a94(self):
        # corrections and dissipation, written by the A94 converter, back to the sections by method
        d, log = self.conv('{"method": "meCCE", "hamiltonian": {"corrections": {"knight": true}, '
                           '"dissipation": {"jump_operators": []}}}')
        self.assertEqual(d['hamiltonian'], {'cCCE': {'knight': True}, 'meCCE': {'jump_operators': []}})
        self.assertIn('hamiltonian.corrections -> hamiltonian.cCCE', log)
        d, _ = self.conv('{"method": "megCCE", "hamiltonian": {"dissipation": {"jump_operators": []}}}')
        self.assertEqual(d['hamiltonian'], {'megCCE': {'jump_operators': []}})
        d, _ = self.conv('{"hamiltonian": {"corrections": {"hfmedi": true}}}')          # no method: cCCE
        self.assertEqual(d['hamiltonian'], {'cCCE': {'hfmedi': True}})
        with self.assertRaisesRegex(ValueError, 'both written'):
            convert_text('{"hamiltonian": {"cCCE": {}, "corrections": {"knight": true}}}')
        # the comments and the column of the values stay
        new, _ = convert_text('{\n    "method"      : "meCCE",\n    "hamiltonian" : {\n'
                              '        "dissipation" : {"jump_operators": []}   // none yet\n    }\n}')
        self.assertEqual(new, '{\n    "method"      : "meCCE",\n    "hamiltonian" : {\n'
                              '        "meCCE"       : {"jump_operators": []}   // none yet\n    }\n}')

    def test_adaptive(self):
        d, _ = self.conv('{"cluster": {"mode": "bundle", "adaptive": {"output": "a", "max-rounds": 3, "policy": "merge"}}}')
        self.assertEqual(d['cluster']['adaptive'], {'workdir': 'a', 'max_rounds': 3})

    def test_auto_convergence(self):
        d, log = self.conv('{"auto_convergence": {"workdir": "w", '
                           '"time_window": {"converge": false, "value": 2, "unit": "us", "threshold": 0.1}, '
                           '"rdip": {"converge": true, "candidates": [4, 6]}, '
                           '"sampling": {"converge": false, "number": 10, "states": 3, "test_number": 4}}, '
                           '"output": {"path": "out", "file": "x"}}')
        self.assertEqual(d['auto_convergence'], {'workdir': 'w',
                                                 'time_window': {'value': {'value': 2, 'unit': 'us'}, 'decay_level': 0.1},
                                                 'rdip': {'candidates': [4, 6]},
                                                 'sampling': {'configurations': 10, 'states': 3, 'test_configurations': 4}})
        self.assertNotIn('output', d)
        with self.assertRaisesRegex(ValueError, 'converge true with a fixed value'):
            convert_text('{"auto_convergence": {"rbath": {"converge": true, "value": 10}}}')
        # auto_convergence off: the ordinary run reads output
        d, _ = self.conv('{"auto_convergence": {"enabled": false, "workdir": "w"}, "output": {"path": "out"}}')
        self.assertEqual(d['output'], {'path': 'out'})


if __name__ == '__main__':
    unittest.main()
