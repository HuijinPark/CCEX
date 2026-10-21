#!/usr/bin/env python3
"""Rewrite CCEX inputs written with older input keys (v1.2.2 and before, and the development versions before 2.0.0) in
the keys of 2.0.0: the names and values of the keys; the configuration, execution, cluster, output and frame blocks; one
state or axis file per bath file; frame.rotation with qubit_axis and enabled only, and one key "frame" ("file" |
"qubit") beside the numbers it describes, written out where the rotation is on and CCEX used to take a hidden default;
the qubit states in Qubit.alpha_state/beta_state only ({"ms": [...]} for the old alphams/betams), a qubit from a file as
a Qubit.qubit entry and its ZFS as a Qubit.intmap self entry; configuration.s_range for nstate and the seed of a k_range
run; pulse.bath.levels, cluster.pcce.print_trials as a number, Defect zfs.axis, cluster.adaptive.workdir, the
auto_convergence keys without converge; the hamiltonian sections by method (cCCE, meCCE/megCCE).

    python3 convert_input_keys.py FILE...            rewrite each file in place and list what changed
    python3 convert_input_keys.py --check FILE...    list only; exit 1 when a file still has old keys
    python3 convert_input_keys.py --stdout FILE      print the new text, leave the file alone

The table of old and new names is docs/user/InputFile.md; (internal notes) has the decisions.
Comments (#, ! or // to the end of a line, outside strings; .json and .jsonc files alike) stay: a key that moves into
a block takes the comment lines above it and the comment after it on its line. A value is copied as written unless its
form changes (method, a 0/1 switch, hf_readmode, ...). An object nothing changed in is left byte for byte. A file in
the new form is left alone.
A key written twice (an old key beside its new name, or two spellings that differ in case only) is kept once when the
values are the same; different values stop the file before anything is written. Tests: tests/ (python3 -m unittest).
"""
import json
import os
import shutil
import sys
import tempfile

WS = ' \t\r\n'


def comment_at(text, i):
    """A comment starts at i (the caller knows i is outside a string): #, ! or //, to the end of the line."""
    return i < len(text) and (text[i] in '#!' or text.startswith('//', i))


class Node:
    def __init__(self, text, start, end):
        self.text, self.start, self.end, self.dirty = text, start, end, False

    def raw(self):
        return self.text[self.start:self.end]


class Scalar(Node):
    def __init__(self, raw_text):
        super().__init__(raw_text, 0, len(raw_text))

    def render(self):
        return self.text


class Member:
    """pre: blank lines and comments before the key; tail: what follows the value on its line (a comment)."""
    def __init__(self, pre, keyraw, mid, value, tail='', shift=0):
        self.pre, self.keyraw, self.mid, self.value, self.tail, self.shift = pre, keyraw, mid, value, tail, shift

    @property
    def key(self):
        return json.loads(self.keyraw)

    def indent(self):
        return self.pre.rsplit('\n', 1)[-1] if '\n' in self.pre else None


class Container(Node):
    open_ch = close_ch = ''

    def __init__(self, text, start, end, items, close):
        super().__init__(text, start, end)
        self.items, self.close = items, close

    def render(self):
        if not self.dirty and not any(changed(m.value) for m in self.items):
            return self.raw()
        out = [self.open_ch]
        for i, m in enumerate(self.items):
            body = render_value(m.value, m.shift)
            last = i == len(self.items) - 1
            out.append(m.pre + (m.keyraw + m.mid if m.keyraw is not None else '') + body + ('' if last else ',') + m.tail)
        out.append(self.close + self.close_ch)
        return ''.join(out)


class Obj(Container):
    open_ch, close_ch = '{', '}'

    def find(self, key):
        for i, m in enumerate(self.items):
            if m.key.lower() == key.lower():
                return i, m
        return None, None

    def get(self, key):
        return self.find(key)[1]

    def take(self, key):
        """Remove the member. The comment lines right above it go with it; what is above them (blank lines, comments
        that belong to the lines before) stays where it was, in front of the member that comes next."""
        i, m = self.find(key)
        if m is not None:
            del self.items[i]
            self.dirty = True
            lines = m.pre.split('\n')
            k = len(lines) - 1
            while k > 0 and lines[k - 1].strip():
                k -= 1
            keep, m.pre = '\n'.join(lines[:k]), '\n' + '\n'.join(lines[k:])
            if keep.strip():
                if i < len(self.items):
                    nxt = self.items[i]
                    nxt.pre = keep + nxt.pre
                    nxt.left = keep + getattr(nxt, 'left', '')
                else:
                    self.close = keep + self.close
            elif keep and i < len(self.items):
                self.items[i].blank = keep          # blank lines that a member put in its place takes over
        return i, m

    def drop(self, key):
        """Remove the member for good; its comment lines stay where they were."""
        i, m = self.find(key)
        if m is None:
            return None
        return self.drop_at(i)

    def drop_at(self, i):
        m = self.items[i]
        del self.items[i]
        self.dirty = True
        body = m.pre.rsplit('\n', 1)[0] if '\n' in m.pre else ''
        if body.strip() or (body and i < len(self.items) and not self.items[i].pre.startswith('\n\n')):
            if i < len(self.items):
                self.items[i].pre = body + self.items[i].pre
            else:
                self.close = body + self.close
        return m

    def insert(self, i, member):
        """A new member at i: the text another removal left in front of the member there goes in front of it."""
        if i < len(self.items):
            nxt = self.items[i]
            left = getattr(nxt, 'left', '')
            if left and nxt.pre.startswith(left):
                nxt.pre, nxt.left = nxt.pre[len(left):], ''
                member.pre = left + member.pre
            elif getattr(nxt, 'blank', '') and '\n' in member.pre:
                member.pre, nxt.blank = nxt.blank + member.pre, ''
        self.items.insert(i, member)
        self.dirty = True


class Arr(Container):
    open_ch, close_ch = '[', ']'


def changed(node):
    return isinstance(node, Container) and (node.dirty or any(changed(m.value) for m in node.items))


def render_value(node, shift=0):
    text = node.render()
    if shift and '\n' in text:
        lines = text.split('\n')
        if shift > 0:
            lines = [lines[0]] + [(' ' * shift + l) if l.strip() else l for l in lines[1:]]
        else:
            lines = [lines[0]] + [l[min(-shift, len(l) - len(l.lstrip(' '))):] for l in lines[1:]]
        text = '\n'.join(lines)
    return text


class Parser:
    def __init__(self, text):
        self.t = text

    def ws(self, i):
        t = self.t
        while i < len(t):
            if t[i] in WS:
                i += 1
            elif comment_at(t, i):
                while i < len(t) and t[i] != '\n':
                    i += 1
            else:
                break
        return i

    def string_end(self, i):
        t = self.t
        assert t[i] == '"', f'expected a string at {i}'
        i += 1
        while t[i] != '"':
            i += 2 if t[i] == '\\' else 1
        return i + 1

    def tail(self, i, close):
        """After a value: (tail text, index where the next member's text starts)."""
        t = self.t
        j = self.ws(i)
        if t[j] == ',':
            k = j + 1                       # a comment after the comma on its line belongs to this member;
            while k < len(t) and t[k] in ' \t':   # the spaces before the next key on the same line, to that key
                k += 1
            if comment_at(t, k):
                while k < len(t) and t[k] != '\n':
                    k += 1
            else:
                k = j + 1
            return t[i:j] + t[j + 1:k], k, False
        assert t[j] == close, f'expected , or {close} at {j}'
        nl = t.find('\n', i, j)
        if nl < 0:
            return t[i:j], j, True
        return t[i:nl], nl, True

    def value(self, i):
        t = self.t
        c = t[i]
        if c == '{' or c == '[':
            close = '}' if c == '{' else ']'
            items, k = [], i + 1
            while True:
                j = self.ws(k)
                if t[j] == close:
                    node = (Obj if c == '{' else Arr)(t, i, j + 1, items, t[k:j])
                    return node, j + 1
                if c == '{':
                    ke = self.string_end(j)
                    colon = self.ws(ke)
                    assert t[colon] == ':', f'expected : at {colon}'
                    v0 = self.ws(colon + 1)
                    val, e = self.value(v0)
                    tail, nxt, last = self.tail(e, close)
                    items.append(Member(t[k:j], t[j:ke], t[ke:v0], val, tail))
                else:
                    val, e = self.value(j)
                    tail, nxt, last = self.tail(e, close)
                    items.append(Member(t[k:j], None, '', val, tail))
                k = nxt
                if last:
                    j = self.ws(k)
                    node = (Obj if c == '{' else Arr)(t, i, j + 1, items, t[k:j])
                    return node, j + 1
        if c == '"':
            e = self.string_end(i)
        else:
            e = i
            while e < len(t) and t[e] not in WS + ',}]' and not comment_at(t, e):
                e += 1
        node = Node(t, i, e)
        node.render = node.raw
        return node, e


def plain(node):
    """The Python value of a node (comments stripped)."""
    return json.loads(strip_comments(node.render()))


def strip_comments(text):
    out, in_str, esc, i = [], False, False, 0
    while i < len(text):
        c = text[i]
        if in_str:
            out.append(c)
            if esc:
                esc = False
            elif c == '\\':
                esc = True
            elif c == '"':
                in_str = False
        elif comment_at(text, i):
            while i < len(text) and text[i] != '\n':
                i += 1
            continue
        else:
            if c == '"':
                in_str = True
            out.append(c)
        i += 1
    return ''.join(out)


class Converter:
    def __init__(self, text):
        self.text = text
        p = Parser(text)
        s = p.ws(0)
        self.root, e = p.value(s)
        if not isinstance(self.root, Obj):
            raise ValueError('the input is not a JSON object')
        self.prefix, self.suffix = text[:s], text[e:]
        self.log = []
        self.unit = self.detect_unit()
        first = self.root.items[0] if self.root.items else None
        self.mid = first.mid if first is not None else ': '
        self.width = len(first.keyraw) if first is not None else 0      # keys padded to a column, as the file pads them

    def detect_unit(self):
        top = self.member_indent(self.root) or ''
        steps = []                              # how much deeper the entries of a block or a list are written
        for m in self.root.items:
            if isinstance(m.value, Container) and m.value.items:
                inner = m.value.items[0].indent()
                if inner is not None and len(inner) > len(top):
                    steps.append(len(inner) - len(top))
        return ' ' * min(steps) if steps and min(steps) <= 8 else '    '

    @staticmethod
    def member_indent(obj):
        for m in obj.items:
            ind = m.indent()
            if ind is not None:
                return ind
        return None

    def indent_of(self, obj, parent_indent):
        ind = self.member_indent(obj)
        return ind if ind is not None else parent_indent + self.unit

    def note(self, text):
        self.log.append(text)

    # ---- building ---------------------------------------------------------------------------------------------
    @staticmethod
    def aligned(mid, width, keyraw):
        """mid for keyraw, keeping the colon in the column it was in (key + spaces = width)."""
        lead = len(mid) - len(mid.lstrip(' '))
        if lead == 0:
            return mid
        return ' ' * max(1, width + lead - len(keyraw)) + mid.lstrip(' ')

    def new_member(self, key, value, indent, moved=None, from_indent=None):
        """A member for key at indent: a new value (text or node), or the value of a moved member with its comments."""
        if moved is not None:
            comments = [l.strip() for l in moved.pre.split('\n') if l.strip()]
            pre = ''.join('\n' + indent + c for c in comments) + '\n' + indent
            # the moved value keeps its own lines, shifted to the new depth; a new value is built at its depth
            shift = len(indent) - len(from_indent) if (from_indent is not None and value is None) else 0
            node = moved.value if value is None else value
            if isinstance(node, str):
                node = Scalar(node)
            keyraw = json.dumps(key)
            return Member(pre, keyraw, self.aligned(moved.mid, len(moved.keyraw), keyraw), node, moved.tail, shift)
        node = Scalar(value) if isinstance(value, str) else value
        keyraw = json.dumps(key)
        return Member('\n' + indent, keyraw, self.aligned(self.mid, self.width, keyraw), node, '')

    def rename(self, member, key):
        keyraw = json.dumps(key)
        member.mid = self.aligned(member.mid, len(member.keyraw), keyraw)
        member.keyraw = keyraw

    def new_block(self, indent):
        """An empty object whose members will sit at indent + unit and whose } sits at indent."""
        o = Obj('', 0, 0, [], '\n' + indent)
        o.dirty = o.new = True
        return o

    def child_block(self, parent, parent_indent, key, position=None):
        """parent[key] as an object, made (at position, or at the end) when absent."""
        m = parent.get(key)
        if m is not None:
            if not isinstance(m.value, Obj):
                raise ValueError(f'"{key}" is not an object')
            return m.value
        o = self.new_block(parent_indent)
        mem = self.new_member(key, o, parent_indent)
        if position is None:
            parent.items.append(mem)
        else:
            parent.insert(position, mem)
        parent.dirty = True
        return o

    def put(self, block, block_indent, key, moved=None, value=None, from_indent=None):
        ind = self.indent_of(block, block_indent)
        block.items.append(self.new_member(key, value, ind, moved, from_indent))
        block.dirty = True

    # ---- the conversion -----------------------------------------------------------------------------------------
    def convert(self):
        R = self.root
        top = self.member_indent(R) or ''
        self.method_shape(R)
        self.radii(R)
        self.renames(R)
        self.savemode(R, top)
        self.single_keys(R, '')         # a renamed key beside its new name, before a block takes one of the two
        self.configuration(R, top)
        self.blocks(R, top)
        self.qubit(R, top)
        self.hamiltonian(R, top)
        self.defect(R, top)
        self.defect_files(R, top)
        self.a65(R, top)
        self.a90(R, top)
        self.a93(R, top)
        self.a94(R, top)
        self.a94_sections(R)
        self.single_keys(R, '')
        self.compact(R)
        return self.prefix + R.render() + self.suffix

    def a94_sections(self, R):
        """2026-10-17 (A94 G1 reverted): the hamiltonian sections stay named by method. hamiltonian.corrections -> cCCE (meCCE
        reads it too), hamiltonian.dissipation -> megCCE in a megCCE run, meCCE otherwise."""
        H = self.obj(R.get('hamiltonian'))
        if H is None:
            return
        _, m = R.find('method')
        method = str(plain(m.value)).lower() if m is not None else ''
        me = 'megCCE' if method == 'megcce' else 'meCCE'
        for old, new in (('corrections', 'cCCE'), ('dissipation', me)):
            o = H.get(old)
            if o is None:
                continue
            if H.get(new) is not None:
                raise ValueError(f'hamiltonian.{old} and hamiltonian.{new} are both written')
            self.rename(o, new)
            H.dirty = True
            self.note(f'hamiltonian.{old} -> hamiltonian.{new}')

    def single_keys(self, node, path):
        """One member per key, case ignored as CCEX ignores it (review 20261013-180049 C1). An old key whose new name is
        already written (npulse beside pulse.count) gave the key twice, and CCEX reads the first, Python json the last:
        the same value is kept once, different values refuse the file before anything is written."""
        if isinstance(node, Obj):
            i = 0
            while i < len(node.items):
                m = node.items[i]
                first = next(x for x in node.items if x.key.lower() == m.key.lower())
                if first is m:
                    i += 1
                    continue
                a, b = plain(first.value), plain(m.value)
                if json.dumps(a, sort_keys=True) != json.dumps(b, sort_keys=True):     # true is not 1 here
                    raise ValueError(f'{path}{first.key} is written twice with different values ({json.dumps(a)} and '
                                     f'{json.dumps(b)}): an old key and its new name were both in the input, or a key in '
                                     f'two cases; keep one')
                node.drop_at(i)                 # its comment lines stay; the comment after it goes above the one kept
                if m.tail.strip():
                    ind = first.indent()
                    if ind is not None:
                        first.pre = first.pre[:len(first.pre) - len(ind)] + ind + m.tail.strip() + '\n' + ind
                    else:
                        self.note(f'{path}{m.key}: its comment {m.tail.strip()!r} was not kept (the block is on one line)')
                self.note(f'{path}{m.key}: the same value was written twice, kept once')
        for n, m in enumerate(node.items):
            if isinstance(m.value, Container):
                self.single_keys(m.value, f'{path}{m.key}.' if m.keyraw is not None else f'{path[:-1]}[{n}].')

    def compact(self, obj):
        """A new block of a few plain values on one line: "pulse": {"count": 1}. The comment lines that came with its
        values go above it; a block holding a block, or a value with a comment after it, stays as it is."""
        for m in obj.items:
            v = m.value
            if not isinstance(v, Container):
                continue
            self.compact(v) if isinstance(v, Obj) else [self.compact(x.value) for x in v.items if isinstance(x.value, Obj)]
            if not getattr(v, 'new', False) or not v.items:
                continue
            if any(isinstance(x.value, Container) or getattr(x.value, 'inline', False) or '\n' in x.value.render()
                   or x.tail.strip() for x in v.items):
                continue
            inline = '{' + ', '.join(f'{x.keyraw}: {x.value.render()}' for x in v.items) + '}'
            if len(inline) > 72:
                continue
            ind = m.indent() or ''
            # line breaks that followed a moved value (an input with the comma at the start of the next line) stay
            # after the block written on one line, where the value stood
            breaks = ''.join(x.tail for x in v.items if '\n' in x.tail and not x.tail.strip())
            if breaks and not m.tail:
                m.tail = breaks
            comments = [l.strip() for x in v.items for l in x.pre.split('\n') if l.strip()]
            if comments:
                m.pre = m.pre[:len(m.pre) - len(ind)] + ''.join(ind + c + '\n' for c in comments) + ind
            m.value = Scalar(inline)
            m.value.inline = True
            obj.dirty = True

    def switch(self, member, name):
        """0/1 -> false/true."""
        v = plain(member.value)
        if isinstance(v, bool):
            return member.value
        if v in (0, 1):
            self.note(f'{name}: {v} -> {"true" if v else "false"}')
            return Scalar('true' if v else 'false')
        raise ValueError(f'{name} must be true or false (or 0 or 1), not {v!r}')

    def method_shape(self, R):
        _, m = R.find('method')
        shape = None
        if m is not None:
            v = plain(m.value)
            low = str(v).lower()
            if low == 'cce':
                m.value = Scalar('"cCCE"')
                R.dirty = True
                self.note(f'method "{v}" -> "cCCE"')
            elif low in ('pcce', 'dsj', 'itb', 'dsjitb'):
                m.value = Scalar('"cCCE"')
                shape = low
                R.dirty = True
                self.note(f'method "{v}" -> "cCCE" with cluster_mode "{low}"')
        i, cm = R.find('cluster_mode')
        if cm is not None and str(plain(cm.value)).lower() == 'adaptive':
            cm.value = Scalar('"bundle"')
            R.dirty = True
            self.note('cluster_mode "adaptive" -> "bundle"')
        if shape is not None:
            if cm is not None:
                cm.value = Scalar(json.dumps(shape))
            else:
                i, _ = R.find('method')
                R.insert(i + 1, self.new_member('cluster_mode', json.dumps(shape), R.items[i].indent() or ''))
        self.shape = shape or (str(plain(cm.value)).lower() if cm is not None else 'standard')
        self.gcce = m is not None and str(plain(m.value)).lower() == 'gcce'

    def radii(self, R):
        for outer, inner in (('rbath', 'rbathcut'), ('rdip', 'rdipcut')):
            cut = R.get(inner)
            if cut is None:
                continue
            o = R.get(outer)
            v = plain(cut.value)
            if v == 0:
                R.drop(inner)
                self.note(f'{inner}: 0 removed (the default)')
                continue
            R.take(inner)
            if o is None:
                raise ValueError(f'{inner} without {outer}')
            o.value = Scalar(f'[{cut.value.render()}, {o.value.render()}]')
            if cut.tail.strip() and not o.tail.strip():
                o.tail = cut.tail
            self.note(f'{inner} -> {outer}: [inner, outer]')

    def renames(self, R):
        for old, new in (('bathadjust', 'bath_shift'), ('statefile', 'bath_state_file'),
                         ('nk', 'max_clusters_per_order'), ('addsubclus', 'add_subclusters')):
            _, m = R.find(old)
            if m is None:
                continue
            self.rename(m, new)
            if new == 'add_subclusters':
                m.value = self.switch(m, old)
            R.dirty = True
            self.note(f'{old} -> {new}')
        _, a = R.find('Adaptive')
        if a is not None and a.key != 'adaptive':
            self.rename(a, 'adaptive')
            R.dirty = True
            self.note('Adaptive -> adaptive')
        for name in ('DefectTotSpin', 'CorrTotSpin'):
            m = R.drop(name)
            if m is None:
                continue
            v = plain(m.value)
            if name == 'CorrTotSpin' and v != 0:
                raise ValueError('CorrTotSpin is not 0: it served old-format VASP tensor files, which have to be remade '
                                 'with example/tools/SpinTensorReader (format v3); the input cannot be converted as it is')
            self.note(f'{name} removed ({v}; ' + ('the qubit spin is used)' if name == 'DefectTotSpin' else 'it was 0)'))

    SAVEMODES = {'normal': 'coherence', 'info': 'coherence_tagged', 'all': 'cluster_coherence',
                 'allfull': 'cluster_density_matrix'}

    def savemode(self, R, top):
        """The savemode names before 2026-10-13 (normal, info, all, allfull) and the top-level save_nodiv, which CCEX
        refuses with their new form: "savemode": "<mode>", or {"mode": ..., "nodiv": ...} when save_nodiv was given."""
        i, m = R.find('savemode')
        if m is not None and not isinstance(m.value, Obj):
            v = plain(m.value)
            if isinstance(v, str) and v.lower() in self.SAVEMODES:
                m.value = Scalar('"' + self.SAVEMODES[v.lower()] + '"')
                R.dirty = True
                self.note(f'savemode "{v}" -> "{self.SAVEMODES[v.lower()]}"')
        j, n = R.find('save_nodiv')
        if n is None:
            return
        if m is not None and isinstance(m.value, Obj):
            raise ValueError('save_nodiv beside a savemode block: write nodiv inside the block')
        if not isinstance(plain(n.value), bool):     # 0/1 too: CCEX took only a JSON boolean, before and after A47
            raise ValueError('save_nodiv must be true or false')
        mode = m.value.render() if m is not None else '"coherence"'
        pos = i if i is not None else j
        if m is not None:
            R.take('savemode')
        R.take('save_nodiv')
        block = self.new_block(top)
        self.put(block, top, 'mode', value=mode)
        self.put(block, top, 'nodiv', value=n.value.render())
        mem = self.new_member('savemode', block, top, moved=m) if m is not None else self.new_member('savemode', block, top)
        mem.value, mem.shift = block, 0
        R.insert(min(pos, len(R.items)), mem)
        R.dirty = True
        self.note('save_nodiv -> savemode {mode, nodiv}')

    def configuration(self, R, top):
        """A64 (2026-10-15): the bath configuration in one block, "configuration": {"gyrofile", "bathfile", "bath_shift",
        "bath_state_file", "nstate", "seed"} (path and k_range are new keys; nothing to move into them)."""
        self.move(R, top, 'configuration', [('gyrofile', 'gyrofile'), ('bathfile', 'bathfile'),
                                            ('bath_shift', 'bath_shift'), ('bath_state_file', 'bath_state_file'),
                                            ('nstate', 'nstate'), ('seed', 'seed')])
        c = R.get('configuration')
        if c is not None and isinstance(c.value, Obj):
            self.per_bath_file(c.value, 'configuration.', 'bath_state_file')

    def defect_files(self, R, top):
        """A64: Defect.axis_file and Defect.state_file are the configuration's defect_axis_file and defect_state_file (the
        Defect block keeps the species alone), one name per bath file."""
        d = R.get('Defect')
        if d is None or not isinstance(d.value, Obj):
            return
        D = d.value
        found = [(old, new) for old, new in (('axis_file', 'defect_axis_file'), ('state_file', 'defect_state_file'))
                 if D.get(old) is not None]
        if not found:
            return
        C = self.child_block(R, top, 'configuration')
        cind = self.indent_of(C, top)
        for old, new in found:
            _, m = D.take(old)
            dind = self.indent_of(D, top)
            if C.get(new) is not None:
                raise ValueError(f'Defect.{old} beside configuration.{new}: keep one')
            self.put(C, top, new, m, None, dind)
            self.note(f'Defect.{old} -> configuration.{new}')
            self.per_bath_file(C, 'configuration.', new)

    def a65(self, R, top):
        """A65 (2026-10-15): execution {backend, gpu_options, propagator, evolution, evolution_update, gcce_refinement,
        column_evolution, profile}; cluster {mode, selection {add_subclusters, max_per_order}, pcce, bundle, adaptive
        (with probe_memory_gb)}; output {file, savemode, cluster_contributions, cluster_only}; frame {rotation}, on when
        given (an old block with "enabled" false or absent was off: removed); defect_axis_reference removed; the Defect
        species' n_onsite_spins and n_axis_states; the auto_convergence block without launcher, its configurations
        (without the file patterns, which go to the configuration and Defect blocks) as sampling."""
        self.move(R, top, 'execution', [(k, k) for k in ('backend', 'gpu_options', 'propagator', 'evolution',
                                                          'evolution_update', 'gcce_refinement', 'column_evolution',
                                                          'profile')])
        probe = R.get('probe_memory_gb')
        if probe is not None and R.get('adaptive') is None:
            R.drop('probe_memory_gb')
            self.note('probe_memory_gb removed (it is read with an adaptive block only)')
        self.move(R, top, 'cluster', [('cluster_mode', 'mode'), ('add_subclusters', ('selection', 'add_subclusters')),
                                      ('max_clusters_per_order', ('selection', 'max_per_order')), ('pcce', 'pcce'),
                                      ('bundle', 'bundle'), ('adaptive', 'adaptive')])
        if R.get('probe_memory_gb') is not None:
            c = R.get('cluster').value
            a = c.get('adaptive').value
            _, m = R.take('probe_memory_gb')
            self.put(a, self.indent_of(c, top), 'probe_memory_gb', m, None, top)
            self.note('probe_memory_gb -> cluster.adaptive.probe_memory_gb')
        self.move(R, top, 'output', [('outfile', 'file'), ('savemode', 'savemode'),
                                     ('cluster_contributions', 'cluster_contributions'), ('cluster_only', 'cluster_only')])
        for name in ('coordinate_frame_rotation', 'bath_coordinate_rotation'):
            m = R.get(name)
            if m is None:
                continue
            v = plain(m.value)
            if not isinstance(v, dict) or v.get('enabled', False) is not True:
                R.drop(name)
                self.note(f'{name} removed (the rotation was off; an absent frame.rotation is no rotation)')
                continue
            self.move(R, top, 'frame', [(name, 'rotation')])
        if R.drop('defect_axis_reference') is not None:
            self.note('defect_axis_reference removed (it only named the rotation\'s qubit axis in the log)')
        d = R.get('Defect')
        species = d.value.get('species') if d is not None and isinstance(d.value, Obj) else None
        if species is not None and isinstance(species.value, Arr):
            for n, item in enumerate(species.value.items):
                e = item.value
                if not isinstance(e, Obj):
                    continue
                for old, new in (('n_nuclear_spins', 'n_onsite_spins'), ('n_configurations', 'n_axis_states')):
                    m = e.get(old)
                    if m is not None:
                        self.rename(m, new)
                        e.dirty = True
                        self.note(f'Defect.species[{n}].{old} -> {new}')
        self.auto_convergence(R, top)

    # ---- A90 (2026-10-16): one key "frame" ("file" | "qubit") beside each set of directional numbers -------------------
    def set_frame(self, obj, value, ind):
        """obj.frame = value: written over an old value, or added at the end (on the line of an inline object)."""
        m = obj.get('frame')
        if m is not None:
            m.value = Scalar(value)
        elif self.member_indent(obj) is None and obj.items:            # {"a": 1, "b": 2}: stays on its line
            last = obj.items[-1]
            mid = last.mid.lstrip(' ') if last.mid.strip() == ':' else ': '
            obj.items.append(Member(' ', '"frame"', mid, Scalar(value), ''))
        else:
            self.put(obj, ind, 'frame', value=value)
        obj.dirty = True

    def rename_frame(self, obj, old, path):
        """obj[old] -> obj.frame, "bath" -> "file"; True when obj has a frame afterwards."""
        m = obj.get(old) if old != 'frame' else None
        if m is not None:
            self.rename(m, 'frame')
            obj.dirty = True
            self.note(f'{path}.{old} -> {path}.frame')
        m = obj.get('frame')
        if m is None:
            return False
        v = plain(m.value)
        if isinstance(v, str) and v.lower() == 'bath':
            m.value = Scalar('"file"')
            obj.dirty = True
            self.note(f'{path}.frame "bath" -> "file"')
        return True

    def a90(self, R, top):
        """A90 (2026-10-16): frame.rotation keeps qubit_axis and enabled. bath_axis, reference_qubit and
        qubit_position_frame are removed; hf_tensor_frame -> hamiltonian.common.hyperfine.frame, qd_tensor_frame ->
        hamiltonian.common.quadrupole.frame; Qubit.intmap[].tensor_frame and Defect.species[].coordinate_frame -> frame;
        the value "bath" -> "file". With the rotation on every place that held directional numbers has to state its frame;
        where CCEX took a hidden default the converter writes it: "qubit" for a Defect species without coordinate_frame, a
        single qubit's Qubit.intmap entry without tensor_frame and Qubit.zfs (all read as written, never turned)."""
        f = R.get('frame')
        rot = f.value.get('rotation') if f is not None and isinstance(f.value, Obj) else None
        rot = rot.value if rot is not None and isinstance(rot.value, Obj) else None
        en = rot.get('enabled') if rot is not None else None
        on = rot is not None and (en is None or plain(en.value) is not False)      # on when the block is given (A65)
        h = R.get('hamiltonian')
        common = h.value.get('common') if h is not None and isinstance(h.value, Obj) else None
        common = common.value if common is not None and isinstance(common.value, Obj) else None
        if rot is not None:
            first = rot.items[0] if rot.items else None
            first_pre = first.pre if first is not None else ''
            m = rot.get('bath_axis')
            if m is not None:
                v = plain(m.value)
                rot.drop('bath_axis')
                if isinstance(v, list) and len(v) == 3 and v[0] == 0 and v[1] == 0 and v[2] > 0:
                    self.note('frame.rotation.bath_axis removed (it was the default [0,0,1])')
                else:
                    self.note(f'frame.rotation.bath_axis {json.dumps(v)} removed: the computational x and y axes now turn '
                              'about qubit_axis to qubit_axis x [0,0,1] (the old default): cCCE and systems symmetric '
                              'about qubit_axis give the same results, a gCCE run with a qubit tensor that is not can '
                              'change (A90)')
            m = rot.get('reference_qubit')
            if m is not None:
                v = plain(m.value)
                rot.drop('reference_qubit')
                self.note(f'frame.rotation.reference_qubit {json.dumps(v)} removed (the rotation origin is the first qubit; '
                          'another origin is a translation, which changes no result)')
            if rot.drop('qubit_position_frame') is not None:
                self.note('frame.rotation.qubit_position_frame removed (its one value was "bath")')
            for old, block in (('hf_tensor_frame', 'hyperfine'), ('qd_tensor_frame', 'quadrupole')):
                m = rot.get(old)
                if m is None:
                    continue
                rot.drop(old)
                b = common.get(block) if common is not None else None
                if b is None or not isinstance(b.value, Obj):
                    self.note(f'frame.rotation.{old} removed (there is no hamiltonian.common.{block} block, so it was '
                              'not used)')
                    continue
                v = plain(m.value)
                value = '"file"' if isinstance(v, str) and v.lower() == 'bath' else m.value.render()
                self.set_frame(b.value, value, self.indent_of(common, top))
                self.note(f'frame.rotation.{old} -> hamiltonian.common.{block}.frame'
                          + (' ("bath" -> "file")' if value == '"file"' else ''))
            if rot.items and rot.items[0] is not first and '\n' not in first_pre:
                rot.items[0].pre = first_pre          # {"qubit_axis": ...} keeps the spacing the first key had
        if common is not None:
            for block in ('hyperfine', 'quadrupole'):
                b = common.get(block)
                if b is not None and isinstance(b.value, Obj):
                    self.rename_frame(b.value, 'frame', f'hamiltonian.common.{block}')     # "bath" -> "file"
        q = R.get('Qubit')
        Q = q.value if q is not None and isinstance(q.value, Obj) else None
        if Q is not None:
            qind = self.indent_of(Q, top)
            if Q.get('file') is not None:
                has = self.rename_frame(Q, 'frame', 'Qubit')
                if on and not has and Q.get('zfs') is not None:
                    self.set_frame(Q, '"qubit"', top)
                    self.note('Qubit.frame "qubit" written: Qubit.zfs was read as written (never turned) with the rotation on')
            qubits = Q.get('qubit')
            nqubit = len(qubits.value.items) if qubits is not None and isinstance(qubits.value, Arr) else None
            im = Q.get('intmap')
            if im is not None and isinstance(im.value, Arr):
                for n, item in enumerate(im.value.items):
                    e = item.value
                    if not isinstance(e, Obj):
                        continue
                    has = self.rename_frame(e, 'tensor_frame', f'Qubit.intmap[{n}]')
                    if on and not has and nqubit == 1:
                        self.set_frame(e, '"qubit"', qind)
                        self.note(f'Qubit.intmap[{n}].frame "qubit" written: a single qubit\'s entry without tensor_frame '
                                  'was read as written with the rotation on')
        d = R.get('Defect')
        species = d.value.get('species') if d is not None and isinstance(d.value, Obj) else None
        if species is not None and isinstance(species.value, Arr):
            sind = self.indent_of(d.value, top)
            for n, item in enumerate(species.value.items):
                e = item.value
                if not isinstance(e, Obj):
                    continue
                has = self.rename_frame(e, 'coordinate_frame', f'Defect.species[{n}]')
                if on and not has:
                    self.set_frame(e, '"qubit"', sind)
                    self.note(f'Defect.species[{n}].frame "qubit" written: a species without coordinate_frame was read '
                              'as written with the rotation on')

    # ---- A93 (2026-10-16): the qubit states in Qubit.alpha_state/beta_state only; the Qubit.file ZFS in Qubit.intmap ------
    @staticmethod
    def inline(node):
        """The value as one line: as written when it is on one line, else compact JSON."""
        text = node.render()
        return text if '\n' not in text else json.dumps(plain(node), separators=(', ', ': '))

    def state_member(self, Q, key, text, comments):
        """Q[key] = text after "qubit" (alpha_state before an existing beta_state), with the comments of what it replaces."""
        ind = self.member_indent(Q)
        first = Q.items[0]
        keyraw = json.dumps(key)
        if ind is None:
            mem = Member(' ', keyraw, first.mid.lstrip(' ') if first.mid.strip() == ':' else ': ', Scalar(text), '')
        else:
            pre = '\n' + ''.join(ind + c + '\n' for c in comments) + ind
            mem = Member(pre, keyraw, self.aligned(first.mid, len(first.keyraw), keyraw), Scalar(text), '')
        if key == 'alpha_state' and Q.get('beta_state') is not None:
            pos = Q.find('beta_state')[0]
        elif key == 'beta_state' and Q.get('alpha_state') is not None:
            pos = Q.find('alpha_state')[0] + 1
        else:
            pos = Q.find('qubit')[0] + 1
        Q.insert(pos, mem)

    def a93(self, R, top):
        """A93: Qubit.qubit[].alphams/betams -> Qubit.alpha_state/beta_state {"ms": [one per qubit]} (a qubit that wrote none
        takes the default, +1 or 0); Qubit.alphams/betams of the Qubit.file form -> {"ms": [m]}. Beside an alpha_state
        (beta_state) already written, which CCEX used, the old keys are removed. Qubit.zfs of the Qubit.file form (kHz)
        -> a Qubit.intmap self entry {"between": [name, name], "tensor": ..., "unit": "kHz"} with Qubit.frame as its
        "frame" (A90 wrote "qubit" there with the rotation on); a Qubit.frame without zfs is removed."""
        q = R.get('Qubit')
        Q = q.value if q is not None and isinstance(q.value, Obj) else None
        if Q is None or not Q.items:
            return
        if Q.get('file') is not None:
            for side in ('alpha', 'beta'):
                m = Q.get(side + 'ms')
                if m is None:
                    continue
                if Q.get(side + '_state') is not None:
                    Q.drop(side + 'ms')
                    self.note(f'Qubit.{side}ms removed (Qubit.{side}_state is written, and CCEX used that)')
                    continue
                self.rename(m, side + '_state')
                m.value = Scalar('{"ms": [' + self.inline(m.value) + ']}')
                Q.dirty = True
                self.note(f'Qubit.{side}ms -> Qubit.{side}_state {{"ms": [...]}}')
            z, f = Q.get('zfs'), Q.get('frame')
            if z is not None:
                if Q.get('intmap') is not None:
                    raise ValueError('the Qubit block has both "file" and "intmap" (CCEX refused it before too)')
                n = Q.get('name')
                name = json.dumps(plain(n.value)) if n is not None else '"q0"'
                entry = (f'[{{"between": [{name}, {name}], "tensor": {self.inline(z.value)}, "unit": "kHz"'
                         + (f', "frame": {f.value.render()}' if f is not None else '') + '}]')
                self.rename(z, 'intmap')
                z.value = Scalar(entry)
                if f is not None:
                    Q.drop('frame')
                Q.dirty = True
                self.note('Qubit.zfs' + (' and Qubit.frame' if f is not None else '') + ' -> Qubit.intmap self entry '
                          f'(between [{name}, {name}], unit kHz)')
            elif f is not None:
                Q.drop('frame')
                self.note('Qubit.frame removed (there is no Qubit.zfs, so it was not used)')
            return
        qubits = Q.get('qubit')
        if qubits is None or not isinstance(qubits.value, Arr):
            return
        entries = [x.value for x in qubits.value.items]
        for side, default in (('alpha', '1'), ('beta', '0')):
            key = side + 'ms'
            written = [e.get(key) if isinstance(e, Obj) else None for e in entries]
            if all(m is None for m in written):
                continue
            keep = Q.get(side + '_state') is not None
            values, comments = [], []
            for e, m in zip(entries, written):
                if m is None:
                    values.append(default)
                    continue
                values.append(self.inline(m.value))
                if m.tail.strip():
                    comments.append(m.tail.strip())
                e.drop(key)
            qubits.value.dirty = True
            if keep:
                self.note(f'Qubit.qubit[].{key} removed: Qubit.{side}_state is written and CCEX used it for the state '
                          f'(with pulse.per_qubit the level of each qubit is now the default m_s {default})')
                continue
            self.state_member(Q, side + '_state', '{"ms": [' + ', '.join(values) + ']}', comments)
            self.note(f'Qubit.qubit[].{key} -> Qubit.{side}_state {{"ms": [{", ".join(values)}]}}')

    # ---- A94 (2026-10-16): the input keys review ((internal review)) ----------------------------------
    FREQ_UNITS = ('Hz', 'kHz', 'MHz', 'GHz')
    GYRO_UNITS = ('radkHz/G', 'rad/ms/G', 'kHz/G', 'MHz/G', 'MHz/T', 'Hz/T', 'radMHz/T', 'rad/us/T', 'rad/s/T')
    EQ_UNITS = ('1e-30 m^2', '1e-30m^2', '10^-30 m^2', '10^-30m^2', 'fm^2', 'fm2', '10 mbarn', 'm^2', 'm2', 'cm^2', 'cm2',
                'barn', 'b', 'mbarn', 'millibarn', 'mb')
    EFG_UNITS = ('Hartree/Bohr^2', 'Ha/Bohr^2', 'Ha/a0^2', 'a.u.', 'au', 'atomic_unit', 'V/m^2', 'V/m2', 'V/angstrom^2',
                 'V/angstrom2', 'V/A^2', 'V/A2')

    @staticmethod
    def obj(member):
        return member.value if member is not None and isinstance(member.value, Obj) else None

    def unit_case(self, block, path, units):
        """A94 (R3): units are case-sensitive now; a spelling that differs in case only is written as in the table."""
        m = block.get('unit') if block is not None else None
        if m is None:
            return
        v = plain(m.value)
        if not isinstance(v, str) or v in units:
            return
        same = [u for u in units if u.lower() == v.lower()]
        if same:
            m.value = Scalar(json.dumps(same[0]))
            block.dirty = True
            self.note(f'{path}.unit "{v}" -> "{same[0]}" (units are case-sensitive)')

    def a94(self, R, top):
        """A94 (2026-10-16): configuration.nstate N -> s_range [1, N] (0: removed, the ensemble); a seed number of a k_range
        run that draws -> a list of the same number (one per configuration: the seed number is now the first configuration's,
        k taking seed + k - first); the Qubit.file form -> a Qubit.qubit entry {"file", "name", "spin"}; Qubit.qubit[].detuning
        [x] -> x; pulse.bath.spin removed, alphams/betams -> levels; cluster.pcce.print_trials a number (print_limit gone);
        Defect.species[].axis -> zfs.axis, zfs.E 0 removed, units in the table's case; cluster.adaptive.output -> workdir,
        - in its keys -> _,
        policy "merge" removed; auto_convergence: converge removed, time_window.unit into value, time_window.threshold ->
        decay_level, sampling.number -> configurations, test_number -> test_configurations, output.path and output.file of the
        base input removed (the runner never used them)."""
        C = self.obj(R.get('configuration'))
        A = self.obj(R.get('auto_convergence'))
        autoconv = A is not None and (A.get('enabled') is None or plain(A.get('enabled').value) is not False)
        d = self.obj(R.get('Defect'))
        species = d.get('species') if d is not None else None
        nspecies = len(species.value.items) if species is not None and isinstance(species.value, Arr) else 0
        if C is not None:
            m = C.get('nstate')
            if m is not None:
                v = plain(m.value)
                if not isinstance(v, int) or isinstance(v, bool) or v < 0:
                    raise ValueError(f'configuration.nstate {json.dumps(v)} is not a nonnegative integer')
                if v == 0:
                    C.drop('nstate')
                    self.note('configuration.nstate 0 removed (the ensemble: no configuration.s_range)')
                else:
                    self.rename(m, 's_range')
                    m.value = Scalar(f'[1, {v}]')
                    C.dirty = True
                    self.note(f'configuration.nstate {v} -> configuration.s_range [1, {v}]')
            k, s = C.get('k_range'), C.get('seed')
            kr = plain(k.value) if k is not None else None
            sv = plain(s.value) if s is not None else None
            cm = self.obj(R.get('cluster'))
            mode = cm.get('mode') if cm is not None else None
            pcce = self.shape == 'pcce' or (mode is not None and str(plain(mode.value)).lower() == 'pcce')
            draws = C.get('s_range') is not None or nspecies > 0 or pcce
            if (isinstance(kr, list) and len(kr) == 2 and kr[1] > kr[0] and isinstance(sv, int) and not isinstance(sv, bool)
                    and draws and not autoconv):
                n = kr[1] - kr[0] + 1
                s.value = Scalar('[' + ', '.join([s.value.render()] * n) + ']')
                C.dirty = True
                self.note(f'configuration.seed {sv} -> a list of {n} times {sv}: until A94 one seed number was every '
                          f'configuration\'s, now it is the first one\'s (k takes seed + k - first)')
        q = R.get('Qubit')
        Q = self.obj(q)
        if Q is not None and Q.get('file') is not None:
            f = Q.get('file')
            parts = [f'"file": {f.value.render()}']
            for key in ('name', 'spin'):
                m = Q.get(key)
                if m is not None:
                    parts.append(f'{json.dumps(key)}: {self.inline(m.value)}')
            self.rename(f, 'qubit')
            f.value = Scalar('[{' + ', '.join(parts) + '}]')
            for key in ('name', 'spin'):
                if Q.get(key) is not None:
                    Q.drop(key)
            Q.dirty = True
            self.note('Qubit.file' + ''.join(f' and Qubit.{k}' for k in ('name', 'spin') if f'"{k}"' in ''.join(parts))
                      + ' -> Qubit.qubit [{"file": ...}] (an entry with the position from the file)')
        if Q is not None:
            qs = Q.get('qubit')
            if qs is not None and isinstance(qs.value, Arr):
                for n, item in enumerate(qs.value.items):
                    e = item.value
                    if not isinstance(e, Obj):
                        continue
                    m = e.get('detuning')
                    if m is not None and isinstance(m.value, Arr):
                        if len(m.value.items) != 1:
                            raise ValueError(f'Qubit.qubit[{n}].detuning must be a number or [number]')
                        m.value = m.value.items[0].value
                        e.dirty = True
                        self.note(f'Qubit.qubit[{n}].detuning [x] -> x')
            im = Q.get('intmap')
            if im is not None and isinstance(im.value, Arr):
                for n, item in enumerate(im.value.items):
                    if isinstance(item.value, Obj):
                        self.unit_case(item.value, f'Qubit.intmap[{n}]', self.FREQ_UNITS)
        p = self.obj(R.get('pulse'))
        b = self.obj(p.get('bath')) if p is not None else None
        if b is not None:
            if b.drop('spin') is not None:
                self.note('pulse.bath.spin removed (the spin of the driven bath spins is used; it had to be the same)')
            a, z = b.get('alphams'), b.get('betams')
            if a is not None or z is not None:
                first = a if a is not None else z
                text = (f'[{a.value.render() if a is not None else "0.5"}, '
                        f'{z.value.render() if z is not None else "-0.5"}]')
                self.rename(first, 'levels')
                first.value = Scalar(text)
                if a is not None and z is not None:
                    b.drop('betams')
                b.dirty = True
                self.note(f'pulse.bath.alphams/betams -> pulse.bath.levels {text}')
        c = self.obj(R.get('cluster'))
        pc = self.obj(c.get('pcce')) if c is not None else None
        if pc is not None:
            t, lim = pc.get('print_trials'), pc.get('print_limit')
            tv = plain(t.value) if t is not None else None
            if t is not None and isinstance(tv, bool):
                on = bool(tv)
                if on:
                    t.value = Scalar(lim.value.render() if lim is not None else '32')
                    pc.dirty = True
                    self.note(f'cluster.pcce.print_trials true -> {t.value.render()} (the number of trials printed, '
                              'print_limit before)')
                else:
                    pc.drop('print_trials')
                    self.note('cluster.pcce.print_trials false removed (0, the default)')
            if lim is not None:
                pc.drop('print_limit')
                self.note('cluster.pcce.print_limit removed (it is the number of print_trials now)')
        ad = self.obj(c.get('adaptive')) if c is not None else None
        if ad is not None:
            for m in list(ad.items):
                if '-' in m.key:
                    new = m.key.replace('-', '_')
                    self.note(f'cluster.adaptive.{m.key} -> {new}')
                    self.rename(m, new)
                    ad.dirty = True
            m = ad.get('output')
            if m is not None:
                self.rename(m, 'workdir')
                ad.dirty = True
                self.note('cluster.adaptive.output -> cluster.adaptive.workdir')
            m = ad.get('policy')
            if m is not None and str(plain(m.value)).lower() == 'merge':
                ad.drop('policy')
                self.note('cluster.adaptive.policy "merge" removed (the only policy, the default)')
        if species is not None and isinstance(species.value, Arr):
            for n, item in enumerate(species.value.items):
                e = item.value
                if not isinstance(e, Obj):
                    continue
                zm = e.get('zfs')
                z = zm.value if zm is not None and isinstance(zm.value, Obj) else None
                if z is not None and z.get('E') is not None:
                    ev = plain(z.get('E').value)
                    if ev != 0:
                        raise ValueError(f'Defect.species[{n}].zfs.E is {ev}: CCEX refused an E other than 0 (write zfs.tensor)')
                    z.drop('E')
                    self.note(f'Defect.species[{n}].zfs.E 0 removed')
                ax = e.get('axis')
                if ax is not None:
                    if z is not None and z.get('D') is not None:
                        text = self.inline(ax.value)
                        if self.member_indent(z) is None:
                            mid = z.items[-1].mid.lstrip(' ') if z.items and z.items[-1].mid.strip() == ':' else ': '
                            z.items.insert(1 if z.items else 0, Member(' ', '"axis"', mid, Scalar(text), ''))
                            z.dirty = True
                        else:
                            self.put(z, self.indent_of(e, top), 'axis', value=text)
                        e.drop('axis')
                        self.note(f'Defect.species[{n}].axis -> Defect.species[{n}].zfs.axis')
                    else:
                        e.drop('axis')
                        self.note(f'Defect.species[{n}].axis removed (only the D form of zfs used it)')
                    e.dirty = True
                for key, units in (('gyros', self.GYRO_UNITS), ('quadrupole_moments', self.EQ_UNITS),
                                   ('hyperfine', self.FREQ_UNITS), ('efg', self.EFG_UNITS), ('zfs', self.FREQ_UNITS),
                                   ('detuning', self.FREQ_UNITS)):
                    self.unit_case(self.obj(e.get(key)), f'Defect.species[{n}].{key}', units)
        if A is not None:
            for block in ('time_window', 'rdip', 'rbath', 'sampling'):
                B = self.obj(A.get(block))
                if B is None:
                    continue
                cv = B.get('converge')
                if cv is not None:
                    fixed = 'configurations' if block == 'sampling' else 'value'
                    has = B.get(fixed) is not None or (block == 'sampling' and B.get('number') is not None)
                    if plain(cv.value) is True and has:
                        raise ValueError(f'auto_convergence.{block}: converge true with a fixed value was refused')
                    B.drop('converge')
                    self.note(f'auto_convergence.{block}.converge removed (a value given is fixed, none is found)')
                if block == 'time_window':
                    u = B.get('unit')
                    if u is not None:
                        v = B.get('value')
                        if v is None or isinstance(v.value, Obj):
                            raise ValueError('auto_convergence.time_window.unit without a plain value (or beside a value '
                                             'object) was refused')
                        v.value = Scalar('{"value": ' + v.value.render() + ', "unit": ' + u.value.render() + '}')
                        B.drop('unit')
                        self.note('auto_convergence.time_window.unit -> time_window.value {"value", "unit"}')
                    m = B.get('threshold')
                    if m is not None:
                        self.rename(m, 'decay_level')
                        self.note('auto_convergence.time_window.threshold -> decay_level')
                if block == 'sampling':
                    for old, new in (('number', 'configurations'), ('test_number', 'test_configurations')):
                        m = B.get(old)
                        if m is not None:
                            self.rename(m, new)
                            self.note(f'auto_convergence.sampling.{old} -> {new}')
                B.dirty = True
            out = self.obj(R.get('output'))
            if autoconv and out is not None:
                for key in ('path', 'file'):
                    if out.drop(key) is not None:
                        self.note(f'output.{key} removed (the automatic convergence put its runs in auto_convergence.workdir '
                                  'and never used it)')
                if not out.items:
                    R.drop('output')

    def auto_convergence(self, R, top):
        """A65: launcher removed (the runs use the mpirun -n N main.out was started with); configurations -> sampling; its
        file patterns go to the configuration block: bathfile -> bathfile, statefile -> bath_state_file, avaaxfile ->
        defect_axis_file, exstatefile -> defect_state_file (lists; without k_range the tool finds the k that exist)."""
        a = R.get('auto_convergence')
        if a is None or not isinstance(a.value, Obj):
            return
        A = a.value
        if A.drop('launcher') is not None:
            self.note('auto_convergence.launcher removed (the runs use the mpirun -n N main.out was started with)')
        c = A.get('configurations')
        if c is None or not isinstance(c.value, Obj):
            return
        C = c.value
        block = 'configuration'
        for old, new in (('bathfile', 'bathfile'), ('statefile', 'bath_state_file'), ('avaaxfile', 'defect_axis_file'),
                         ('exstatefile', 'defect_state_file')):
            m = C.get(old)
            if m is None:
                continue
            C.take(old)
            B = self.child_block(R, top, block)
            if B.get(new) is not None:
                raise ValueError(f'auto_convergence.configurations.{old} beside {block}.{new}: keep one')
            value = m.value if isinstance(m.value, Container) else Scalar('[' + m.value.render() + ']')
            self.put(B, top, new, m, value, top)
            self.note(f'auto_convergence.configurations.{old} -> {block}.{new}')
        if C.items:                    # a line left blank by the keys taken out of the block's first line goes
            first = C.items[0]
            lines = first.pre.split('\n')
            if len(lines) > 2 and lines[0] == '' and not lines[1].strip():
                first.pre = '\n'.join([''] + lines[2:])
                C.dirty = True
        self.rename(c, 'sampling')
        A.dirty = True
        self.note('auto_convergence.configurations -> auto_convergence.sampling')

    def bath_files(self):
        c = self.root.get('configuration')
        b = c.value.get('bathfile') if c is not None and isinstance(c.value, Obj) else None
        v = plain(b.value) if b is not None else None
        return len(v) if isinstance(v, list) else None

    def per_bath_file(self, block, path, key):
        """A64: bath_state_file, Defect.axis_file and Defect.state_file are lists of one name per bath file. A single name
        becomes ["name"] with one bath file; with more, CCEX did not read it (the states and axes were drawn at random),
        so it is removed."""
        m = block.get(key)
        if m is None or isinstance(m.value, Container):
            return
        n = self.bath_files()
        if n is None or n == 1:
            m.value = Scalar('[' + m.value.render() + ']')
            block.dirty = True
            self.note(f'{path}{key}: "name" -> ["name"] (one per bath file)')
        else:
            block.drop(key)
            self.note(f'{path}{key} removed: with {n} bath files it was not read (drawn at random); give one per bath file')

    def move(self, R, top, block_key, pairs, position_keys=None, values=None):
        """Move R[old] to R[block_key][new] for (old, new) in pairs; returns the block (None if nothing moved)."""
        found = [(R.find(old)[0], old, new) for old, new in pairs if R.find(old)[1] is not None]
        if not found:
            return None
        pos = min(i for i, _, _ in found)
        first_pre = R.items[pos].pre                        # a block in place of a key written inline stays inline
        moved = []
        for _, old, new in sorted(found, key=lambda x: [p[0] for p in pairs].index(x[1])):
            _, m = R.take(old)
            moved.append((old, new, m))
        pos = min(pos, len(R.items))
        made = R.get(block_key) is None
        block = self.child_block(R, top, block_key, pos)
        if made and '\n' not in first_pre:
            R.items[pos].pre = first_pre
        for old, new, m in moved:
            value = values(old, m) if values else None
            if isinstance(new, tuple):                      # a block inside the block
                sub = self.child_block(block, self.indent_of(block, top), new[0])
                self.put(sub, self.indent_of(block, top), new[1], m, value, top)
                self.note(f'{old} -> {block_key}.{new[0]}.{new[1]}')
            else:
                self.put(block, top, new, m, value, top)
                self.note(f'{old} -> {block_key}.{new}')
        return block

    def blocks(self, R, top):
        pcce = [('sK', 'partition_size'), ('max_trial', 'trials'), ('max_iter', 'max_iterations'),
                ('kmeans_pp', 'kmeans_pp'), ('iter_detail', 'print_trials'), ('pcce_detail_limit', 'print_limit'),
                ('pcce_partition', 'partition'), ('pcce_connections', 'connections'),
                ('pcce_cache_mode', ('cache', 'mode')), ('pcce_cache_file', ('cache', 'file'))]
        if self.shape != 'pcce':
            for old, _ in pcce:
                m = R.get(old)
                if m is None or (old == 'pcce_cache_mode' and str(plain(m.value)).lower() != 'off'):
                    continue
                R.drop(old)
                self.note(f'{old} removed (it is read with cluster_mode pcce only)')
        self.move(R, top, 'pcce', pcce,
                  values=lambda old, m: self.switch(m, old) if old in ('kmeans_pp', 'iter_detail') else None)
        self.move(R, top, 'bundle', [('bundle_max_dim', 'max_dim'), ('bundle_merge_cutoff', 'merge_cutoff'),
                                     ('bundle_connections', 'connections'),
                                     ('bundle_connection_budget', 'connection_budget'),
                                     ('expansion_max_clusters', 'max_clusters'), ('expansion_max_work', 'max_work')])
        if R.get('adaptive_work_file') is not None:
            self.move(R, top, 'adaptive', [('adaptive_work_file', 'work_file')])
        self.move(R, top, 'pulse', [('npulse', 'count'), ('pulsename', 'name'), ('sequence', 'sequence'),
                                    ('pulseiter', 'per_qubit'), ('pulse_time', 'duration'),
                                    ('detuning_factor', 'detuning_factor')],
                  values=lambda old, m: self.switch(m, old) if old == 'pulseiter' else None)
        for name, inner, sub in (('gcce_refinement', 'gcce_refinement_max_sweeps', 'max_sweeps'),
                                 ('cluster_contributions', 'cluster_contribution_max_mb', 'max_mb')):
            i, m = R.find(name)
            _, n = R.find(inner)
            if m is not None and isinstance(m.value, Obj):
                continue
            if m is None and n is None:
                continue
            on = m.value.render() if m is not None else 'false'
            if m is not None and not isinstance(plain(m.value), bool):   # 0/1 too, as CCEX (readStrictBool)
                raise ValueError(f'{name} must be true or false')
            pos = i if i is not None else R.find(inner)[0]
            if m is not None:
                R.take(name)
            block = self.new_block(top)
            self.put(block, top, 'enabled', value=on)
            if n is not None:
                R.take(inner)
                self.put(block, top, sub, n, None, top)
            mem = self.new_member(name, block, top, moved=m) if m is not None else self.new_member(name, block, top)
            mem.value, mem.shift = block, 0
            R.insert(min(pos, len(R.items)), mem)
            R.dirty = True
            self.note(f'{name}' + (f' and {inner}' if n is not None else '') + f' -> {name} {{enabled, {sub}}}')

    def qubit(self, R, top):
        i, qf = R.find('qubitfile')
        if qf is not None:
            if R.get('Qubit') is not None:
                raise ValueError('the input has both qubitfile and a Qubit block (CCEX refused it before too)')
            self.move(R, top, 'Qubit', [('qubitfile', 'file'), ('qspin', 'spin'), ('qzfs', 'zfs'),
                                        ('qalphams', 'alphams'), ('qbetams', 'betams')])
            return
        _, q = R.find('Qubit')
        if q is None or not isinstance(q.value, Obj):
            for k in ('qspin', 'qzfs', 'qalphams', 'qbetams'):
                if R.get(k) is not None:
                    R.drop(k)
                    self.note(f'{k} removed (it was read with qubitfile only)')
            return
        Q = q.value
        qind = self.indent_of(Q, top)
        n = Q.drop('nqubit')
        if n is not None:
            qubits = Q.get('qubit')
            count = len(plain(qubits.value)) if qubits is not None else None
            if count != plain(n.value):
                raise ValueError(f'Qubit.nqubit is {plain(n.value)} but Qubit.qubit has {count} entries')
            self.note('Qubit.nqubit removed (the length of Qubit.qubit)')
        # the states: psia/psib were used when both were given, alphaidx/betaidx when both were given, else neither
        names = {'psia': None, 'psib': None, 'alphaidx': None, 'betaidx': None}
        first = None
        for k in names:
            j, m = Q.find(k)
            if m is not None:
                first = j if first is None else min(first, j)
        for k in names:
            names[k] = Q.take(k)[1]
        if first is not None:
            first = min(first, len(Q.items))
            if names['psia'] is not None and names['psib'] is not None:
                use = (('psia', 'alpha_state', 'vector'), ('psib', 'beta_state', 'vector'))
                dropped = [k for k in ('alphaidx', 'betaidx') if names[k] is not None]
            elif names['alphaidx'] is not None and names['betaidx'] is not None:
                use = (('alphaidx', 'alpha_state', 'index'), ('betaidx', 'beta_state', 'index'))
                dropped = [k for k in ('psia', 'psib') if names[k] is not None]
            else:
                use = ()
                dropped = [k for k in names if names[k] is not None]
            for k in dropped:
                self.note(f'Qubit.{k} removed (it was not used: the run took ' +
                          ('psia/psib' if use and use[0][2] == 'vector' else
                           'alphaidx/betaidx' if use else 'alphams/betams') + ')')
            for old, new, form in reversed(use):
                m = names[old]
                inner = self.new_block(qind)
                self.put(inner, qind, form, m, None, qind)
                mem = self.new_member(new, inner, qind)
                Q.insert(first, mem)
                Q.dirty = True
            for old, new, form in use:
                self.note(f'Qubit.{old} -> Qubit.{new}.{form}')
        self.qubit_overhaus = Q.take('overhaus')[1]

    def hamiltonian(self, R, top):
        common, ccce = [], []
        positions = []

        def take(key):
            i, m = R.find(key)
            if m is not None:
                positions.append(i)
                R.take(key)
            return m

        ovh = getattr(self, 'qubit_overhaus', None)
        if ovh is not None:
            common.append(('qubit_bath_mean_field', ovh, self.switch(ovh, 'Qubit.overhaus'), 'Qubit.overhaus'))
        vmf = take('vector_mean_field')
        if vmf is not None:
            v = plain(vmf.value)
            common.append(('bath_bath_mean_field', vmf, '"xyz"' if v in (True, 1) else '"z"', 'vector_mean_field'))
        m = R.get('hf_readmode')
        mode = plain(m.value) if m is not None else 0
        if mode not in (0, 1, 2, 3):
            raise ValueError(f'hf_readmode {mode} is not 0, 1, 2 or 3')
        if mode == 0:
            for k in ('hf_tensorfile', 'hf_ignore_oor'):
                if R.drop(k) is not None:
                    self.note(f'{k} removed (it is not read with hf_readmode 0)')
            if R.drop('hf_readmode') is not None:
                self.note('hf_readmode 0 removed (the default: hyperfine dipole "point", contact false)')
        else:
            hf = {k: take(k) for k in ('hf_readmode', 'hf_tensorfile', 'hf_ignore_oor')}
            block = self.new_block('')
            common.append(('hyperfine', hf['hf_readmode'], block, 'hf_readmode'))
            self.hf_block = (block, mode, hf)
        m = R.get('qd_readmode')
        mode = plain(m.value) if m is not None else 0
        if mode not in (0, 2, 3, 4):
            raise ValueError(f'qd_readmode {mode} is not 0, 2, 3 or 4')
        unread = ('qd_tensorfile', 'qd_tensorfile_woqubit', 'qd_cellpara') if mode == 0 else \
                 ('qd_tensorfile_woqubit', 'qd_cellpara') if mode == 2 else ()
        for k in unread:
            if R.drop(k) is not None:
                self.note(f'{k} removed (it is not read with qd_readmode {mode})')
        if mode == 0:
            if R.drop('qd_readmode') is not None:
                self.note('qd_readmode 0 removed (the default: quadrupole "off")')
        else:
            qd = {k: take(k) for k in ('qd_readmode', 'qd_tensorfile', 'qd_tensorfile_woqubit', 'qd_cellpara')}
            block = self.new_block('')
            common.append(('quadrupole', qd['qd_readmode'], block, 'qd_readmode'))
            self.qd_block = (block, mode, qd)
        for k in ('knight', 'hfmedi'):
            if self.gcce:
                if R.drop(k) is not None:
                    self.note(f'{k} removed (it has no effect with method gCCE)')
                continue
            m = take(k)
            if m is not None:
                ccce.append((k, m, self.switch(m, k), k))
        if not common and not ccce:
            return
        pos = min(positions) if positions else len(R.items)
        H = self.child_block(R, top, 'hamiltonian', min(pos, len(R.items)))
        hind = self.indent_of(H, top)
        for name, items in (('common', common), ('cCCE', ccce)):
            if not items:
                continue
            B = self.child_block(H, hind, name)
            bind = self.indent_of(B, hind)
            for key, moved, value, old in items:
                if isinstance(value, Obj):
                    value.close = '\n' + bind
                    self.fill(value, bind, key, top)
                self.put(B, hind, key, moved, value, top)
                self.note(f'{old} -> hamiltonian.{name}.{key}')

    def fill(self, block, bind, key, top):
        if key == 'hyperfine':
            _, mode, hf = self.hf_block
            self.put(block, bind, 'dipole', value='"dft"' if mode in (2, 3) else '"point"')
            self.put(block, bind, 'contact', value='true' if mode in (1, 3) else 'false')
            if mode != 0:
                if hf['hf_tensorfile'] is not None:
                    self.put(block, bind, 'file', hf['hf_tensorfile'], None, top)
                    self.note('hf_tensorfile -> hamiltonian.common.hyperfine.file')
                if hf['hf_ignore_oor'] is not None:
                    v = plain(hf['hf_ignore_oor'].value)
                    if v not in (0, 1):
                        raise ValueError(f'hf_ignore_oor {v} is not 0 or 1')
                    self.put(block, bind, 'missing_row', hf['hf_ignore_oor'], '"point_dipole"' if v else '"stop"', top)
                    self.note('hf_ignore_oor -> hamiltonian.common.hyperfine.missing_row')
        else:
            _, mode, qd = self.qd_block
            self.put(block, bind, 'structure', value=json.dumps({2: 'plain', 3: 'strained', 4: 'bubble'}[mode]))
            for old, new in (('qd_tensorfile', 'file'), ('qd_tensorfile_woqubit', 'file_without_defect'),
                             ('qd_cellpara', 'cell')):
                if qd[old] is not None and (mode != 2 or old == 'qd_tensorfile'):
                    self.put(block, bind, new, qd[old], None, top)
                    self.note(f'{old} -> hamiltonian.common.quadrupole.{new}')

    def defect(self, R, top):
        i, d = R.find('Defect')
        ex, ax = R.get('exstatefile'), R.get('avaaxfile')
        if d is None:
            for k in ('exstatefile', 'avaaxfile'):
                if R.get(k) is not None:
                    R.drop(k)
                    self.note(f'{k} removed (it is read with Defect species only)')
            return
        species = d.value
        if isinstance(species, Arr):
            renames = {'dfname': 'name', 'naddspin': 'n_nuclear_spins', 'navaax': 'n_configurations',
                       'eqs': 'quadrupole_moments', 'rxyzs': 'nuclear_positions', 'hypf': 'hyperfine'}
            for n, item in enumerate(species.items):
                e = item.value
                if not isinstance(e, Obj):
                    continue
                for m in list(e.items):
                    if m.key in renames:
                        self.rename(m, renames[m.key])
                        e.dirty = True
                if e.drop('apprx') is not None:
                    self.note(f'Defect[{n}].apprx removed (it changed nothing)')
            if any(isinstance(it.value, Obj) and it.value.dirty for it in species.items):
                self.note('Defect entries: dfname -> name, naddspin -> n_nuclear_spins, navaax -> n_configurations, '
                          'eqs -> quadrupole_moments, rxyzs -> nuclear_positions, hypf -> hyperfine')
            block = self.new_block(top)
            for key, old in (('state_file', 'exstatefile'), ('axis_file', 'avaaxfile')):
                if R.get(old) is not None:
                    _, m = R.take(old)
                    self.put(block, top, key, m, None, top)
                    self.note(f'{old} -> Defect.{key}')
            ind = self.indent_of(block, top)
            mem = Member('\n' + ind, '"species"', self.aligned(self.mid, self.width, '"species"'), species, '',
                         len(ind) - len(top))
            block.items.append(mem)
            d.value = block
            R.dirty = True
            self.note('Defect [...] -> Defect {"species": [...]}')
        elif isinstance(species, Obj):
            for key, old in (('state_file', 'exstatefile'), ('axis_file', 'avaaxfile')):
                if R.get(old) is not None:
                    _, m = R.take(old)
                    self.put(species, top, key, m, None, top)
                    self.note(f'{old} -> Defect.{key}')


def convert_text(text):
    c = Converter(text)
    new = c.convert()
    json.loads(strip_comments(new))          # what CCEX will read has to be JSON
    return new, c.log


def main(argv):
    mode = 'write'
    if argv and argv[0] in ('--check', '--stdout'):
        mode, argv = argv[0][2:], argv[1:]
    if not argv or argv[0] in ('-h', '--help'):
        print(__doc__)
        return 0 if argv else 2
    status = 0
    for path in argv:
        try:                   # a file that cannot be read (missing, not UTF-8, ...) is reported like one that cannot be converted
            with open(path, newline='') as f:
                raw = f.read()
            crlf = '\r\n' in raw          # converted with \n and written back with \r\n (lines added by the converter too)
            text = raw.replace('\r\n', '\n') if crlf else raw
            new, log = convert_text(text)
        except Exception as error:
            print(f'{path}: not converted: {error}', file=sys.stderr)
            status = 2
            continue
        out = new.replace('\n', '\r\n') if crlf else new
        if mode == 'stdout':
            sys.stdout.write(out)
            continue
        if new == text:
            print(f'{path}: already in the new form')
            continue
        print(f'{path}:')
        for line in log:
            print(f'    {line}')
        if mode == 'check':
            status = max(status, 1)
            continue
        try:
            write_atomically(path, out)
        except Exception as error:
            print(f'{path}: not written (the file is unchanged): {error}', file=sys.stderr)
            status = 2
    return status


def write_atomically(path, text):
    """A temporary file beside the target, then os.replace: an interrupted write leaves the original whole."""
    target = os.path.realpath(path)       # through a symbolic link, as open(path, 'w') did
    fd, tmp = tempfile.mkstemp(dir=os.path.dirname(target), prefix='.' + os.path.basename(target) + '.')
    try:
        with os.fdopen(fd, 'w', newline='') as f:
            f.write(text)
        shutil.copymode(target, tmp)
        os.replace(tmp, target)
    except BaseException:
        os.unlink(tmp)
        raise


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
