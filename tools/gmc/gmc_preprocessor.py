"""
GSPICE Model Compiler (GMC) - Verilog-A text preprocessor.

Resolves the Verilog-AMS preprocessor directives that appear in real compact
model sources (BSIM3/BSIM4 and friends) before the GMC lexer/parser runs:

  `define  NAME <body>            object-like macro
  `define  NAME(a, b) <body>      function-like macro (body may span lines
                                  using trailing-backslash continuations)
  `undef   NAME
  `ifdef   NAME / `ifndef NAME / `else / `endif
  `include "file"                 resolved relative to the including file,
                                  then against the configured include dirs
  `resetall                       clears all defined macros

Macros are expanded lazily at the point of use: an expanded body is rescanned
for further macro references, and function-like macro arguments are expanded
before substitution. Directives must start a line (ignoring leading
whitespace); anything else is ordinary text that is passed through macro
expansion.

The output keeps one line per input line so diagnostics stay line-aligned.
"""

import re
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

_DIRECTIVES = {
    "define", "undef", "ifdef", "ifndef", "else", "endif",
    "include", "resetall",
}

_MACRO_TOKEN = re.compile(r"`\s*([A-Za-z_][A-Za-z0-9_]*)")

_MAX_EXPANSIONS = 100000


class _Macro:
    __slots__ = ("kind", "params", "body")

    def __init__(self, kind: str, params: Sequence[str], body: str):
        self.kind = kind          # "object" | "function"
        self.params = list(params)
        self.body = body


class PreprocessedResult:
    """Preprocessed source text plus provenance for diagnostics."""

    __slots__ = ("text", "line_map", "files")

    def __init__(self, text: str, line_map: List[Tuple[str, int]],
                 files: List[str]):
        self.text = text
        self.line_map = line_map          # per output line: (file, physical line)
        self.files = files                # canonical paths of all files read


class VerilogAPreprocessor:
    def __init__(self, include_dirs: Sequence[Path] = (), max_depth: int = 64):
        self.include_dirs = [Path(d) for d in include_dirs]
        self.max_depth = max_depth
        self.macros = {}                  # name -> _Macro
        self.stack: List[List] = []       # conditional stack: [parent_active, taken]
        self.active = True
        self._output_lines: List[str] = []
        self._line_map: List[Tuple[str, int]] = []
        self._file_order: List[str] = []
        self._file_seen = set()
        self._include_stack: List[str] = []

    # ------------------------------------------------------------------ #
    # Public entry
    # ------------------------------------------------------------------ #
    def run(self, text: str, filename: Optional[str] = None) -> None:
        self._run_file(text, filename, 0)
        if self.stack:
            raise ValueError(
                "unterminated `ifdef (missing `endif) at top level")

    def result(self) -> PreprocessedResult:
        return PreprocessedResult(
            "\n".join(self._output_lines), self._line_map, list(self._file_order))

    # ------------------------------------------------------------------ #
    # File / line processing
    # ------------------------------------------------------------------ #
    def _run_file(self, text: str, filename: Optional[str], depth: int) -> None:
        if depth > self.max_depth:
            raise ValueError("Verilog-A `include nesting too deep "
                             f"(limit {self.max_depth})")
        fname = filename or "<source>"
        if filename is not None:
            canon = str(Path(filename))
            if canon not in self._file_seen:
                self._file_seen.add(canon)
                self._file_order.append(canon)
        entry_depth = len(self.stack)
        in_block = False
        lines = text.splitlines()
        i = 0
        while i < len(lines):
            # Build one logical line (joining trailing-backslash continuations).
            logical = []
            while i < len(lines):
                stripped, in_block = self._strip_line_comments(lines[i], in_block)
                if not in_block and stripped.rstrip().endswith("\\"):
                    logical.append(stripped.rstrip()[:-1])
                    i += 1
                    continue
                logical.append(stripped)
                i += 1
                break
            line = (" ".join(p for p in logical if p.strip())
                    if len(logical) > 1 else logical[0])
            stripped_line = line.lstrip()
            is_directive = False
            if stripped_line.startswith("`"):
                m = _MACRO_TOKEN.match(stripped_line)
                keyword = m.group(1) if m else ""
                is_directive = keyword in _DIRECTIVES
            if is_directive:
                self._handle_directive(stripped_line, fname, i, depth)
                continue
            if not self.active:
                continue
            # A macro invocation (or plain expression) may span physical
            # lines without a backslash continuation: pull more lines while
            # an opening parenthesis stays unbalanced.
            while self._unclosed_parens(line) > 0 and i < len(lines):
                extra, in_block = self._strip_line_comments(lines[i], in_block)
                if not in_block and extra.rstrip().endswith("\\"):
                    extra = extra.rstrip()[:-1]
                i += 1
                line = line + " " + extra
            self._emit(self._expand(line, fname, i), fname, i)
        if in_block:
            raise ValueError(f"unterminated /* comment in {fname}")
        if len(self.stack) != entry_depth:
            raise ValueError(f"unterminated or unbalanced `ifdef/`endif "
                             f"across file {fname}")

    def _emit(self, text: str, fname: str, lineno: int) -> None:
        self._output_lines.append(text)
        self._line_map.append((fname, lineno))

    # ------------------------------------------------------------------ #
    # Comments
    # ------------------------------------------------------------------ #
    @staticmethod
    def _strip_line_comments(line: str, in_block: bool) -> Tuple[str, bool]:
        out = []
        i = 0
        n = len(line)
        while i < n:
            if in_block:
                j = line.find("*/", i)
                if j < 0:
                    return "".join(out), True
                in_block = False
                i = j + 2
                continue
            if line[i:i + 2] == "/*":
                in_block = True
                i += 2
                continue
            if line[i:i + 2] == "//":
                break
            out.append(line[i])
            i += 1
        return "".join(out), in_block

    # ------------------------------------------------------------------ #
    # Directives
    # ------------------------------------------------------------------ #
    def _handle_directive(self, s: str, fname: str, lineno: int,
                          depth: int) -> None:
        m = _MACRO_TOKEN.match(s)
        directive = m.group(1)
        rest = s[m.end():]

        if directive in ("ifdef", "ifndef"):
            cond_name = rest.strip()
            parent = self.active
            defined = cond_name in self.macros
            if directive == "ifndef":
                defined = not defined
            taken = parent and defined
            self.stack.append([parent, taken])
            self.active = taken
            return

        if directive == "else":
            if not self.stack:
                raise ValueError(f"unexpected `else at {fname}:{lineno}")
            parent, taken = self.stack[-1]
            if taken:
                self.stack[-1] = [parent, False]
                self.active = False
            else:
                self.stack[-1] = [parent, True]
                self.active = parent
            return

        if directive == "endif":
            if not self.stack:
                raise ValueError(f"unexpected `endif at {fname}:{lineno}")
            parent, _ = self.stack.pop()
            self.active = parent
            return

        if directive == "include":
            if not self.active:
                return
            self._handle_include(rest, fname, lineno, depth)
            return

        if not self.active:
            return

        if directive == "resetall":
            self.macros.clear()
            return

        if directive == "undef":
            name = rest.strip().split(None, 1)[0] if rest.strip() else ""
            self.macros.pop(name, None)
            return

        if directive == "define":
            self._handle_define(rest, fname, lineno)
            return

        raise ValueError(f"unknown Verilog-A directive `{directive} "
                         f"at {fname}:{lineno}")

    def _handle_define(self, rest: str, fname: str, lineno: int) -> None:
        body_start = rest.lstrip()
        m = re.match(r"([A-Za-z_][A-Za-z0-9_]*)", body_start)
        if not m:
            raise ValueError(f"malformed `define at {fname}:{lineno}")
        name = m.group(1)
        tail = body_start[m.end():]
        if tail and tail[0] == "(":
            end = self._match_paren(tail, 0, fname, lineno)
            arg_list = self._split_args(tail[1:end])
            macro_body = tail[end + 1:].strip()
            self.macros[name] = _Macro("function", arg_list, macro_body)
        else:
            self.macros[name] = _Macro("object", [], tail.strip())

    # ------------------------------------------------------------------ #
    # Include
    # ------------------------------------------------------------------ #
    def _handle_include(self, rest: str, fname: str, lineno: int,
                        depth: int) -> None:
        m = re.search(r'["<]([^">]+)[">]', rest.strip())
        if not m:
            raise ValueError(f"malformed `include at {fname}:{lineno}")
        name = m.group(1)
        resolved = self._resolve_include(name, fname)
        if resolved is None:
            raise ValueError(
                f"`include file not found: {name} "
                f"(from {fname}:{lineno}, search: "
                f"{[str(d) for d in self.include_dirs]})")
        if resolved in self._include_stack:
            raise ValueError(f"recursive `include of {name} "
                             f"at {fname}:{lineno}")
        source = Path(resolved).read_text(encoding="utf-8", errors="replace")
        self._include_stack.append(resolved)
        try:
            self._run_file(source, resolved, depth + 1)
        finally:
            self._include_stack.pop()

    def _resolve_include(self, name: str, fname: Optional[str]) -> Optional[str]:
        candidates = []
        if fname and fname != "<source>":
            candidates.append(str(Path(fname).resolve().parent / name))
        for d in self.include_dirs:
            candidates.append(str(Path(d).resolve() / name))
        for cand in candidates:
            if Path(cand).is_file():
                return str(Path(cand).resolve())
        return None

    # ------------------------------------------------------------------ #
    # Macro expansion
    # ------------------------------------------------------------------ #
    def _expand(self, text: str, fname: str, lineno: int, depth: int = 0) -> str:
        if depth > 64:
            raise ValueError(f"macro expansion too deep (recursive macro?) "
                             f"at {fname}:{lineno}")
        iterations = 0
        while True:
            idx = text.find("`")
            if idx < 0:
                return text
            iterations += 1
            if iterations > _MAX_EXPANSIONS:
                raise ValueError(f"macro expansion limit exceeded "
                                 f"(recursive macro?) at {fname}:{lineno}")
            m = _MACRO_TOKEN.match(text[idx:])
            if m is None:
                raise ValueError(f"invalid macro token at {fname}:{lineno}")
            name = m.group(1)
            after = idx + m.end()
            macro = self.macros.get(name)
            if macro is None:
                raise ValueError(f"undefined Verilog-A macro `{name}` "
                                 f"at {fname}:{lineno}")

            if after < len(text) and text[after] == "(":
                if macro.kind != "function":
                    raise ValueError(
                        f"object-like macro `{name}` used with arguments "
                        f"at {fname}:{lineno}")
                end = self._match_paren(text, after, fname, lineno)
                raw_args = self._split_args(text[after + 1:end]) \
                    if text[after + 1:end].strip() else []
                if len(raw_args) != len(macro.params):
                    raise ValueError(
                        f"macro `{name}` expects {len(macro.params)} "
                        f"argument(s), got {len(raw_args)} "
                        f"at {fname}:{lineno}")
                expanded_args = [
                    self._expand(a, fname, lineno, depth + 1)
                    for a in raw_args
                ]
                body = macro.body
                for param, value in zip(macro.params, expanded_args):
                    body = re.sub(
                        r"(?<![A-Za-z0-9_])" + re.escape(param)
                        + r"(?![A-Za-z0-9_])", value, body)
                text = text[:idx] + body + text[end + 1:]
            else:
                if macro.kind == "function":
                    raise ValueError(
                        f"function-like macro `{name}` invoked without "
                        f"arguments at {fname}:{lineno}")
                text = text[:idx] + macro.body + text[after:]

    @staticmethod
    def _match_paren(text: str, open_idx: int, fname: str,
                     lineno: int) -> int:
        depth = 1
        i = open_idx + 1
        while i < len(text):
            ch = text[i]
            if ch == '"':
                i += 1
                while i < len(text) and text[i] != '"':
                    i += 1
                i += 1
                continue
            if ch == "(":
                depth += 1
            elif ch == ")":
                depth -= 1
                if depth == 0:
                    return i
            i += 1
        raise ValueError(f"unbalanced '(' in macro at {fname}:{lineno}")

    @staticmethod
    def _split_args(raw: str) -> List[str]:
        args = []
        depth = 0
        start = 0
        i = 0
        n = len(raw)
        while i < n:
            ch = raw[i]
            if ch == '"':
                i += 1
                while i < n and raw[i] != '"':
                    i += 1
                i += 1
                continue
            if ch == "(":
                depth += 1
            elif ch == ")":
                depth -= 1
            elif ch == "," and depth == 0:
                args.append(raw[start:i].strip())
                start = i + 1
            i += 1
        args.append(raw[start:].strip())
        return args

    @staticmethod
    def _unclosed_parens(text: str) -> int:
        """Count of unmatched '(' outside string literals in the line."""
        depth = 0
        i = 0
        n = len(text)
        while i < n:
            ch = text[i]
            if ch == '"':
                i += 1
                while i < n and text[i] != '"':
                    i += 1
                i += 1
                continue
            if ch == "(":
                depth += 1
            elif ch == ")":
                depth -= 1
                if depth < 0:
                    return 0
            i += 1
        return depth


def preprocess_source(text: str, *, filename: Optional[str] = None,
                      include_dirs: Sequence[Path] = ()) -> PreprocessedResult:
    """Preprocess Verilog-A source text into plain GMC-readable text."""
    pp = VerilogAPreprocessor(include_dirs)
    pp.run(text, filename)
    return pp.result()