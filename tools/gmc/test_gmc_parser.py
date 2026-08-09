import json
import tempfile
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))

from gmc_compile import build_gsdi_artifact, module_to_dict
from gmc_cpp_emitter import GMCCppEmitter
from gmc_parser import GMCParser
from gmc_preprocessor import preprocess_source


def demo():
    source = """
    module va_probe(p, n);
        inout p, n;
        electrical p, n;
        parameter real Is = 1e-14;
        parameter real Cj = 2p;
        real vd, id;
        analog begin
            vd = V(p, n);
            id = Is * (exp(vd / 0.02585) - 1);
            I(p, n) <+ id + ddt(Cj * vd);
            if (vd > 0)
                I(p, n) <+ 1u;
            else
                I(p, n) <+ -1u;
        end
    endmodule
    """
    module = GMCParser(source).parse_module()
    ir = module_to_dict(module)
    assert ir["schema"] == "gmc-ir-v0"
    assert ir["name"] == "va_probe"
    assert ir["ports"] == ["p", "n"]
    assert [p["name"] for p in ir["parameters"]] == ["Is", "Cj"]
    assert ir["parameters"][1]["default"] == 2e-12
    assert ir["analog"][2]["kind"] == "branch_contribution"
    assert ir["analog"][2]["expr"]["right"]["name"] == "ddt"
    assert ir["analog"][3]["kind"] == "conditional"
    assert ir["analog"][3]["condition"]["op"] == ">"
    header = GMCCppEmitter(module).emit_cpp_header()
    assert "class GmcVaProbeInstance" in header
    assert "eval_model" in header
    assert "result.static_jacobian" in header
    assert "result.dynamic_residual" in header
    assert "result.dynamic_jacobian" in header
    assert "charges[0] +=" in header
    assert "V_anode" not in header
    assert "GsdiModelDescriptor" in header
    assert "terminalCount() const override { return 2; }" in header
    assert "GsdiNoCollapse" in header
    assert "GmcDual<2>" in header
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        source = td / "va_probe.va"
        header_path = td / "va_probe.hpp"
        source.write_text(source_text := """
module va_probe(p, n);
    electrical p, n;
    parameter real G = 2m;
    analog begin
        I(p, n) <+ G * V(p, n);
    end
endmodule
""", encoding="utf-8")
        artifact_module = GMCParser(source_text, source_name=str(source)).parse_module()
        artifact = build_gsdi_artifact(artifact_module, source, header_path)
        assert artifact["schema"] == "gsdi-artifact-v0"
        assert artifact["producer"] == "gmc"
        assert artifact["model_type"] == "va_probe"
        assert artifact["source_sha256"]
        assert artifact["generated_header"].endswith("va_probe.hpp")
    print(json.dumps(ir, indent=2, sort_keys=True))


def internal_nodes():
    source = """
    module psp_like(d, g, s, b, Td, Tint);
        electrical d, g, s, b, Td, Tint;
        terminal 4;
        collapsible Tint with s;
        parameter real Rd = 100;
        parameter real KP = 2m;
        real vch;
        analog begin
            vch = KP * (V(g) - V(s));
            I(d, Td) <+ (V(d) - V(Td)) / Rd;
            I(Td, s) <+ vch;
            I(Tint, s) <+ 0.0;
        end
    endmodule
    """
    module = GMCParser(source).parse_module()
    assert module.terminal_count == 4
    assert module.collapsible_pairs == {"Tint": "s"}

    ir = module_to_dict(module)
    assert ir["terminal_count"] == 4
    assert ir["collapsible_pairs"] == {"Tint": "s"}
    assert ir["ports"] == ["d", "g", "s", "b", "Td", "Tint"]

    emitter = GMCCppEmitter(module)
    terminal_count, roles = emitter.local_node_info()
    assert terminal_count == 4
    assert roles == [("Terminal", -1), ("Terminal", -1), ("Terminal", -1),
                     ("Terminal", -1), ("Internal", -1), ("Collapsible", 2)]

    header = emitter.emit_cpp_header()
    assert "GsdiNodeRole::Internal" in header
    assert "GsdiNodeRole::Collapsible, 2" in header
    assert "terminalCount() const override { return 4; }" in header
    assert "internalNodeCount() const override { return 1; }" in header
    assert "GmcDual<6>" in header

    # Fail closed: collapsible node must merge into a declared external terminal.
    bad_source = source.replace("collapsible Tint with s;", "collapsible Tint with Td;")
    bad_module = GMCParser(bad_source).parse_module()
    try:
        GMCCppEmitter(bad_module).local_node_info()
        raise AssertionError("expected ValueError for non-terminal collapse partner")
    except ValueError:
        pass

    # Fail closed: terminal count cannot exceed the port list.
    bad_source2 = source.replace("terminal 4;", "terminal 7;")
    bad_module2 = GMCParser(bad_source2).parse_module()
    try:
        GMCCppEmitter(bad_module2).local_node_info()
        raise AssertionError("expected ValueError for terminal count > ports")
    except ValueError:
        pass

    print(json.dumps(ir, indent=2, sort_keys=True))


def deep():
    source = """
    module deep_probe(p, n);
        inout p, n;
        electrical p, n;
        parameter real G = 2m;
        parameter real Cj = 3u;
        real v, w;
        analog begin
            v = V(p, n);
            w = exp(v) + ln(v + 1.0) + sqrt(abs(v) + 1.0)
                + tanh(v) + sin(v) + cosh(v) + sinh(v) / 4.0
                + (v > 0.5 ? v ** 2 : 0.5)
                + max(v - 0.5, 0.0) + min(v, 1.5)
                + log10(v * 10.0) + $ln(v + 2.0)
                + limexp(v - 2.0) + $temperature / 300.0 - 1.0
                + pow(v, 2);
            I(p, n) <+ G * w;
            if (v > 0.0)
                I(p, n) <+ 1e-12 * v;
            else
                I(p, n) <+ 0.0;
            I(p, n) <+ ddt(Cj * v);
        end
    endmodule
    """
    module = GMCParser(source).parse_module()
    ir = module_to_dict(module)
    assert ir["name"] == "deep_probe"
    assign = next(s for s in ir["analog"]
                  if s["kind"] == "assignment" and s["name"] == "w")
    # Ternary parses as its own node kind.
    ternary = next(e for e in walk(assign["expr"]) if e["kind"] == "ternary")
    assert ternary["condition"]["op"] == ">"
    assert ternary["then"]["op"] == "**"
    # $--prefixed and positional intrinsics land as calls/idents in the IR.
    calls = {e["name"] for e in walk(assign["expr"])
             if e["kind"] == "call"}
    identifiers = {e["name"] for e in walk(assign["expr"])
                   if e["kind"] == "identifier"}
    assert "$ln" in calls and "$temperature" in identifiers
    assert {"exp", "ln", "sqrt", "abs", "tanh", "sin", "cosh", "sinh",
            "log10", "max", "min", "pow", "limexp"} <= calls

    header = GMCCppEmitter(module).emit_cpp_header()
    for fragment in ("gmc_pow2(", "gmc_log10(", "gmc_sin(", "gmc_tanh(",
                     "gmc_min(", "gmc_max(", "temperature_", "vt_",
                     "8.617333262145179e-5"):
        assert fragment in header, fragment

    # Fail closed: unknown functions are rejected instead of emitting raw calls.
    bad_source = source.replace("+ pow(v, 2);", "+ foo(v);")
    bad_module = GMCParser(bad_source).parse_module()
    try:
        GMCCppEmitter(bad_module).emit_cpp_header()
        raise AssertionError("expected ValueError for unsupported function")
    except ValueError as exc:
        assert "unsupported function 'foo'" in str(exc)

    # Fail closed: ddt() inside an assignment cannot preserve its charge.
    bad_ddt = source.replace("I(p, n) <+ ddt(Cj * v);", "w = ddt(Cj * v);")
    bad_ddt_module = GMCParser(bad_ddt).parse_module()
    try:
        GMCCppEmitter(bad_ddt_module).emit_cpp_header()
        raise AssertionError("expected ValueError for ddt in assignment")
    except ValueError as exc:
        assert "ddt() inside an assignment" in str(exc)

    print(json.dumps(ir, indent=2, sort_keys=True))


def functions():
    source = """
    module fn_probe(p, n);
        inout p, n;
        electrical p, n;
        parameter real G = 2m;
        analog function real gclamp;
            input x, lo, hi;
            real x, lo, hi;
            begin
                gclamp = min(max(x, lo), hi);
            end
        endfunction
        analog function real gsoft;
            input x;
            real x;
            begin
                gsoft = 1.0 + 0.5 * gclamp(x, 0.0, 1.0);
            end
        endfunction
        real w;
        analog begin
            w = gsoft(0.75);
            I(p, n) <+ G * w;
        end
    endmodule
    """
    module = GMCParser(source).parse_module()
    assert list(module.functions) == ["gclamp", "gsoft"]
    gclamp = module.functions["gclamp"]
    assert gclamp.return_type == "real"
    assert [(a.name, a.vtype) for a in gclamp.args] == [
        ("x", "real"), ("lo", "real"), ("hi", "real")]
    assert len(gclamp.body) == 1

    ir = module_to_dict(module)
    assert [f["name"] for f in ir["functions"]] == ["gclamp", "gsoft"]
    assert ir["functions"][0]["args"][0]["name"] == "x"

    header = GMCCppEmitter(module).emit_cpp_header()
    # callee-first emission: gclamp lambda appears before gsoft, and gsoft
    # calls gclamp by name.
    assert header.index("auto gclamp =") < header.index("auto gsoft =")
    assert "auto gclamp = [&](const auto& x, const auto& lo, const auto& hi)" in header
    assert "Scalar gclamp_ret = 0.0;" in header
    assert "return gclamp_ret;" in header
    assert "gclamp(x, 0.0, 1.0)" in header
    assert "gsoft((v - 0.25))" not in header  # call sites are in the body

    # Fail closed: recursion (direct or indirect) is rejected.
    recursive_source = source.replace(
        "gsoft = 1.0 + 0.5 * gclamp(x, 0.0, 1.0);",
        "gsoft = 1.0 + 0.5 * gsoft(x);")
    recursive_module = GMCParser(recursive_source).parse_module()
    try:
        GMCCppEmitter(recursive_module).emit_cpp_header()
        raise AssertionError("expected ValueError for recursive analog function")
    except ValueError as exc:
        assert "recursive analog function" in str(exc)

    # Fail closed: a function referencing a module-scoped variable (not passed
    # as an argument) is out of scope.
    out_of_scope_source = source.replace(
        "gsoft = 1.0 + 0.5 * gclamp(x, 0.0, 1.0);",
        "gsoft = 1.0 + 0.5 * G;")
    out_of_scope_module = GMCParser(out_of_scope_source).parse_module()
    try:
        GMCCppEmitter(out_of_scope_module).emit_cpp_header()
        raise AssertionError("expected ValueError for out-of-scope identifier")
    except ValueError as exc:
        assert "out-of-scope identifier 'G'" in str(exc)

    print(json.dumps(ir, indent=2, sort_keys=True))


def case_statement():
    source = """
    module sel_probe(p, n);
        inout p, n;
        electrical p, n;
        parameter real G = 2m;
        integer sel;
        real v, w;
        analog begin
            sel = 2;
            v = V(p, n);
            w = sel * v;
            I(p, n) <+ G * w;
            case (sel)
                1: I(p, n) <+ 0.0;
                2, 3: I(p, n) <+ 1e-9 * v;
                default: I(p, n) <+ 2e-9 * v;
            endcase
        end
    endmodule
    """
    module = GMCParser(source).parse_module()
    # integer variables parse as module variables (type carried).
    assert module.variables["sel"].vtype == "integer"
    assert module.variables["v"].vtype == "real"

    ir = module_to_dict(module)
    case = next(s for s in ir["analog"] if s["kind"] == "case")
    assert case["expr"]["name"] == "sel"
    assert len(case["items"]) == 3
    assert [len(i["values"]) for i in case["items"]] == [1, 2, 0]
    assert case["items"][1]["values"][1]["value"] == 3

    header = GMCCppEmitter(module).emit_cpp_header()
    assert "if (gmc_nz(sel == 1.0)) {" in header
    assert "else if (gmc_nz(sel == 2.0) || gmc_nz(sel == 3.0)) {" in header
    assert "else {" in header
    assert "Scalar sel = 0.0;" in header
    # case in a function body is out of scope if it reads module state -> the
    # non-argument identifier is rejected for functions, but module-body case
    # emits fine (verified above).

    # Fail closed: a case item calling an unknown function is rejected.
    bad_source = source.replace("I(p, n) <+ 0.0;", "I(p, n) <+ foo(v);")
    bad_module = GMCParser(bad_source).parse_module()
    try:
        GMCCppEmitter(bad_module).emit_cpp_header()
        raise AssertionError("expected ValueError for unsupported function in case")
    except ValueError as exc:
        assert "unsupported function 'foo'" in str(exc)

    print(json.dumps(ir, indent=2, sort_keys=True))


def preprocessor():
    # Object-like macro substitution.
    out = preprocess_source("`define SCALE 2.5\nreal x;\nx = `SCALE * 3;\n").text
    assert out.count("`") == 0
    assert "x = 2.5 * 3;" in out

    # Macros referencing other macros expand lazily at the use site.
    out = preprocess_source("""
`define A 2
`define B (`A + 3)
y = `B * `A;
""").text
    assert "y = (2 + 3) * 2;" in out

    # Multi-line macro invocation without backslash continuations (BSIM
    # statement style: args span physical lines).
    out = preprocess_source(
        "`define DEVAL(nvtm, ijth, satc, xexpbv, vjm) \\\n"
        "   vjm = nvtm * log(ijth / satc + xexpbv);\n"
        "`DEVAL(Nvtm, BSIM3ijth, BSIM3SourceSatCurrent,\n"
        "       BSIM3XExpBVS, BSIM3vjsmFwd)\n").text
    assert "BSIM3vjsmFwd = Nvtm * log(BSIM3ijth / BSIM3SourceSatCurrent + BSIM3XExpBVS);" in out

    # String arguments containing commas and empty strings (VBIC style).
    out = preprocess_source(
        "`define MPRnb(nam, def, uni, des) parameter real nam=def;\n"
        "`MPRnb(xii, 3.0, \"\", \"temperature exponent of ibei, ibci, ibeip\")\n"
        "`MPRnb(eaie, 1.12, \"V\", \"activation energy for ibei\")\n").text
    assert "parameter real xii=3.0;" in out
    assert "parameter real eaie=1.12;" in out

    # Function-like macro with multiple args, nested parens and nested calls.
    out = preprocess_source(
        "`define MAND(a, b, c) (a * (b + c))\n"
        "z = `MAND(g(1, 2), h(3), 4);\n").text
    assert "z = (g(1, 2) * (h(3) + 4));" in out

    # Multiline `define via trailing-backslash continuations.
    out = preprocess_source(
        "`define DIO(nvtm, isb, xexpbv, vjm)   \\\n"
        "   if (isb != 0.0)                    \\\n"
        "      vjm = nvtm * log(xexpbv);       \\\n"
        "   else                               \\\n"
        "      vjm = 0.0;\n"
        "`DIO(2.0, 1.0, 3.0, t)\n").text
    expanded = next(ln for ln in out.splitlines() if "t = 2.0" in ln)
    assert " ".join(expanded.split()) == \
        "if (1.0 != 0.0) t = 2.0 * log(3.0); else t = 0.0;"

    # ifdef/else/endif: defined branch kept, undefined branch dropped.
    out = preprocess_source("""
`define ENABLED
`ifdef ENABLED
keep_a = 1;
`else
drop_a = 1;
`endif
`ifdef MISSING
drop_b = 1;
`else
keep_b = 2;
`endif
""").text
    assert "keep_a = 1;" in out and "drop_a" not in out
    assert "keep_b = 2;" in out and "drop_b" not in out

    # Nested conditionals (including inside an inactive branch).
    out = preprocess_source("""
`ifdef OUTER
`ifdef INNER
inner_kept = 1;
`endif
outer_kept = 1;
`endif
""").text
    assert out.strip() == "" or "outer_kept" not in out

    # ifndef and undef.
    out = preprocess_source(
        "`ifndef FIRST\nfirst = 1;\n`endif\n"
        "`define FLAG 1\n`undef FLAG\n"
        "`ifdef FLAG\nbad = 1;\n`endif\n").text
    assert "first = 1;" in out and "bad" not in out

    # Empty-body macro is defined (BSIM style `define NOISE).
    out = preprocess_source(
        "`define NOISE\n`ifdef NOISE\nn = 1;\n`endif\n"
        "`ifdef EMPTY_USE\nskip = 1;\n`endif\n").text
    assert "n = 1;" in out and "skip" not in out

    # Macros inside comments and disabled regions are inert.
    out = preprocess_source("""
// `define INSIDE_COMMENT 1
/* `define INSIDE_BLOCK 1
   `ifdef NEVER
*/
`ifdef NOT_SET
`define SKIP_ME 1
`endif
v = 3;
""").text
    assert "v = 3;" in out and "SKIP_ME" not in out

    # Fail closed: undefined macro.
    try:
        preprocess_source("x = `NOT_DEFINED;\n")
        raise AssertionError("expected ValueError for undefined macro")
    except ValueError as exc:
        assert "NOT_DEFINED" in str(exc)

    # Fail closed: unterminated / unbalanced conditionals.
    try:
        preprocess_source("`ifdef A\nx = 1;\n")
        raise AssertionError("expected ValueError for unterminated ifdef")
    except ValueError as exc:
        assert "unterminated" in str(exc).lower()
    try:
        preprocess_source("`else\n")
        raise AssertionError("expected ValueError for stray else")
    except ValueError as exc:
        assert "else" in str(exc)
    try:
        preprocess_source("`endif\n")
        raise AssertionError("expected ValueError for stray endif")
    except ValueError as exc:
        assert "endif" in str(exc)

    # Fail closed: function-like macro used without arguments.
    try:
        preprocess_source("`define F(x) (x + 1)\nq = `F;\n")
        raise AssertionError("expected ValueError for bare function-like use")
    except ValueError as exc:
        assert "without arguments" in str(exc)

    # Fail closed: wrong argument count.
    try:
        preprocess_source("`define F2(a, b) (a + b)\nq = `F2(1);\n")
        raise AssertionError("expected ValueError for arg count mismatch")
    except ValueError as exc:
        assert "expects 2" in str(exc)

    # Includes: standard .vams guard pattern, relative + include-dir search.
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        consts = td / "constants.vams"
        consts.write_text(
            "`ifdef CONSTANTS_VAMS\n`else\n"
            "`define CONSTANTS_VAMS 1\n"
            "`define P_Q 1.602176462e-19\n"
            "`define P_CELSIUS0 273.15\n"
            "`endif\n", encoding="utf-8")
        include_dir = td / "inc"
        include_dir.mkdir()
        (include_dir / "extra.vams").write_text(
            "`define EXTRA 7.0\n", encoding="utf-8")
        main = td / "model.va"
        main.write_text(
            "`include \"constants.vams\"\n"
            "`include \"constants.vams\"\n"
            "`include \"extra.vams\"\n"
            "k = `P_Q / `P_CELSIUS0 + `EXTRA;\n", encoding="utf-8")
        result = preprocess_source(
            main.read_text(encoding="utf-8"), filename=str(main),
            include_dirs=[str(include_dir)])
        assert "k = 1.602176462e-19 / 273.15 + 7.0;" in result.text
        assert len([f for f in result.files
                    if f.endswith("constants.vams")]) == 1

        # Missing include file is a hard error.
        try:
            preprocess_source("`include \"nope.vams\"\n",
                              filename=str(main), include_dirs=[str(include_dir)])
            raise AssertionError("expected ValueError for missing include")
        except ValueError as exc:
            assert "not found" in str(exc)

        # Cyclic include is a hard error.
        a = td / "a.va"
        b = td / "b.va"
        a.write_text("`include \"b.va\"\n", encoding="utf-8")
        b.write_text("`include \"a.va\"\n", encoding="utf-8")
        try:
            preprocess_source(a.read_text(encoding="utf-8"),
                              filename=str(a))
            raise AssertionError("expected ValueError for cyclic include")
        except ValueError as exc:
            assert "recursive" in str(exc)

        # End to end: include artifacts with nature/discipline sections parse
        # into a module; macros expand into IR and the emitted header.
        disc = td / "disciplines.vams"
        disc.write_text(
            "`ifdef DISCIPLINES_VAMS\n`else\n"
            "`define DISCIPLINES_VAMS 1\n"
            "nature Voltage;\n  units = \"V\";\n  access = V;\nendnature\n"
            "discipline electrical;\n  potential Voltage;\n  flow Current;\n"
            "enddiscipline\n"
            "`endif\n", encoding="utf-8")
        model = td / "macro_model.va"
        model.write_text(
            "`include \"disciplines.vams\"\n"
            "`define GM_SCALE 1.0e-9\n"
            "module macro_probe(p, n);\n"
            "  inout p, n;\n  electrical p, n;\n"
            "  parameter real G = 2m;\n  real v;\n"
            "  analog begin\n"
            "    v = V(p, n);\n"
            "`ifdef GM_EXTRA\n"
            "    v = v + `GM_EXTRA;\n"
            "`endif\n"
            "    I(p, n) <+ G * v + `GM_SCALE * v;\n"
            "  end\n"
            "endmodule\n", encoding="utf-8")
        module = GMCParser(model.read_text(encoding="utf-8"),
                           source_name=str(model)).parse_module()
        ir = module_to_dict(module)
        contrib = next(s for s in ir["analog"]
                       if s["kind"] == "branch_contribution")
        assert "GM_SCALE" not in json.dumps(contrib)
        header = GMCCppEmitter(module).emit_cpp_header()
        assert "1e-09 * v" in header
        assert "GmcDual<2>" in header

    print(json.dumps({"preprocessor": "ok"}, indent=2, sort_keys=True))


def real_models():
    # Constructs required to parse real BSIM3/BSIM4 sources: logical
    # operators, no-op system tasks ($strobe/$discontinuity/...), for/while
    # loops and @(initial_step) event blocks, and noise contributions.
    source = """
    module bsim_like(d, g, s, b);
        inout d, g, s, b;
        electrical d, g, s, b;
        parameter real vth0 = 0.5;
        real vgs, ids, niter, toxpf;
        analog begin
            vgs = V(g, s);
            if ((vgs <= 0.0) && (vgs > -1.0)) begin
                ids = 0.0;
            end else if ((vgs > 0.0) || (vgs < -2.0)) begin
                ids = vgs * vgs;
            end else begin
                ids = !ids;
            end
            $discontinuity(4.0);
            $strobe("vgs=%g temp exp of ..., ...", vgs, ids);
            niter = 0;
            while ((niter <= 4) && (abs(toxpf - 300.0) > 1e-12)) begin
                niter = niter + 1;
            end
            for (i = 0; i < 4; i = i + 1) begin
                ids = ids + i;
            end
            @(initial_step) begin
                toxpf = 0.0;
            end
            I(d, s) <+ ids;
        end
    endmodule
    """
    module = GMCParser(source).parse_module()

    ir = module_to_dict(module)
    kinds = [s["kind"] for s in ir["analog"]]
    assert kinds == ["assignment", "conditional", "assignment",
                     "while", "for", "event", "branch_contribution"]
    # else-if is the first conditional's else_body holding a nested
    # conditional whose condition is the || branch.
    else_if = ir["analog"][1]["else"][0]["condition"]
    assert else_if["kind"] == "binary" and else_if["op"] == "||"
    # No-op system tasks must not appear in the IR at all.
    assert "$strobe" not in json.dumps(ir)
    assert "$discontinuity" not in json.dumps(ir)
    # Logical && / || land as binary nodes inside the conditions.
    cond = ir["analog"][1]["condition"]
    assert cond["kind"] == "binary" and cond["op"] == "&&"
    assert cond["right"]["kind"] == "binary"
    wh = ir["analog"][3]
    assert wh["condition"]["op"] == "&&"
    assert wh["body"][0]["kind"] == "assignment"
    fo = ir["analog"][4]
    assert fo["init"]["name"] == "i" and fo["condition"]["op"] == "<"
    assert fo["increment"]["expr"]["op"] == "+"
    ev = ir["analog"][5]
    assert ev["event"] == "initial_step"

    header = GMCCppEmitter(module).emit_cpp_header()
    assert "gmc_nz" in header
    assert "while (gmc_nz((gmc_nz((niter <= 4.0)) &&" in header
    assert "for (i = 0.0; gmc_nz((i < 4.0)); i = (i + 1.0)) {" in header
    assert "@(initial_step) event block" in header

    # Bitwise operators parse but fail closed at emission (no continuous
    # analog semantics in GMC); modulo emits via gmc_fmod.
    for op, emit_fragment in (("&", "bitwise operator '&'"),):
        bad = GMCParser("""
module bad_op; inout p, n; electrical p, n;
real w;
analog begin
    w = 5 %s 2;
    I(p, n) <+ w;
end
endmodule
""" % op).parse_module()
        try:
            GMCCppEmitter(bad).emit_cpp_header()
            raise AssertionError("expected ValueError for '%s'" % op)
        except ValueError as exc:
            assert emit_fragment in str(exc)

    # Noise contributions parse and emit through the GSDI noise-source path.
    noise_module = GMCParser("""
module noisy; inout p, n; electrical p, n;
analog begin
    I(p, n) <+ white_noise(1e-9, "wgn");
end
endmodule
""").parse_module()
    header = GMCCppEmitter(noise_module).emit_cpp_header()
    assert "d.supports_noise = true" in header
    assert "add_noise(0, 1, 1e-09, \"wgn\")" in header

    # Internal nodes declared via "electrical" but not in the module port
    # list (BSIM's di/si pattern) become hidden nodes after the terminals.
    hidden_source = """
module bsim_hidden(d, g, s, b);
    inout d, g, s, b;
    electrical d, g, s, b, di, si;
    parameter real rd = 100;
    analog begin
        I(di, si) <+ (V(d) - V(s)) / rd;
        I(d, di) <+ (V(d) - V(di)) / rd;
    end
endmodule
"""
    hidden_module = GMCParser(hidden_source).parse_module()
    assert hidden_module.ports == ["d", "g", "s", "b", "di", "si"]
    assert hidden_module.terminal_count == 4
    hidden_header = GMCCppEmitter(hidden_module).emit_cpp_header()
    assert "terminalCount() const override { return 4; }" in hidden_header
    assert "internalNodeCount() const override { return 2; }" in hidden_header
    assert "GmcDual<6>" in hidden_header
    assert "currents[4] +=" in hidden_header
    assert "currents[5] -=" in hidden_header

    # Bare begin/end grouping blocks (BSIM "if (x) begin begin ... end end")
    # must not truncate the enclosing analog block at a stray 'end', and the
    # grouped statements must survive into the emitted header.
    begin_source = """
module bsim_begins(d, g, s, b);
    electrical d, g, s, b;
    real ids;
    analog begin
        if (1.0) begin
            begin
                $strobe("warning");
            end
            ids = 1.0;
        end
        I(d, s) <+ ids;
    end
endmodule
"""
    begin_module = GMCParser(begin_source).parse_module()
    assert len(begin_module.analog_body) == 2
    assert begin_module.analog_body[0].then_body[0].__class__.__name__ == "GroupNode"
    assert begin_module.analog_body[0].then_body[1].__class__.__name__ == "AssignmentNode"
    assert begin_module.analog_body[1].__class__.__name__ == "BranchContributionNode"
    begin_header = GMCCppEmitter(begin_module).emit_cpp_header()
    assert "currents[0] += ids;" in begin_header
    assert "currents[2] -= ids;" in begin_header

    # Labeled grouping blocks ("begin : label ... end", VBIC's
    # evaluateStatic/loadStatic) carry no control flow; their statements must
    # be flattened into the emitted body, including branch contributions.
    label_source = """
module bsim_labeled(d, g, s, b);
    electrical d, g, s, b;
    parameter real gm = 1m;
    real ids;
    analog begin
        ids = gm * (V(g) - V(s));
        begin : loadStatic
            I(d, s) <+ ids;
        end
    end
endmodule
"""
    label_module = GMCParser(label_source).parse_module()
    assert len(label_module.analog_body) == 2
    assert label_module.analog_body[1].__class__.__name__ == "GroupNode"
    assert label_module.analog_body[1].label == "loadStatic"
    assert label_module.analog_body[1].body[0].__class__.__name__ == "BranchContributionNode"
    label_header = GMCCppEmitter(label_module).emit_cpp_header()
    assert "currents[0] += ids;" in label_header
    assert "currents[2] -= ids;" in label_header

    # Return-type-less analog functions (BSIM3 "analog function limvgs;").
    fn_source = """
module fn_no_type; inout p, n; electrical p, n;
real w;
analog function limvgs;
    input x, y;
    real x, y, limited;
    begin
        if (x > 0.0) begin
            limited = x + y;
        end else begin
            limited = -x;
        end
        limvgs = limited;
    end
endfunction
analog begin
    w = limvgs(1.0, 2.0);
    I(p, n) <+ w;
end
endmodule
"""
    fn_module = GMCParser(fn_source).parse_module()
    assert fn_module.functions["limvgs"].return_type == "real"
    assert [a.name for a in fn_module.functions["limvgs"].args] == ["x", "y"]
    header = GMCCppEmitter(fn_module).emit_cpp_header()
    assert "limvgs_ret" in header

    # Non-literal parameter defaults (BSIM3 "parameter real dsub=drout;")
    # parse and appear in the IR, but fail closed at emission.
    param_source = """
module param_expr; inout p, n; electrical p, n;
parameter real drout = 0.56;
parameter real dsub = drout;
analog begin
    I(p, n) <+ dsub;
end
endmodule
"""
    param_module = GMCParser(param_source).parse_module()
    assert param_module.parameters["drout"].default_expr is None
    assert param_module.parameters["dsub"].default_expr is not None
    param_ir = module_to_dict(param_module)
    dsub = next(x for x in param_ir["parameters"] if x["name"] == "dsub")
    assert dsub["default_expr"]["name"] == "drout"
    param_header = GMCCppEmitter(param_module).emit_cpp_header()
    assert "instance->dsub_ = instance->drout_;" in param_header

    print(json.dumps(ir, indent=2, sort_keys=True))


def walk(expr):
    yield expr
    if expr["kind"] == "binary":
        yield from walk(expr["left"])
        yield from walk(expr["right"])
    elif expr["kind"] == "unary":
        yield from walk(expr["expr"])
    elif expr["kind"] == "call":
        for arg in expr["args"]:
            yield from walk(arg)
    elif expr["kind"] == "ternary":
        yield from walk(expr["condition"])
        yield from walk(expr["then"])
        yield from walk(expr["else"])


if __name__ == "__main__":
    demo()
    internal_nodes()
    deep()
    functions()
    case_statement()
    preprocessor()
    real_models()
