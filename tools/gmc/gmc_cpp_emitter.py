"""GMC C++ emitter for the restricted BSIM/PSP-required Verilog-A subset.

Generated code is GSDI-native: a GsdiModelDescriptor + GsdiModel /
GsdiInstance pair, evaluated with GmcDual<N> dual-number arithmetic where N is
the number of local nodes (ports for now; internal/collapsible nodes are
reserved for the next emitter slice).
"""

import json
import re
from gmc_ast import (
    AssignmentNode, BinaryOpNode, BranchContributionNode, BranchProbeNode,
    CaseNode, ConditionalNode, EventNode, ForNode, FunctionCallNode,
    GroupNode, IdentifierNode, ModuleNode, NumberNode, TernaryNode, UnaryOpNode, WhileNode
)

# Verilog-A intrinsic -> generated evaluation lambda in the emitted header.
# The lambdas dispatch between plain double evaluation and GmcDual<N> AD.
MATH_INTRINSICS = {
    "exp": "gmc_exp",
    "ln": "gmc_log",
    "log": "gmc_log",
    "log10": "gmc_log10",
    "sqrt": "gmc_sqrt",
    "abs": "gmc_abs",
    "tanh": "gmc_tanh",
    "sin": "gmc_sin",
    "cos": "gmc_cos",
    "tan": "gmc_tan",
    "atan": "gmc_atan",
    "sinh": "gmc_sinh",
    "cosh": "gmc_cosh",
    "min": "gmc_min",
    "max": "gmc_max",
    "pow": "gmc_pow2",
}

class GMCCppEmitter:
    def __init__(self, module: ModuleNode):
        self.module = module
        self.current_function = None
        self.current_function_arg_names = set()

    def cpp_name(self, name: str) -> str:
        clean = re.sub(r'[^a-zA-Z0-9_]', '_', name)
        if clean and clean[0].isdigit():
            clean = "_" + clean
        return clean

    def class_stem(self) -> str:
        return "".join(part.capitalize() for part in self.cpp_name(self.module.name).split("_") if part)

    def fold_constants(self, expr):
        if isinstance(expr, UnaryOpNode):
            sub = self.fold_constants(expr.expr)
            if isinstance(sub, NumberNode):
                if expr.op == '-':
                    return NumberNode(-sub.value)
                if expr.op == '+':
                    return NumberNode(+sub.value)
            return UnaryOpNode(expr.op, sub)
        if isinstance(expr, BinaryOpNode):
            left = self.fold_constants(expr.left)
            right = self.fold_constants(expr.right)
            if isinstance(left, NumberNode) and isinstance(right, NumberNode):
                if expr.op == '+':
                    return NumberNode(left.value + right.value)
                if expr.op == '-':
                    return NumberNode(left.value - right.value)
                if expr.op == '*':
                    return NumberNode(left.value * right.value)
                if expr.op == '/' and right.value != 0:
                    return NumberNode(left.value / right.value)
                if expr.op == '**':
                    return NumberNode(left.value ** right.value)
            if expr.op == '*':
                if isinstance(left, NumberNode):
                    if left.value == 0:
                        return NumberNode(0.0)
                    if left.value == 1:
                        return right
                if isinstance(right, NumberNode):
                    if right.value == 0:
                        return NumberNode(0.0)
                    if right.value == 1:
                        return left
            if expr.op == '+':
                if isinstance(left, NumberNode) and left.value == 0:
                    return right
                if isinstance(right, NumberNode) and right.value == 0:
                    return left
            if expr.op == '-':
                if isinstance(right, NumberNode) and right.value == 0:
                    return left
            return BinaryOpNode(expr.op, left, right)
        if isinstance(expr, TernaryNode):
            cond = self.fold_constants(expr.condition)
            then_e = self.fold_constants(expr.then_expr)
            else_e = self.fold_constants(expr.else_expr)
            if isinstance(cond, NumberNode):
                return then_e if cond.value != 0 else else_e
            return TernaryNode(cond, then_e, else_e)
        return expr

    def expr_cpp(self, expr) -> str:
        expr = self.fold_constants(expr)
        if isinstance(expr, NumberNode):
            return repr(expr.value)
        if isinstance(expr, IdentifierNode):
            if expr.name == "$temperature":
                return "temperature_"
            if expr.name == "$vt":
                return "vt_"
            if self.current_function is not None:
                if expr.name == self.current_function.name:
                    return f"{self.cpp_name(expr.name)}_ret"
                if expr.name in self.current_function_arg_names:
                    return self.cpp_name(expr.name)
                if any(l.name == expr.name for l in self.module.functions[self.current_function.name].locals):
                    return self.cpp_name(expr.name)
                raise ValueError(
                    f"module '{self.module.name}': function '{self.current_function.name}' "
                    f"references out-of-scope identifier '{expr.name}' (pass it as an "
                    f"argument or declare it locally)")
            if expr.name in self.module.parameters:
                return f"{self.cpp_name(expr.name)}_"
            return self.cpp_name(expr.name)
        if isinstance(expr, UnaryOpNode):
            if expr.op == '!':
                return f"(!gmc_nz({self.expr_cpp(expr.expr)}))"
            return f"(-{self.expr_cpp(expr.expr)})"
        if isinstance(expr, BinaryOpNode):
            if expr.op in ('&&', '||'):
                return (f"(gmc_nz({self.expr_cpp(expr.left)}) {expr.op} "
                        f"gmc_nz({self.expr_cpp(expr.right)}))")
            if expr.op in ('&', '|'):
                raise ValueError(
                    f"module '{self.module.name}': bitwise operator '{expr.op}' "
                    f"is not supported in GMC emission")
            if expr.op == '%':
                return f"gmc_fmod({self.expr_cpp(expr.left)}, {self.expr_cpp(expr.right)})"
            if expr.op == "**":
                return f"gmc_pow2({self.expr_cpp(expr.left)}, {self.expr_cpp(expr.right)})"
            return f"({self.expr_cpp(expr.left)} {expr.op} {self.expr_cpp(expr.right)})"
        if isinstance(expr, TernaryNode):
            return (f"((gmc_nz({self.expr_cpp(expr.condition)})) ? ({self.expr_cpp(expr.then_expr)})"
                    f" : ({self.expr_cpp(expr.else_expr)}))")
        if isinstance(expr, BranchProbeNode):
            if expr.quantity != "V":
                return "0.0"
            n1 = self.module.ports.index(expr.node1) if expr.node1 in self.module.ports else -1
            if expr.node2 and expr.node2 in self.module.ports:
                n2 = self.module.ports.index(expr.node2)
                return f"(vn[{n1}] - vn[{n2}])"
            return f"vn[{n1}]"
        if isinstance(expr, FunctionCallNode):
            name = expr.func_name
            if name == "ddt":
                return "0.0"
            if name in ("limexp", "$limexp"):
                args = ", ".join(self.expr_cpp(arg) for arg in expr.args)
                return f"limexp({args})"
            if name == "$simparam":
                # Simulator parameter query: the generated model has no
                # simulator context, so a two-argument call resolves to its
                # default (the standard Verilog-A behavior when the parameter
                # is not set) and a one-argument call to 0.0.
                if len(expr.args) == 2:
                    return f"({self.expr_cpp(expr.args[1])})"
                if len(expr.args) == 1:
                    return "0.0"
                raise ValueError(
                    f"module '{self.module.name}': $simparam requires 1 or 2 "
                    f"arguments, got {len(expr.args)}")
            if name == "analysis":
                # Analysis-kind query resolved against the instance's
                # analysis_ member (default "dc").
                if len(expr.args) != 1 or not isinstance(expr.args[0], IdentifierNode) \
                        or not expr.args[0].name.startswith('"'):
                    raise ValueError(
                        f"module '{self.module.name}': analysis() requires a "
                        f"single string argument")
                return f"(analysis_ == {expr.args[0].name})"
            if name == "$limit":
                # $limit(x, ...) clamps an intermediate quantity to keep the
                # evaluation inside the physically valid domain; all other
                # arguments are domain bounds/hints, so only x survives.
                if len(expr.args) < 1:
                    raise ValueError(
                        f"module '{self.module.name}': $limit requires at least "
                        f"one argument")
                return self.expr_cpp(expr.args[0])
            if name == "$param_given":
                # $param_given(p) reports whether the caller supplied p
                # explicitly (vs. it falling back to its default). Static
                # codegen can only know this at card-binding time, so each
                # parameter carries a given_<name> flag on the instance.
                if len(expr.args) != 1:
                    raise ValueError(
                        f"module '{self.module.name}': $param_given requires "
                        f"exactly one argument")
                arg = expr.args[0]
                if isinstance(arg, IdentifierNode):
                    pname = arg.name
                elif isinstance(arg, LiteralNode) and isinstance(arg.value, str):
                    pname = arg.value
                else:
                    raise ValueError(
                        f"module '{self.module.name}': $param_given argument "
                        f"must name a parameter")
                if pname not in self.module.parameters:
                    raise ValueError(
                        f"module '{self.module.name}': $param_given('{pname}') "
                        f"does not name a declared parameter")
                return f"given_{self.cpp_name(pname)}"
            args = ", ".join(self.expr_cpp(arg) for arg in expr.args)
            raw = name[1:] if name.startswith('$') else name
            if raw == "temperature":
                return "temperature_"
            if raw == "vt":
                return "vt_"
            if raw in self.module.functions:
                return f"{self.cpp_name(raw)}({args})"
            if raw in ('white_noise', 'flicker_noise'):
                # Noise calls are emitted by the separate request.noise pass;
                # they contribute no DC/transient branch current.
                return "0.0 /* noise handled separately */"
            if raw == 'Temp':
                # Temp(branch) probes the temperature of a thermal-discipline
                # (self-heating) branch. GMC does not MNA-solve thermal nodes:
                # with the thermal outputs (Pwr contributions) collapsed the
                # consistent value is a zero temperature rise.
                return "0.0 /* thermal node unsolved */"
            target = MATH_INTRINSICS.get(raw)
            if target is None:
                raise ValueError(
                    f"module '{self.module.name}': unsupported function '{name}'")
            return f"{target}({args})"
        raise TypeError(f"Unsupported expression node: {type(expr).__name__}")

    def expr_has_noise(self, expr) -> bool:
        if isinstance(expr, FunctionCallNode):
            raw = expr.func_name[1:] if expr.func_name.startswith('$') else expr.func_name
            return raw in ('white_noise', 'flicker_noise') or any(
                self.expr_has_noise(arg) for arg in expr.args)
        if isinstance(expr, UnaryOpNode):
            return self.expr_has_noise(expr.expr)
        if isinstance(expr, BinaryOpNode):
            return self.expr_has_noise(expr.left) or self.expr_has_noise(expr.right)
        if isinstance(expr, TernaryNode):
            return (self.expr_has_noise(expr.condition)
                    or self.expr_has_noise(expr.then_expr)
                    or self.expr_has_noise(expr.else_expr))
        return False

    def has_noise(self, statements) -> bool:
        for stmt in statements:
            if isinstance(stmt, ConditionalNode):
                if self.has_noise(stmt.then_body) or self.has_noise(stmt.else_body):
                    return True
            elif isinstance(stmt, CaseNode):
                if self.has_noise([child for item in stmt.items for child in item.body]):
                    return True
            elif isinstance(stmt, (ForNode, WhileNode, EventNode, GroupNode)):
                if self.has_noise(stmt.body):
                    return True
            elif isinstance(stmt, BranchContributionNode):
                if self.expr_has_noise(stmt.expr):
                    return True
            elif isinstance(stmt, AssignmentNode):
                if self.expr_has_noise(stmt.expr):
                    return True
        return False

    def noise_name_cpp(self, expr, fallback: str) -> str:
        if isinstance(expr, IdentifierNode) and expr.name.startswith('"') and expr.name.endswith('"'):
            return expr.name
        return json.dumps(fallback)

    def emit_noise_terms(self, n1: int, n2: int, expr, indent: str) -> str:
        if isinstance(expr, FunctionCallNode):
            raw = expr.func_name[1:] if expr.func_name.startswith('$') else expr.func_name
            if raw == 'white_noise' and expr.args:
                name = self.noise_name_cpp(expr.args[1], "white_noise") if len(expr.args) > 1 else '"white_noise"'
                return f"{indent}add_noise({n1}, {n2}, {self.expr_cpp(expr.args[0])}, {name});\n"
            if raw == 'flicker_noise' and expr.args:
                exponent = self.expr_cpp(expr.args[1]) if len(expr.args) > 1 else "1.0"
                name_arg = expr.args[2] if len(expr.args) > 2 else None
                name = self.noise_name_cpp(name_arg, "flicker_noise") if name_arg else '"flicker_noise"'
                return (
                    f"{indent}add_noise({n1}, {n2}, ({self.expr_cpp(expr.args[0])}) / "
                    f"std::pow(std::max(std::abs(request.omega) / 6.28318530717958647692, 1.0), "
                    f"{exponent}), {name});\n"
                )
        if isinstance(expr, BinaryOpNode) and expr.op in ('+', '-'):
            code = self.emit_noise_terms(n1, n2, expr.left, indent)
            code += self.emit_noise_terms(n1, n2, expr.right, indent)
            return code
        if isinstance(expr, UnaryOpNode):
            return self.emit_noise_terms(n1, n2, expr.expr, indent)
        if isinstance(expr, TernaryNode):
            code = f"{indent}if (gmc_nz({self.expr_cpp(expr.condition)})) {{\n"
            code += self.emit_noise_terms(n1, n2, expr.then_expr, indent + "    ")
            code += f"{indent}}} else {{\n"
            code += self.emit_noise_terms(n1, n2, expr.else_expr, indent + "    ")
            code += f"{indent}}}\n"
            return code
        return ""

    def emit_noise_statement_cpp(self, stmt, indent: str = "            ") -> str:
        if isinstance(stmt, AssignmentNode):
            lhs = self.cpp_name(stmt.var_name)
            return f"{indent}{lhs} = {self.expr_cpp(stmt.expr)};\n"
        if isinstance(stmt, BranchContributionNode):
            if stmt.quantity != "I" or not self.expr_has_noise(stmt.expr):
                return ""
            n1 = self.module.ports.index(stmt.node1) if stmt.node1 in self.module.ports else -1
            n2 = self.module.ports.index(stmt.node2) if stmt.node2 in self.module.ports else -1
            return self.emit_noise_terms(n1, n2, stmt.expr, indent)
        if isinstance(stmt, ConditionalNode):
            code = f"{indent}if (gmc_nz({self.expr_cpp(stmt.condition)})) {{\n"
            for child in stmt.then_body:
                code += self.emit_noise_statement_cpp(child, indent + "    ")
            code += f"{indent}}} else {{\n"
            for child in stmt.else_body:
                code += self.emit_noise_statement_cpp(child, indent + "    ")
            code += f"{indent}}}\n"
            return code
        if isinstance(stmt, GroupNode):
            return "".join(self.emit_noise_statement_cpp(child, indent) for child in stmt.body)
        return ""

    def split_ddt(self, expr, sign: int = 1):
        if isinstance(expr, BinaryOpNode) and expr.op in ("+", "-"):
            static_l, dynamic_l = self.split_ddt(expr.left, sign)
            static_r, dynamic_r = self.split_ddt(expr.right, sign if expr.op == "+" else -sign)
            return static_l + static_r, dynamic_l + dynamic_r
        if isinstance(expr, UnaryOpNode) and expr.op == "-":
            return self.split_ddt(expr.expr, -sign)
        if isinstance(expr, FunctionCallNode) and expr.func_name == "ddt" and len(expr.args) == 1:
            return [], [(sign, expr.args[0])]
        if isinstance(expr, BinaryOpNode) and expr.op == "*":
            # Charge partition idiom "C * ddt(q)": when the multiplier C is
            # time-invariant (no ddt, no branch probes) the derivative factor
            # is constant and ddt can be hoisted out: C * ddt(q) -> ddt(C * q).
            if (isinstance(expr.right, FunctionCallNode)
                    and expr.right.func_name == "ddt"
                    and len(expr.right.args) == 1
                    and self.expr_is_time_invariant(expr.left)):
                return [], [(sign, BinaryOpNode("*", expr.left, expr.right.args[0]))]
            if (isinstance(expr.left, FunctionCallNode)
                    and expr.left.func_name == "ddt"
                    and len(expr.left.args) == 1
                    and self.expr_is_time_invariant(expr.right)):
                return [], [(sign, BinaryOpNode("*", expr.right, expr.left.args[0]))]
        return [(sign, expr)], []

    def expr_is_time_invariant(self, expr) -> bool:
        if isinstance(expr, BranchProbeNode):
            return False
        if isinstance(expr, (FunctionCallNode, UnaryOpNode, BinaryOpNode, TernaryNode)):
            for node in expr.__dict__.values():
                if isinstance(node, ExpressionNode):
                    if not self.expr_is_time_invariant(node):
                        return False
                elif isinstance(node, list):
                    for item in node:
                        if isinstance(item, ExpressionNode) \
                                and not self.expr_is_time_invariant(item):
                            return False
        return True

    def sum_cpp(self, terms) -> str:
        if not terms:
            return "0.0"
        chunks = []
        for sign, expr in terms:
            text = self.expr_cpp(expr)
            chunks.append(text if sign > 0 else f"-({text})")
        return " + ".join(chunks)

    def emit_branch_add(self, n1: int, n2: int, expr: str, target: str, indent: str) -> str:
        code = f"{indent}{target}[{n1}] += {expr};\n"
        if n2 >= 0:
            code += f"{indent}{target}[{n2}] -= {expr};\n"
        return code

    def emit_statement_cpp(self, stmt, indent: str = "            ") -> str:
        if isinstance(stmt, AssignmentNode):
            if self.expr_has_ddt(stmt.expr):
                raise ValueError(
                    f"module '{self.module.name}': ddt() inside an assignment is not "
                    f"supported; apply ddt directly in a branch contribution")
            lhs = self.cpp_name(stmt.var_name)
            if self.current_function is not None and stmt.var_name == self.current_function.name:
                lhs = f"{self.cpp_name(stmt.var_name)}_ret"
            return f"{indent}{lhs} = {self.expr_cpp(stmt.expr)};\n"
        if isinstance(stmt, BranchContributionNode):
            if stmt.quantity != "I":
                return ""
            n1 = self.module.ports.index(stmt.node1) if stmt.node1 in self.module.ports else -1
            n2 = self.module.ports.index(stmt.node2) if stmt.node2 in self.module.ports else -1
            static_terms, dynamic_terms = self.split_ddt(stmt.expr)
            for sign, term_expr in static_terms:
                if self.expr_has_ddt(term_expr):
                    raise ValueError(
                        f"module '{self.module.name}': ddt() nested inside a function "
                        f"call or product cannot be split; restructure the contribution")
            code = ""
            if static_terms:
                code += self.emit_branch_add(n1, n2, self.sum_cpp(static_terms), "currents", indent)
            if dynamic_terms:
                code += self.emit_branch_add(n1, n2, self.sum_cpp(dynamic_terms), "charges", indent)
            return code
        if isinstance(stmt, ConditionalNode):
            code = f"{indent}if (gmc_nz({self.expr_cpp(stmt.condition)})) {{\n"
            for child in stmt.then_body:
                code += self.emit_statement_cpp(child, indent + "    ")
            code += f"{indent}}} else {{\n"
            for child in stmt.else_body:
                code += self.emit_statement_cpp(child, indent + "    ")
            code += f"{indent}}}\n"
            return code
        if isinstance(stmt, CaseNode):
            # Verilog-A case emits as a chained if / else-if ladder; a
            # "default" item with no values becomes the trailing else.
            sel = self.expr_cpp(stmt.expr)
            generated = False
            for item in stmt.items:
                if item.values:
                    conditions = " || ".join(
                        f"gmc_nz({sel} == {self.expr_cpp(v)})" for v in item.values)
                    if generated:
                        code += f"{indent}else if ({conditions}) {{\n"
                    else:
                        code = f"{indent}if ({conditions}) {{\n"
                        generated = True
                else:
                    if generated:
                        code += f"{indent}else {{\n"
                    else:
                        code = f"{indent}if (true) {{\n"
                        generated = True
                for child in item.body:
                    code += self.emit_statement_cpp(child, indent + "    ")
                code += f"{indent}}}\n"
            if not generated:
                code = ""
            return code
        if isinstance(stmt, ForNode):
            init = self._loop_part_cpp(stmt.init)
            cond = (f"gmc_nz({self.expr_cpp(stmt.condition)})"
                    if stmt.condition else "")
            incr = self._loop_part_cpp(stmt.increment)
            code = f"{indent}for ({init}; {cond}; {incr}) {{\n"
            for child in stmt.body:
                code += self.emit_statement_cpp(child, indent + "    ")
            code += f"{indent}}}\n"
            return code
        if isinstance(stmt, WhileNode):
            code = f"{indent}while (gmc_nz({self.expr_cpp(stmt.condition)})) {{\n"
            for child in stmt.body:
                code += self.emit_statement_cpp(child, indent + "    ")
            code += f"{indent}}}\n"
            return code
        if isinstance(stmt, EventNode):
            # GMC evaluates the model body once per solver evaluation; event
            # scoping (@initial_step, @(tran), ...) is currently collapsed
            # into that single evaluation pass.
            code = f"{indent}// @({stmt.event_name}) event block\n"
            for child in stmt.body:
                code += self.emit_statement_cpp(child, indent)
            return code
        if isinstance(stmt, GroupNode):
            # Plain begin...end grouping/scope block: no control flow, so its
            # statements execute inline in program order (VBIC's labeled
            # evaluateStatic/loadStatic blocks, BSIM's bare grouped if bodies).
            code = ""
            for child in stmt.body:
                code += self.emit_statement_cpp(child, indent)
            return code
        return ""

    def _loop_part_cpp(self, part: Optional[AssignmentNode]) -> str:
        if part is None:
            return ""
        return f"{self.cpp_name(part.var_name)} = {self.expr_cpp(part.expr)}"

    def expr_has_ddt(self, expr) -> bool:
        if isinstance(expr, FunctionCallNode):
            if expr.func_name == "ddt":
                return True
            return any(self.expr_has_ddt(arg) for arg in expr.args)
        if isinstance(expr, UnaryOpNode):
            return self.expr_has_ddt(expr.expr)
        if isinstance(expr, BinaryOpNode):
            return self.expr_has_ddt(expr.left) or self.expr_has_ddt(expr.right)
        if isinstance(expr, TernaryNode):
            return (self.expr_has_ddt(expr.condition)
                    or self.expr_has_ddt(expr.then_expr)
                    or self.expr_has_ddt(expr.else_expr))
        return False

    def has_ddt(self, statements) -> bool:
        for stmt in statements:
            if isinstance(stmt, ConditionalNode):
                if self.has_ddt(stmt.then_body) or self.has_ddt(stmt.else_body):
                    return True
            elif isinstance(stmt, CaseNode):
                if self.has_ddt([child for item in stmt.items for child in item.body]):
                    return True
            elif isinstance(stmt, ForNode):
                if (self.has_ddt(stmt.body)
                        or (stmt.condition is not None and self.expr_has_ddt(stmt.condition))):
                    return True
            elif isinstance(stmt, WhileNode):
                if self.has_ddt(stmt.body) or self.expr_has_ddt(stmt.condition):
                    return True
            elif isinstance(stmt, EventNode):
                if self.has_ddt(stmt.body):
                    return True
            elif isinstance(stmt, GroupNode):
                if self.has_ddt(stmt.body):
                    return True
            elif isinstance(stmt, BranchContributionNode):
                if self.expr_has_ddt(stmt.expr):
                    return True
            elif isinstance(stmt, AssignmentNode):
                if self.expr_has_ddt(stmt.expr):
                    return True
        return False

    def expr_calls_functions(self, expr, names):
        if isinstance(expr, FunctionCallNode):
            if expr.func_name in names:
                return True
            return any(self.expr_calls_functions(arg, names) for arg in expr.args)
        if isinstance(expr, UnaryOpNode):
            return self.expr_calls_functions(expr.expr, names)
        if isinstance(expr, BinaryOpNode):
            return self.expr_calls_functions(expr.left, names) or self.expr_calls_functions(expr.right, names)
        if isinstance(expr, TernaryNode):
            return (self.expr_calls_functions(expr.condition, names)
                    or self.expr_calls_functions(expr.then_expr, names)
                    or self.expr_calls_functions(expr.else_expr, names))
        return False

    def stmts_call_functions(self, statements, names):
        for stmt in statements:
            if isinstance(stmt, ConditionalNode):
                if self.stmts_call_functions(stmt.then_body, names) or self.stmts_call_functions(stmt.else_body, names):
                    return True
            elif isinstance(stmt, CaseNode):
                if self.stmts_call_functions([child for item in stmt.items for child in item.body], names):
                    return True
            elif isinstance(stmt, ForNode):
                if self.stmts_call_functions(stmt.body, names):
                    return True
            elif isinstance(stmt, WhileNode):
                if self.stmts_call_functions(stmt.body, names):
                    return True
            elif isinstance(stmt, EventNode):
                if self.stmts_call_functions(stmt.body, names):
                    return True
            elif isinstance(stmt, GroupNode):
                if self.stmts_call_functions(stmt.body, names):
                    return True
            elif isinstance(stmt, AssignmentNode):
                if self.expr_calls_functions(stmt.expr, names):
                    return True
            elif isinstance(stmt, BranchContributionNode):
                if self.expr_calls_functions(stmt.expr, names):
                    return True
        return False

    # Returns user functions in definition order such that every callee
    # appears before its callers (fallback: source order). Raises on cycles.
    def function_order(self):
        names = set(self.module.functions)
        order = []
        visited = set()
        visiting = set()

        def visit(name):
            if name in visited:
                return
            if name in visiting:
                raise ValueError(
                    f"module '{self.module.name}': recursive analog function "
                    f"'{name}' is not supported")
            visiting.add(name)
            fn = self.module.functions[name]
            for dep in names:
                if not self.stmts_call_functions(fn.body, {dep}):
                    continue
                if dep == name:
                    raise ValueError(
                        f"module '{self.module.name}': recursive analog function "
                        f"'{name}' is not supported")
                if dep not in visited:
                    visit(dep)
            visiting.discard(name)
            visited.add(name)
            order.append(name)

        for name in names:
            visit(name)
        return order

    def emit_function(self, name: str):
        fn = self.module.functions[name]
        indent = "            "
        if fn.return_type not in ("real", "integer"):
            raise ValueError(
                f"module '{self.module.name}': analog function '{name}' return "
                f"type '{fn.return_type}' is not supported (only 'real'/'integer')")
        args = ", ".join(f"const auto& {self.cpp_name(a.name)}" for a in fn.args)
        arg_names = {self.cpp_name(a.name) for a in fn.args}
        first_arg = self.cpp_name(fn.args[0].name) if fn.args else None
        scalar_using = (
            f"using Scalar = std::decay_t<decltype({first_arg})>;"
            if first_arg else "using Scalar = double;"
        )
        old_fn = self.current_function
        old_arg_names = self.current_function_arg_names
        self.current_function = fn
        self.current_function_arg_names = arg_names
        try:
            body = "".join(self.emit_statement_cpp(stmt, indent + "    ") for stmt in fn.body)
        finally:
            self.current_function = old_fn
            self.current_function_arg_names = old_arg_names
        return (
            f"        auto {self.cpp_name(fn.name)} = [&]({args}) {{\n"
            f"{indent}    {scalar_using}\n"
            f"{indent}    Scalar {self.cpp_name(fn.name)}_ret = 0.0;\n"
            + "".join(f"{indent}    Scalar {self.cpp_name(l.name)} = 0.0;\n" for l in fn.locals)
            + body
            + f"{indent}    return {self.cpp_name(fn.name)}_ret;\n"
            + f"{indent}}};\n"
        )

    # Returns (terminal_count, [(role, partner_index)] per port).
    # role: "Terminal" | "Internal" | "Collapsible"; partner_index is the
    # local index of the terminal the collapsible node merges into. Fail
    # closed: invalid node declarations raise instead of generating broken
    # models.
    def local_node_info(self):
        ports = self.module.ports
        terminal_count = self.module.terminal_count or len(ports)
        if terminal_count < 0 or terminal_count > len(ports):
            raise ValueError(
                f"module '{self.module.name}': terminal {terminal_count} exceeds "
                f"{len(ports)} declared ports")
        roles = [("Terminal", -1)] * len(ports)
        for hidden, partner in self.module.collapsible_pairs.items():
            if hidden not in ports:
                raise ValueError(
                    f"module '{self.module.name}': collapsible node '{hidden}' is "
                    f"not a declared port")
            if partner not in ports[:terminal_count]:
                raise ValueError(
                    f"module '{self.module.name}': collapsible node '{hidden}' must "
                    f"merge into an external terminal, got '{partner}'")
            roles[ports.index(hidden)] = ("Collapsible", ports.index(partner))
        for i in range(terminal_count, len(ports)):
            if roles[i][0] != "Collapsible":
                roles[i] = ("Internal", -1)
        return terminal_count, roles

    def emit_cpp_header(self) -> str:
        class_name = f"Gmc{self.class_stem()}Instance"
        model_name = f"Gmc{self.class_stem()}Model"

        ports = self.module.ports
        num_ports = len(ports)
        terminal_count, roles = self.local_node_info()
        internal_count = sum(1 for role, _ in roles if role == "Internal")

        params_decl = ""
        params_init = ""
        params_info = ""
        aliases_init = ""
        for name, pdecl in self.module.parameters.items():
            cpp = self.cpp_name(name)
            if pdecl.default_expr is not None:
                # Only bare aliases ("parameter real dsub = drout;") are
                # supported: resolved per instance after explicit bindings.
                target = pdecl.default_expr
                if not (isinstance(target, IdentifierNode)
                        and target.name in self.module.parameters):
                    raise ValueError(
                        f"module '{self.module.name}': parameter '{name}' has a "
                        f"non-literal default expression, which GMC emission "
                        f"supports only for plain parameter aliases")
                target_cpp = self.cpp_name(target.name)
                params_decl += f"    double {cpp}_ = {pdecl.default_val};\n"
                params_decl += f"    bool given_{cpp} = false;\n"
                aliases_init += (f"        if (card.parameters.count(\"{name}\") == 0) "
                                 f"instance->{cpp}_ = instance->{target_cpp}_;\n")
                aliases_init += (f"        instance->given_{cpp} = "
                                 f"(card.parameters.count(\"{name}\") != 0);\n")
                params_info += f'            {{"{name}", {pdecl.default_val}, "", "", true}},\n'
                continue
            params_decl += f"    double {cpp}_ = {pdecl.default_val};\n"
            params_decl += f"    bool given_{cpp} = false;\n"
            params_init += f"        auto it_{cpp} = card.parameters.find(\"{name}\"); if (it_{cpp} != card.parameters.end()) instance->{cpp}_ = it_{cpp}->second;\n"
            params_init += f"        instance->given_{cpp} = (card.parameters.count(\"{name}\") != 0);\n"
            params_info += f'            {{"{name}", {pdecl.default_val}, "", "", true}},\n'

        nodes_info = ""
        for i, port in enumerate(ports):
            role, partner = roles[i]
            if role == "Collapsible":
                nodes_info += '                {"%s", %d, GsdiNodeRole::Collapsible, %d},\n' % (port, i, partner)
            elif role == "Internal":
                nodes_info += '                {"%s", %d, GsdiNodeRole::Internal, GsdiNoCollapse},\n' % (port, i)
            else:
                nodes_info += '                {"%s", %d, GsdiNodeRole::Terminal, GsdiNoCollapse},\n' % (port, i)

        jacobian_pattern = ""
        for row in range(num_ports):
            for col in range(num_ports):
                jacobian_pattern += f"                {{{row}, {col}}},\n"

        vars_decl = ""
        noise_vars_decl = ""
        for name in self.module.variables:
            vars_decl += f"            Scalar {self.cpp_name(name)} = 0.0;\n"
            noise_vars_decl += f"            double {self.cpp_name(name)} = 0.0;\n"

        body = ""
        for stmt in self.module.analog_body:
            body += self.emit_statement_cpp(stmt)

        supports_transient = "true" if self.has_ddt(self.module.analog_body) else "false"
        supports_noise = "true" if self.has_noise(self.module.analog_body) else "false"
        noise_body = "".join(self.emit_noise_statement_cpp(stmt) for stmt in self.module.analog_body)

        # User analog functions are emitted as generic lambdas before
        # eval_model so the model body can call them. Definitions are ordered
        # callee-first; recursion/cycles fail closed.
        functions_cpp = ""
        if self.module.functions:
            funcs = "\n".join(self.emit_function(name)
                              for name in self.function_order())
            functions_cpp = "\n" + funcs

        code = f"""#ifndef GSPICE_GMC_{self.module.name.upper()}_HPP
#define GSPICE_GMC_{self.module.name.upper()}_HPP

#include "gmc_dual.hpp"
#include "gsdi.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

namespace gspice {{

class {class_name} final : public GsdiInstance {{
public:
    // System functions: $temperature / $vt resolve to the instance's ambient
    // temperature (300.0 K unless a model-specific value is injected later).
    double temperature_ = 300.0;
    double vt_ = temperature_ * 8.617333262145179e-5;
    // analysis("...") resolves against this; defaults to DC.
    std::string analysis_ = "dc";
{params_decl}
    bool evaluate(
        const GsdiEvalRequest& request,
        GsdiEvalResult& result) override {{
        result.clear();
        if (!request.solution || request.solution_size < {num_ports}) return false;

        auto gmc_exp = [](const auto& x) {{
            using Scalar = std::decay_t<decltype(x)>;
            if constexpr (is_gmc_dual_v<Scalar>) return gmcExp(x);
            else return std::exp(x);
        }};
        auto gmc_log = [](const auto& x) {{
            using Scalar = std::decay_t<decltype(x)>;
            if constexpr (is_gmc_dual_v<Scalar>) return gmcLog(x);
            else return std::log(x);
        }};
        auto gmc_sqrt = [](const auto& x) {{
            using Scalar = std::decay_t<decltype(x)>;
            if constexpr (is_gmc_dual_v<Scalar>) return gmcSqrt(x);
            else return std::sqrt(x);
        }};
        auto gmc_pow = [](const auto& x, double p) {{
            using Scalar = std::decay_t<decltype(x)>;
            if constexpr (is_gmc_dual_v<Scalar>) return gmcPow(x, p);
            else return std::pow(x, p);
        }};
        auto gmc_abs = [](const auto& x) {{
            return x < 0.0 ? -x : x;
        }};
        auto gmc_log10 = [](const auto& x) {{
            using Scalar = std::decay_t<decltype(x)>;
            if constexpr (is_gmc_dual_v<Scalar>) return gmcLog(x) / 2.302585092994046;
            else return std::log10(x);
        }};
        auto gmc_tanh = [](const auto& x) {{
            using Scalar = std::decay_t<decltype(x)>;
            if constexpr (is_gmc_dual_v<Scalar>) return gmcTanh(x);
            else return std::tanh(x);
        }};
        auto gmc_sin = [](const auto& x) {{
            using Scalar = std::decay_t<decltype(x)>;
            if constexpr (is_gmc_dual_v<Scalar>) return gmcSin(x);
            else return std::sin(x);
        }};
        auto gmc_cos = [](const auto& x) {{
            using Scalar = std::decay_t<decltype(x)>;
            if constexpr (is_gmc_dual_v<Scalar>) return gmcCos(x);
            else return std::cos(x);
        }};
        auto gmc_tan = [](const auto& x) {{
            using Scalar = std::decay_t<decltype(x)>;
            if constexpr (is_gmc_dual_v<Scalar>) return gmcTan(x);
            else return std::tan(x);
        }};
        auto gmc_atan = [](const auto& x) {{
            using Scalar = std::decay_t<decltype(x)>;
            if constexpr (is_gmc_dual_v<Scalar>) return gmcAtan(x);
            else return std::atan(x);
        }};
        auto gmc_fmod = [](const auto& a, const auto& b) {{
            using A = std::decay_t<decltype(a)>;
            using B = std::decay_t<decltype(b)>;
            if constexpr (is_gmc_dual_v<A>) {{
                if constexpr (is_gmc_dual_v<B>) return gmcFmod(a, b);
                else return gmcFmod(a, b);
            }} else {{
                return std::fmod(a, b);
            }}
        }};
        auto gmc_sinh = [](const auto& x) {{
            using Scalar = std::decay_t<decltype(x)>;
            if constexpr (is_gmc_dual_v<Scalar>) return gmcSinh(x);
            else return std::sinh(x);
        }};
        auto gmc_cosh = [](const auto& x) {{
            using Scalar = std::decay_t<decltype(x)>;
            if constexpr (is_gmc_dual_v<Scalar>) return gmcCosh(x);
            else return std::cosh(x);
        }};
        auto gmc_min = [](const auto& a, const auto& b) {{
            using A = std::decay_t<decltype(a)>;
            using B = std::decay_t<decltype(b)>;
            if constexpr (is_gmc_dual_v<A>) return gmcMin(a, b);
            else if constexpr (is_gmc_dual_v<B>) {{ return gmcMin(b, a); }}
            else return std::min(a, b);
        }};
        auto gmc_max = [](const auto& a, const auto& b) {{
            using A = std::decay_t<decltype(a)>;
            using B = std::decay_t<decltype(b)>;
            if constexpr (is_gmc_dual_v<A>) return gmcMax(a, b);
            else if constexpr (is_gmc_dual_v<B>) {{ return gmcMax(b, a); }}
            else return std::max(a, b);
        }};
        auto gmc_pow2 = [](const auto& a, const auto& b) {{
            using A = std::decay_t<decltype(a)>;
            using B = std::decay_t<decltype(b)>;
            if constexpr (is_gmc_dual_v<A>) {{
                if constexpr (is_gmc_dual_v<B>) return gmcExp(gmcLog(a) * b);
                else return gmcPow(a, b);
            }} else {{
                if constexpr (is_gmc_dual_v<B>) return gmcExp(b * std::log(a));
                else return std::pow(a, b);
            }}
        }};
        auto limexp = [](const auto& x) {{
            using Scalar = std::decay_t<decltype(x)>;
            if constexpr (is_gmc_dual_v<Scalar>) {{
                if (x.value < -40.0) return Scalar::constant(std::exp(-40.0));
                if (x.value > 40.0) return Scalar::constant(std::exp(40.0));
                return gmcExp(x);
            }} else {{
                return std::exp(std::clamp(x, -40.0, 40.0));
            }}
        }};
        auto gmc_nz = [](const auto& x) {{
            using X = std::decay_t<decltype(x)>;
            if constexpr (is_gmc_dual_v<X>) return x.value != 0.0;
            else return x != 0.0;
        }};{functions_cpp}
        auto eval_model = [&](const auto& vn) {{
            using Scalar = std::decay_t<decltype(vn[0])>;
            std::array<Scalar, {num_ports}> currents{{}};
            std::array<Scalar, {num_ports}> charges{{}};
"""
        code += vars_decl
        code += body
        code += f"""            return std::make_pair(currents, charges);
        }};

        std::array<double, {num_ports}> v{{}};
        for (std::size_t i = 0; i < {num_ports}; ++i) v[i] = request.solution[i];
        const auto evaluated = eval_model(v);

        if (request.residual) {{
            for (int row = 0; row < {num_ports}; ++row) {{
                result.static_residual.push_back({{row, evaluated.first[row], 0}});
            }}
        }}

        if (request.jacobian) {{
            std::array<GmcDual<{num_ports}>, {num_ports}> vd{{}};
            for (std::size_t i = 0; i < {num_ports}; ++i) vd[i] = GmcDual<{num_ports}>::variable(request.solution[i], i);
            const auto differentiated = eval_model(vd);
            for (int col = 0; col < {num_ports}; ++col) {{
                for (int row = 0; row < {num_ports}; ++row) {{
                    result.static_jacobian.push_back({{row, col, differentiated.first[row].derivative[col], 0}});
                }}
            }}
        }}

        if (request.dynamic_residual) {{
            for (int row = 0; row < {num_ports}; ++row) {{
                result.dynamic_residual.push_back({{row, evaluated.second[row], 0}});
            }}
        }}

        if (request.dynamic_jacobian) {{
            std::array<GmcDual<{num_ports}>, {num_ports}> vd{{}};
            for (std::size_t i = 0; i < {num_ports}; ++i) vd[i] = GmcDual<{num_ports}>::variable(request.solution[i], i);
            const auto differentiated = eval_model(vd);
            for (int col = 0; col < {num_ports}; ++col) {{
                for (int row = 0; row < {num_ports}; ++row) {{
                    result.dynamic_jacobian.push_back({{row, col, differentiated.second[row].derivative[col], 0}});
                }}
            }}
        }}

        if (request.noise) {{
            std::array<double, {num_ports}> vn{{}};
            for (std::size_t i = 0; i < {num_ports}; ++i) vn[i] = request.solution[i];
            auto add_noise = [&](int pos, int neg, double psd, const std::string& name) {{
                if (psd > 0.0 && std::isfinite(psd)) result.noise.push_back({{pos, neg, psd, name}});
            }};
{noise_vars_decl}{noise_body}        }}

        return true;
    }}

std::size_t terminalCount() const override {{ return {terminal_count}; }}
        std::size_t internalNodeCount() const override {{ return {internal_count}; }}
    std::size_t stateBytes() const override {{ return 0; }}
}};

class {model_name} final : public GsdiModel {{
public:
    const GsdiModelDescriptor& descriptor() const override {{
        static const GsdiModelDescriptor d = [] {{
            GsdiModelDescriptor d;
            d.model_type = "{self.module.name}";
            d.version = "1.0";
            d.terminal_count = {terminal_count};
            d.nodes = {{
{nodes_info}            }};
            d.parameters = {{
{params_info}            }};
            d.jacobian_pattern = {{
{jacobian_pattern}            }};
            d.supports_op = true;
            d.supports_transient = {supports_transient};
            d.supports_ac = {supports_transient};
            d.supports_noise = {supports_noise};
            return d;
        }}();
        return d;
    }}

    std::unique_ptr<GsdiInstance> createInstance(
        const GsdiModelCard& card,
        const std::vector<int>& terminal_nodes) const override {{
        (void)terminal_nodes;
        auto instance = std::make_unique<{class_name}>();
{params_init}{aliases_init}        return instance;
    }}
}};

}} // namespace gspice

#endif // GSPICE_GMC_{self.module.name.upper()}_HPP
"""
        return code
