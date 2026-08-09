"""
GSPICE Model Compiler (GMC) - native model-source lexer and parser.
Parses GMC model source into the internal AST representation.
"""

import re
import sys
from pathlib import Path
from typing import List, Optional
from gmc_ast import (
    ModuleNode, ParameterDecl, VariableDecl, ExpressionNode, NumberNode,
    IdentifierNode, BinaryOpNode, UnaryOpNode, FunctionCallNode, BranchProbeNode,
    BranchContributionNode, AssignmentNode, StatementNode, ConditionalNode,
    TernaryNode, FunctionDeclNode, CaseNode, CaseItem, ForNode, WhileNode,
    EventNode, GroupNode
)
from gmc_preprocessor import preprocess_source

NOOP_SYSTEM_TASKS = {
    '$strobe', '$display', '$fatal', '$finish', '$error', '$write',
    '$monitor', '$warning', '$discontinuity', '$bound_step',
}

class GMCParser:
    def __init__(self, source_code: str, *, include_dirs=(), source_name=None):
        self._preprocessed = preprocess_source(
            source_code, filename=source_name,
            include_dirs=[Path(d) for d in include_dirs])
        self.source = self.strip_comments(self._preprocessed.text)
        self.line_map = self._preprocessed.line_map
        self.preprocessed_files = self._preprocessed.files
        self.pos = 0
        self.tokens = self.tokenize(self.source)
        self.token_idx = 0
        self._block_locals = {}

    def strip_comments(self, code: str) -> str:
        # Strip line comments //... and block comments /*...*/
        code = re.sub(r'//.*', '', code)
        code = re.sub(r'/\*.*?\*/', '', code, flags=re.DOTALL)
        return code

    def tokenize(self, code: str) -> List[str]:
        # Restricted Verilog-A lexer for GMC model-source input.
        token_spec = [
            ('STRING',   r'"(?:[^"\\]|\\.)*"'),
            ('LOGAND',   r'&&'),
            ('LOGOR',    r'\|\|'),
            ('NUMBER',   r'(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?(?:[uUnNpPfFkKmMgGtT])?'),
            ('ASSIGN',   r'<\+|==|!=|<=|>=|=|<|>'),
            ('STARSTAR', r'\*\*'),
            ('OP',       r'\+|\-|\*|\/|\?|\:'),
            ('AMP',      r'&'),
            ('PIPE',     r'\|'),
            ('BANG',     r'!'),
            ('PERCENT',  r'%'),
            ('AT',       r'@'),
            ('PAREN',    r'\(|\)'),
            ('BRACE',    r'\{|\}'),
            ('SEMI',     r';'),
            ('COMMA',    r','),
            ('ID',       r'\$?[a-zA-Z_][a-zA-Z0-9_]*'),
            ('SKIP',     r'[ \t\r\n]+'),
        ]
        tok_regex = '|'.join(f'(?P<{pair[0]}>{pair[1]})' for pair in token_spec)
        tokens = []
        for mo in re.finditer(tok_regex, code):
            kind = mo.lastgroup
            value = mo.group()
            if kind == 'SKIP':
                continue
            tokens.append(value)
        return tokens

    def peek(self) -> Optional[str]:
        if self.token_idx >= len(self.tokens):
            return None
        return self.tokens[self.token_idx]

    def take(self, expected: Optional[str] = None) -> str:
        tok = self.peek()
        if tok is None:
            raise ValueError("Unexpected end of GMC source")
        if expected is not None and tok != expected:
            raise ValueError(f"Expected '{expected}', got '{tok}'")
        self.token_idx += 1
        return tok

    def parse_module(self) -> ModuleNode:
        while self.token_idx < len(self.tokens):
            tok = self.tokens[self.token_idx]
            if tok == 'module':
                self.token_idx += 1
                name = self.tokens[self.token_idx]
                self.token_idx += 1
                # Parse ports (p1, p2, p3)
                ports = []
                if self.peek() == '(':
                    self.take('(')
                    while self.peek() != ')':
                        if self.peek() != ',':
                            ports.append(self.take())
                        else:
                            self.take(',')
                    self.take(')')
                if self.peek() == ';':
                    self.take(';')
                
                module = ModuleNode(name, ports)
                self.parse_module_body(module)
                return module
            self.token_idx += 1
        raise ValueError("No GMC module declaration found")

    def parse_module_body(self, module: ModuleNode):
        while self.token_idx < len(self.tokens):
            tok = self.tokens[self.token_idx]
            if tok == 'endmodule':
                break
            elif tok == 'parameter':
                self.parse_parameter(module)
            elif tok == 'terminal':
                self.parse_terminal_decl(module)
            elif tok == 'collapsible':
                self.parse_collapsible_decl(module)
            elif tok in ('electrical', 'inout', 'input', 'output'):
                self.parse_node_decl(module, tok)
            elif tok in ('real', 'integer'):
                self.parse_variable(module)
            elif tok in ('nature', 'discipline', 'connectrules'):
                self.skip_block('end' + tok)
            elif tok == 'analog':
                if self.token_idx + 1 < len(self.tokens) and self.tokens[self.token_idx + 1] == 'function':
                    self.token_idx += 1
                    func = self.parse_function()
                    module.functions[func.name] = func
                else:
                    self.token_idx += 1
                    self.parse_analog_block(module)
            else:
                if tok.startswith('$') and self.token_idx + 1 < len(self.tokens) \
                        and self.tokens[self.token_idx + 1] == '(':
                    self.skip_paren_call()
                else:
                    self.token_idx += 1

    def parse_parameter(self, module: ModuleNode):
        self.take('parameter')
        ptype = "real"
        if self.peek() in ('real', 'integer'):
            ptype = self.take()
        pname = self.take()
        default_val = 0.0
        default_expr = None
        if self.peek() == '=':
            self.take('=')
            if self.peek() is not None and self.tokens[self.token_idx] not in (';', ','):
                nxt = self.tokens[self.token_idx]
                nxt2 = (self.tokens[self.token_idx + 1]
                        if self.token_idx + 1 < len(self.tokens) else None)
                # Literal unless it is a lone sign (BSIM macro defaults such as
                # "= -9999999" tokenize as '-' NUMBER).
                is_literal = (re.match(r'[+-]?(?:\d|\.)', nxt) is not None) or (
                    nxt in ('+', '-') and nxt2 is not None
                    and re.match(r'(?:\d|\.)', nxt2) is not None)
                if is_literal:
                    default_val = self.parse_const_number()
                else:
                    # Non-literal default (e.g. "parameter real dsub=drout;").
                    default_expr = self.parse_expression()
        while self.peek() not in (';', None):
            self.token_idx += 1
        if self.peek() == ';':
            self.take(';')
        module.parameters[pname] = ParameterDecl(pname, default_val, ptype, default_expr)

    def parse_terminal_decl(self, module: ModuleNode):
        self.take('terminal')
        tok = self.take()
        module.terminal_count = int(tok)
        if self.peek() == ';':
            self.take(';')

    def parse_node_decl(self, module: ModuleNode, kind: str):
        # "electrical d, g, s, b, di, si;" — names not declared as module
        # ports are hidden (internal) nodes. "inout/input/output" are mere
        # port/type aliases and never add local nodes.
        self.take(kind)
        names = []
        while self.peek() not in (';', None):
            name = self.take()
            if name != ',':
                names.append(name)
        if self.peek() == ';':
            self.take(';')
        if kind == 'electrical':
            for name in names:
                if name not in module.ports:
                    module.ports.append(name)
                    if module.terminal_count == 0:
                        module.terminal_count = len(module.ports) - 1

    def parse_collapsible_decl(self, module: ModuleNode):
        self.take('collapsible')
        name = self.take()
        if self.peek() == 'with':
            self.take('with')
        partner = self.take()
        if self.peek() == ';':
            self.take(';')
        module.collapsible_pairs[name] = partner

    def parse_variable(self, module: ModuleNode):
        vtype = self.take()
        while self.peek() not in (';', None):
            vname = self.take()
            if vname != ',':
                module.variables[vname] = VariableDecl(vname, vtype)
        if self.peek() == ';':
            self.take(';')

    def parse_function(self) -> FunctionDeclNode:
        self.take('function')
        # Verilog-A allows the return type to be omitted; it defaults to real.
        if self.peek() in ('real', 'integer'):
            return_type = self.take()
        else:
            return_type = 'real'
        name = self.take()
        if self.peek() == ';':
            self.take(';')
        # Port declarations ("input x, y;") followed by type declarations
        # ("real x, y, vlimited;"). Names declared as input are arguments;
        # any other typed name is a local variable.
        port_names = []
        typed = {}
        while self.peek() in ('input', 'output', 'inout', 'real', 'integer'):
            tok = self.take()
            if tok in ('input', 'output', 'inout'):
                while self.peek() not in (';', None):
                    pn = self.take()
                    if pn != ',' and pn not in port_names:
                        port_names.append(pn)
                self.take(';')
            elif tok in ('real', 'integer'):
                while self.peek() not in (';', None):
                    pn = self.take()
                    if pn != ',':
                        typed[pn] = tok
                self.take(';')
            else:
                raise ValueError(f"Unexpected token '{tok}' in analog function '{name}'")
        args = [VariableDecl(pn, typed.get(pn, 'real')) for pn in port_names]
        locals_ = [VariableDecl(pn, ty) for pn, ty in typed.items() if pn not in port_names]
        body = []
        if self.peek() == 'begin':
            self.take('begin')
            while self.peek() not in ('end', 'endfunction', None):
                stmt = self.parse_statement()
                if stmt:
                    body.append(stmt)
            if self.peek() == 'end':
                self.take('end')
        else:
            while self.peek() not in ('endfunction', None):
                stmt = self.parse_statement()
                if stmt:
                    body.append(stmt)
        if self.peek() == 'endfunction':
            self.take('endfunction')
        return FunctionDeclNode(name, return_type, args, locals_, body)

    def skip_decl(self):
        while self.peek() not in (';', None):
            self.token_idx += 1
        if self.peek() == ';':
            self.take(';')

    def skip_block(self, end_token: str):
        # Ignore a whole region between two keywords (e.g. nature/endnature,
        # discipline/enddiscipline, connectrules/endconnectrules) so that
        # `include "disciplines.vams"` artifacts parse without error.
        while self.token_idx < len(self.tokens) and self.peek() != end_token:
            self.token_idx += 1
        if self.peek() == end_token:
            self.token_idx += 1

    def parse_analog_block(self, module: ModuleNode):
        if self.peek() == 'begin':
            self.take('begin')
        while self.peek() not in ('end', 'endmodule', None):
            stmt = self.parse_statement()
            if stmt:
                module.analog_body.append(stmt)
        if self.peek() == 'end':
            self.take('end')
        # Block-scoped "real x, y; ..." locals (VBIC declares them inside
        # labeled begin : ... end scopes) are hoisted to module scope: the
        # emitter declares every module variable once per evaluation, so a
        # straight-line scope is all they need.
        for name, decl in self._block_locals.items():
            if name not in module.variables:
                module.variables[name] = decl

    def parse_statement(self) -> Optional[StatementNode]:
        tok = self.peek()
        if tok == 'if':
            return self.parse_conditional()
        if tok == 'case':
            return self.parse_case()
        if tok in ('for', 'while'):
            return self.parse_loop()
        if tok == '@':
            self.take('@')
            self.take('(')
            event_name = self.take()
            self.take(')')
            body = self.parse_statement_block()
            return EventNode(event_name, body)
        if tok in ('real', 'integer'):
            # Block-scoped real/integer local declaration (VBIC declares
            # temporaries inside begin : label scopes). These carry no
            # statement semantics; the names are collected and hoisted to
            # module scope by parse_analog_block.
            vtype = self.take()
            while self.peek() not in (';', None):
                vname = self.take()
                if vname != ',':
                    self._block_locals[vname] = VariableDecl(vname, vtype)
            if self.peek() == ';':
                self.take(';')
            return None
        if tok in NOOP_SYSTEM_TASKS and self.token_idx + 1 < len(self.tokens) and self.tokens[self.token_idx + 1] == '(':
            # $strobe/$display/$discontinuity/... are diagnostics or solver
            # hints; GMC compiles them to nothing.
            self.skip_paren_call()
            return None
        if tok == 'begin':
            # A bare begin...end block with no leading control keyword is a
            # plain grouping/scope construct: either "begin : label ... end"
            # (VBIC's evaluateStatic/loadStatic blocks) or bare "begin ... end"
            # nested inside an if body (BSIM "if (x) begin begin ... end end").
            # It carries no control flow, so it is preserved as a transparent
            # GroupNode and the emitter flattens its statements in program order.
            self.take('begin')
            label = ""
            if self.peek() == ':':
                self.take(':')
                label = self.take()
            body = []
            while self.peek() not in ('end', None):
                stmt = self.parse_statement()
                if stmt:
                    body.append(stmt)
            if self.peek() == 'end':
                self.take('end')
            return GroupNode(body, label)
        if tok in ('V', 'I') and self.token_idx + 1 < len(self.tokens) and self.tokens[self.token_idx + 1] == '(':
            quantity, n1, n2 = self.parse_branch_head()
            if self.peek() == '<+':
                self.take('<+')
                expr = self.parse_expression()
                if self.peek() == ';':
                    self.take(';')
                return BranchContributionNode(quantity, n1, n2, expr)
        elif self.token_idx + 1 < len(self.tokens) and self.tokens[self.token_idx + 1] == '=':
            var_name = self.take()
            self.take('=')
            expr = self.parse_expression()
            if self.peek() == ';':
                self.take(';')
            return AssignmentNode(var_name, expr)
        self.token_idx += 1
        return None

    def skip_paren_call(self):
        # Consume "name ( balanced-tokens )" and an optional trailing ';'.
        self.take()
        if self.peek() == '(':
            self.take('(')
            depth = 1
            while depth > 0:
                t = self.take()
                if t == '(':
                    depth += 1
                elif t == ')':
                    depth -= 1
        if self.peek() == ';':
            self.take(';')

    def parse_loop(self) -> StatementNode:
        keyword = self.take()
        self.take('(')
        if keyword == 'for':
            init = None
            if self.peek() != ';':
                var_name = self.take()
                self.take('=')
                init = AssignmentNode(var_name, self.parse_expression())
            self.take(';')
            condition = None
            if self.peek() != ';':
                condition = self.parse_expression()
            self.take(';')
            increment = None
            if self.peek() != ')':
                var_name = self.take()
                self.take('=')
                increment = AssignmentNode(var_name, self.parse_expression())
            self.take(')')
            body = self.parse_statement_block()
            return ForNode(init, condition, increment, body)
        condition = self.parse_expression()
        self.take(')')
        body = self.parse_statement_block()
        return WhileNode(condition, body)

    def parse_conditional(self) -> ConditionalNode:
        self.take('if')
        self.take('(')
        condition = self.parse_expression()
        self.take(')')
        then_body = self.parse_statement_block()
        else_body = []
        if self.peek() == 'else':
            self.take('else')
            else_body = self.parse_statement_block()
        return ConditionalNode(condition, then_body, else_body)

    def parse_case(self) -> CaseNode:
        self.take('case')
        self.take('(')
        expr = self.parse_expression()
        self.take(')')
        items = []
        while self.peek() not in ('endcase', None):
            values = []
            if self.peek() == 'default':
                self.take('default')
            else:
                while self.peek() not in (':', None):
                    values.append(self.parse_expression())
                    if self.peek() == ',':
                        self.take(',')
            self.take(':')
            body = self.parse_statement_block()
            items.append(CaseItem(values, body))
        if self.peek() == 'endcase':
            self.take('endcase')
        return CaseNode(expr, items)

    def parse_statement_block(self) -> List[StatementNode]:
        body = []
        if self.peek() == 'begin':
            self.take('begin')
            while self.peek() not in ('end', None):
                stmt = self.parse_statement()
                if stmt:
                    body.append(stmt)
            self.take('end')
        else:
            stmt = self.parse_statement()
            if stmt:
                body.append(stmt)
        return body

    def parse_branch_head(self):
        quantity = self.take()
        self.take('(')
        n1 = self.take()
        n2 = None
        if self.peek() == ',':
            self.take(',')
            n2 = self.take()
        self.take(')')
        return quantity, n1, n2

    def parse_expression(self) -> ExpressionNode:
        expr = self.parse_logical_or()
        if self.peek() == '?':
            self.take('?')
            then_expr = self.parse_expression()
            self.take(':')
            else_expr = self.parse_expression()
            expr = TernaryNode(expr, then_expr, else_expr)
        return expr

    def parse_logical_or(self) -> ExpressionNode:
        left = self.parse_logical_and()
        while self.peek() in ('||', 'or'):
            self.take()
            left = BinaryOpNode('||', left, self.parse_logical_and())
        return left

    def parse_logical_and(self) -> ExpressionNode:
        left = self.parse_bitor()
        while self.peek() in ('&&', 'and'):
            self.take()
            left = BinaryOpNode('&&', left, self.parse_bitor())
        return left

    def parse_bitor(self) -> ExpressionNode:
        left = self.parse_bitand()
        while self.peek() == '|':
            self.take()
            left = BinaryOpNode('|', left, self.parse_bitand())
        return left

    def parse_bitand(self) -> ExpressionNode:
        left = self.parse_equality()
        while self.peek() == '&':
            self.take()
            left = BinaryOpNode('&', left, self.parse_equality())
        return left

    def parse_equality(self) -> ExpressionNode:
        left = self.parse_relational()
        while self.peek() in ('==', '!='):
            op = self.take()
            left = BinaryOpNode(op, left, self.parse_relational())
        return left

    def parse_relational(self) -> ExpressionNode:
        left = self.parse_additive()
        while self.peek() in ('<=', '>=', '<', '>'):
            op = self.take()
            left = BinaryOpNode(op, left, self.parse_additive())
        return left

    def parse_additive(self) -> ExpressionNode:
        left = self.parse_multiplicative()
        while self.peek() in ('+', '-'):
            op = self.take()
            left = BinaryOpNode(op, left, self.parse_multiplicative())
        return left

    def parse_multiplicative(self) -> ExpressionNode:
        left = self.parse_unary()
        while self.peek() in ('*', '/', '**', '%'):
            op = self.take()
            left = BinaryOpNode(op, left, self.parse_unary())
        return left

    def parse_unary(self) -> ExpressionNode:
        if self.peek() in ('+', '-', '!'):
            op = self.take()
            expr = self.parse_unary()
            if op == '+':
                return expr
            return UnaryOpNode(op, expr)
        return self.parse_primary()

    def parse_primary(self) -> ExpressionNode:
        tok = self.peek()
        if tok in ('V', 'I') and self.token_idx + 1 < len(self.tokens) and self.tokens[self.token_idx + 1] == '(':
            quantity, n1, n2 = self.parse_branch_head()
            return BranchProbeNode(quantity, n1, n2)
        if tok == '(':
            self.take('(')
            expr = self.parse_expression()
            self.take(')')
            return expr
        if tok and re.match(r'(?:\d|\.)', tok):
            self.take()
            return NumberNode(self.parse_number(tok))
        name = self.take()
        if self.peek() == '(':
            self.take('(')
            args = []
            while self.peek() != ')':
                args.append(self.parse_expression())
                if self.peek() == ',':
                    self.take(',')
            self.take(')')
            return FunctionCallNode(name, args)
        return IdentifierNode(name)

    def parse_const_number(self) -> float:
        sign = 1.0
        if self.peek() == '-':
            self.take('-')
            sign = -1.0
        tok = self.take()
        return sign * self.parse_number(tok)

    def parse_number(self, val_str: str) -> float:
        units = {'u': 1e-6, 'n': 1e-9, 'p': 1e-12, 'f': 1e-15, 'k': 1e3, 'm': 1e-3, 'g': 1e9}
        if val_str[-1].lower() in units:
            mult = units[val_str[-1].lower()]
            return float(val_str[:-1]) * mult
        return float(val_str)

if __name__ == "__main__":
    sample_va = """
    module diode(anode, cathode);
        electrical anode, cathode;
        parameter real Is = 1e-14;
        parameter real N = 1.0;
        real vd, id;
        analog begin
            vd = V(anode, cathode);
            id = Is * (exp(vd / (N * 0.02585)) - 1.0);
            I(anode, cathode) <+ id;
        end
    endmodule
    """
    parser = GMCParser(sample_va)
    mod = parser.parse_module()
    print("[+] GMC Parsed Module:", mod)
    print("    Parameters:", mod.parameters)
    print("    Analog Body:", mod.analog_body)
