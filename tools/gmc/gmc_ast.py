"""
GSPICE Model Compiler (GMC) - AST Definitions
Defines node classes for GMC native model-source AST construction.
"""

from typing import List, Dict, Any, Optional

class ASTNode:
    pass

class ExpressionNode(ASTNode):
    pass

class NumberNode(ExpressionNode):
    def __init__(self, value: float):
        self.value = value

    def __repr__(self):
        return f"Number({self.value})"

class IdentifierNode(ExpressionNode):
    def __init__(self, name: str):
        self.name = name

    def __repr__(self):
        return f"Var({self.name})"

class BinaryOpNode(ExpressionNode):
    def __init__(self, op: str, left: ExpressionNode, right: ExpressionNode):
        self.op = op
        self.left = left
        self.right = right

    def __repr__(self):
        return f"({self.left} {self.op} {self.right})"

class UnaryOpNode(ExpressionNode):
    def __init__(self, op: str, expr: ExpressionNode):
        self.op = op
        self.expr = expr

    def __repr__(self):
        return f"({self.op}{self.expr})"

class TernaryNode(ExpressionNode):
    def __init__(self, condition: ExpressionNode, then_expr: ExpressionNode, else_expr: ExpressionNode):
        self.condition = condition
        self.then_expr = then_expr
        self.else_expr = else_expr

    def __repr__(self):
        return f"({self.condition} ? {self.then_expr} : {self.else_expr})"

class FunctionCallNode(ExpressionNode):
    def __init__(self, func_name: str, args: List[ExpressionNode]):
        self.func_name = func_name
        self.args = args

    def __repr__(self):
        return f"{self.func_name}({', '.join(map(str, self.args))})"

class BranchProbeNode(ExpressionNode):
    def __init__(self, quantity: str, node1: str, node2: Optional[str] = None):
        self.quantity = quantity  # 'V' or 'I'
        self.node1 = node1
        self.node2 = node2

    def __repr__(self):
        n2 = f", {self.node2}" if self.node2 else ""
        return f"{self.quantity}({self.node1}{n2})"

class StatementNode(ASTNode):
    pass

class GroupNode(StatementNode):
    """A labeled or plain `begin ... end` grouping block. Verilog-A treats
    `begin : name ... end` at statement level as a plain sequential scope
    (no control flow), so the emitted C++ simply flattens its statements
    in program order."""

    def __init__(self, body: List[StatementNode], label: str = ""):
        self.body = body
        self.label = label

    def __repr__(self):
        return f"begin{(' : ' + self.label) if self.label else ''} {{{len(self.body)} stmts}}"

class BranchContributionNode(StatementNode):
    def __init__(self, quantity: str, node1: str, node2: Optional[str], expr: ExpressionNode):
        self.quantity = quantity  # 'V' or 'I'
        self.node1 = node1
        self.node2 = node2
        self.expr = expr

    def __repr__(self):
        n2 = f", {self.node2}" if self.node2 else ""
        return f"{self.quantity}({self.node1}{n2}) <+ {self.expr};"

class AssignmentNode(StatementNode):
    def __init__(self, var_name: str, expr: ExpressionNode):
        self.var_name = var_name
        self.expr = expr

    def __repr__(self):
        return f"{self.var_name} = {self.expr};"

class ConditionalNode(StatementNode):
    def __init__(self, condition: ExpressionNode, then_body: List[StatementNode], else_body: List[StatementNode] = None):
        self.condition = condition
        self.then_body = then_body
        self.else_body = else_body or []

class CaseItem:
    def __init__(self, values: List[ExpressionNode], body: List[StatementNode]):
        self.values = values  # empty list == default item
        self.body = body

    def __repr__(self):
        label = "default" if not self.values else ", ".join(str(v) for v in self.values)
        return f"[{label}: {len(self.body)} stmts]"

class CaseNode(StatementNode):
    def __init__(self, expr: ExpressionNode, items: List[CaseItem]):
        self.expr = expr
        self.items = items

    def __repr__(self):
        return f"case({self.expr}) {self.items}"

class ForNode(StatementNode):
    def __init__(self, init: Optional[AssignmentNode], condition: Optional[ExpressionNode],
                 increment: Optional[AssignmentNode], body: List[StatementNode]):
        self.init = init
        self.condition = condition
        self.increment = increment
        self.body = body

    def __repr__(self):
        return f"for(...) {{{len(self.body)} stmts}}"

class WhileNode(StatementNode):
    def __init__(self, condition: ExpressionNode, body: List[StatementNode]):
        self.condition = condition
        self.body = body

    def __repr__(self):
        return f"while(...) {{{len(self.body)} stmts}}"

class EventNode(StatementNode):
    def __init__(self, event_name: str, body: List[StatementNode]):
        self.event_name = event_name  # e.g. "initial_step"
        self.body = body

    def __repr__(self):
        return f"@({self.event_name}) {{{len(self.body)} stmts}}"

class ParameterDecl:
    def __init__(self, name: str, default_val: float, ptype: str = "real",
                 default_expr: Optional[ExpressionNode] = None):
        self.name = name
        self.default_val = default_val
        self.ptype = ptype
        # Set when the default is an expression (e.g. "parameter real dsub=drout;")
        # instead of a literal; emission fails closed on this.
        self.default_expr = default_expr

class VariableDecl:
    def __init__(self, name: str, vtype: str = "real"):
        self.name = name
        self.vtype = vtype

class FunctionDeclNode(ASTNode):
    def __init__(self, name: str, return_type: str,
                 args: List[VariableDecl], locals_: List[VariableDecl],
                 body: List[StatementNode]):
        self.name = name
        self.return_type = return_type
        self.args = args
        self.locals = locals_
        self.body = body

    def __repr__(self):
        return f"Function({self.name}({', '.join(a.name for a in self.args)}) -> {self.return_type})"

class ModuleNode(ASTNode):
    def __init__(self, name: str, ports: List[str]):
        self.name = name
        self.ports = ports
        # First terminal_count ports are external terminals; the remaining
        # ports are hidden nodes (internal or collapsible). 0 = all ports.
        self.terminal_count: int = 0
        # hidden port name -> external terminal port name it collapses into
        self.collapsible_pairs: Dict[str, str] = {}
        self.parameters: Dict[str, ParameterDecl] = {}
        self.variables: Dict[str, VariableDecl] = {}
        # user-defined analog functions (return type, args, body)
        self.functions: Dict[str, FunctionDeclNode] = {}
        self.analog_body: List[StatementNode] = []

    def __repr__(self):
        return f"VerilogAModule({self.name}, ports={self.ports}, params={list(self.parameters.keys())})"
