"""
GSPICE Model Compiler (GMC) CLI Tool
Usage: python tools/gmc/gmc_compile.py --input <model.gmc> --output <include/devices/gmc_model.hpp>
"""

import sys
import argparse
import json
import hashlib
from pathlib import Path
from gmc_parser import GMCParser
from gmc_cpp_emitter import GMCCppEmitter
from gmc_ast import (
    AssignmentNode, BinaryOpNode, BranchContributionNode, BranchProbeNode,
    CaseNode, ConditionalNode, EventNode, ForNode, FunctionCallNode,
    GroupNode, IdentifierNode, NumberNode, TernaryNode, UnaryOpNode, WhileNode
)

def expr_to_dict(expr):
    if isinstance(expr, NumberNode):
        return {"kind": "number", "value": expr.value}
    if isinstance(expr, IdentifierNode):
        return {"kind": "identifier", "name": expr.name}
    if isinstance(expr, UnaryOpNode):
        return {"kind": "unary", "op": expr.op, "expr": expr_to_dict(expr.expr)}
    if isinstance(expr, BinaryOpNode):
        return {"kind": "binary", "op": expr.op, "left": expr_to_dict(expr.left), "right": expr_to_dict(expr.right)}
    if isinstance(expr, FunctionCallNode):
        return {"kind": "call", "name": expr.func_name, "args": [expr_to_dict(arg) for arg in expr.args]}
    if isinstance(expr, BranchProbeNode):
        return {"kind": "branch_probe", "quantity": expr.quantity, "node1": expr.node1, "node2": expr.node2}
    if isinstance(expr, TernaryNode):
        return {
            "kind": "ternary",
            "condition": expr_to_dict(expr.condition),
            "then": expr_to_dict(expr.then_expr),
            "else": expr_to_dict(expr.else_expr),
        }
    raise TypeError(f"Unsupported expression node: {type(expr).__name__}")

def stmt_to_dict(stmt):
    if isinstance(stmt, AssignmentNode):
        return {"kind": "assignment", "name": stmt.var_name, "expr": expr_to_dict(stmt.expr)}
    if isinstance(stmt, BranchContributionNode):
        return {
            "kind": "branch_contribution",
            "quantity": stmt.quantity,
            "node1": stmt.node1,
            "node2": stmt.node2,
            "expr": expr_to_dict(stmt.expr),
        }
    if isinstance(stmt, ConditionalNode):
        return {
            "kind": "conditional",
            "condition": expr_to_dict(stmt.condition),
            "then": [stmt_to_dict(s) for s in stmt.then_body],
            "else": [stmt_to_dict(s) for s in stmt.else_body],
        }
    if isinstance(stmt, CaseNode):
        return {
            "kind": "case",
            "expr": expr_to_dict(stmt.expr),
            "items": [
                {
                    "values": [expr_to_dict(v) for v in item.values],
                    "body": [stmt_to_dict(s) for s in item.body],
                }
                for item in stmt.items
            ],
        }
    if isinstance(stmt, ForNode):
        return {
            "kind": "for",
            "init": stmt_to_dict(stmt.init) if stmt.init else None,
            "condition": expr_to_dict(stmt.condition) if stmt.condition else None,
            "increment": stmt_to_dict(stmt.increment) if stmt.increment else None,
            "body": [stmt_to_dict(s) for s in stmt.body],
        }
    if isinstance(stmt, WhileNode):
        return {
            "kind": "while",
            "condition": expr_to_dict(stmt.condition),
            "body": [stmt_to_dict(s) for s in stmt.body],
        }
    if isinstance(stmt, EventNode):
        return {
            "kind": "event",
            "event": stmt.event_name,
            "body": [stmt_to_dict(s) for s in stmt.body],
        }
    if isinstance(stmt, GroupNode):
        return {
            "kind": "group",
            "label": stmt.label,
            "body": [stmt_to_dict(s) for s in stmt.body],
        }
    raise TypeError(f"Unsupported statement node: {type(stmt).__name__}")

def module_to_dict(module):
    return {
        "schema": "gmc-ir-v0",
        "name": module.name,
        "ports": module.ports,
        "terminal_count": module.terminal_count,
        "collapsible_pairs": dict(module.collapsible_pairs),
        "parameters": [
            {"name": p.name, "type": p.ptype, "default": p.default_val,
             "default_expr": expr_to_dict(p.default_expr) if p.default_expr else None}
            for p in module.parameters.values()
        ],
        "variables": [
            {"name": v.name, "type": v.vtype}
            for v in module.variables.values()
        ],
        "functions": [
            {
                "name": f.name,
                "return_type": f.return_type,
                "args": [{"name": a.name, "type": a.vtype} for a in f.args],
                "locals": [{"name": v.name, "type": v.vtype} for v in f.locals],
                "body": [stmt_to_dict(stmt) for stmt in f.body],
            }
            for f in module.functions.values()
        ],
        "analog": [stmt_to_dict(stmt) for stmt in module.analog_body],
    }

def build_gsdi_artifact(module, source_path, header_path):
    emitter = GMCCppEmitter(module)
    return {
        "schema": "gsdi-artifact-v0",
        "producer": "gmc",
        "model_type": module.name,
        "source": str(source_path),
        "source_sha256": hashlib.sha256(source_path.read_bytes()).hexdigest(),
        "generated_header": str(header_path) if header_path else "",
        "ports": module.ports,
        "terminal_count": module.terminal_count or len(module.ports),
        "collapsible_pairs": dict(module.collapsible_pairs),
        "parameters": [
            {"name": p.name, "type": p.ptype, "default": p.default_val}
            for p in module.parameters.values()
        ],
        "features": {
            "op": True,
            "tran": emitter.has_ddt(module.analog_body),
            "ac": emitter.has_ddt(module.analog_body),
            "noise": False,
        },
    }

def main():
    parser = argparse.ArgumentParser(description="GSPICE Model Compiler (GMC)")
    parser.add_argument("--input", required=True, help="Path to input GMC model source (.gmc)")
    parser.add_argument("--output", help="Path to output compiled C++ header (.hpp)")
    parser.add_argument("--gsdi-output",
                        help="Path to output .gsdi artifact manifest (default: output with .gsdi suffix)")
    parser.add_argument("--include", action="append", default=[],
                        help="Additional Verilog-A include directory (repeatable)")
    parser.add_argument("--dump-ir", action="store_true", help="Print parsed GMC IR JSON and exit")
    args = parser.parse_args()

    input_path = Path(args.input)

    if not input_path.exists():
        print(f"[-] Input file not found: {input_path}")
        sys.exit(1)

    target_name = Path(args.output).name if args.output else "IR"
    print(f"[*] GMC compiling native model: {input_path.name} -> {target_name}", file=sys.stderr)
    va_code = input_path.read_text(encoding="utf-8")

    va_parser = GMCParser(
        va_code,
        source_name=str(input_path),
        include_dirs=[str(input_path.parent)] + args.include,
    )
    for dep in va_parser.preprocessed_files:
        print(f"[*] GMC dependency: {dep}", file=sys.stderr)
    module_ast = va_parser.parse_module()

    if args.dump_ir:
        print(json.dumps(module_to_dict(module_ast), indent=2, sort_keys=True))
        return

    if not args.output:
        print("[-] --output is required unless --dump-ir is used")
        sys.exit(1)

    output_path = Path(args.output)
    emitter = GMCCppEmitter(module_ast)
    cpp_code = emitter.emit_cpp_header()

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(cpp_code, encoding="utf-8")
    gsdi_path = Path(args.gsdi_output) if args.gsdi_output else output_path.with_suffix(".gsdi")
    gsdi_path.parent.mkdir(parents=True, exist_ok=True)
    gsdi_path.write_text(
        json.dumps(build_gsdi_artifact(module_ast, input_path, output_path), indent=2, sort_keys=True) + "\n",
        encoding="utf-8")
    print(f"[+] GMC successfully compiled {module_ast.name} model header to {output_path}")
    print(f"[+] GMC wrote GSDI artifact manifest to {gsdi_path}")

if __name__ == "__main__":
    main()
