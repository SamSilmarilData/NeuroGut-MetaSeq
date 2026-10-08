#!/usr/bin/env python3
"""
scripts/02b_qc_audit.py
Entrypoint wrapper for 01b_audit_data_landscape.py.
Provides quality control, microglial purity, and ex vivo dissociation stress audits.
"""

import sys
import os
import importlib.util

def main():
    script_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "01b_audit_data_landscape.py")
    spec = importlib.util.spec_from_file_location("audit_data_landscape", script_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    if hasattr(module, "main"):
        module.main()

if __name__ == "__main__":
    main()
