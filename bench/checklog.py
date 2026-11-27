"""Pass/fail bookkeeping for the check-style test scripts (boundary_cases, cli_smoke, edge_cases, spill_check, validate_tags)."""


class CheckLog:
    """Prints one PASS or FAIL line per check and keeps the names of the failed checks in `failed`.

    cli_smoke prints a check's detail only when it fails, so its output stays one short line per passing check.
    """

    def __init__(self, detail_on_pass=True):
        self.failed = []
        self.detail_on_pass = detail_on_pass

    def __call__(self, name, passed, detail=""):
        shown = detail and (self.detail_on_pass or not passed)
        print(f"  [{'PASS' if passed else 'FAIL'}] {name}" + (f"  {detail}" if shown else ""), flush=True)
        if not passed:
            self.failed.append(name)
