"""Desktop shell for CMM. This is the only package that depends on Qt."""

_DESKTOP_MODULES = frozenset({"matplotlib", "qtpy", "PyQt5"})

try:
    from cmm.app.main_window import CmmMainWindow
except ModuleNotFoundError as error:  # pragma: no cover - depends on the install extras
    # Name the extra rather than leaving a bare "No module named 'qtpy'": the desktop shell
    # is optional, and every analysis it shows is reachable from the Python API and the CLI.
    if error.name not in _DESKTOP_MODULES:
        raise
    raise ModuleNotFoundError(
        f"cmm.app requires {error.name}, which ships in the 'desktop' extra: install it "
        "with `pip install 'cmm[desktop]'`. The same analyses run headlessly through the "
        "Python API and the `cmm` command without it."
    ) from error

__all__ = ["CmmMainWindow"]
