"""Ejecuta los planes de exploración y deja sus tablas en el workspace Codex."""
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "configs/exploracion_planes_central.jsonc"
SALIDA = Path.home() / "Documents/Codex/2026-10-15/h/exploracion_central_M1_SS"


def main() -> int:
    comando = [
        sys.executable, "-u", str(ROOT / "scripts/evaluar_dossier_final.py"),
        str(CONFIG), "--salida", str(SALIDA),
    ]
    return subprocess.run(comando, cwd=ROOT, check=False).returncode


if __name__ == "__main__":
    raise SystemExit(main())
