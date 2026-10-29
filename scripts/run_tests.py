"""Ejecuta todas las pruebas del proyecto desde el botón Play de VS Code."""

from pathlib import Path
import sys
import unittest


RAIZ_PROYECTO = Path(__file__).resolve().parents[1]
RUTA_SRC = RAIZ_PROYECTO / "src"
RUTA_TESTS = RAIZ_PROYECTO / "tests"

# Permite importar mcplan al ejecutar este archivo directamente desde VS Code.
sys.path.insert(0, str(RUTA_SRC))

suite = unittest.defaultTestLoader.discover(
    start_dir=str(RUTA_TESTS),
    pattern="test_*.py",
)

ejecutor = unittest.TextTestRunner(verbosity=2)
resultado = ejecutor.run(suite)
raise SystemExit(not resultado.wasSuccessful())
