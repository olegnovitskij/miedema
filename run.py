import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent / 'src'))

from miedema.app import MiedemaApp

if __name__ == "__main__":
    MiedemaApp().run()

