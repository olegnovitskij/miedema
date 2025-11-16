from pathlib import Path
from kivy.app import App
from kivy.lang import Builder
from kivy.uix.screenmanager import ScreenManager

from miedema.ui.screens import MainWindow, VisualizationWindow


class WindowManager(ScreenManager):
    pass


class MiedemaApp(App):
    def build(self):
        project_root = Path(__file__).parent.parent.parent
        kv_file = project_root / "assets" / "ui.kv"
        kv = Builder.load_file(str(kv_file))
        return kv

