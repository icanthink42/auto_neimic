import tkinter as tk
from tkinter import ttk
from typing import Callable


class BoundaryForm(ttk.LabelFrame):
    def __init__(self, master: tk.Widget, on_change: Callable[[str, str], None]):
        super().__init__(master, text="Boundary")
        self.on_change = on_change
        self.left_var = tk.StringVar(value="fixed")
        self.right_var = tk.StringVar(value="free")

        ttk.Label(self, text="Left").grid(row=0, column=0, padx=4, pady=2, sticky="w")
        ttk.Radiobutton(self, text="Fixed", variable=self.left_var, value="fixed", command=self._emit).grid(
            row=0, column=1, padx=2, pady=2, sticky="w"
        )
        ttk.Radiobutton(self, text="Free", variable=self.left_var, value="free", command=self._emit).grid(
            row=0, column=2, padx=2, pady=2, sticky="w"
        )

        ttk.Label(self, text="Right").grid(row=1, column=0, padx=4, pady=2, sticky="w")
        ttk.Radiobutton(self, text="Fixed", variable=self.right_var, value="fixed", command=self._emit).grid(
            row=1, column=1, padx=2, pady=2, sticky="w"
        )
        ttk.Radiobutton(self, text="Free", variable=self.right_var, value="free", command=self._emit).grid(
            row=1, column=2, padx=2, pady=2, sticky="w"
        )

    def set_values(self, left: bool, right: bool):
        self.left_var.set("fixed" if left else "free")
        self.right_var.set("fixed" if right else "free")

    def _emit(self):
        self.on_change(self.left_var.get(), self.right_var.get())

