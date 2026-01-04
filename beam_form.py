import tkinter as tk
from tkinter import ttk
from typing import Callable


class BeamForm(ttk.LabelFrame):
    def __init__(self, master: tk.Widget, on_change: Callable[[tuple], None]):
        super().__init__(master, text="Beam")
        self.on_change = on_change
        self.vars = {
            "length": tk.StringVar(value="2.0"),
            "radius": tk.StringVar(value="0.05"),
            "elastic_modulus": tk.StringVar(value="200e9"),
            "shear_modulus": tk.StringVar(value="79.3e9"),
            "density": tk.StringVar(value="7850"),
            "elements": tk.StringVar(value="20"),
        }
        labels = [
            ("Length [m]", "length"),
            ("Radius [m]", "radius"),
            ("E [Pa]", "elastic_modulus"),
            ("G [Pa]", "shear_modulus"),
            ("Density [kg/m3]", "density"),
            ("Elements", "elements"),
        ]
        for row, (text, key) in enumerate(labels):
            ttk.Label(self, text=text).grid(row=row, column=0, sticky="w", padx=4, pady=2)
            entry = ttk.Entry(self, textvariable=self.vars[key], width=14)
            entry.grid(row=row, column=1, sticky="ew", padx=4, pady=2)
            entry.bind("<FocusOut>", self._emit_change)
            entry.bind("<Return>", self._emit_change)
        self.columnconfigure(1, weight=1)

    def set_values(self, values: tuple):
        keys = ["length", "radius", "elastic_modulus", "shear_modulus", "density", "elements"]
        for k, v in zip(keys, values):
            self.vars[k].set(str(v))

    def get_values(self):
        try:
            return (
                float(self.vars["length"].get()),
                float(self.vars["radius"].get()),
                float(self.vars["elastic_modulus"].get()),
                float(self.vars["shear_modulus"].get()),
                float(self.vars["density"].get()),
                int(float(self.vars["elements"].get())),
            )
        except ValueError:
            return None

    def _emit_change(self, _event=None):
        values = self.get_values()
        if values is None:
            return
        self.on_change(values)

