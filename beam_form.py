import tkinter as tk
from tkinter import ttk
from typing import Callable

from si_prefix import prefix_labels, prefix_multiplier


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
        self.prefix_vars = {
            "length": tk.StringVar(value="m"),
            "radius": tk.StringVar(value="m"),
            "elastic_modulus": tk.StringVar(value="Pa"),
            "shear_modulus": tk.StringVar(value="Pa"),
            "density": tk.StringVar(value="kg/m3"),
        }
        labels = [
            ("Length [m]", "length", "m"),
            ("Radius [m]", "radius", "m"),
            ("E [Pa]", "elastic_modulus", "Pa"),
            ("G [Pa]", "shear_modulus", "Pa"),
            ("Density [kg/m3]", "density", "kg/m3"),
            ("Elements", "elements"),
        ]
        for row, (text, key, *unit) in enumerate(labels):
            ttk.Label(self, text=text).grid(row=row, column=0, sticky="w", padx=4, pady=2)
            entry = ttk.Entry(self, textvariable=self.vars[key], width=14)
            entry.grid(row=row, column=1, sticky="ew", padx=4, pady=2)
            entry.bind("<FocusOut>", self._emit_change)
            entry.bind("<Return>", self._emit_change)
            if unit:
                unit_value = unit[0]
                ttk.OptionMenu(
                    self,
                    self.prefix_vars[key],
                    self.prefix_vars[key].get(),
                    *prefix_labels(unit_value),
                    command=lambda _value, field=key: self._emit_change(),
                ).grid(row=row, column=2, sticky="w", padx=4, pady=2)
        self.columnconfigure(1, weight=1)

    def set_values(self, values: tuple):
        keys = ["length", "radius", "elastic_modulus", "shear_modulus", "density", "elements"]
        for k, v in zip(keys, values):
            self.vars[k].set(str(v))
            if k in self.prefix_vars:
                unit = self._unit_for_key(k)
                self.prefix_vars[k].set(unit)

    def get_values(self):
        try:
            return (
                float(self.vars["length"].get()) * self._prefix_multiplier("length"),
                float(self.vars["radius"].get()) * self._prefix_multiplier("radius"),
                float(self.vars["elastic_modulus"].get()) * self._prefix_multiplier("elastic_modulus"),
                float(self.vars["shear_modulus"].get()) * self._prefix_multiplier("shear_modulus"),
                float(self.vars["density"].get()) * self._prefix_multiplier("density"),
                int(float(self.vars["elements"].get())),
            )
        except ValueError:
            return None

    def _emit_change(self, _event=None):
        values = self.get_values()
        if values is None:
            return
        self.on_change(values)

    def _unit_for_key(self, key: str) -> str:
        return {
            "length": "m",
            "radius": "m",
            "elastic_modulus": "Pa",
            "shear_modulus": "Pa",
            "density": "kg/m3",
        }[key]

    def _prefix_multiplier(self, key: str) -> float:
        unit = self._unit_for_key(key)
        return prefix_multiplier(self.prefix_vars[key].get(), unit)
