import tkinter as tk
from tkinter import ttk
from typing import Sequence


class FrequencyPanel(ttk.Frame):
    def __init__(self, master: tk.Widget):
        super().__init__(master)
        ttk.Label(self, text="Bending [Hz]").grid(row=0, column=0, sticky="w", padx=4, pady=(2, 0))
        ttk.Label(self, text="Torsion [Hz]").grid(row=0, column=1, sticky="w", padx=4, pady=(2, 0))

        self.bend_text = tk.Text(self, height=5, width=30)
        self.tors_text = tk.Text(self, height=5, width=30)
        self.bend_text.grid(row=1, column=0, padx=4, pady=2, sticky="ew")
        self.tors_text.grid(row=1, column=1, padx=4, pady=2, sticky="ew")
        for widget in (self.bend_text, self.tors_text):
            widget.configure(state="disabled")

        stats = ttk.Frame(self)
        stats.grid(row=1, column=2, sticky="nsew", padx=6, pady=2)
        ttk.Label(stats, text="Max shear [N]").grid(row=0, column=0, sticky="w")
        ttk.Label(stats, text="Avg shear [N]").grid(row=1, column=0, sticky="w")
        ttk.Label(stats, text="Max moment [N·m]").grid(row=2, column=0, sticky="w")
        ttk.Label(stats, text="Avg moment [N·m]").grid(row=3, column=0, sticky="w")
        self.max_shear = ttk.Label(stats, text="0")
        self.avg_shear = ttk.Label(stats, text="0")
        self.max_moment = ttk.Label(stats, text="0")
        self.avg_moment = ttk.Label(stats, text="0")
        self.max_shear.grid(row=0, column=1, sticky="e")
        self.avg_shear.grid(row=1, column=1, sticky="e")
        self.max_moment.grid(row=2, column=1, sticky="e")
        self.avg_moment.grid(row=3, column=1, sticky="e")

    def update_values(self, bend: Sequence[float], tors: Sequence[float]):
        self._fill(self.bend_text, bend)
        self._fill(self.tors_text, tors)

    def update_shear_stats(self, shear: Sequence[float], moment: Sequence[float]):
        if shear is None or moment is None or len(shear) == 0 or len(moment) == 0:
            for label in (self.max_shear, self.avg_shear, self.max_moment, self.avg_moment):
                label.config(text="—")
            self.max_moment.config(text="No fixed points")
            return
        max_s = max(max(shear), -min(shear))
        avg_s = sum(abs(v) for v in shear) / len(shear)
        max_m = max(max(moment), -min(moment))
        avg_m = sum(abs(v) for v in moment) / len(moment)
        self.max_shear.config(text=f"{max_s:.2f}")
        self.avg_shear.config(text=f"{avg_s:.2f}")
        self.max_moment.config(text=f"{max_m:.2f}")
        self.avg_moment.config(text=f"{avg_m:.2f}")

    def _fill(self, widget: tk.Text, data: Sequence[float]):
        widget.configure(state="normal")
        widget.delete("1.0", tk.END)
        for i, val in enumerate(data, 1):
            widget.insert(tk.END, f"{i}: {val:.3f}\n")
        widget.configure(state="disabled")

