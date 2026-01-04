import tkinter as tk
from tkinter import ttk
from typing import Sequence


class ShearMomentView(ttk.Frame):
    def __init__(self, master: tk.Widget, width: int = 800, height: int = 320):
        super().__init__(master)
        self._last_x = []
        self._last_shear = []
        self._last_moment = []
        self.shear_canvas = tk.Canvas(
            self, width=width, height=height // 2, bg="white", highlightthickness=1, highlightbackground="#ccc"
        )
        self.moment_canvas = tk.Canvas(
            self, width=width, height=height // 2, bg="white", highlightthickness=1, highlightbackground="#ccc"
        )
        ttk.Label(self, text="Shear force").pack(anchor="w")
        self.shear_canvas.pack(fill="both", expand=True, pady=(0, 4))
        ttk.Label(self, text="Bending moment").pack(anchor="w")
        self.moment_canvas.pack(fill="both", expand=True)

        self.shear_canvas.bind("<Configure>", lambda _e: self._redraw())
        self.moment_canvas.bind("<Configure>", lambda _e: self._redraw())

    def _draw_curve(self, canvas: tk.Canvas, x: Sequence[float], y: Sequence[float], color: str):
        canvas.delete("all")
        if len(x) == 0 or y is None:
            return
        w = canvas.winfo_width() or int(canvas["width"])
        h = canvas.winfo_height() or int(canvas["height"])
        if w <= 0 or h <= 0:
            return
        margin = 30
        x0, x1 = margin, w - margin
        y0, y1 = margin / 2, h - margin / 2
        canvas.create_line(x0, (y0 + y1) / 2, x1, (y0 + y1) / 2, fill="#888", width=2)
        max_abs = max(max(y), -min(y), 1e-9)

        def to_canvas(ix):
            cx = x0 + (x[ix] - x[0]) / max(x[-1] - x[0], 1e-9) * (x1 - x0)
            cy = (y0 + y1) / 2 - (y[ix] / max_abs) * (y1 - y0) / 2
            return cx, cy

        pts = []
        for i in range(len(x)):
            pts.extend(to_canvas(i))
        canvas.create_line(*pts, fill=color, width=3)

    def update_view(self, x: Sequence[float], shear: Sequence[float], moment: Sequence[float]):
        self._last_x = list(x) if x is not None else []
        self._last_shear = list(shear) if shear is not None else None
        self._last_moment = list(moment) if moment is not None else None
        self._redraw()

    def _redraw(self):
        if not self._last_x:
            self.shear_canvas.delete("all")
            self.moment_canvas.delete("all")
            return
        if self._last_shear is None or self._last_moment is None:
            for canvas in (self.shear_canvas, self.moment_canvas):
                canvas.delete("all")
                canvas.create_text(
                    canvas.winfo_width() / 2,
                    canvas.winfo_height() / 2,
                    text="No fixed points. Nonphysical results suppressed.",
                    fill="#444",
                )
            return

        self._draw_curve(self.shear_canvas, self._last_x, self._last_shear, "#2c3e50")
        self._draw_curve(self.moment_canvas, self._last_x, self._last_moment, "#8e44ad")

