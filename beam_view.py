import tkinter as tk
from tkinter import ttk

from constraint import Constraint
from gui_state import BeamUIState
from point_mass import PointMass
from torsional_spring import TorsionalSpring
from translational_spring import TranslationalSpring


class BeamView(ttk.LabelFrame):
    def __init__(self, master: tk.Widget, state: BeamUIState, on_change):
        super().__init__(master, text="Beam")
        self.state = state
        self.on_change = on_change
        self.canvas = tk.Canvas(
            self, width=800, height=220, bg="white", highlightthickness=1, highlightbackground="#ccc"
        )
        self.canvas.pack(fill="both", expand=True)
        self.margin = 40

        self.canvas.bind("<Configure>", lambda _e: self.update_view())
        self.canvas.bind("<Button-3>", self._on_right_click)

    def _dims(self):
        w = max(300, self.canvas.winfo_width())
        h = max(220, self.canvas.winfo_height())
        return w, h

    def _to_canvas_x(self, x: float, length: float, x0: float, x1: float) -> float:
        x = max(0.0, min(length, x))
        return x0 + (x / max(length, 1e-9)) * (x1 - x0)

    def _from_canvas_x(self, cx: float, length: float, x0: float, x1: float) -> float:
        cx = max(x0, min(x1, cx))
        return (cx - x0) / max((x1 - x0), 1e-9) * length

    def update_view(self):
        self.canvas.delete("all")
        length = self.state.length
        if length <= 0:
            return
        w, h = self._dims()
        m = self.margin
        x0, x1 = m, w - m
        y = h // 2

        model = self.state.to_model()
        base_radius = max(self.state.radius, 1e-6)
        base_thickness = 10.0
        min_thickness = 4.0
        samples = max(40, model.elements * 2)
        top_points = []
        bottom_points = []
        for i in range(samples + 1):
            x = length * i / samples
            radius = model.radius_at(x)
            thickness = max(min_thickness, base_thickness * (radius / base_radius))
            half = thickness / 2
            cx = self._to_canvas_x(x, length, x0, x1)
            top_points.append((cx, y - half))
            bottom_points.append((cx, y + half))
        beam_points = top_points + bottom_points[::-1]
        self.canvas.create_polygon(beam_points, fill="#666", outline="#555", width=1, tags=("beam",))
        self.canvas.create_text(x0, y + 18, text="0 m", anchor="w", fill="#444", tags=("beam",))
        self.canvas.create_text(x1, y + 18, text=f"{length:g} m", anchor="e", fill="#444", tags=("beam",))

        for i, pm in enumerate(self.state.point_masses):
            cx = self._to_canvas_x(pm.position, length, x0, x1)
            r = 9
            tags = ("mass", f"mass-{i}")
            self.canvas.create_oval(cx - r, y - r, cx + r, y + r, fill="#e74c3c", outline="", tags=tags)
            self.canvas.create_text(cx, y - 16, text=f"m={pm.mass:g}", fill="#c0392b", tags=tags)

        for i, sp in enumerate(self.state.trans_springs):
            cx = self._to_canvas_x(sp.position, length, x0, x1)
            tags = ("spring", "trans", f"trans-{i}")
            self.canvas.create_rectangle(cx - 8, y - 18, cx + 8, y + 18, outline="#27ae60", width=2, tags=tags)
            self.canvas.create_text(cx, y - 24, text=f"k={sp.k:g}", fill="#27ae60", tags=tags)

        for i, sp in enumerate(self.state.tors_springs):
            cx = self._to_canvas_x(sp.position, length, x0, x1)
            tags = ("spring", "tors", f"tors-{i}")
            self.canvas.create_arc(
                cx - 12, y - 20, cx + 12, y + 20, start=0, extent=300, style="arc", outline="#8e44ad", width=2, tags=tags
            )
            self.canvas.create_text(cx, y + 24, text=f"k={sp.k:g}", fill="#8e44ad", tags=tags)

        for i, constraint in enumerate(self.state.constraints):
            cx = self._to_canvas_x(constraint.position, length, x0, x1)
            tags = ("constraint", f"constraint-{i}")
            # Draw triangle support
            if constraint.fix_translation:
                pts = [cx, y + 5, cx - 10, y + 20, cx + 10, y + 20]
                self.canvas.create_polygon(pts, fill="#f39c12", outline="#e67e22", width=2, tags=tags)
                # Draw ground lines
                for j in range(3):
                    self.canvas.create_line(cx - 12 + j * 8, y + 20, cx - 16 + j * 8, y + 26, fill="#e67e22", width=2, tags=tags)
            if constraint.fix_rotation:
                # Draw a small clamp rectangle to indicate rotation is fixed
                self.canvas.create_rectangle(cx - 6, y - 8, cx + 6, y + 8, fill="#d35400", outline="#e67e22", width=2, tags=tags)
            label_parts = []
            if constraint.fix_translation:
                label_parts.append("T")
            if constraint.fix_rotation:
                label_parts.append("R")
            label = "+".join(label_parts) if label_parts else "?"
            self.canvas.create_text(cx, y + 32, text=label, fill="#e67e22", tags=tags)

    def _on_right_click(self, event):
        item = self.canvas.find_withtag("current")
        if not item:
            return
        tags = self.canvas.gettags(item[0])
        for tag in tags:
            if tag.startswith("mass-"):
                idx = int(tag.split("-")[1])
                if 0 <= idx < len(self.state.point_masses):
                    self.state.point_masses.pop(idx)
                    self.on_change()
                    self.update_view()
                    return
            if tag.startswith("trans-"):
                idx = int(tag.split("-")[1])
                if 0 <= idx < len(self.state.trans_springs):
                    self.state.trans_springs.pop(idx)
                    self.on_change()
                    self.update_view()
                    return
            if tag.startswith("tors-"):
                idx = int(tag.split("-")[1])
                if 0 <= idx < len(self.state.tors_springs):
                    self.state.tors_springs.pop(idx)
                    self.on_change()
                    self.update_view()
                    return
            if tag.startswith("constraint-"):
                idx = int(tag.split("-")[1])
                if 0 <= idx < len(self.state.constraints):
                    self.state.constraints.pop(idx)
                    self.on_change()
                    self.update_view()
                    return
