import tkinter as tk
from tkinter import ttk, simpledialog

from beam_form import BeamForm
from boundary_form import BoundaryForm
from beam_view import BeamView
from constraint import Constraint
from cross_section import CrossSection
from fem import natural_frequencies
from frequency_panel import FrequencyPanel
from gui_state import BeamUIState
from point_mass import PointMass
from torsional_spring import TorsionalSpring
from translational_spring import TranslationalSpring
from shear_moment import shear_moment
from shear_moment_view import ShearMomentView


class BeamApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Auto Niemiec Beam")
        self.state = BeamUIState()
        self.status_var = tk.StringVar(value="Ready")

        toolbar = ttk.Frame(self, padding=6)
        toolbar.grid(row=0, column=0, sticky="ew")
        ttk.Button(toolbar, text="Beam Params", command=self._open_beam_params).pack(side="left", padx=4)
        ttk.Button(toolbar, text="Cross Sections", command=self._open_cross_section_dialog).pack(side="left", padx=4)
        ttk.Button(toolbar, text="Add Mass", command=self._add_mass_dialog).pack(side="left", padx=4)
        ttk.Button(toolbar, text="Add Spring", command=self._add_spring_dialog).pack(side="left", padx=4)
        ttk.Button(toolbar, text="Add Constraint", command=self._add_constraint_dialog).pack(side="left", padx=4)
        ttk.Button(toolbar, text="Boundary Conds", command=self._open_boundary_dialog).pack(side="left", padx=4)
        ttk.Label(toolbar, textvariable=self.status_var).pack(side="right")

        main = ttk.Frame(self, padding=8)
        main.grid(row=1, column=0, sticky="nsew")
        self.rowconfigure(1, weight=1)
        self.columnconfigure(0, weight=1)
        main.columnconfigure(0, weight=1)
        main.rowconfigure(0, weight=0)
        main.rowconfigure(1, weight=1)
        main.rowconfigure(2, weight=1)
        main.rowconfigure(3, weight=0)

        top_bar = ttk.Frame(main)
        top_bar.grid(row=0, column=0, sticky="ew", padx=4, pady=2)
        self.boundary_form = BoundaryForm(top_bar, on_change=self._on_boundary_change)
        self.boundary_form.pack(side="left", padx=4)

        self.beam_view = BeamView(main, self.state, on_change=self._refresh_model)
        self.beam_view.grid(row=1, column=0, sticky="nsew", padx=4, pady=4)

        self.shear_view = ShearMomentView(main)
        self.shear_view.grid(row=2, column=0, sticky="nsew", padx=4, pady=4)

        self.freq_panel = FrequencyPanel(main)
        self.freq_panel.grid(row=3, column=0, sticky="ew", padx=4, pady=4)

        self.after(50, self._refresh_model)

    def _add_mass_dialog(self):
        length = self.state.length
        pos = simpledialog.askfloat("Position", f"x along beam [0, {length} m]:", parent=self, initialvalue=length / 2)
        if pos is None:
            return
        pos = max(0.0, min(length, pos))
        mass = simpledialog.askfloat("Add mass", "Mass [kg]:", parent=self)
        if mass is None:
            return
        inertia = simpledialog.askfloat("Add mass", "Rotary inertia J [kg*m^2]:", parent=self, initialvalue=0.0)
        if inertia is None:
            inertia = 0.0
        self.state.point_masses.append(PointMass(position=pos, mass=mass, rotary_inertia=inertia))
        self.status_var.set("Mass added")
        self._refresh_model()

    def _add_spring_dialog(self):
        length = self.state.length
        pos = simpledialog.askfloat("Position", f"x along beam [0, {length} m]:", parent=self, initialvalue=length / 2)
        if pos is None:
            return
        pos = max(0.0, min(length, pos))
        k = simpledialog.askfloat("Add spring", "Stiffness k:", parent=self, initialvalue=1e5)
        if k is None:
            return

        win = tk.Toplevel(self)
        win.title("Spring type")
        choice = tk.StringVar(value="trans")
        ttk.Radiobutton(win, text="Translational", variable=choice, value="trans").pack(anchor="w", padx=8, pady=4)
        ttk.Radiobutton(win, text="Torsional", variable=choice, value="tors").pack(anchor="w", padx=8, pady=4)
        ttk.Button(win, text="OK", command=win.destroy).pack(pady=6)
        win.transient(self)
        win.grab_set()
        self.wait_window(win)

        spring_type = choice.get()

        if spring_type == "tors":
            self.state.tors_springs.append(TorsionalSpring(position=pos, k=k))
        else:
            self.state.trans_springs.append(TranslationalSpring(position=pos, k=k))
        self.status_var.set("Spring added")
        self._refresh_model()

    def _open_beam_params(self):
        win = tk.Toplevel(self)
        win.title("Beam parameters")
        form = BeamForm(win, on_change=lambda _vals: None)
        form.pack(fill="both", expand=True, padx=8, pady=8)
        form.set_values(
            (
                self.state.length,
                self.state.radius,
                self.state.elastic_modulus,
                self.state.shear_modulus,
                self.state.density,
                self.state.elements,
            )
        )
        btns = ttk.Frame(win)
        btns.pack(fill="x", padx=8, pady=4)
        ttk.Button(
            btns,
            text="OK",
            command=lambda: self._apply_beam_params_from_form(form, win),
        ).pack(side="right", padx=4)
        ttk.Button(btns, text="Cancel", command=win.destroy).pack(side="right")

    def _apply_beam_params_from_form(self, form: BeamForm, win: tk.Toplevel):
        params = form.get_values()
        if params is None:
            self.status_var.set("Invalid beam parameters")
            return
        win.destroy()
        self.state.set_beam(params)
        self.status_var.set("Beam updated")
        self._refresh_model()

    def _open_cross_section_dialog(self):
        dialog = CrossSectionDialog(self, self.state.length, self.state.cross_sections)
        self.wait_window(dialog)
        if dialog.result is None:
            return
        self.state.cross_sections = dialog.result
        self.status_var.set("Cross sections updated")
        self._refresh_model()

    def _refresh_model(self):
        try:
            model = self.state.to_model()
            bending_fixed, torsion_fixed = self._build_bc(model)
            if not bending_fixed and not torsion_fixed:
                bend, tors = [], []
            else:
                bend, tors = natural_frequencies(
                    model, bending_fixed=bending_fixed, torsion_fixed=torsion_fixed, n_modes=5
                )
            x, shear, moment = shear_moment(
                model, left_fixed=self.state.left_fixed, right_fixed=self.state.right_fixed,
                constraints=self.state.constraints
            )
        except Exception:
            bend, tors = [], []
            x, shear, moment = [], None, None
        self.freq_panel.update_values(bend, tors)
        self.freq_panel.update_shear_stats(shear, moment)
        self.beam_view.update_view()
        self.shear_view.update_view(x, shear, moment)
        self.status_var.set("Ready")

    def _build_bc(self, model):
        import numpy as np
        bending_fixed = []
        torsion_fixed = []
        if self.state.left_fixed:
            bending_fixed.extend([0, 1])
            torsion_fixed.append(0)
        if self.state.right_fixed:
            last_bend = 2 * (model.elements)
            bending_fixed.extend([last_bend, last_bend + 1])
            torsion_fixed.append(model.elements)

        # Add constraints at arbitrary positions
        x_nodes = np.linspace(0.0, model.length, model.elements + 1)
        for constraint in self.state.constraints:
            # Find nearest node
            node_idx = np.argmin(np.abs(x_nodes - constraint.position))
            if constraint.fix_translation:
                # Fix vertical displacement (w DOF)
                bending_fixed.append(2 * node_idx)
            if constraint.fix_rotation:
                # Fix rotation (θ DOF for bending, φ DOF for torsion)
                bending_fixed.append(2 * node_idx + 1)
                torsion_fixed.append(node_idx)

        return bending_fixed, torsion_fixed

    def _on_boundary_change(self, left: str, right: str):
        self.state.left_fixed = left == "fixed"
        self.state.right_fixed = right == "fixed"
        self.status_var.set("Boundary updated")
        self._refresh_model()

    def _add_constraint_dialog(self):
        length = self.state.length
        pos = simpledialog.askfloat("Position", f"x along beam [0, {length} m]:", parent=self, initialvalue=length / 2)
        if pos is None:
            return
        pos = max(0.0, min(length, pos))

        win = tk.Toplevel(self)
        win.title("Constraint type")
        fix_trans = tk.BooleanVar(value=True)
        fix_rot = tk.BooleanVar(value=False)
        ttk.Checkbutton(win, text="Fix translation (vertical)", variable=fix_trans).pack(anchor="w", padx=8, pady=4)
        ttk.Checkbutton(win, text="Fix rotation (slope)", variable=fix_rot).pack(anchor="w", padx=8, pady=4)
        ttk.Button(win, text="OK", command=win.destroy).pack(pady=6)
        win.transient(self)
        win.grab_set()
        self.wait_window(win)

        self.state.constraints.append(Constraint(position=pos, fix_translation=fix_trans.get(), fix_rotation=fix_rot.get()))
        self.status_var.set("Constraint added")
        self._refresh_model()

    def _open_boundary_dialog(self):
        win = tk.Toplevel(self)
        win.title("Boundary Conditions")
        form = BoundaryForm(win, on_change=lambda _l, _r: None)
        form.pack(fill="both", expand=True, padx=8, pady=8)
        form.set_values(self.state.left_fixed, self.state.right_fixed)
        btns = ttk.Frame(win)
        btns.pack(fill="x", padx=8, pady=4)
        ttk.Button(
            btns,
            text="OK",
            command=lambda: self._apply_boundary_from_form(form, win),
        ).pack(side="right", padx=4)
        ttk.Button(btns, text="Cancel", command=win.destroy).pack(side="right")

    def _apply_boundary_from_form(self, form: BoundaryForm, win: tk.Toplevel):
        left_val = form.left_var.get()
        right_val = form.right_var.get()
        self.state.left_fixed = left_val == "fixed"
        self.state.right_fixed = right_val == "fixed"
        win.destroy()
        self.status_var.set("Boundary conditions updated")
        self._refresh_model()


def run():
    app = BeamApp()
    app.mainloop()


class CrossSectionDialog(tk.Toplevel):
    def __init__(self, master: tk.Widget, length: float, sections):
        super().__init__(master)
        self.title("Cross Sections")
        self.result = None
        self.length = length
        self.sections = list(sections)

        self.start_var = tk.StringVar(value="0.0")
        self.end_var = tk.StringVar(value=f"{length:g}")
        self.radius_var = tk.StringVar(value="0.05")

        form = ttk.Frame(self, padding=8)
        form.pack(fill="both", expand=True)
        ttk.Label(form, text="Start x [m]").grid(row=0, column=0, sticky="w", padx=4, pady=2)
        ttk.Entry(form, textvariable=self.start_var, width=12).grid(row=0, column=1, sticky="ew", padx=4, pady=2)
        ttk.Label(form, text="End x [m]").grid(row=0, column=2, sticky="w", padx=4, pady=2)
        ttk.Entry(form, textvariable=self.end_var, width=12).grid(row=0, column=3, sticky="ew", padx=4, pady=2)
        ttk.Label(form, text="Radius [m]").grid(row=0, column=4, sticky="w", padx=4, pady=2)
        ttk.Entry(form, textvariable=self.radius_var, width=12).grid(row=0, column=5, sticky="ew", padx=4, pady=2)
        ttk.Button(form, text="Add", command=self._add_section).grid(row=0, column=6, padx=4, pady=2)

        self.tree = ttk.Treeview(form, columns=("start", "end", "radius"), show="headings", height=6)
        self.tree.heading("start", text="Start [m]")
        self.tree.heading("end", text="End [m]")
        self.tree.heading("radius", text="Radius [m]")
        self.tree.column("start", width=90, anchor="center")
        self.tree.column("end", width=90, anchor="center")
        self.tree.column("radius", width=90, anchor="center")
        self.tree.grid(row=1, column=0, columnspan=7, sticky="nsew", padx=4, pady=6)

        btns = ttk.Frame(form)
        btns.grid(row=2, column=0, columnspan=7, sticky="ew")
        ttk.Button(btns, text="Remove Selected", command=self._remove_selected).pack(side="left", padx=4)
        ttk.Button(btns, text="Clear", command=self._clear_sections).pack(side="left", padx=4)
        ttk.Button(btns, text="OK", command=self._on_ok).pack(side="right", padx=4)
        ttk.Button(btns, text="Cancel", command=self.destroy).pack(side="right", padx=4)

        for section in self.sections:
            self._insert_row(section)

        for col in range(7):
            form.columnconfigure(col, weight=1)

        self.transient(master)
        self.grab_set()

    def _insert_row(self, section: CrossSection):
        self.tree.insert(
            "", "end",
            values=(f"{section.start:g}", f"{section.end:g}", f"{section.radius:g}")
        )

    def _add_section(self):
        try:
            start = float(self.start_var.get())
            end = float(self.end_var.get())
            radius = float(self.radius_var.get())
        except ValueError:
            return

        start = max(0.0, min(self.length, start))
        end = max(0.0, min(self.length, end))
        if end <= start or radius <= 0:
            return
        section = CrossSection(start=start, end=end, radius=radius)
        self.sections.append(section)
        self._insert_row(section)

    def _remove_selected(self):
        selected = list(self.tree.selection())
        if not selected:
            return
        indices = [self.tree.index(item) for item in selected]
        for item in selected:
            self.tree.delete(item)
        for idx in sorted(indices, reverse=True):
            if 0 <= idx < len(self.sections):
                self.sections.pop(idx)

    def _clear_sections(self):
        for item in self.tree.get_children():
            self.tree.delete(item)
        self.sections = []

    def _on_ok(self):
        self.result = list(self.sections)
        self.destroy()
