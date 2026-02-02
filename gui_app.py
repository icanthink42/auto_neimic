import pickle
import threading
import tkinter as tk
from tkinter import ttk, simpledialog, filedialog, messagebox

from beam_form import BeamForm
from boundary_form import BoundaryForm
from beam_view import BeamView
from constraint import Constraint
from cross_section import CrossSection
from distributed_load import DistributedLoad
from fem import natural_frequencies, natural_frequencies_and_modes
from frequency_panel import FrequencyPanel
from gui_state import BeamUIState
from mode_shape_view import ModeShapeView
from point_mass import PointMass
from torsional_spring import TorsionalSpring
from translational_spring import TranslationalSpring
from shear_moment import shear_moment, GRAVITY
from shear_moment_view import ShearMomentView
from si_prefix import prefix_labels, prefix_multiplier


class BeamApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Auto Niemiec Beam")
        self.state = BeamUIState()
        self.status_var = tk.StringVar(value="Ready")
        self._run_token = 0
        self._active_run_token = None
        self._cancel_requested = False
        self._is_running = False
        self._cancel_event = threading.Event()
        self._progress_var = tk.DoubleVar(value=0.0)
        self._progress_text = tk.StringVar(value="0%")
        self.current_file = None

        # Create menu bar
        menubar = tk.Menu(self)
        self.config(menu=menubar)

        file_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="File", menu=file_menu)
        file_menu.add_command(label="Save", command=self._save_setup)
        file_menu.add_command(label="Save As...", command=self._save_setup_as)
        file_menu.add_command(label="Load", command=self._load_setup)
        file_menu.add_separator()
        file_menu.add_command(label="Exit", command=self.quit)

        toolbar = ttk.Frame(self, padding=6)
        toolbar.grid(row=0, column=0, sticky="ew")

        # File dropdown
        file_btn = ttk.Menubutton(toolbar, text="File")
        file_btn.pack(side="left", padx=4)
        file_dropdown = tk.Menu(file_btn, tearoff=0)
        file_btn.config(menu=file_dropdown)
        file_dropdown.add_command(label="Save", command=self._save_setup)
        file_dropdown.add_command(label="Save As...", command=self._save_setup_as)
        file_dropdown.add_command(label="Load...", command=self._load_setup)

        # Add dropdown
        add_btn = ttk.Menubutton(toolbar, text="Add")
        add_btn.pack(side="left", padx=4)
        add_dropdown = tk.Menu(add_btn, tearoff=0)
        add_btn.config(menu=add_dropdown)
        add_dropdown.add_command(label="Add Mass", command=self._add_mass_dialog)
        add_dropdown.add_command(label="Add Spring", command=self._add_spring_dialog)
        add_dropdown.add_command(label="Add Dist. Mass", command=self._add_distributed_load_dialog)
        add_dropdown.add_command(label="Add Constraint", command=self._add_constraint_dialog)

        ttk.Button(toolbar, text="Beam Params", command=self._open_beam_params).pack(side="left", padx=4)
        ttk.Button(toolbar, text="Cross Sections", command=self._open_cross_section_dialog).pack(side="left", padx=4)
        ttk.Button(toolbar, text="Boundary Conds", command=self._open_boundary_dialog).pack(side="left", padx=4)
        ttk.Label(toolbar, textvariable=self.status_var).pack(side="right")
        self.progress = ttk.Progressbar(
            toolbar,
            mode="determinate",
            length=120,
            maximum=100,
            variable=self._progress_var,
        )
        self.progress_label = ttk.Label(toolbar, textvariable=self._progress_text)
        self.run_button = ttk.Button(toolbar, text="Run", command=self._on_run_click)
        self.run_button.pack(side="right", padx=4)
        self.progress_label.pack(side="right", padx=(0, 4))
        self.progress.pack(side="right", padx=4)
        self.progress.pack_forget()
        self.progress_label.pack_forget()

        # Create a canvas with scrollbar for the entire main area
        main_canvas = tk.Canvas(self, highlightthickness=0)
        main_canvas.grid(row=1, column=0, sticky="nsew")
        self.rowconfigure(1, weight=1)
        self.columnconfigure(0, weight=1)

        main_scrollbar = ttk.Scrollbar(self, orient="vertical", command=main_canvas.yview)
        main_scrollbar.grid(row=1, column=1, sticky="ns")
        main_canvas.configure(yscrollcommand=main_scrollbar.set)

        # Create the main frame inside the canvas
        main = ttk.Frame(main_canvas, padding=8)
        main_window = main_canvas.create_window((0, 0), window=main, anchor="nw")

        # Configure the main frame
        main.columnconfigure(0, weight=1)
        main.rowconfigure(0, weight=0)
        main.rowconfigure(1, weight=0)
        main.rowconfigure(2, weight=0)
        main.rowconfigure(3, weight=0)
        main.rowconfigure(4, weight=0)

        # Set equal heights for the three display rows
        display_height = 300
        main.rowconfigure(1, minsize=display_height)
        main.rowconfigure(2, minsize=display_height)
        main.rowconfigure(3, minsize=display_height)

        top_bar = ttk.Frame(main)
        top_bar.grid(row=0, column=0, sticky="ew", padx=4, pady=2)
        self.boundary_form = BoundaryForm(top_bar, on_change=self._on_boundary_change)
        self.boundary_form.pack(side="left", padx=4)

        self.beam_view = BeamView(main, self.state, on_change=self._mark_dirty)
        self.beam_view.grid(row=1, column=0, sticky="nsew", padx=4, pady=4)

        self.shear_view = ShearMomentView(main)
        self.shear_view.grid(row=2, column=0, sticky="nsew", padx=4, pady=4)

        self.mode_view = ModeShapeView(main)
        self.mode_view.grid(row=3, column=0, sticky="nsew", padx=4, pady=4)

        self.freq_panel = FrequencyPanel(main)
        self.freq_panel.grid(row=4, column=0, sticky="ew", padx=4, pady=4)

        # Update scroll region when content changes
        def _configure_scroll(event):
            main_canvas.configure(scrollregion=main_canvas.bbox("all"))
            # Make sure the canvas window is wide enough
            canvas_width = event.width
            main_canvas.itemconfig(main_window, width=canvas_width)

        main_canvas.bind("<Configure>", _configure_scroll)
        main.bind("<Configure>", lambda e: main_canvas.configure(scrollregion=main_canvas.bbox("all")))

        # Enable mousewheel scrolling
        def _on_mousewheel(event):
            main_canvas.yview_scroll(int(-1 * (event.delta / 120)), "units")

        main_canvas.bind_all("<MouseWheel>", _on_mousewheel)  # Windows/MacOS
        main_canvas.bind_all("<Button-4>", lambda e: main_canvas.yview_scroll(-1, "units"))  # Linux scroll up
        main_canvas.bind_all("<Button-5>", lambda e: main_canvas.yview_scroll(1, "units"))  # Linux scroll down

        self.status_var.set("Ready (press Run)")

    def _collect_support_positions(self, model, left_fixed, right_fixed, constraints):
        supports = []
        if left_fixed:
            supports.append(0.0)
        if right_fixed:
            supports.append(model.length)
        for constraint in constraints:
            if constraint.fix_translation and not any(abs(s - constraint.position) < 1e-9 for s in supports):
                supports.append(constraint.position)
        supports.sort()
        return supports

    def _estimate_progress_steps(self, model, supports, include_frequencies, shear_points):
        steps = 0
        if include_frequencies:
            steps += model.elements * 2
            steps += len(model.point_masses) * 2
            steps += len(model.trans_springs)
            steps += len(model.tors_springs)
            steps += 2
        if not supports:
            steps += 1
        else:
            steps += shear_points
            if len(supports) == 2:
                steps += 1
            elif len(supports) > 2:
                steps += len(supports) + len(supports) * len(supports) + 1
        return max(1, steps)

    def _build_progress_tracker(self, total_steps, token):
        class _ProgressTracker:
            def __init__(self, total, update):
                self.total = max(1, total)
                self.completed = 0
                self.update = update

            def step(self, count=1):
                self.completed += count
                percent = min(100.0, (self.completed / self.total) * 100.0)
                self.update(percent)

        return _ProgressTracker(total_steps, lambda percent: self._queue_progress(token, percent))

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
        self.status_var.set("Mass added (press Run)")
        self._mark_dirty()

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
        self.status_var.set("Spring added (press Run)")
        self._mark_dirty()

    def _add_distributed_load_dialog(self):
        length = self.state.length
        start = simpledialog.askfloat("Start position", f"Start x along beam [0, {length} m]:", parent=self, initialvalue=0.0)
        if start is None:
            return
        start = max(0.0, min(length, start))

        end = simpledialog.askfloat("End position", f"End x along beam [{start}, {length} m]:", parent=self, initialvalue=length)
        if end is None:
            return
        end = max(start, min(length, end))

        if end <= start:
            return

        mass = simpledialog.askfloat("Add distributed mass", "Mass per length [kg/m]:", parent=self, initialvalue=100.0)
        if mass is None:
            return

        self.state.distributed_loads.append(DistributedLoad(start=start, end=end, mass_per_length=mass))
        self.status_var.set("Distributed mass added (press Run)")
        self._mark_dirty()

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
        self.status_var.set("Beam updated (press Run)")
        self._mark_dirty()

    def _open_cross_section_dialog(self):
        dialog = CrossSectionDialog(self, self.state.length, self.state.cross_sections)
        self.wait_window(dialog)
        if dialog.result is None:
            return
        self.state.cross_sections = dialog.result
        self.status_var.set("Cross sections updated (press Run)")
        self._mark_dirty()

    def _refresh_model(self, bend, tors, x, shear, moment, bend_modes=None, tors_modes=None, bending_fixed=None, torsion_fixed=None):
        self.freq_panel.update_values(bend, tors)
        self.freq_panel.update_shear_stats(shear, moment)
        self.beam_view.update_view()
        self.shear_view.update_view(x, shear, moment)

        # Update mode shape view if we have mode data
        if bend_modes is not None and tors_modes is not None:
            model = self.state.to_model()
            n_nodes = model.elements + 1
            self.mode_view.update_modes(
                beam_length=model.length,
                n_nodes=n_nodes,
                bend_freqs=bend,
                tors_freqs=tors,
                bend_modes=bend_modes,
                tors_modes=tors_modes,
                bending_fixed=bending_fixed or [],
                torsion_fixed=torsion_fixed or [],
            )

        self.status_var.set("Ready")

    def _build_bc(self, model, left_fixed, right_fixed, constraints):
        import numpy as np
        bending_fixed = []
        torsion_fixed = []
        if left_fixed:
            bending_fixed.extend([0, 1])
            torsion_fixed.append(0)
        if right_fixed:
            last_bend = 2 * (model.elements)
            bending_fixed.extend([last_bend, last_bend + 1])
            torsion_fixed.append(model.elements)

        # Add constraints at arbitrary positions
        x_nodes = np.linspace(0.0, model.length, model.elements + 1)
        for constraint in constraints:
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

    def _mark_dirty(self):
        if not self._is_running:
            self.beam_view.update_view()
            if "press Run" not in self.status_var.get():
                self.status_var.set("Changes pending (press Run)")

    def _on_run_click(self):
        if self._is_running:
            self._cancel_run()
            return
        self._start_run()

    def _start_run(self):
        if self._is_running:
            return
        self._is_running = True
        self._cancel_requested = False
        self._cancel_event = threading.Event()
        self._run_token += 1
        token = self._run_token
        self._active_run_token = token
        self.status_var.set("Running...")
        self.run_button.config(text="Cancel")
        self._set_progress(0)
        self.progress.pack(side="right", padx=4)
        self.progress_label.pack(side="right", padx=(0, 4))
        snapshot = {
            "model": self.state.to_model(),
            "left_fixed": self.state.left_fixed,
            "right_fixed": self.state.right_fixed,
            "constraints": list(self.state.constraints),
        }
        thread = threading.Thread(target=self._run_model_thread, args=(token, snapshot), daemon=True)
        thread.start()

    def _cancel_run(self):
        if not self._is_running:
            return
        self._cancel_requested = True
        self._cancel_event.set()
        self.status_var.set("Canceling...")

    def _run_model_thread(self, token, snapshot):
        bend, tors = [], []
        bend_modes, tors_modes = None, None
        bending_fixed_out, torsion_fixed_out = [], []
        x, shear, moment = [], None, None
        error = None
        try:
            model = snapshot["model"]
            bending_fixed, torsion_fixed = self._build_bc(
                model,
                snapshot["left_fixed"],
                snapshot["right_fixed"],
                snapshot["constraints"],
            )
            supports = self._collect_support_positions(
                model,
                snapshot["left_fixed"],
                snapshot["right_fixed"],
                snapshot["constraints"],
            )
            include_frequencies = bool(bending_fixed or torsion_fixed)
            total_steps = self._estimate_progress_steps(model, supports, include_frequencies, 200)
            tracker = self._build_progress_tracker(total_steps, token)
            progress_step = tracker.step
            if include_frequencies:
                bend, tors, bend_modes, tors_modes, bending_fixed_out, torsion_fixed_out = natural_frequencies_and_modes(
                    model,
                    bending_fixed=bending_fixed,
                    torsion_fixed=torsion_fixed,
                    n_modes=5,
                    progress=progress_step,
                    cancel_event=self._cancel_event,
                )
            x, shear, moment = shear_moment(
                model,
                g=GRAVITY,
                left_fixed=snapshot["left_fixed"],
                right_fixed=snapshot["right_fixed"],
                constraints=snapshot["constraints"],
                n=200,
                progress=progress_step,
                cancel_event=self._cancel_event,
            )
        except Exception as exc:
            error = str(exc)
        self.after(0, lambda: self._on_run_complete(token, bend, tors, x, shear, moment, error, bend_modes, tors_modes, bending_fixed_out, torsion_fixed_out))

    def _on_run_complete(self, token, bend, tors, x, shear, moment, error, bend_modes=None, tors_modes=None, bending_fixed=None, torsion_fixed=None):
        if token != self._active_run_token:
            return
        self._is_running = False
        self.progress.pack_forget()
        self.progress_label.pack_forget()
        self.run_button.config(text="Run")
        self._set_progress(0)
        if self._cancel_requested or error == "Canceled":
            self.status_var.set("Canceled")
            self._cancel_requested = False
            return
        if error:
            self.status_var.set("Run failed")
            self._cancel_requested = False
            return
        self._cancel_requested = False
        self._refresh_model(bend, tors, x, shear, moment, bend_modes, tors_modes, bending_fixed, torsion_fixed)

    def _queue_progress(self, token, value: float):
        self.after(0, lambda: self._set_progress_for_token(token, value))

    def _set_progress_for_token(self, token, value: float):
        if token != self._active_run_token:
            return
        self._set_progress(value)

    def _set_progress(self, value: float):
        safe_value = max(0.0, min(100.0, value))
        self._progress_var.set(safe_value)
        self._progress_text.set(f"{int(round(safe_value))}%")

    def _on_boundary_change(self, left: str, right: str):
        self.state.left_fixed = left == "fixed"
        self.state.right_fixed = right == "fixed"
        self.status_var.set("Boundary updated (press Run)")
        self._mark_dirty()

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
        self.status_var.set("Constraint added (press Run)")
        self._mark_dirty()

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
        self.status_var.set("Boundary conditions updated (press Run)")
        self._mark_dirty()

    def _save_setup(self):
        """Save the current setup to the current file, or prompt for a filename if none exists."""
        if self.current_file is None:
            self._save_setup_as()
        else:
            self._save_to_file(self.current_file)

    def _save_setup_as(self):
        """Prompt for a filename and save the current setup."""
        filename = filedialog.asksaveasfilename(
            defaultextension=".niemiec",
            filetypes=[("Niemiec files", "*.niemiec"), ("All files", "*.*")],
            title="Save Setup As"
        )
        if filename:
            self._save_to_file(filename)
            self.current_file = filename

    def _save_to_file(self, filename: str):
        """Save the current state to a file using pickle."""
        try:
            # Create a dictionary with all the state data
            save_data = {
                'length': self.state.length,
                'radius': self.state.radius,
                'elastic_modulus': self.state.elastic_modulus,
                'shear_modulus': self.state.shear_modulus,
                'density': self.state.density,
                'elements': self.state.elements,
                'cross_sections': self.state.cross_sections,
                'point_masses': self.state.point_masses,
                'trans_springs': self.state.trans_springs,
                'tors_springs': self.state.tors_springs,
                'distributed_loads': self.state.distributed_loads,
                'constraints': self.state.constraints,
                'left_fixed': self.state.left_fixed,
                'right_fixed': self.state.right_fixed,
            }

            with open(filename, 'wb') as f:
                pickle.dump(save_data, f)

            self.status_var.set(f"Saved to {filename}")
            messagebox.showinfo("Success", f"Setup saved successfully to {filename}")
        except Exception as e:
            self.status_var.set("Save failed")
            messagebox.showerror("Error", f"Failed to save setup: {str(e)}")

    def _load_setup(self):
        """Prompt for a filename and load a setup."""
        filename = filedialog.askopenfilename(
            defaultextension=".niemiec",
            filetypes=[("Niemiec files", "*.niemiec"), ("All files", "*.*")],
            title="Load Setup"
        )
        if filename:
            self._load_from_file(filename)
            self.current_file = filename

    def _load_from_file(self, filename: str):
        """Load state from a file using pickle."""
        try:
            with open(filename, 'rb') as f:
                save_data = pickle.load(f)

            # Restore all state data
            self.state.length = save_data.get('length', 2.0)
            self.state.radius = save_data.get('radius', 0.05)
            self.state.elastic_modulus = save_data.get('elastic_modulus', 200e9)
            self.state.shear_modulus = save_data.get('shear_modulus', 79.3e9)
            self.state.density = save_data.get('density', 7850)
            self.state.elements = save_data.get('elements', 20)
            self.state.cross_sections = save_data.get('cross_sections', [])
            self.state.point_masses = save_data.get('point_masses', [])
            self.state.trans_springs = save_data.get('trans_springs', [])
            self.state.tors_springs = save_data.get('tors_springs', [])
            self.state.distributed_loads = save_data.get('distributed_loads', [])
            self.state.constraints = save_data.get('constraints', [])
            self.state.left_fixed = save_data.get('left_fixed', True)
            self.state.right_fixed = save_data.get('right_fixed', False)

            # Update the boundary form to reflect the loaded state
            self.boundary_form.set_values(self.state.left_fixed, self.state.right_fixed)

            # Mark as dirty to trigger a refresh
            self._mark_dirty()
            self.status_var.set(f"Loaded from {filename} (press Run)")
            messagebox.showinfo("Success", f"Setup loaded successfully from {filename}")
        except Exception as e:
            self.status_var.set("Load failed")
            messagebox.showerror("Error", f"Failed to load setup: {str(e)}")


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
        self.start_prefix_var = tk.StringVar(value="m")
        self.end_prefix_var = tk.StringVar(value="m")
        self.radius_prefix_var = tk.StringVar(value="m")

        form = ttk.Frame(self, padding=8)
        form.pack(fill="both", expand=True)
        ttk.Label(form, text="Start x [m]").grid(row=0, column=0, sticky="w", padx=4, pady=2)
        ttk.Entry(form, textvariable=self.start_var, width=12).grid(row=0, column=1, sticky="ew", padx=4, pady=2)
        ttk.OptionMenu(
            form,
            self.start_prefix_var,
            self.start_prefix_var.get(),
            *prefix_labels("m"),
        ).grid(row=0, column=2, sticky="w", padx=4, pady=2)
        ttk.Label(form, text="End x [m]").grid(row=0, column=3, sticky="w", padx=4, pady=2)
        ttk.Entry(form, textvariable=self.end_var, width=12).grid(row=0, column=4, sticky="ew", padx=4, pady=2)
        ttk.OptionMenu(
            form,
            self.end_prefix_var,
            self.end_prefix_var.get(),
            *prefix_labels("m"),
        ).grid(row=0, column=5, sticky="w", padx=4, pady=2)
        ttk.Label(form, text="Radius [m]").grid(row=0, column=6, sticky="w", padx=4, pady=2)
        ttk.Entry(form, textvariable=self.radius_var, width=12).grid(row=0, column=7, sticky="ew", padx=4, pady=2)
        ttk.OptionMenu(
            form,
            self.radius_prefix_var,
            self.radius_prefix_var.get(),
            *prefix_labels("m"),
        ).grid(row=0, column=8, sticky="w", padx=4, pady=2)
        ttk.Button(form, text="Add", command=self._add_section).grid(row=0, column=9, padx=4, pady=2)

        self.tree = ttk.Treeview(form, columns=("start", "end", "radius"), show="headings", height=6)
        self.tree.heading("start", text="Start [m]")
        self.tree.heading("end", text="End [m]")
        self.tree.heading("radius", text="Radius [m]")
        self.tree.column("start", width=90, anchor="center")
        self.tree.column("end", width=90, anchor="center")
        self.tree.column("radius", width=90, anchor="center")
        self.tree.grid(row=1, column=0, columnspan=10, sticky="nsew", padx=4, pady=6)

        btns = ttk.Frame(form)
        btns.grid(row=2, column=0, columnspan=10, sticky="ew")
        ttk.Button(btns, text="Remove Selected", command=self._remove_selected).pack(side="left", padx=4)
        ttk.Button(btns, text="Clear", command=self._clear_sections).pack(side="left", padx=4)
        ttk.Button(btns, text="OK", command=self._on_ok).pack(side="right", padx=4)
        ttk.Button(btns, text="Cancel", command=self.destroy).pack(side="right", padx=4)

        for section in self.sections:
            self._insert_row(section)

        for col in range(10):
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
            start = float(self.start_var.get()) * prefix_multiplier(self.start_prefix_var.get(), "m")
            end = float(self.end_var.get()) * prefix_multiplier(self.end_prefix_var.get(), "m")
            radius = float(self.radius_var.get()) * prefix_multiplier(self.radius_prefix_var.get(), "m")
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
